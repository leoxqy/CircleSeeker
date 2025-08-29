#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_split_regions_final_v10.py
环状DNA候选检测脚本 - SR-events修复版（基于v9最小改动）

变更要点（v10）：
- SR-events 只统计 SR reads（has SA 或 is_supplementary）的对齐条数，限定在候选区间内
- sr_events_per_read = SR-events / unique_split_count（唯一SR-read数）
- 环状方向保持 v9 的 L-R 判定（左端L/B，右端R/B）
- 保持表头与v9一致，新增统计逻辑不影响其他列语义

示例：
python filter_split_regions_final_v10.py \
  -i step11.win50.per-base.bed \
  -b step10.bam \
  -o step12.split_regions_refined.tsv \
  -q 20 -s 20 --breakpoint-distance 50 \
  --min-breakpoint-support 1 \
  --min-bridge-reads 1 \
  --min-interval-size 100 \
  --require-sa \
  --normalize-ids \
  -v
"""

import sys, re, argparse, logging, shutil, subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional
from collections import defaultdict

import pysam

# ---------------------- Logging ----------------------
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("split_regions_v10")

# ---------------------- Dataclasses ----------------------
@dataclass
class Region:
    chrom: str
    start: int  # 0-based half-open
    end: int
    depth_sum: float = 0.0

@dataclass
class BPCluster:
    chrom: str
    pos: int
    start_pos: int
    end_pos: int
    support_reads: Set[str]
    left_clip_reads: Set[str]   # reads with left softclip at this position
    right_clip_reads: Set[str]  # reads with right softclip at this position
    near_start_boundary: bool
    near_end_boundary: bool
    
    @property
    def clip_pattern(self) -> str:
        """返回软剪切模式：L（左）、R（右）、B（双向）、N（无）"""
        has_left = len(self.left_clip_reads) > 0
        has_right = len(self.right_clip_reads) > 0
        if has_left and has_right:
            return "B"
        elif has_left:
            return "L"
        elif has_right:
            return "R"
        else:
            return "N"

@dataclass
class CandidateInterval:
    chrom: str
    start: int
    end: int
    start_bp: int
    end_bp: int
    start_support: int
    end_support: int
    bridge_reads: int
    read_count: int
    unique_read_count: int
    split_read_count: int
    unique_split_count: int
    sa_unique_count: int
    sc_support_reads: Set[str]
    all_reads: Set[str]
    split_reads: Set[str]
    sa_reads: Set[str]
    # 软剪切模式
    start_clip_pattern: str = "N"
    end_clip_pattern: str = "N"
    start_left_clips: int = 0
    start_right_clips: int = 0
    end_left_clips: int = 0
    end_right_clips: int = 0
    # 环状模式（L-R）
    is_circular_pattern: bool = False
    # SR事件统计（仅SR reads）
    sr_events: int = 0               # 区间内 SR 对齐条数总计
    sr_events_per_read: float = 0.0  # SR 对齐条数 / 唯一SR-read数
    # 其他统计
    mean_mapq_support: float = 0.0
    mean_mapq_region: float = 0.0
    kind: str = "normal"
    depth_sum: float = 0.0

# ---------------------- Helpers ----------------------
def check_requirements():
    if not shutil.which("bedtools"):
        logger.warning("bedtools 未找到，将无法计算 depth_sum。建议安装：conda install -c bioconda bedtools")
        return False
    return True

def parse_region_string(s: str) -> Region:
    if ":" in s and "-" in s:
        chrom, rest = s.split(":", 1)
        start, end = rest.replace(",", "").split("-", 1)
        return Region(chrom=chrom, start=int(start), end=int(end))
    parts = re.split(r"[\s,]+", s.strip())
    if len(parts) >= 3:
        return Region(chrom=parts[0], start=int(parts[1]), end=int(parts[2]))
    raise ValueError(f"无法解析 region: {s}")

def _normalize_read_id(qname: str, normalize: bool) -> str:
    rid = qname.split()[0]
    if normalize and (rid.endswith("/1") or rid.endswith("/2")):
        rid = rid[:-2]
    return rid

def cigar_softclips(aln) -> Tuple[int, int]:
    if aln.cigartuples is None: 
        return (0, 0)
    left = aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0
    right = aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0
    return (left, right)

# ---------------------- Core ----------------------
class SplitRegionRefiner:
    def __init__(
        self,
        min_depth: int = 1,
        merge_distance: int = 5,
        min_mapq: int = 0,
        min_softclip: int = 20,
        breakpoint_distance: int = 50,
        min_breakpoint_support: int = 1,
        min_bridge_reads: int = 1,
        min_interval_size: int = 50,  # 最小区间大小
        max_nm: Optional[int] = None,
        min_as: Optional[int] = None,
        normalize_ids: bool = False,
        require_sa: bool = True,
        require_bridge: bool = True,
        require_circular: bool = False,  # 是否要求L-R模式
        verbose: bool = False
    ):
        self.min_depth = min_depth
        self.merge_distance = merge_distance
        self.min_mapq = min_mapq
        self.min_softclip = min_softclip
        self.breakpoint_distance = breakpoint_distance
        self.min_breakpoint_support = min_breakpoint_support
        self.min_bridge_reads = min_bridge_reads
        self.min_interval_size = min_interval_size
        self.max_nm = max_nm
        self.min_as = min_as
        self.normalize_ids = normalize_ids
        self.require_sa = require_sa
        self.require_bridge = require_bridge
        self.require_circular = require_circular
        self.has_bedtools = check_requirements()

        if verbose: 
            logger.setLevel(logging.DEBUG)

    def _rid(self, qname: str) -> str:
        return _normalize_read_id(qname, self.normalize_ids)

    # ---------- region acquisition ----------
    def merge_depth_regions(self, per_base_bed: Path) -> List[Region]:
        if not self.has_bedtools:
            logger.warning("bedtools 不可用，跳过深度区域合并")
            return []
            
        tmp = per_base_bed.parent
        f1 = tmp / "tmp.perbase.filtered.bed"
        f2 = tmp / "tmp.perbase.sorted.bed"
        f3 = tmp / "tmp.perbase.merged.bed"
        
        try:
            # 过滤低深度区域
            with per_base_bed.open() as fin, f1.open("w") as fout:
                for line in fin:
                    if not line.strip(): 
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4: 
                        continue
                    try: 
                        dep = float(p[3])
                    except: 
                        continue
                    if dep >= self.min_depth: 
                        fout.write(line)
            
            # 排序
            with f2.open("w") as fout:
                subprocess.run(["bedtools", "sort", "-i", str(f1)], stdout=fout, check=True)
            
            # 合并
            with f3.open("w") as fout:
                subprocess.run([
                    "bedtools", "merge", "-i", str(f2), 
                    "-d", str(self.merge_distance),
                    "-c", "4", "-o", "sum"
                ], stdout=fout, check=True)
            
            # 读取结果
            regions: List[Region] = []
            with f3.open() as fin:
                for line in fin:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4: 
                        continue
                    regions.append(Region(
                        chrom=p[0], 
                        start=int(p[1]), 
                        end=int(p[2]), 
                        depth_sum=float(p[3])
                    ))
            
            # 确保排序
            regions.sort(key=lambda r: (r.chrom, r.start, r.end))
            logger.info(f"合并得到 {len(regions)} 个候选区间")
            return regions
            
        finally:
            for f in (f1, f2, f3):
                try:
                    if f.exists(): 
                        f.unlink()
                except: 
                    pass

    # ---------- breakpoint discovery ----------
    def collect_breakpoints(
        self,
        bamf: pysam.AlignmentFile,
        region: Region
    ) -> Tuple[List[BPCluster], Set[str], int, int, Dict[str, int]]:
        """
        收集断点信息，记录软剪切方向
        返回：断点簇列表、区间内split reads（唯一）、读数统计、以及（DEBUG用）每read的SR对齐条数
        """
        margin = self.breakpoint_distance
        fetch_start = max(0, region.start - margin)
        fetch_end = region.end + margin

        # 分别记录左/右软剪切
        pos2reads_left: Dict[int, Set[str]] = {}
        pos2reads_right: Dict[int, Set[str]] = {}
        pos2reads_all: Dict[int, Set[str]] = {}
        
        # DEBUG：SR 对齐条数（仅 SR reads）
        read_sr_events: Dict[str, int] = defaultdict(int)
        
        region_all_reads: Set[str] = set()
        region_split_reads: Set[str] = set()
        region_read_count = 0

        for aln in bamf.fetch(region.chrom, fetch_start, fetch_end):
            if aln.is_unmapped: 
                continue
            if aln.mapping_quality < self.min_mapq: 
                continue
            if self.max_nm is not None and aln.has_tag("NM") and aln.get_tag("NM") > self.max_nm: 
                continue
            if self.min_as is not None and aln.has_tag("AS") and aln.get_tag("AS") < self.min_as: 
                continue

            rid = self._rid(aln.query_name)

            # 区间内统计（与 SR 判定分开）
            if not (aln.reference_end <= region.start or aln.reference_start >= region.end):
                region_read_count += 1
                region_all_reads.add(rid)
                if aln.has_tag("SA") or aln.is_supplementary:
                    region_split_reads.add(rid)

            # SR 对齐条数（DEBUG用，不参与最终 sr_events 计算）
            if aln.is_supplementary or aln.has_tag("SA"):
                # 此处不强制要求软剪切；只要是 SR 对齐就记 1 条
                read_sr_events[rid] += 1

            # 收集软剪切断点（用于端点模式/方向判定）
            ls, rs = cigar_softclips(aln)
            if ls >= self.min_softclip:
                bp = aln.reference_start
                if fetch_start <= bp <= fetch_end:
                    pos2reads_left.setdefault(bp, set()).add(rid)
                    pos2reads_all.setdefault(bp, set()).add(rid)
            if rs >= self.min_softclip:
                bp = aln.reference_end
                if fetch_start <= bp <= fetch_end:
                    pos2reads_right.setdefault(bp, set()).add(rid)
                    pos2reads_all.setdefault(bp, set()).add(rid)

        # 调试输出
        if logger.level <= logging.DEBUG and read_sr_events:
            logger.debug(f"Region {region.chrom}:{region.start}-{region.end} SR-read 及其对齐条数（前5条）：")
            for rid, cnt in list(sorted(read_sr_events.items()))[:5]:
                logger.debug(f"  {rid}: {cnt} SR alignments")

        region_unique_reads = len(region_all_reads)
        
        # 聚类断点
        clusters = self._cluster_breakpoints(
            region, pos2reads_all, pos2reads_left, pos2reads_right
        )
        
        return clusters, region_split_reads, region_read_count, region_unique_reads, dict(read_sr_events)

    def _cluster_breakpoints(
        self, 
        region: Region, 
        pos2reads_all: Dict[int, Set[str]],
        pos2reads_left: Dict[int, Set[str]],
        pos2reads_right: Dict[int, Set[str]]
    ) -> List[BPCluster]:
        """聚类断点，记录软剪切方向"""
        if not pos2reads_all: 
            return []
            
        bpd = self.breakpoint_distance
        sorted_pos = sorted(pos2reads_all.keys())
        
        # 找出聚类范围
        idx_ranges: List[Tuple[int,int]] = []
        cur = 0
        for i in range(1, len(sorted_pos)):
            if sorted_pos[i] - sorted_pos[i-1] > bpd:
                idx_ranges.append((cur, i-1))
                cur = i
        idx_ranges.append((cur, len(sorted_pos)-1))

        clusters: List[BPCluster] = []
        for i0, i1 in idx_ranges:
            sub = sorted_pos[i0:i1+1]
            
            # 计算加权中心
            weights = [len(pos2reads_all[p]) for p in sub]
            den = sum(weights)
            num = sum(p*w for p,w in zip(sub, weights))
            center = int(round(num/den)) if den else sub[0]
            
            # 合并支持reads
            reads_all: Set[str] = set()
            reads_left: Set[str] = set()
            reads_right: Set[str] = set()
            
            for p in sub:
                reads_all |= pos2reads_all[p]
                if p in pos2reads_left:
                    reads_left |= pos2reads_left[p]
                if p in pos2reads_right:
                    reads_right |= pos2reads_right[p]
            
            clusters.append(BPCluster(
                chrom=region.chrom,
                pos=center,
                start_pos=sub[0],
                end_pos=sub[-1],
                support_reads=reads_all,
                left_clip_reads=reads_left,
                right_clip_reads=reads_right,
                near_start_boundary=(abs(center - region.start) <= bpd),
                near_end_boundary=(abs(center - region.end) <= bpd)
            ))
        
        # 确保排序
        clusters.sort(key=lambda c: c.pos)
        logger.debug(f"Region {region.chrom}:{region.start}-{region.end} 聚类得到 {len(clusters)} 个断点簇")
        return clusters

    # ---------- derive candidate intervals ----------
    def derive_candidates(
        self,
        region: Region,
        clusters: List[BPCluster],
        region_split_reads: Set[str],
        read_sr_events: Dict[str, int],  # 仅用于DEBUG，不参与最终 sr_events 计算
        bamf: pysam.AlignmentFile
    ) -> List[CandidateInterval]:
        """派生候选区间，过滤掉起止断点相同或过近的情况"""
        clusters = sorted(clusters, key=lambda c: c.pos)
        bpd = self.breakpoint_distance
        min_sup = self.min_breakpoint_support

        # 分类断点
        start_cands = [c for c in clusters if c.near_start_boundary]
        end_cands = [c for c in clusters if c.near_end_boundary]
        internals = [c for c in clusters if not (c.near_start_boundary or c.near_end_boundary)]

        def pick_best(cands: List[BPCluster], is_start: bool) -> Optional[BPCluster]:
            if not cands: 
                return None
            ref_pos = region.start if is_start else region.end
            # 确定性排序
            return sorted(cands, key=lambda x: (
                -len(x.support_reads),
                abs(x.pos - ref_pos),
                x.pos,
                min(sorted(x.support_reads)) if x.support_reads else ""
            ))[0]

        start_bp = pick_best(start_cands, True)
        end_bp = pick_best(end_cands, False)

        # 构建锚点序列
        anchors: List[BPCluster] = []
        if start_bp and len(start_bp.support_reads) >= min_sup: 
            anchors.append(start_bp)
        for c in internals:
            if len(c.support_reads) >= min_sup: 
                anchors.append(c)
        if end_bp and len(end_bp.support_reads) >= min_sup: 
            anchors.append(end_bp)

        # 去重并排序
        pos_seen = set()
        uniq_anchors: List[BPCluster] = []
        for c in sorted(anchors, key=lambda x: x.pos):
            if c.pos not in pos_seen:
                uniq_anchors.append(c)
                pos_seen.add(c.pos)
        anchors = uniq_anchors

        cands: List[CandidateInterval] = []
        made_normal = False

        # normal类型：边界到边界
        if start_bp and end_bp and \
           len(start_bp.support_reads) >= min_sup and \
           len(end_bp.support_reads) >= min_sup:
            # 关键检查：起止断点不能相同或太近
            if abs(end_bp.pos - start_bp.pos) >= self.min_interval_size:
                cands.append(self._build_candidate(
                    region, start_bp, end_bp, region_split_reads, read_sr_events, bamf, "normal"
                ))
                made_normal = True
            else:
                logger.debug(f"跳过过小的候选区间: {start_bp.pos}-{end_bp.pos} (长度={end_bp.pos-start_bp.pos})")

        # subregion类型：相邻锚点
        for i in range(len(anchors) - 1):
            L, R = anchors[i], anchors[i+1]
            # 避免重复
            if made_normal and start_bp and end_bp and \
               L.pos == start_bp.pos and R.pos == end_bp.pos:
                continue
            # 检查间距
            if R.pos - L.pos < self.min_interval_size:
                logger.debug(f"跳过过小的子区间: {L.pos}-{R.pos} (长度={R.pos-L.pos})")
                continue
            cands.append(self._build_candidate(
                region, L, R, region_split_reads, read_sr_events, bamf, "subregion"
            ))

        # 去重（同坐标只留最佳）
        kind_rank = {"normal": 2, "subregion": 1}
        def score(ci: CandidateInterval):
            return (
                kind_rank.get(ci.kind, 0),
                ci.start_support + ci.end_support,
                ci.bridge_reads,
                ci.sa_unique_count,
                1 if ci.is_circular_pattern else 0  # 优先环状模式
            )
        
        chosen: Dict[Tuple[str,int,int], CandidateInterval] = {}
        for c in sorted(cands, key=lambda x: (x.chrom, x.start, x.end, score(x))):
            key = (c.chrom, c.start, c.end)
            if key not in chosen or score(c) > score(chosen[key]):
                chosen[key] = c
        
        return sorted(chosen.values(), key=lambda x: (x.chrom, x.start, x.end))

    def _build_candidate(
        self,
        region: Region,
        left: BPCluster,
        right: BPCluster,
        region_split_reads: Set[str],
        read_sr_events: Dict[str, int],  # 未使用，仅为接口兼容
        bamf: pysam.AlignmentFile,
        kind: str
    ) -> CandidateInterval:
        """构建候选区间，判定环状模式（L-R），并按 SR reads 计算 SR-events"""
        s = max(region.start, min(left.pos, right.pos))
        e = min(region.end, max(left.pos, right.pos))
        if e <= s: 
            e = s + 1

        left_reads = set(left.support_reads)
        right_reads = set(right.support_reads)
        sc_union = left_reads | right_reads
        bridge = len(left_reads & right_reads)

        # 环状模式（L-R）：左端 L/B，右端 R/B
        is_circular = (left.clip_pattern in ["L", "B"]) and (right.clip_pattern in ["R", "B"])
        
        if logger.level <= logging.DEBUG:
            logger.debug(f"候选区间 {region.chrom}:{s}-{e}:")
            logger.debug(f"  左断点@{left.pos}: 模式={left.clip_pattern}, L={len(left.left_clip_reads)}, R={len(left.right_clip_reads)}")
            logger.debug(f"  右断点@{right.pos}: 模式={right.clip_pattern}, L={len(right.left_clip_reads)}, R={len(right.right_clip_reads)}")
            logger.debug(f"  环状判定(L-R): {is_circular}")

        # —— 区间统计（含 SR 对齐条数）——
        (read_count, uniq_reads, split_cnt, uniq_split,
         all_reads_set, split_reads_set, mean_mapq_region,
         sa_unique_cnt, sa_reads_set, sr_event_count) = self._region_stats(bamf, region.chrom, s, e)

        # 支持reads的平均 MAPQ（对 sc_union）
        mapqs_support = []
        for aln in bamf.fetch(region.chrom, s, e):
            if aln.is_unmapped or aln.mapping_quality < self.min_mapq:
                continue
            rid = self._rid(aln.query_name)
            if rid in sc_union:
                mapqs_support.append(int(aln.mapping_quality))
        mean_mapq_support = (sum(mapqs_support) / float(len(mapqs_support))) if mapqs_support else 0.0

        # —— SR 事件统计：只用 _region_stats 的 sr_event_count（仅 SR reads 的对齐条数）——
        sr_events_total = sr_event_count
        sr_events_per_read = (sr_events_total / uniq_split) if uniq_split > 0 else 0.0

        return CandidateInterval(
            chrom=region.chrom,
            start=s,
            end=e,
            start_bp=left.pos,
            end_bp=right.pos,
            start_support=len(left_reads),
            end_support=len(right_reads),
            bridge_reads=bridge,
            read_count=read_count,
            unique_read_count=uniq_reads,
            split_read_count=split_cnt,
            unique_split_count=uniq_split,
            sa_unique_count=sa_unique_cnt,
            sc_support_reads=sc_union,
            all_reads=all_reads_set,
            split_reads=split_reads_set,
            sa_reads=sa_reads_set,
            # 软剪切模式
            start_clip_pattern=left.clip_pattern,
            end_clip_pattern=right.clip_pattern,
            start_left_clips=len(left.left_clip_reads),
            start_right_clips=len(left.right_clip_reads),
            end_left_clips=len(right.left_clip_reads),
            end_right_clips=len(right.right_clip_reads),
            # 环状模式
            is_circular_pattern=is_circular,
            # SR 事件
            sr_events=sr_events_total,
            sr_events_per_read=sr_events_per_read,
            # 质量统计
            mean_mapq_support=mean_mapq_support,
            mean_mapq_region=mean_mapq_region,
            kind=kind,
            depth_sum=0.0
        )

    def _region_stats(
        self, 
        bamf: pysam.AlignmentFile, 
        chrom: str, 
        start: int, 
        end: int
    ) -> Tuple[int,int,int,int,Set[str],Set[str],float,int,Set[str],int]:
        """计算区间统计
        返回：
          read_count                - 对齐记录条数（所有）
          uniq_reads                - 唯一 read 数（所有）
          split_cnt                 - 唯一 SR-read 数（保持与既有语义一致）
          uniq_split                - 同上（保持两列一致，避免下游破坏）
          all_reads                 - 所有唯一 read
          split_reads               - 唯一 SR-read 集合
          mean_mapq                 - 区间平均 MAPQ
          sa_unique_cnt             - 唯一 SA-read 数
          sa_reads                  - 唯一 SA-read 集合
          sr_event_count            - **SR 对齐条数总计（仅 SR reads：is_supplementary 或 has SA）**
        """
        read_count = 0
        all_reads: Set[str] = set()
        split_reads: Set[str] = set()
        sa_reads: Set[str] = set()
        mapqs: List[int] = []

        # SR 对齐条数（仅 SR reads）
        sr_event_count = 0
        
        for aln in bamf.fetch(chrom, start, end):
            if aln.is_unmapped: 
                continue
            if aln.mapping_quality < self.min_mapq: 
                continue
            if self.max_nm is not None and aln.has_tag("NM") and aln.get_tag("NM") > self.max_nm: 
                continue
            if self.min_as is not None and aln.has_tag("AS") and aln.get_tag("AS") < self.min_as: 
                continue
            if aln.reference_end <= start or aln.reference_start >= end: 
                continue
                
            read_count += 1
            rid = self._rid(aln.query_name)
            all_reads.add(rid)
            mapqs.append(int(aln.mapping_quality))
            
            is_sr_alignment = (aln.is_supplementary or aln.has_tag("SA"))
            if is_sr_alignment:
                split_reads.add(rid)     # 唯一 SR-read（集）
                sr_event_count += 1      # 该 SR 对齐计 1 条事件

            if aln.has_tag("SA"):
                sa_reads.add(rid)
                
        uniq_reads = len(all_reads)
        uniq_split = len(split_reads)
        # 历史兼容：split_cnt 与 uniq_split 保持一致（若想改为“SR对齐条数”，把下一行改成 sr_event_count）
        split_cnt = len(split_reads)
        sa_unique_cnt = len(sa_reads)
        mean_mapq = (sum(mapqs)/float(len(mapqs))) if mapqs else 0.0
        
        return (read_count, uniq_reads, split_cnt, uniq_split,
                all_reads, split_reads, mean_mapq, sa_unique_cnt, sa_reads,
                sr_event_count)

    # ---------- depth map ----------
    def fill_depth_sum_for_candidates(self, per_base_bed: Path, cands: List[CandidateInterval]) -> None:
        """使用bedtools计算深度总和"""
        if not cands or not self.has_bedtools: 
            return
            
        tmp = per_base_bed.parent
        fa_u = tmp / "tmp.cands.a.bed"
        fa_s = tmp / "tmp.cands.a.sorted.bed"
        fb_s = tmp / "tmp.perbase.sorted.bed"
        fm = tmp / "tmp.cands.mapped.bed"
        
        # 写入候选区间
        with fa_u.open("w") as fa:
            for idx, c in enumerate(sorted(cands, key=lambda x: (x.chrom, x.start, x.end))):
                fa.write(f"{c.chrom}\t{c.start}\t{c.end}\tcand_{idx}\n")
        
        try:
            # 排序
            with fa_s.open("w") as f: 
                subprocess.run(["bedtools", "sort", "-i", str(fa_u)], stdout=f, check=True)
            with fb_s.open("w") as f: 
                subprocess.run(["bedtools", "sort", "-i", str(per_base_bed)], stdout=f, check=True)
            
            # 映射
            with fm.open("w") as f:
                subprocess.run([
                    "bedtools", "map", "-a", str(fa_s), "-b", str(fb_s),
                    "-c", "4", "-o", "sum", "-sorted"
                ], stdout=f, check=True)
            
            # 读取结果
            id2sum: Dict[str,float] = {}
            with fm.open() as fin:
                for line in fin:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5: 
                        continue
                    dep = 0.0 if p[4] in (".", "") else float(p[4])
                    id2sum[p[3]] = dep
            
            # 回填
            for idx, c in enumerate(sorted(cands, key=lambda x: (x.chrom, x.start, x.end))):
                c.depth_sum = float(id2sum.get(f"cand_{idx}", 0.0))
                
        finally:
            for f in (fa_u, fa_s, fb_s, fm):
                try:
                    if f.exists(): 
                        f.unlink()
                except: 
                    pass

    # ---------- run ----------
    def run(
        self, 
        per_base_bed: Path, 
        bam_path: Path, 
        output_tsv: Path, 
        direct_regions: Optional[List[Region]] = None
    ) -> Path:
        """主运行函数"""
        if not per_base_bed.exists(): 
            raise FileNotFoundError(f"per-base bed 不存在: {per_base_bed}")
        if not bam_path.exists(): 
            raise FileNotFoundError(f"BAM 不存在: {bam_path}")
            
        output_tsv.parent.mkdir(parents=True, exist_ok=True)

        # 获取区间
        if direct_regions:
            regions = sorted(direct_regions, key=lambda r: (r.chrom, r.start, r.end))
            logger.info(f"使用用户提供的 {len(regions)} 个区间")
        else:
            regions = self.merge_depth_regions(per_base_bed)

        bamf = pysam.AlignmentFile(str(bam_path), "rb")

        all_candidates: List[CandidateInterval] = []
        for idx, reg in enumerate(regions, 1):
            clusters, region_split_reads, _, _, read_sr_events = self.collect_breakpoints(bamf, reg)
            cands = self.derive_candidates(reg, clusters, region_split_reads, read_sr_events, bamf)
            all_candidates.extend(cands)
            
            if idx % 100 == 0: 
                logger.info(f"已处理 {idx}/{len(regions)} 个区间")
                
        bamf.close()

        # 过滤
        logger.info(f"开始过滤 {len(all_candidates)} 个候选区间...")
        prefiltered: List[CandidateInterval] = []
        
        stats = {
            "total": len(all_candidates),
            "pass_size": 0,
            "pass_sa": 0,
            "pass_bridge": 0,
            "pass_support": 0,
            "pass_circular": 0,
            "final": 0
        }
        
        for c in all_candidates:
            # 大小过滤
            if (c.end - c.start) < self.min_interval_size:
                continue
            stats["pass_size"] += 1
            
            # SA过滤
            if self.require_sa and c.sa_unique_count < 1:
                continue
            stats["pass_sa"] += 1
            
            # 桥接读数过滤
            if self.require_bridge and c.bridge_reads < self.min_bridge_reads:
                continue
            stats["pass_bridge"] += 1
            
            # 断点支持过滤
            if c.start_support < self.min_breakpoint_support or \
               c.end_support < self.min_breakpoint_support:
                continue
            stats["pass_support"] += 1
            
            # 环状模式过滤（L-R模式）
            if self.require_circular and not c.is_circular_pattern:
                continue
            stats["pass_circular"] += 1
            
            prefiltered.append(c)
        
        stats["final"] = len(prefiltered)
        
        # 打印过滤统计
        logger.info("过滤统计:")
        logger.info(f"  初始候选: {stats['total']}")
        logger.info(f"  大小合格(>={self.min_interval_size}bp): {stats['pass_size']}")
        if self.require_sa:
            logger.info(f"  含SA标签: {stats['pass_sa']}")
        if self.require_bridge:
            logger.info(f"  满足桥接读数(>={self.min_bridge_reads}): {stats['pass_bridge']}")
        logger.info(f"  满足断点支持(>={self.min_breakpoint_support}): {stats['pass_support']}")
        if self.require_circular:
            logger.info(f"  L-R环状模式: {stats['pass_circular']}")
        logger.info(f"  最终输出: {stats['final']}")

        # 统计环状模式
        circular_count = sum(1 for c in prefiltered if c.is_circular_pattern)
        logger.info(f"环状模式统计: {circular_count}/{len(prefiltered)} 为L-R模式")

        # 计算深度
        self.fill_depth_sum_for_candidates(per_base_bed, prefiltered)
        
        # 输出
        self.write_output(prefiltered, output_tsv)
        logger.info(f"结果已保存至: {output_tsv}")
        
        return output_tsv

    def write_output(self, cands: List[CandidateInterval], out_tsv: Path) -> None:
        """输出结果"""
        headers = [
            "chrom", "start", "end", "length",
            "start_bp", "end_bp",
            "start_support", "end_support",
            "bridge_reads",
            "read_count", "unique_read_count",
            "split_read_count", "unique_split_count",
            "sa_unique_count",
            # 软剪切信息
            "start_clip_pattern", "end_clip_pattern",
            "start_left_clips", "start_right_clips",
            "end_left_clips", "end_right_clips",
            # 环状模式
            "is_circular_pattern",
            # SR事件（仅SR reads）
            "sr_events", "sr_events_per_read",
            # 质量和深度
            "mean_mapq_support", "mean_mapq_region",
            "depth_sum",
            "kind",
            # reads列表
            "support_read_ids",
            "all_read_ids",
            "split_read_ids",
            "sa_read_ids"
        ]
        
        with out_tsv.open("w") as f:
            f.write("\t".join(headers) + "\n")
            
            # 确保输出顺序一致
            for c in sorted(cands, key=lambda x: (x.chrom, x.start, x.end, x.kind)):
                row = [
                    c.chrom,
                    str(c.start),
                    str(c.end),
                    str(c.end - c.start),
                    str(c.start_bp),
                    str(c.end_bp),
                    str(c.start_support),
                    str(c.end_support),
                    str(c.bridge_reads),
                    str(c.read_count),
                    str(c.unique_read_count),
                    str(c.split_read_count),
                    str(c.unique_split_count),
                    str(c.sa_unique_count),
                    # 软剪切模式
                    c.start_clip_pattern,
                    c.end_clip_pattern,
                    str(c.start_left_clips),
                    str(c.start_right_clips),
                    str(c.end_left_clips),
                    str(c.end_right_clips),
                    # 环状模式
                    "Yes" if c.is_circular_pattern else "No",
                    # SR事件（仅SR reads）
                    str(c.sr_events),
                    f"{c.sr_events_per_read:.2f}",
                    # 质量
                    f"{c.mean_mapq_support:.3f}",
                    f"{c.mean_mapq_region:.3f}",
                    f"{c.depth_sum:.1f}",
                    c.kind,
                    # reads
                    ";".join(sorted(c.sc_support_reads)),
                    ";".join(sorted(c.all_reads)),
                    ";".join(sorted(c.split_reads)),
                    ";".join(sorted(c.sa_reads)),
                ]
                f.write("\t".join(row) + "\n")

# ---------------------- CLI ----------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="环状DNA候选检测脚本 v10 - 修复SR-events（仅SR reads），保持L-R判定"
    )
    p.add_argument("-i", "--input-bed", type=Path, required=True,
                   help="mosdepth per-base BED（chrom start end depth）")
    p.add_argument("-b", "--bam", type=Path, required=True,
                   help="坐标排序且已索引的 BAM")
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="输出 TSV")
    p.add_argument("-d", "--merge-distance", type=int, default=5,
                   help="合并小区间的最大间距（默认5）")
    p.add_argument("-m", "--min-depth", type=int, default=1,
                   help="最小深度阈值（默认1）")
    p.add_argument("-q", "--min-mapq", type=int, default=0,
                   help="最小 MAPQ（默认0）")
    p.add_argument("-s", "--min-softclip", type=int, default=20,
                   help="最小软剪切长度（默认20）")
    p.add_argument("--breakpoint-distance", type=int, default=50,
                   help="断点聚类距离（默认50）")
    p.add_argument("--min-breakpoint-support", type=int, default=1,
                   help="断点最小支持读数（默认1）")
    p.add_argument("--min-bridge-reads", type=int, default=1,
                   help="最小桥接读数（默认1）")
    p.add_argument("--min-interval-size", type=int, default=50,
                   help="候选区间最小长度（默认50bp）")
    p.add_argument("--max-nm", type=int, default=None,
                   help="最大编辑距离 NM（可选）")
    p.add_argument("--min-as", type=int, default=None,
                   help="最小比对分数 AS（可选）")
    p.add_argument("--normalize-ids", action="store_true",
                   help="去除 read ID 的 /1 /2 后缀")
    
    # 过滤开关
    p.add_argument("--no-require-sa", dest="require_sa", action="store_false",
                   help="不要求SA标签")
    p.add_argument("--no-require-bridge", dest="require_bridge", action="store_false",
                   help="不要求桥接读数阈值")
    p.add_argument("--require-circular", action="store_true",
                   help="只输出L-R环状模式的候选")
    p.set_defaults(require_sa=True, require_bridge=True)
    
    # 直接指定区间
    p.add_argument("--region", action="append", default=[],
                   help="chr:start-end 或 'chr start end'（可多次）")
    p.add_argument("--regions-bed", type=Path, default=None,
                   help="三列 BED（chrom start end）")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="DEBUG 日志")
    
    return p

def load_direct_regions(args) -> Optional[List[Region]]:
    """加载用户指定的区间"""
    regs: List[Region] = []
    if args.region:
        for s in args.region:
            regs.append(parse_region_string(s))
    if args.regions_bed and args.regions_bed.exists():
        with args.regions_bed.open() as fin:
            for line in fin:
                if not line.strip() or line.startswith("#"): 
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 3: 
                    continue
                regs.append(Region(chrom=p[0], start=int(p[1]), end=int(p[2])))
    return regs if regs else None

def main():
    args = build_argparser().parse_args()
    if args.verbose: 
        logger.setLevel(logging.DEBUG)
        
    # 检查输入文件
    if not args.input_bed.exists(): 
        logger.error(f"输入 BED 不存在：{args.input_bed}")
        sys.exit(1)
    if not args.bam.exists(): 
        logger.error(f"BAM 不存在：{args.bam}")
        sys.exit(1)
        
    direct_regions = load_direct_regions(args)
    
    try:
        refiner = SplitRegionRefiner(
            min_depth=args.min_depth,
            merge_distance=args.merge_distance,
            min_mapq=args.min_mapq,
            min_softclip=args.min_softclip,
            breakpoint_distance=args.breakpoint_distance,
            min_breakpoint_support=args.min_breakpoint_support,
            min_bridge_reads=args.min_bridge_reads,
            min_interval_size=args.min_interval_size,
            max_nm=args.max_nm,
            min_as=args.min_as,
            normalize_ids=args.normalize_ids,
            require_sa=args.require_sa,
            require_bridge=args.require_bridge,
            require_circular=args.require_circular,
            verbose=args.verbose
        )
        
        refiner.run(
            per_base_bed=args.input_bed,
            bam_path=args.bam,
            output_tsv=args.output,
            direct_regions=direct_regions
        )
        
        sys.exit(0)
        
    except Exception as e:
        logger.exception(f"执行失败：{e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

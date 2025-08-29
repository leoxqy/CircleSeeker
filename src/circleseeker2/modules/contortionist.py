"""
Contortionist - Split Regions Detection and Analysis

This module detects circular DNA candidates by analyzing split regions in alignment data.
It processes depth information and BAM files to identify breakpoints and bridge reads
that indicate circular structures.

Key features:
- SR-events counting: Only counts SR reads (has SA or is_supplementary) alignments
- sr_events_per_read = SR-events / unique_split_count
- Circular pattern detection: L-R pattern (left-end L/B, right-end R/B)
- Maintains compatibility with v10 headers

Migrated from step12_Contortionist.py to the new CircleSeeker2 architecture.
"""

from __future__ import annotations

import sys
import re
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict, Set, Optional
from collections import defaultdict

import pysam

from circleseeker2.exceptions import PipelineError


# ---------------------- Dataclasses ----------------------
@dataclass
class Region:
    """Represents a genomic region with depth information."""
    chrom: str
    start: int  # 0-based half-open
    end: int
    depth_sum: float = 0.0


@dataclass
class BPCluster:
    """Represents a breakpoint cluster with clip direction tracking."""
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
        """Return soft-clip pattern: L(left), R(right), B(bidirectional), N(none)"""
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
    """Represents a candidate circular DNA interval."""
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
    # Soft-clip patterns
    start_clip_pattern: str = "N"
    end_clip_pattern: str = "N"
    start_left_clips: int = 0
    start_right_clips: int = 0
    end_left_clips: int = 0
    end_right_clips: int = 0
    # Circular pattern (L-R)
    is_circular_pattern: bool = False
    # SR event statistics (SR reads only)
    sr_events: int = 0               # Total SR alignments in interval
    sr_events_per_read: float = 0.0  # SR alignments / unique SR-reads
    # Other statistics
    mean_mapq_support: float = 0.0
    mean_mapq_region: float = 0.0
    kind: str = "normal"
    depth_sum: float = 0.0


# ---------------------- Helpers ----------------------
def check_requirements():
    """Check if required tools are available."""
    if not shutil.which("bedtools"):
        return False
    return True


def parse_region_string(s: str) -> Region:
    """Parse region string in various formats."""
    if ":" in s and "-" in s:
        chrom, rest = s.split(":", 1)
        start, end = rest.replace(",", "").split("-", 1)
        return Region(chrom=chrom, start=int(start), end=int(end))
    parts = re.split(r"[\s,]+", s.strip())
    if len(parts) >= 3:
        return Region(chrom=parts[0], start=int(parts[1]), end=int(parts[2]))
    raise ValueError(f"Cannot parse region: {s}")


def cigar_softclips(aln) -> Tuple[int, int]:
    """Extract left and right soft-clip sizes from alignment."""
    if aln.cigartuples is None: 
        return (0, 0)
    left = aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0
    right = aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0
    return (left, right)


# ---------------------- Main Class ----------------------
class Contortionist:
    """Split regions detector for circular DNA analysis."""
    
    def __init__(
        self,
        min_depth: int = 1,
        merge_distance: int = 5,
        min_mapq: int = 0,
        min_softclip: int = 20,
        breakpoint_distance: int = 50,
        min_breakpoint_support: int = 1,
        min_bridge_reads: int = 1,
        min_interval_size: int = 50,
        max_nm: Optional[int] = None,
        min_as: Optional[int] = None,
        normalize_ids: bool = False,
        require_sa: bool = True,
        require_bridge: bool = True,
        require_circular: bool = False,
        logger: Optional[logging.Logger] = None
    ):
        """Initialize Contortionist with parameters."""
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
        self.logger = logger or logging.getLogger(self.__class__.__name__)
    
    def _normalize_read_id(self, qname: str) -> str:
        """Normalize read ID by removing /1 or /2 suffix if requested."""
        rid = qname.split()[0]
        if self.normalize_ids and (rid.endswith("/1") or rid.endswith("/2")):
            rid = rid[:-2]
        return rid
    
    def _rid(self, qname: str) -> str:
        """Shorthand for normalize_read_id."""
        return self._normalize_read_id(qname)
    
    # ---------- Region Acquisition ----------
    def merge_depth_regions(self, per_base_bed: Path) -> List[Region]:
        """Merge depth regions using bedtools."""
        if not self.has_bedtools:
            self.logger.warning("bedtools not available, skipping depth region merge")
            return []
            
        tmp = per_base_bed.parent
        f1 = tmp / "tmp.perbase.filtered.bed"
        f2 = tmp / "tmp.perbase.sorted.bed"
        f3 = tmp / "tmp.perbase.merged.bed"
        
        try:
            # Filter low-depth regions
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
            
            # Sort
            with f2.open("w") as fout:
                subprocess.run(["bedtools", "sort", "-i", str(f1)], stdout=fout, check=True)
            
            # Merge
            with f3.open("w") as fout:
                subprocess.run([
                    "bedtools", "merge", "-i", str(f2), 
                    "-d", str(self.merge_distance),
                    "-c", "4", "-o", "sum"
                ], stdout=fout, check=True)
            
            # Read results
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
            
            # Ensure sorted
            regions.sort(key=lambda r: (r.chrom, r.start, r.end))
            self.logger.info(f"Merged to {len(regions)} candidate intervals")
            return regions
            
        finally:
            for f in (f1, f2, f3):
                try:
                    if f.exists(): 
                        f.unlink()
                except: 
                    pass
    
    # ---------- Breakpoint Discovery ----------
    def collect_breakpoints(
        self,
        bamf: pysam.AlignmentFile,
        region: Region
    ) -> Tuple[List[BPCluster], Set[str], int, int, Dict[str, int]]:
        """
        Collect breakpoint information, record soft-clip directions.
        Returns: breakpoint clusters, unique split reads in region, read count stats, 
                 and (for DEBUG) per-read SR alignment counts
        """
        margin = self.breakpoint_distance
        fetch_start = max(0, region.start - margin)
        fetch_end = region.end + margin

        # Separate left/right soft-clips
        pos2reads_left: Dict[int, Set[str]] = {}
        pos2reads_right: Dict[int, Set[str]] = {}
        pos2reads_all: Dict[int, Set[str]] = {}
        
        # DEBUG: SR alignment counts (SR reads only)
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

            # In-region statistics (separate from SR determination)
            if not (aln.reference_end <= region.start or aln.reference_start >= region.end):
                region_read_count += 1
                region_all_reads.add(rid)
                if aln.has_tag("SA") or aln.is_supplementary:
                    region_split_reads.add(rid)

            # SR alignment counts (for DEBUG, not used in final sr_events calculation)
            if aln.is_supplementary or aln.has_tag("SA"):
                # Count each SR alignment as 1 event
                read_sr_events[rid] += 1

            # Collect soft-clip breakpoints (for endpoint pattern/direction determination)
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

        # Debug output
        if self.logger.level <= logging.DEBUG and read_sr_events:
            self.logger.debug(f"Region {region.chrom}:{region.start}-{region.end} SR-read alignments (first 5):")
            for rid, cnt in list(sorted(read_sr_events.items()))[:5]:
                self.logger.debug(f"  {rid}: {cnt} SR alignments")

        region_unique_reads = len(region_all_reads)
        
        # Cluster breakpoints
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
        """Cluster breakpoints, recording soft-clip directions."""
        if not pos2reads_all: 
            return []
            
        bpd = self.breakpoint_distance
        sorted_pos = sorted(pos2reads_all.keys())
        
        # Find cluster ranges
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
            
            # Calculate weighted center
            weights = [len(pos2reads_all[p]) for p in sub]
            den = sum(weights)
            num = sum(p*w for p,w in zip(sub, weights))
            center = int(round(num/den)) if den else sub[0]
            
            # Merge support reads
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
        
        # Ensure sorted
        clusters.sort(key=lambda c: c.pos)
        self.logger.debug(f"Region {region.chrom}:{region.start}-{region.end} clustered to {len(clusters)} breakpoints")
        return clusters
    
    # ---------- Derive Candidate Intervals ----------
    def derive_candidates(
        self,
        region: Region,
        clusters: List[BPCluster],
        region_split_reads: Set[str],
        read_sr_events: Dict[str, int],  # Only for DEBUG, not used in final sr_events calculation
        bamf: pysam.AlignmentFile
    ) -> List[CandidateInterval]:
        """Derive candidate intervals, filtering out same or too-close start/end breakpoints."""
        clusters = sorted(clusters, key=lambda c: c.pos)
        bpd = self.breakpoint_distance
        min_sup = self.min_breakpoint_support

        # Categorize breakpoints
        start_cands = [c for c in clusters if c.near_start_boundary]
        end_cands = [c for c in clusters if c.near_end_boundary]
        internals = [c for c in clusters if not (c.near_start_boundary or c.near_end_boundary)]

        def pick_best(cands: List[BPCluster], is_start: bool) -> Optional[BPCluster]:
            if not cands: 
                return None
            ref_pos = region.start if is_start else region.end
            # Deterministic sorting
            return sorted(cands, key=lambda x: (
                -len(x.support_reads),
                abs(x.pos - ref_pos),
                x.pos,
                min(sorted(x.support_reads)) if x.support_reads else ""
            ))[0]

        start_bp = pick_best(start_cands, True)
        end_bp = pick_best(end_cands, False)

        # Build anchor sequence
        anchors: List[BPCluster] = []
        if start_bp and len(start_bp.support_reads) >= min_sup: 
            anchors.append(start_bp)
        for c in internals:
            if len(c.support_reads) >= min_sup: 
                anchors.append(c)
        if end_bp and len(end_bp.support_reads) >= min_sup: 
            anchors.append(end_bp)

        # Deduplicate and sort
        pos_seen = set()
        uniq_anchors: List[BPCluster] = []
        for c in sorted(anchors, key=lambda x: x.pos):
            if c.pos not in pos_seen:
                uniq_anchors.append(c)
                pos_seen.add(c.pos)
        anchors = uniq_anchors

        cands: List[CandidateInterval] = []
        made_normal = False

        # Normal type: boundary to boundary
        if start_bp and end_bp and \
           len(start_bp.support_reads) >= min_sup and \
           len(end_bp.support_reads) >= min_sup:
            # Key check: start and end breakpoints cannot be same or too close
            if abs(end_bp.pos - start_bp.pos) >= self.min_interval_size:
                cands.append(self._build_candidate(
                    region, start_bp, end_bp, region_split_reads, read_sr_events, bamf, "normal"
                ))
                made_normal = True
            else:
                self.logger.debug(f"Skipping too-small candidate interval: {start_bp.pos}-{end_bp.pos} (length={end_bp.pos-start_bp.pos})")

        # Subregion type: adjacent anchors
        for i in range(len(anchors) - 1):
            L, R = anchors[i], anchors[i+1]
            # Avoid duplicates
            if made_normal and start_bp and end_bp and \
               L.pos == start_bp.pos and R.pos == end_bp.pos:
                continue
            # Check distance
            if R.pos - L.pos < self.min_interval_size:
                self.logger.debug(f"Skipping too-small subregion: {L.pos}-{R.pos} (length={R.pos-L.pos})")
                continue
            cands.append(self._build_candidate(
                region, L, R, region_split_reads, read_sr_events, bamf, "subregion"
            ))

        # Deduplicate (keep best for same coordinates)
        kind_rank = {"normal": 2, "subregion": 1}
        def score(ci: CandidateInterval):
            return (
                kind_rank.get(ci.kind, 0),
                ci.start_support + ci.end_support,
                ci.bridge_reads,
                ci.sa_unique_count,
                1 if ci.is_circular_pattern else 0  # Prioritize circular pattern
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
        read_sr_events: Dict[str, int],  # Not used, only for interface compatibility
        bamf: pysam.AlignmentFile,
        kind: str
    ) -> CandidateInterval:
        """Build candidate interval, determine circular pattern (L-R), and count SR-events by SR reads."""
        s = max(region.start, min(left.pos, right.pos))
        e = min(region.end, max(left.pos, right.pos))
        if e <= s: 
            e = s + 1

        left_reads = set(left.support_reads)
        right_reads = set(right.support_reads)
        sc_union = left_reads | right_reads
        bridge = len(left_reads & right_reads)

        # Circular pattern (L-R): left end L/B, right end R/B
        is_circular = (left.clip_pattern in ["L", "B"]) and (right.clip_pattern in ["R", "B"])
        
        if self.logger.level <= logging.DEBUG:
            self.logger.debug(f"Candidate interval {region.chrom}:{s}-{e}:")
            self.logger.debug(f"  Left breakpoint@{left.pos}: pattern={left.clip_pattern}, L={len(left.left_clip_reads)}, R={len(left.right_clip_reads)}")
            self.logger.debug(f"  Right breakpoint@{right.pos}: pattern={right.clip_pattern}, L={len(right.left_clip_reads)}, R={len(right.right_clip_reads)}")
            self.logger.debug(f"  Circular determination (L-R): {is_circular}")

        # Region statistics (including SR alignment counts)
        (read_count, uniq_reads, split_cnt, uniq_split,
         all_reads_set, split_reads_set, mean_mapq_region,
         sa_unique_cnt, sa_reads_set, sr_event_count) = self._region_stats(bamf, region.chrom, s, e)

        # Support reads' average MAPQ (for sc_union)
        mapqs_support = []
        for aln in bamf.fetch(region.chrom, s, e):
            if aln.is_unmapped or aln.mapping_quality < self.min_mapq:
                continue
            rid = self._rid(aln.query_name)
            if rid in sc_union:
                mapqs_support.append(int(aln.mapping_quality))
        mean_mapq_support = (sum(mapqs_support) / float(len(mapqs_support))) if mapqs_support else 0.0

        # SR event statistics: use sr_event_count from _region_stats (SR reads' alignment counts)
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
            # Soft-clip patterns
            start_clip_pattern=left.clip_pattern,
            end_clip_pattern=right.clip_pattern,
            start_left_clips=len(left.left_clip_reads),
            start_right_clips=len(left.right_clip_reads),
            end_left_clips=len(right.left_clip_reads),
            end_right_clips=len(right.right_clip_reads),
            # Circular pattern
            is_circular_pattern=is_circular,
            # SR events
            sr_events=sr_events_total,
            sr_events_per_read=sr_events_per_read,
            # Quality statistics
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
        """Calculate region statistics.
        Returns:
          read_count                - Alignment record count (all)
          uniq_reads                - Unique read count (all)
          split_cnt                 - Unique SR-read count (maintain compatibility)
          uniq_split                - Same as above (keep both columns consistent)
          all_reads                 - All unique reads
          split_reads               - Unique SR-read set
          mean_mapq                 - Region average MAPQ
          sa_unique_cnt             - Unique SA-read count
          sa_reads                  - Unique SA-read set
          sr_event_count            - **SR alignment count total (SR reads only: is_supplementary or has SA)**
        """
        read_count = 0
        all_reads: Set[str] = set()
        split_reads: Set[str] = set()
        sa_reads: Set[str] = set()
        mapqs: List[int] = []

        # SR alignment count (SR reads only)
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
                split_reads.add(rid)     # Unique SR-read (set)
                sr_event_count += 1      # This SR alignment counts as 1 event

            if aln.has_tag("SA"):
                sa_reads.add(rid)
                
        uniq_reads = len(all_reads)
        uniq_split = len(split_reads)
        # Historical compatibility: split_cnt equals uniq_split
        split_cnt = len(split_reads)
        sa_unique_cnt = len(sa_reads)
        mean_mapq = (sum(mapqs)/float(len(mapqs))) if mapqs else 0.0
        
        return (read_count, uniq_reads, split_cnt, uniq_split,
                all_reads, split_reads, mean_mapq, sa_unique_cnt, sa_reads,
                sr_event_count)
    
    # ---------- Depth Map ----------
    def fill_depth_sum_for_candidates(self, per_base_bed: Path, cands: List[CandidateInterval]) -> None:
        """Calculate depth sum using bedtools."""
        if not cands or not self.has_bedtools: 
            return
            
        tmp = per_base_bed.parent
        fa_u = tmp / "tmp.cands.a.bed"
        fa_s = tmp / "tmp.cands.a.sorted.bed"
        fb_s = tmp / "tmp.perbase.sorted.bed"
        fm = tmp / "tmp.cands.mapped.bed"
        
        # Write candidate intervals
        with fa_u.open("w") as fa:
            for idx, c in enumerate(sorted(cands, key=lambda x: (x.chrom, x.start, x.end))):
                fa.write(f"{c.chrom}\t{c.start}\t{c.end}\tcand_{idx}\n")
        
        try:
            # Sort
            with fa_s.open("w") as f: 
                subprocess.run(["bedtools", "sort", "-i", str(fa_u)], stdout=f, check=True)
            with fb_s.open("w") as f: 
                subprocess.run(["bedtools", "sort", "-i", str(per_base_bed)], stdout=f, check=True)
            
            # Map
            with fm.open("w") as f:
                subprocess.run([
                    "bedtools", "map", "-a", str(fa_s), "-b", str(fb_s),
                    "-c", "4", "-o", "sum", "-sorted"
                ], stdout=f, check=True)
            
            # Read results
            id2sum: Dict[str,float] = {}
            with fm.open() as fin:
                for line in fin:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5: 
                        continue
                    dep = 0.0 if p[4] in (".", "") else float(p[4])
                    id2sum[p[3]] = dep
            
            # Fill back
            for idx, c in enumerate(sorted(cands, key=lambda x: (x.chrom, x.start, x.end))):
                c.depth_sum = float(id2sum.get(f"cand_{idx}", 0.0))
                
        finally:
            for f in (fa_u, fa_s, fb_s, fm):
                try:
                    if f.exists(): 
                        f.unlink()
                except: 
                    pass
    
    # ---------- Run ----------
    def run(
        self, 
        per_base_bed: Path, 
        bam_path: Path, 
        output_tsv: Path, 
        direct_regions: Optional[List[Region]] = None
    ) -> Path:
        """Main run function."""
        if not per_base_bed.exists(): 
            raise FileNotFoundError(f"per-base bed does not exist: {per_base_bed}")
        if not bam_path.exists(): 
            raise FileNotFoundError(f"BAM does not exist: {bam_path}")
            
        output_tsv.parent.mkdir(parents=True, exist_ok=True)

        # Get regions
        if direct_regions:
            regions = sorted(direct_regions, key=lambda r: (r.chrom, r.start, r.end))
            self.logger.info(f"Using user-provided {len(regions)} regions")
        else:
            regions = self.merge_depth_regions(per_base_bed)

        bamf = pysam.AlignmentFile(str(bam_path), "rb")

        all_candidates: List[CandidateInterval] = []
        for idx, reg in enumerate(regions, 1):
            clusters, region_split_reads, _, _, read_sr_events = self.collect_breakpoints(bamf, reg)
            cands = self.derive_candidates(reg, clusters, region_split_reads, read_sr_events, bamf)
            all_candidates.extend(cands)
            
            if idx % 100 == 0: 
                self.logger.info(f"Processed {idx}/{len(regions)} regions")
                
        bamf.close()

        # Filter
        self.logger.info(f"Starting to filter {len(all_candidates)} candidate intervals...")
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
            # Size filter
            if (c.end - c.start) < self.min_interval_size:
                continue
            stats["pass_size"] += 1
            
            # SA filter
            if self.require_sa and c.sa_unique_count < 1:
                continue
            stats["pass_sa"] += 1
            
            # Bridge reads filter
            if self.require_bridge and c.bridge_reads < self.min_bridge_reads:
                continue
            stats["pass_bridge"] += 1
            
            # Breakpoint support filter
            if c.start_support < self.min_breakpoint_support or \
               c.end_support < self.min_breakpoint_support:
                continue
            stats["pass_support"] += 1
            
            # Circular pattern filter (L-R pattern)
            if self.require_circular and not c.is_circular_pattern:
                continue
            stats["pass_circular"] += 1
            
            prefiltered.append(c)
        
        stats["final"] = len(prefiltered)
        
        # Print filter statistics
        self.logger.info("Filter statistics:")
        self.logger.info(f"  Initial candidates: {stats['total']}")
        self.logger.info(f"  Size qualified (>={self.min_interval_size}bp): {stats['pass_size']}")
        if self.require_sa:
            self.logger.info(f"  Contains SA tag: {stats['pass_sa']}")
        if self.require_bridge:
            self.logger.info(f"  Meets bridge reads (>={self.min_bridge_reads}): {stats['pass_bridge']}")
        self.logger.info(f"  Meets breakpoint support (>={self.min_breakpoint_support}): {stats['pass_support']}")
        if self.require_circular:
            self.logger.info(f"  L-R circular pattern: {stats['pass_circular']}")
        self.logger.info(f"  Final output: {stats['final']}")

        # Count circular patterns
        circular_count = sum(1 for c in prefiltered if c.is_circular_pattern)
        self.logger.info(f"Circular pattern statistics: {circular_count}/{len(prefiltered)} are L-R pattern")

        # Calculate depth
        self.fill_depth_sum_for_candidates(per_base_bed, prefiltered)
        
        # Output
        self.write_output(prefiltered, output_tsv)
        self.logger.info(f"Results saved to: {output_tsv}")
        
        return output_tsv
    
    def write_output(self, cands: List[CandidateInterval], out_tsv: Path) -> None:
        """Write output results."""
        headers = [
            "chrom", "start", "end", "length",
            "start_bp", "end_bp",
            "start_support", "end_support",
            "bridge_reads",
            "read_count", "unique_read_count",
            "split_read_count", "unique_split_count",
            "sa_unique_count",
            # Soft-clip info
            "start_clip_pattern", "end_clip_pattern",
            "start_left_clips", "start_right_clips",
            "end_left_clips", "end_right_clips",
            # Circular pattern
            "is_circular_pattern",
            # SR events (SR reads only)
            "sr_events", "sr_events_per_read",
            # Quality and depth
            "mean_mapq_support", "mean_mapq_region",
            "depth_sum",
            "kind",
            # Read lists
            "support_read_ids",
            "all_read_ids",
            "split_read_ids",
            "sa_read_ids"
        ]
        
        with out_tsv.open("w") as f:
            f.write("\t".join(headers) + "\n")
            
            # Ensure consistent output order
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
                    # Soft-clip patterns
                    c.start_clip_pattern,
                    c.end_clip_pattern,
                    str(c.start_left_clips),
                    str(c.start_right_clips),
                    str(c.end_left_clips),
                    str(c.end_right_clips),
                    # Circular pattern
                    "Yes" if c.is_circular_pattern else "No",
                    # SR events (SR reads only)
                    str(c.sr_events),
                    f"{c.sr_events_per_read:.2f}",
                    # Quality
                    f"{c.mean_mapq_support:.3f}",
                    f"{c.mean_mapq_region:.3f}",
                    f"{c.depth_sum:.1f}",
                    c.kind,
                    # Reads
                    ";".join(sorted(c.sc_support_reads)),
                    ";".join(sorted(c.all_reads)),
                    ";".join(sorted(c.split_reads)),
                    ";".join(sorted(c.sa_reads)),
                ]
                f.write("\t".join(row) + "\n")
    
    # ---------- Public Interface for Module ----------
    def load_direct_regions(self, region_strings: List[str] = None, regions_bed: Path = None) -> Optional[List[Region]]:
        """Load user-specified regions from strings or BED file."""
        regs: List[Region] = []
        
        if region_strings:
            for s in region_strings:
                regs.append(parse_region_string(s))
        
        if regions_bed and regions_bed.exists():
            with regions_bed.open() as fin:
                for line in fin:
                    if not line.strip() or line.startswith("#"): 
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 3: 
                        continue
                    regs.append(Region(chrom=p[0], start=int(p[1]), end=int(p[2])))
        
        return regs if regs else None
    
    def run_analysis(
        self,
        per_base_bed: Path,
        bam_path: Path,
        output_tsv: Path,
        region_strings: List[str] = None,
        regions_bed: Path = None
    ) -> Path:
        """Run complete analysis pipeline."""
        self.logger.info("=" * 60)
        self.logger.info("Contortionist - Circular DNA Candidate Detection")
        self.logger.info("=" * 60)
        
        # Load direct regions if provided
        direct_regions = self.load_direct_regions(region_strings, regions_bed)
        
        # Run main pipeline
        result = self.run(per_base_bed, bam_path, output_tsv, direct_regions)
        
        self.logger.info("=" * 60)
        
        return result
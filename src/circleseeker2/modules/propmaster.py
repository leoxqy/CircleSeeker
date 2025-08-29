#!/usr/bin/env python3
# coding: utf-8

"""
CircleSeeker2 File Organizer (Propmaster)
- Collects & renames CircleSeeker2 outputs
- Optionally detects HTML eccDNA report and renames/copies it
- Handles XeccDNA files controlled by configuration
"""

import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass

@dataclass
class PropmasterConfig:
    """Configuration for Propmaster file organizer."""
    prefix: str
    output_dir: Path
    include_xecc: bool = False
    report_name_template: str = "{prefix}_eccDNA_report.html"
    verbose: bool = False

def human_size(num_bytes: int) -> str:
    """Convert bytes to human-readable format."""
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(num_bytes)
    for u in units:
        if size < 1024.0 or u == units[-1]:
            return f"{size:.2f} {u}"
        size /= 1024.0
    return f"{size:.2f} TB"

def copy_with_rename(src: Path, dst: Path, logger: logging.Logger) -> bool:
    """Copy file with rename and create destination directory if needed."""
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        logger.info(f"Copied: {src} → {dst}")
        return True
    except Exception as e:
        logger.error(f"Failed to copy {src} → {dst}: {e}")
        return False

def newest(paths: List[Path]) -> Optional[Path]:
    """Get the newest file from a list of paths."""
    if not paths:
        return None
    try:
        return max(paths, key=lambda p: p.stat().st_mtime)
    except Exception:
        return None

class PropmasterOrganizer:
    """Organizes CircleSeeker2 output files with standardized naming."""
    
    def __init__(self, config: PropmasterConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Track processed destinations to avoid duplicates
        self.processed_destinations: Set[Path] = set()
        
        # Define standard mappings for each step
        self.step_mappings = {
            'uecc': {
                'fasta': f"{config.prefix}_UeccDNA.fa",
                'bed': f"{config.prefix}_UeccDNA.bed",
                'csv': f"{config.prefix}_UeccDNA.core.csv",
            },
            'mecc': {
                'fasta': f"{config.prefix}_MeccDNA.fa",
                'sites_bed': f"{config.prefix}_MeccDNA.sites.bed",
                'bestsite_bed': f"{config.prefix}_MeccDNA.bestsite.bed",
                'csv': f"{config.prefix}_MeccDNA.core.csv",
            },
            'cecc': {
                'fasta': f"{config.prefix}_CeccDNA.fa",
                'junctions_bedpe': f"{config.prefix}_CeccDNA.junctions.bedpe",
                'segments_bed': f"{config.prefix}_CeccDNA.segments.bed",
                'csv': f"{config.prefix}_CeccDNA.segments.core.csv",
            },
            'inferred': {
                'csv': f"{config.prefix}_InferredUeccDNA.csv",
                'fasta': f"{config.prefix}_InferredUeccDNA.fa",
                'bed': f"{config.prefix}_InferredUeccDNA.bed",
            }
        }
        
        if config.include_xecc:
            self.step_mappings['xecc'] = {
                'fasta': f"{config.prefix}_XeccDNA.fa",
                'sites_bed': f"{config.prefix}_XeccDNA.sites.bed",
                'bestsite_bed': f"{config.prefix}_XeccDNA.bestsite.bed",
                'junctions_bedpe': f"{config.prefix}_XeccDNA.junctions.bedpe",
                'segments_bed': f"{config.prefix}_XeccDNA.segments.bed",
                'csv': f"{config.prefix}_XeccDNA.core.csv",
            }
            self.step_mappings['inferred_xecc'] = {
                'csv': f"{config.prefix}_InferredXeccDNA.csv",
                'fasta': f"{config.prefix}_InferredXeccDNA.fa",
                'bed': f"{config.prefix}_InferredXeccDNA.bed",
            }
    
    def organize_files(self, 
                      base_dir: Path,
                      menagerie_outputs: Optional[Dict[str, Path]] = None,
                      harmonizer_outputs: Optional[Dict[str, Path]] = None,
                      juggler_outputs: Optional[Dict[str, Path]] = None,
                      report_file: Optional[Path] = None) -> Dict[str, int]:
        """
        Organize all CircleSeeker2 output files.
        
        Args:
            base_dir: Base directory containing output files
            menagerie_outputs: Dict mapping eccDNA types to their FASTA files
            harmonizer_outputs: Dict mapping eccDNA types to their harmonized files
            juggler_outputs: Dict mapping to final validation outputs
            report_file: Path to HTML report file
        
        Returns:
            Dict with counts of organized files by category
        """
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        
        organized_counts = {
            'uecc': 0,
            'mecc': 0, 
            'cecc': 0,
            'xecc': 0,
            'inferred': 0,
            'inferred_xecc': 0,
            'reports': 0
        }
        
        # Organize menagerie outputs (FASTA files)
        if menagerie_outputs:
            organized_counts.update(self._organize_menagerie_files(menagerie_outputs))
        
        # Organize harmonizer outputs (harmonized CSV files)
        if harmonizer_outputs:
            organized_counts.update(self._organize_harmonizer_files(harmonizer_outputs))
        
        # Organize juggler outputs (final validation results)
        if juggler_outputs:
            organized_counts.update(self._organize_juggler_files(juggler_outputs))
        
        # Copy report file
        if report_file and report_file.exists():
            self._copy_report(report_file)
            organized_counts['reports'] = 1
        else:
            # Try to detect report automatically
            detected_report = self._detect_report(base_dir)
            if detected_report:
                self._copy_report(detected_report)
                organized_counts['reports'] = 1
        
        # Write summary
        self._write_summary()
        
        self.logger.info("File organization completed")
        return organized_counts
    
    def _organize_menagerie_files(self, menagerie_outputs: Dict[str, Path]) -> Dict[str, int]:
        """Organize FASTA files from menagerie step."""
        counts = {}
        
        for ecc_type, fasta_file in menagerie_outputs.items():
            if fasta_file and fasta_file.exists():
                ecc_type_lower = ecc_type.lower()
                if ecc_type_lower in self.step_mappings:
                    dst_name = self.step_mappings[ecc_type_lower]['fasta']
                    dst_path = self.config.output_dir / dst_name
                    
                    if copy_with_rename(fasta_file, dst_path, self.logger):
                        counts[ecc_type_lower] = counts.get(ecc_type_lower, 0) + 1
        
        return counts
    
    def _organize_harmonizer_files(self, harmonizer_outputs: Dict[str, Path]) -> Dict[str, int]:
        """Organize harmonized CSV files."""
        counts = {}
        
        for ecc_type, csv_file in harmonizer_outputs.items():
            if csv_file and csv_file.exists():
                ecc_type_lower = ecc_type.lower()
                if ecc_type_lower in self.step_mappings:
                    dst_name = self.step_mappings[ecc_type_lower]['csv']
                    dst_path = self.config.output_dir / dst_name
                    
                    if copy_with_rename(csv_file, dst_path, self.logger):
                        counts[ecc_type_lower] = counts.get(ecc_type_lower, 0) + 1
        
        return counts
    
    def _organize_juggler_files(self, juggler_outputs: Dict[str, Path]) -> Dict[str, int]:
        """Organize final validation results from juggler."""
        counts = {}
        
        # Handle different types of juggler outputs
        for output_type, output_file in juggler_outputs.items():
            if output_file and output_file.exists():
                if output_type == 'validation_csv':
                    dst_name = f"{self.config.prefix}_final_validation.csv"
                elif output_type == 'inferred_fasta':
                    dst_name = self.step_mappings['inferred']['fasta']
                elif output_type == 'inferred_bed':
                    dst_name = self.step_mappings['inferred']['bed']
                else:
                    continue
                
                dst_path = self.config.output_dir / dst_name
                if copy_with_rename(output_file, dst_path, self.logger):
                    counts['inferred'] = counts.get('inferred', 0) + 1
        
        return counts
    
    def _detect_report(self, base_dir: Path) -> Optional[Path]:
        """Auto-detect HTML report files."""
        search_patterns = [
            "eccdna_report.html",
            "*_report.html", 
            "step*.html",
            "*eccdna*.html",
            "*.html"
        ]
        
        candidates = []
        for pattern in search_patterns:
            candidates.extend(base_dir.glob(pattern))
        
        # Filter existing files and get newest
        existing_files = [f for f in candidates if f.exists() and f.is_file()]
        
        if existing_files:
            report_file = newest(existing_files)
            self.logger.info(f"Auto-detected report: {report_file}")
            return report_file
        
        self.logger.warning("No HTML report found by auto-detection")
        return None
    
    def _copy_report(self, report_file: Path) -> bool:
        """Copy report file to organized directory."""
        dst_name = self.config.report_name_template.format(prefix=self.config.prefix)
        dst_path = self.config.output_dir / dst_name
        return copy_with_rename(report_file, dst_path, self.logger)
    
    def _write_summary(self) -> Path:
        """Write file organization summary."""
        summary_path = self.config.output_dir / f"{self.config.prefix}_file_summary.txt"
        lines = [f"CircleSeeker2 Organized Files - Prefix: {self.config.prefix}\n"]
        
        # Group files by category
        categories = {
            "UeccDNA": [],
            "MeccDNA": [],
            "CeccDNA": [],
            "XeccDNA": [],
            "Inferred UeccDNA": [],
            "Inferred XeccDNA": [],
            "Reports": [],
            "Other": []
        }
        
        # Scan output directory for all files
        for f in sorted(self.config.output_dir.glob(f"{self.config.prefix}_*")):
            if not f.is_file():
                continue
            
            fname = f.name
            size_str = f"({human_size(f.stat().st_size)})"
            file_entry = f"  {fname}  {size_str}"
            
            if "InferredXecc" in fname:
                categories["Inferred XeccDNA"].append(file_entry)
            elif "InferredUecc" in fname:
                categories["Inferred UeccDNA"].append(file_entry)
            elif "Xecc" in fname:
                categories["XeccDNA"].append(file_entry)
            elif "Uecc" in fname:
                categories["UeccDNA"].append(file_entry)
            elif "Mecc" in fname:
                categories["MeccDNA"].append(file_entry)
            elif "Cecc" in fname:
                categories["CeccDNA"].append(file_entry)
            elif "report" in fname.lower():
                categories["Reports"].append(file_entry)
            else:
                categories["Other"].append(file_entry)
        
        # Write categories
        for cat_name, files in categories.items():
            if files:
                lines.append(f"\n{cat_name}:")
                lines.extend(files)
        
        if not any(categories.values()):
            lines.append("\n(No files organized)")
        
        content = "\n".join(lines) + "\n"
        summary_path.write_text(content, encoding="utf-8")
        self.logger.info(f"Wrote summary: {summary_path}")
        return summary_path
    
    def organize_legacy_files(self, 
                            step6_dir: Path,
                            step8_dir: Path, 
                            step13_dir: Path) -> Dict[str, int]:
        """
        Organize files from legacy step structure.
        This method handles the old step-based file organization.
        """
        organized_counts = {
            'uecc': 0,
            'mecc': 0,
            'cecc': 0,
            'xecc': 0,
            'inferred': 0,
            'reports': 0
        }
        
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Legacy step8 mappings (U/M/C eccDNA)
        step8_mappings = {
            "step8_UeccDNA.fa": self.step_mappings['uecc']['fasta'],
            "step8_UeccDNA.bed": self.step_mappings['uecc']['bed'],
            "step8_UeccDNA.core.csv": self.step_mappings['uecc']['csv'],
            "step8_Mecc.fa": self.step_mappings['mecc']['fasta'],
            "step8_MeccSites.bed": self.step_mappings['mecc']['sites_bed'],
            "step8_MeccBestSite.bed": self.step_mappings['mecc']['bestsite_bed'],
            "step8_MeccSites.core.csv": self.step_mappings['mecc']['csv'],
            "step8_Cecc.fa": self.step_mappings['cecc']['fasta'],
            "step8_CeccJunctions.bedpe": self.step_mappings['cecc']['junctions_bedpe'],
            "step8_CeccSegments.bed": self.step_mappings['cecc']['segments_bed'],
            "step8_CeccSegments.core.csv": self.step_mappings['cecc']['csv'],
        }
        
        # Legacy step13 mappings (inferred Uecc)
        step13_mappings = {
            "step13.clean.csv": self.step_mappings['inferred']['csv'],
            "step13.renamed.fasta": self.step_mappings['inferred']['fasta'],
            "step13.simplified.bed": self.step_mappings['inferred']['bed'],
        }
        
        # Copy step8 files
        self.logger.info("Organizing legacy step8 files...")
        organized_counts['uecc'] += self._copy_mapping(step8_dir, step8_mappings, "step8")
        organized_counts['mecc'] += self._copy_mapping(step8_dir, step8_mappings, "step8")
        organized_counts['cecc'] += self._copy_mapping(step8_dir, step8_mappings, "step8")
        
        # Copy step13 files
        self.logger.info("Organizing legacy step13 files...")
        organized_counts['inferred'] += self._copy_mapping(step13_dir, step13_mappings, "step13")
        
        # Handle Xecc files if enabled
        if self.config.include_xecc:
            organized_counts['xecc'] = self._copy_xecc_files(step6_dir, step8_dir, step13_dir)
        
        # Detect and copy report
        for search_dir in [step13_dir, step8_dir, step6_dir, Path(".")]:
            report = self._detect_report(search_dir)
            if report:
                self._copy_report(report)
                organized_counts['reports'] = 1
                break
        
        # Write summary
        self._write_summary()
        
        return organized_counts
    
    def _copy_mapping(self, base_dir: Path, mapping: Dict[str, str], label: str) -> int:
        """Copy files based on mapping and return count of successful copies."""
        copied_count = 0
        for src_name, dst_name in mapping.items():
            src = base_dir / src_name
            dst = self.config.output_dir / dst_name
            
            if src.exists() and copy_with_rename(src, dst, self.logger):
                copied_count += 1
            elif not src.exists():
                self.logger.debug(f"[{label}] File not found: {src_name}")
        
        return copied_count
    
    def _copy_xecc_files(self, step6_dir: Path, step8_dir: Path, step13_dir: Path) -> int:
        """Handle Xecc files from all steps."""
        if not self.config.include_xecc:
            return 0
        
        self.logger.info("Organizing XeccDNA files...")
        copied_count = 0
        
        # Define all possible Xecc file locations and their destinations
        xecc_mappings = [
            # Step 6 files
            (step6_dir / "step6_XeccDNA.fasta", self.step_mappings['xecc']['fasta']),
            (step6_dir / "step6_XeccDNA.fa", self.step_mappings['xecc']['fasta']),
            (step6_dir / "step6_XeccDNA_pre.fasta", self.step_mappings['xecc']['fasta']),
            
            # Step 8 files
            (step8_dir / "step8_Xecc.fa", self.step_mappings['xecc']['fasta']),
            (step8_dir / "step8_XeccSites.bed", self.step_mappings['xecc']['sites_bed']),
            (step8_dir / "step8_XeccBestSite.bed", self.step_mappings['xecc']['bestsite_bed']),
            (step8_dir / "step8_XeccSites.core.csv", self.step_mappings['xecc']['csv']),
            (step8_dir / "step8_XeccJunctions.bedpe", self.step_mappings['xecc']['junctions_bedpe']),
            (step8_dir / "step8_XeccSegments.bed", self.step_mappings['xecc']['segments_bed']),
            
            # Step 13 files
            (step13_dir / "step13.Xecc.clean.csv", self.step_mappings['inferred_xecc']['csv']),
            (step13_dir / "step13.Xecc.renamed.fasta", self.step_mappings['inferred_xecc']['fasta']),
            (step13_dir / "step13.Xecc.simplified.bed", self.step_mappings['inferred_xecc']['bed']),
        ]
        
        # Process each mapping
        for src, dst_name in xecc_mappings:
            if src.exists():
                dst = self.config.output_dir / dst_name
                # Avoid duplicate copies to the same destination
                if dst not in self.processed_destinations:
                    if copy_with_rename(src, dst, self.logger):
                        copied_count += 1
                        self.processed_destinations.add(dst)
                else:
                    self.logger.debug(f"Skipping duplicate copy to: {dst}")
        
        return copied_count
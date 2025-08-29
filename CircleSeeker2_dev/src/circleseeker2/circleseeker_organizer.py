#!/usr/bin/env python3
# coding: utf-8

"""
CircleSeeker2 File Organizer
Renames and organizes output files from CircleSeeker2 pipeline
"""

import pandas as pd
import shutil
import os
import argparse
import logging
from pathlib import Path
from typing import Dict, List

class CircleSeekerFileOrganizer:
    def __init__(self, 
                 prefix: str,
                 step8_dir: str = ".",
                 step13_dir: str = ".",
                 output_dir: str = "organized_results",
                 verbose: bool = False):
        
        self.prefix = prefix
        self.step8_dir = Path(step8_dir)
        self.step13_dir = Path(step13_dir)
        self.output_dir = Path(output_dir)
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        log_level = logging.DEBUG if verbose else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Define file mappings for step8
        self.step8_mappings = {
            'step8_UeccDNA.fa': f'{prefix}_UeccDNA.fa',
            'step8_UeccDNA.bed': f'{prefix}_UeccDNA.bed',
            'step8_UeccDNA.core.csv': f'{prefix}_UeccDNA.core.csv',
            'step8_Mecc.fa': f'{prefix}_MeccDNA.fa',
            'step8_MeccSites.bed': f'{prefix}_MeccDNA.sites.bed',
            'step8_MeccBestSite.bed': f'{prefix}_MeccDNA.bestsite.bed',
            'step8_MeccSites.core.csv': f'{prefix}_MeccDNA.core.csv',
            'step8_Cecc.fa': f'{prefix}_CeccDNA.fa',
            'step8_CeccJunctions.bedpe': f'{prefix}_CeccDNA.junctions.bedpe',
            'step8_CeccSegments.bed': f'{prefix}_CeccDNA.segments.bed'
        }
        
        # Define file mappings for step13
        self.step13_mappings = {
            'step13.clean.csv': f'{prefix}_InferredUeccDNA.csv',
            'step13.renamed.fasta': f'{prefix}_InferredUeccDNA.fa',
            'step13.simplified.bed': f'{prefix}_InferredUeccDNA.bed'
        }
    
    def copy_and_rename_step8(self):
        """Copy and rename step8 files"""
        self.logger.info("Processing Step 8 files...")
        copied_files = []
        
        for old_name, new_name in self.step8_mappings.items():
            old_path = self.step8_dir / old_name
            new_path = self.output_dir / new_name
            
            if old_path.exists():
                shutil.copy2(old_path, new_path)
                self.logger.info(f"Copied: {old_name} → {new_name}")
                copied_files.append(new_name)
            else:
                self.logger.warning(f"File not found: {old_path}")
        
        return copied_files
    
    def copy_and_rename_step13(self):
        """Copy and rename step13 files"""
        self.logger.info("Processing Step 13 files...")
        copied_files = []
        
        for old_name, new_name in self.step13_mappings.items():
            old_path = self.step13_dir / old_name
            new_path = self.output_dir / new_name
            
            if old_path.exists():
                shutil.copy2(old_path, new_path)
                self.logger.info(f"Copied: {old_name} → {new_name}")
                copied_files.append(new_name)
            else:
                self.logger.warning(f"File not found: {old_path}")
        
        return copied_files
    
    def generate_file_summary(self):
        """Generate a summary of all organized files"""
        summary_path = self.output_dir / f"{self.prefix}_file_summary.txt"
        
        with open(summary_path, 'w') as f:
            f.write("CircleSeeker2 Organized Files Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Sample Prefix: {self.prefix}\n")
            f.write(f"Output Directory: {self.output_dir}\n\n")
            
            f.write("Confirmed eccDNA Files:\n")
            f.write("-" * 30 + "\n")
            
            # List UeccDNA files
            f.write("UeccDNA:\n")
            for file in [f"{self.prefix}_UeccDNA.fa", f"{self.prefix}_UeccDNA.bed", f"{self.prefix}_UeccDNA.core.csv"]:
                if (self.output_dir / file).exists():
                    size = (self.output_dir / file).stat().st_size
                    f.write(f"  - {file} ({self.format_size(size)})\n")
            
            # List MeccDNA files
            f.write("\nMeccDNA:\n")
            for file in [f"{self.prefix}_MeccDNA.fa", f"{self.prefix}_MeccDNA.sites.bed", 
                        f"{self.prefix}_MeccDNA.bestsite.bed", f"{self.prefix}_MeccDNA.core.csv"]:
                if (self.output_dir / file).exists():
                    size = (self.output_dir / file).stat().st_size
                    f.write(f"  - {file} ({self.format_size(size)})\n")
            
            # List CeccDNA files
            f.write("\nCeccDNA:\n")
            for file in [f"{self.prefix}_CeccDNA.fa", f"{self.prefix}_CeccDNA.junctions.bedpe", 
                        f"{self.prefix}_CeccDNA.segments.bed"]:
                if (self.output_dir / file).exists():
                    size = (self.output_dir / file).stat().st_size
                    f.write(f"  - {file} ({self.format_size(size)})\n")
            
            # List Inferred UeccDNA files
            f.write("\nInferred UeccDNA:\n")
            for file in [f"{self.prefix}_InferredUeccDNA.csv", f"{self.prefix}_InferredUeccDNA.fa", 
                        f"{self.prefix}_InferredUeccDNA.bed"]:
                if (self.output_dir / file).exists():
                    size = (self.output_dir / file).stat().st_size
                    f.write(f"  - {file} ({self.format_size(size)})\n")
        
        self.logger.info(f"File summary saved: {summary_path}")
    
    def format_size(self, size: int) -> str:
        """Format file size in human-readable format"""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size < 1024.0:
                return f"{size:.1f} {unit}"
            size /= 1024.0
        return f"{size:.1f} TB"
    
    def run(self):
        """Main workflow"""
        self.logger.info(f"Starting file organization for prefix: {self.prefix}")
        
        # Process Step 8 files
        self.copy_and_rename_step8()
        
        # Process Step 13 files (simple rename, no deduplication)
        self.copy_and_rename_step13()
        
        # Generate summary
        self.generate_file_summary()
        
        self.logger.info(f"File organization complete! Results in: {self.output_dir}")
        print(f"\n✅ Success! All files have been organized in: {self.output_dir}")
        print(f"📊 Check {self.prefix}_file_summary.txt for details")

def main():
    parser = argparse.ArgumentParser(
        description="Organize and rename CircleSeeker2 output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Basic usage with all files in current directory
  python circleseeker_organizer.py --prefix Sample001
  
  # Specify different directories for step8 and step13 files
  python circleseeker_organizer.py \\
    --prefix HeLa_Sample \\
    --step8-dir results/step8 \\
    --step13-dir results/step13 \\
    --output organized_results
  
  # With verbose output
  python circleseeker_organizer.py --prefix Sample001 --verbose
        """
    )
    
    # Required arguments
    parser.add_argument("--prefix", required=True,
                       help="Prefix for renamed files (e.g., Sample001)")
    
    # Optional arguments
    parser.add_argument("--step8-dir", default=".",
                       help="Directory containing step8 files (default: current directory)")
    parser.add_argument("--step13-dir", default=".",
                       help="Directory containing step13 files (default: current directory)")
    parser.add_argument("--output", default="organized_results",
                       help="Output directory for organized files (default: organized_results)")
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    
    args = parser.parse_args()
    
    # Create organizer instance and run
    organizer = CircleSeekerFileOrganizer(
        prefix=args.prefix,
        step8_dir=args.step8_dir,
        step13_dir=args.step13_dir,
        output_dir=args.output,
        verbose=args.verbose
    )
    
    try:
        organizer.run()
    except Exception as e:
        logging.error(f"Organization failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()

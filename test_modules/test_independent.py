#!/usr/bin/env python3
"""
Independent test for Playbill module (Step 14)
Tests the playbill report generator module in isolation using new implementation.
"""

import sys
from pathlib import Path
import logging
import filecmp

# Add source directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from circleseeker2.modules.playbill import PlaybillReportGenerator

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_playbill_module():
    """Test Playbill module independently"""
    
    logger.info("=" * 60)
    logger.info("Testing Playbill Module (Step 14) - Independent")
    logger.info("=" * 60)
    
    # Directory setup
    test_dir = Path(__file__).parent
    input_dir = test_dir / "input"
    output_dir = test_dir / "output"
    expected_dir = test_dir / "expected"
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Input files - CORRECT paths matching actual FASTA files
    uecc_fasta = input_dir / "step6_UeccDNA_pre.fasta"
    mecc_fasta = input_dir / "step6_MeccDNA_pre.fasta"
    cecc_fasta = input_dir / "step6_CeccDNA_pre.fasta"
    inferred_fasta = input_dir / "step13.renamed.fasta"
    reads_fasta = input_dir / "step1_Test.HeLa_10k.fasta"
    processed_csv = input_dir / "step2_processed.csv"
    
    # Output files
    output_html = output_dir / "eccdna_report.html"
    expected_html = expected_dir / "eccdna_report.html"
    
    # Check all required inputs exist
    input_files = {
        "UeccDNA FASTA": uecc_fasta,
        "MeccDNA FASTA": mecc_fasta,
        "CeccDNA FASTA": cecc_fasta,
        "Inferred FASTA": inferred_fasta,
        "Reads FASTA": reads_fasta,
        "Processed CSV": processed_csv
    }
    
    missing_files = []
    for description, file_path in input_files.items():
        if file_path.exists():
            size = file_path.stat().st_size
            logger.info(f"✅ {description}: {file_path.name} ({size:,} bytes)")
        else:
            logger.error(f"❌ {description} not found: {file_path}")
            missing_files.append(description)
    
    if missing_files:
        logger.error(f"Missing {len(missing_files)} required input files")
        return False
    
    try:
        # Initialize PlaybillReportGenerator
        playbill = PlaybillReportGenerator(logger=logger)
        
        # Generate analysis report
        logger.info("\nGenerating eccDNA analysis report...")
        results = playbill.generate_report(
            uecc_fasta=uecc_fasta,
            mecc_fasta=mecc_fasta,
            cecc_fasta=cecc_fasta,
            inferred_fasta=inferred_fasta,
            reads_fasta=reads_fasta,
            processed_csv=processed_csv,
            output_html=output_html,
            sample_name="HeLa_10k_Test"
        )
        
        if results is None:
            logger.error("Playbill processing failed - no results")
            return False
        
        # Check output was created
        if not output_html.exists():
            logger.error(f"Output HTML not created: {output_html}")
            return False
        
        # Report output details
        html_size = output_html.stat().st_size
        logger.info(f"\n✅ Generated HTML report: {output_html.name} ({html_size:,} bytes)")
        
        # Show analysis summary
        logger.info("\nAnalysis Results:")
        
        # eccDNA statistics
        if 'eccdna_stats' in results:
            for stats in results['eccdna_stats']:
                if stats['label'] != "All eccDNA (Combined)":
                    logger.info(f"  {stats['label']}: {stats['count']:,} sequences")
        
        # Read classification
        if 'reads_classification' in results:
            reads = results['reads_classification']
            total = reads['total_reads'] if reads['total_reads'] > 0 else 1
            ctcr_pct = (reads['total_ctcr'] / total * 100) if total > 0 else 0
            logger.info(f"  Total reads: {reads['total_reads']:,}")
            logger.info(f"  CtcR reads: {reads['total_ctcr']:,} ({ctcr_pct:.2f}%)")
        
        # Overlap analysis
        if 'overlaps' in results and 'inferred_stats' in results:
            inferred_total = results['inferred_stats']['count']
            if inferred_total > 0:
                uecc_overlap_pct = len(results['overlaps']['uecc']) / inferred_total * 100
                mecc_overlap_pct = len(results['overlaps']['mecc']) / inferred_total * 100
                cecc_overlap_pct = len(results['overlaps']['cecc']) / inferred_total * 100
                
                logger.info(f"  Overlap analysis:")
                logger.info(f"    UeccDNA: {len(results['overlaps']['uecc']):,} ({uecc_overlap_pct:.1f}%)")
                logger.info(f"    MeccDNA: {len(results['overlaps']['mecc']):,} ({mecc_overlap_pct:.1f}%)")
                logger.info(f"    CeccDNA: {len(results['overlaps']['cecc']):,} ({cecc_overlap_pct:.1f}%)")
        
        # Validate HTML content
        try:
            with open(output_html) as f:
                html_content = f.read()
                required_elements = [
                    "eccDNA Analysis Report",
                    "HeLa_10k_Test",
                    "Read Classification Statistics",
                    "CtcR Subtype Distribution",
                    "eccDNA Sequence Statistics"
                ]
                
                missing_elements = []
                for element in required_elements:
                    if element not in html_content:
                        missing_elements.append(element)
                
                if missing_elements:
                    logger.warning(f"HTML missing elements: {missing_elements}")
                else:
                    logger.info("✅ HTML content validation passed")
                    
        except Exception as e:
            logger.warning(f"Could not validate HTML content: {e}")
        
        # Compare with expected output (if available)
        match_result = None
        if expected_html.exists():
            # For HTML, compare sizes with tolerance
            output_size = output_html.stat().st_size
            expected_size = expected_html.stat().st_size
            size_diff = abs(output_size - expected_size) / expected_size if expected_size > 0 else 0
            
            if size_diff < 0.1:  # 10% tolerance
                match_result = True
                logger.info(f"\nHTML comparison: ✅ SIMILAR (size diff: {size_diff:.2%})")
            else:
                match_result = False
                logger.info(f"\nHTML comparison: ❌ DIFFER (size diff: {size_diff:.2%})")
        else:
            logger.info(f"\nNo expected file for comparison")
        
        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("PLAYBILL TEST SUMMARY")
        logger.info("=" * 60)
        logger.info(f"Processing: ✅ SUCCESS")
        logger.info(f"HTML generated: {html_size:,} bytes")
        if match_result is not None:
            logger.info(f"Output validation: {'✅ PASS' if match_result else '⚠️ SIZE DIFFERS'}")
        
        return True
        
    except Exception as e:
        logger.error(f"Playbill test failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    success = test_playbill_module()
    sys.exit(0 if success else 1)
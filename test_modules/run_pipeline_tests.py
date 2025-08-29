#!/usr/bin/env python3
"""
CircleSeeker2 Pipeline Module Testing Script

This script runs individual modules of the CircleSeeker2 pipeline 
to test each step independently with realistic data.

Usage: python run_pipeline_tests.py [--step STEP_NUMBER] [--all]
"""

import argparse
import subprocess
import sys
from pathlib import Path
import logging
import shutil

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PipelineTester:
    """Test runner for CircleSeeker2 pipeline modules"""
    
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.src_dir = base_dir.parent / "src" / "circleseeker2"
        self.legacy_dir = base_dir.parent / "legacy_archive"
        
        # Test data directories
        self.input_dirs = {}
        self.output_dirs = {}
        self.expected_dirs = {}
        
        # Setup directory mappings for each step
        for i in range(16):
            step_name = f"step{i:02d}"
            step_dirs = {
                'input': base_dir / f"{step_name}_*" / "input",
                'output': base_dir / f"{step_name}_*" / "output", 
                'expected': base_dir / f"{step_name}_*" / "expected"
            }
            self.input_dirs[i] = step_dirs['input']
            self.output_dirs[i] = step_dirs['output']
            self.expected_dirs[i] = step_dirs['expected']
        
    def get_step_dirs(self, step_num):
        """Get input, output, and expected directories for a step"""
        step_mappings = {
            0: "step00_make_blastdb",
            1: "step01_tidehunter", 
            2: "step02_carousel",
            3: "step03_run_blast",
            4: "step04_gatekeeper",
            5: "step05_trapeze",
            6: "step06_menagerie",
            7: "step07_cd_hit",
            8: "step08_harmonizer",
            9: "step09_sieve",
            10: "step10_minimap2",
            11: "step11_mosdepth",
            12: "step12_contortionist",
            13: "step13_juggler",
            14: "step14_playbill",
            15: "step15_propmaster"
        }
        
        if step_num not in step_mappings:
            raise ValueError(f"Invalid step number: {step_num}")
            
        step_dir = self.base_dir / step_mappings[step_num]
        return {
            'input': step_dir / "input",
            'output': step_dir / "output",
            'expected': step_dir / "expected",
            'base': step_dir
        }
        
    def check_dependencies(self):
        """Check if required dependencies are available"""
        required_tools = ['blastn', 'makeblastdb', 'cd-hit-est', 'minimap2', 'samtools', 'mosdepth']
        missing_tools = []
        
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            logger.warning(f"Missing tools: {missing_tools}")
            return False
        return True
    
    def run_step_00_make_blastdb(self):
        """Step 0: Build BLAST database"""
        logger.info("Running Step 0: Make BLAST Database")
        step_dir = self.base_dir / "step00_make_blastdb"
        
        if not Path(self.reference_genome).exists():
            logger.error(f"Reference genome not found: {self.reference_genome}")
            return False
            
        cmd = [
            "python", str(self.legacy_dir / "step0_make_blastdb.py"),
            "-i", self.reference_genome,
            "-o", str(step_dir / "chm13v2_db")
        ]
        
        return self._run_command(cmd, step_dir)
    
    def run_step_01_tidehunter(self):
        """Step 1: TideHunter tandem repeat detection"""
        logger.info("Running Step 1: TideHunter")
        dirs = self.get_step_dirs(1)
        
        # Check input files
        input_fasta = dirs['input'] / "step1_Test.HeLa_10k.fasta"
        if not input_fasta.exists():
            logger.error(f"Input FASTA not found: {input_fasta}")
            return False
            
        cmd = [
            "python", str(self.legacy_dir / "step1_tidehunter.py"),
            "-i", str(input_fasta),
            "-o", str(dirs['output'] / "step1_Test.TH.ecc_candidates.txt"),
            "-t", "8"
        ]
        
        result = self._run_command(cmd, dirs['base'])
        
        # Compare output with expected
        if result:
            result = self._compare_outputs(1)
        return result
    
    def _compare_outputs(self, step_num):
        """Compare step outputs with expected results"""
        dirs = self.get_step_dirs(step_num)
        
        # Get all files in output directory
        if not dirs['output'].exists():
            logger.error(f"Output directory not found: {dirs['output']}")
            return False
            
        if not dirs['expected'].exists():
            logger.error(f"Expected directory not found: {dirs['expected']}")
            return False
            
        output_files = list(dirs['output'].glob("*"))
        expected_files = list(dirs['expected'].glob("*"))
        
        if len(output_files) == 0:
            logger.warning(f"No output files generated in step {step_num}")
            return False
            
        logger.info(f"Comparing {len(output_files)} output files with {len(expected_files)} expected files")
        
        # For now, just check that files exist - could be enhanced to do content comparison
        for output_file in output_files:
            expected_file = dirs['expected'] / output_file.name
            if not expected_file.exists():
                logger.warning(f"No expected file found for: {output_file.name}")
            else:
                logger.info(f"✓ Found expected file for: {output_file.name}")
                
        return True
    
    def run_step_02_carousel(self):
        """Step 2: Carousel processing"""
        logger.info("Running Step 2: Carousel")
        step_dir = self.base_dir / "step02_carousel"
        
        input_file = step_dir / "step1_Test.TH.ecc_candidates.txt"
        if not input_file.exists():
            # Try to get from previous step or use existing test data
            alt_input = self.base_dir / "step01_tidehunter" / "step1_Test.TH.ecc_candidates.txt"
            if alt_input.exists():
                shutil.copy2(alt_input, input_file)
            elif (step_dir / "input_tidehunter_results.txt").exists():
                input_file = step_dir / "input_tidehunter_results.txt"
            else:
                logger.error(f"TideHunter input not found: {input_file}")
                return False
        
        cmd = [
            "python", str(self.legacy_dir / "step2_carousel.py"),
            "-i", str(input_file),
            "-o", str(step_dir / "step2_processed.csv"),
            "-f", str(step_dir / "step2_circular.fasta")
        ]
        
        return self._run_command(cmd, step_dir)
    
    def run_step_03_run_blast(self):
        """Step 3: BLAST alignment"""
        logger.info("Running Step 3: BLAST")
        step_dir = self.base_dir / "step03_run_blast"
        
        # Check for database and query file
        db_path = self.base_dir / "step00_make_blastdb" / "chm13v2_db"
        query_file = step_dir / "step2_circular.fasta"
        
        # Try to get circular fasta from carousel step
        if not query_file.exists():
            alt_query = self.base_dir / "step02_carousel" / "step2_circular.fasta"
            if alt_query.exists():
                shutil.copy2(alt_query, query_file)
            else:
                logger.error(f"Circular FASTA not found: {query_file}")
                return False
        
        if not db_path.exists():
            logger.error(f"BLAST database not found: {db_path}")
            return False
        
        cmd = [
            "python", str(self.legacy_dir / "step3_run_blast.py"),
            "-d", str(db_path),
            "-q", str(query_file),
            "-o", str(step_dir / "step3_blast_results.tsv"),
            "-t", "8"
        ]
        
        return self._run_command(cmd, step_dir)
    
    def run_step_04_gatekeeper(self):
        """Step 4: Gatekeeper classification"""
        logger.info("Running Step 4: Gatekeeper")
        step_dir = self.base_dir / "step04_gatekeeper"
        
        input_file = step_dir / "step3_blast_results.tsv"
        if not input_file.exists():
            # Try to get from previous step
            alt_input = self.base_dir / "step03_run_blast" / "step3_blast_results.tsv"
            if alt_input.exists():
                shutil.copy2(alt_input, input_file)
            else:
                logger.error(f"BLAST results not found: {input_file}")
                return False
        
        cmd = [
            "python", str(self.legacy_dir / "step4_gatekeeper.py"),
            "-b", str(input_file),
            "--prefix", "step4"
        ]
        
        return self._run_command(cmd, step_dir)
    
    def run_step_05_trapeze(self):
        """Step 5: Trapeze processing"""
        logger.info("Running Step 5: Trapeze")
        step_dir = self.base_dir / "step05_trapeze"
        
        input_file = step_dir / "step4_unclassified.csv"
        if not input_file.exists():
            # Try to get from previous step
            alt_input = self.base_dir / "step04_gatekeeper" / "step4_unclassified.csv"
            if alt_input.exists():
                shutil.copy2(alt_input, input_file)
            else:
                logger.warning(f"Unclassified file not found: {input_file}, skipping trapeze")
                return True  # This is optional
        
        cmd = [
            "python", str(self.legacy_dir / "step5_trapeze.py"),
            "-i", str(input_file),
            "-o", str(step_dir / "step5_cecc.csv")
        ]
        
        return self._run_command(cmd, step_dir)
    
    def run_step_06_menagerie(self):
        """Step 6: Menagerie FASTA generation"""
        logger.info("Running Step 6: Menagerie")
        step_dir = self.base_dir / "step06_menagerie"
        
        # Get required files
        circular_fasta = step_dir / "step2_circular.fasta"
        uecc_csv = step_dir / "step4_uecc.csv"
        mecc_csv = step_dir / "step4_mecc.csv"
        cecc_csv = step_dir / "step5_cecc.csv"
        
        # Try to get files from previous steps if not present
        if not circular_fasta.exists():
            alt_fasta = self.base_dir / "step02_carousel" / "step2_circular.fasta"
            if alt_fasta.exists():
                shutil.copy2(alt_fasta, circular_fasta)
        
        if not uecc_csv.exists():
            alt_uecc = self.base_dir / "step04_gatekeeper" / "step4_uecc.csv"
            if alt_uecc.exists():
                shutil.copy2(alt_uecc, uecc_csv)
        
        if not mecc_csv.exists():
            alt_mecc = self.base_dir / "step04_gatekeeper" / "step4_mecc.csv"
            if alt_mecc.exists():
                shutil.copy2(alt_mecc, mecc_csv)
        
        if not cecc_csv.exists():
            alt_cecc = self.base_dir / "step05_trapeze" / "step5_cecc.csv"
            if alt_cecc.exists():
                shutil.copy2(alt_cecc, cecc_csv)
        
        cmd = [
            "python", str(self.legacy_dir / "step6_menagerie.py"),
            "-f", str(circular_fasta),
            "-u", str(uecc_csv),
            "-m", str(mecc_csv),
            "-c", str(cecc_csv),
            "-o", str(step_dir),
            "-p", "step6",
            "-X"
        ]
        
        return self._run_command(cmd, step_dir)
    
    def _run_command(self, cmd, working_dir):
        """Execute a command and return success status"""
        try:
            logger.info(f"Executing: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, 
                cwd=working_dir, 
                check=True, 
                capture_output=True, 
                text=True,
                timeout=300  # 5 minutes timeout
            )
            logger.info(f"Command completed successfully")
            if result.stdout:
                logger.debug(f"STDOUT: {result.stdout}")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}")
            if e.stdout:
                logger.error(f"STDOUT: {e.stdout}")
            if e.stderr:
                logger.error(f"STDERR: {e.stderr}")
            return False
            
        except subprocess.TimeoutExpired:
            logger.error("Command timed out")
            return False
            
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            return False
    
    def run_all_steps(self):
        """Run all pipeline steps in sequence"""
        steps = [
            ("Step 0: Make BLAST DB", self.run_step_00_make_blastdb),
            ("Step 1: TideHunter", self.run_step_01_tidehunter),
            ("Step 2: Carousel", self.run_step_02_carousel),
            ("Step 3: BLAST", self.run_step_03_run_blast),
            ("Step 4: Gatekeeper", self.run_step_04_gatekeeper),
            ("Step 5: Trapeze", self.run_step_05_trapeze),
            ("Step 6: Menagerie", self.run_step_06_menagerie),
        ]
        
        results = {}
        
        for step_name, step_func in steps:
            logger.info(f"=" * 60)
            logger.info(f"Starting {step_name}")
            logger.info(f"=" * 60)
            
            success = step_func()
            results[step_name] = success
            
            if success:
                logger.info(f"✅ {step_name} completed successfully")
            else:
                logger.error(f"❌ {step_name} failed")
                # Continue with other steps even if one fails
        
        # Summary
        logger.info(f"=" * 60)
        logger.info("PIPELINE TEST SUMMARY")
        logger.info(f"=" * 60)
        
        for step_name, success in results.items():
            status = "✅ PASS" if success else "❌ FAIL"
            logger.info(f"{status:10} {step_name}")
        
        total_passed = sum(results.values())
        total_steps = len(results)
        logger.info(f"\nPassed: {total_passed}/{total_steps}")
        
        return results


def main():
    parser = argparse.ArgumentParser(description='Test CircleSeeker2 pipeline modules')
    parser.add_argument('--step', type=int, help='Run specific step (0-6)')
    parser.add_argument('--all', action='store_true', help='Run all steps')
    parser.add_argument('--check-deps', action='store_true', help='Check dependencies')
    
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent
    tester = PipelineTester(base_dir)
    
    if args.check_deps:
        if tester.check_dependencies():
            print("✅ All dependencies are available")
        else:
            print("❌ Some dependencies are missing")
        return
    
    if args.step is not None:
        step_methods = {
            0: tester.run_step_00_make_blastdb,
            1: tester.run_step_01_tidehunter,
            2: tester.run_step_02_carousel,
            3: tester.run_step_03_run_blast,
            4: tester.run_step_04_gatekeeper,
            5: tester.run_step_05_trapeze,
            6: tester.run_step_06_menagerie,
        }
        
        if args.step in step_methods:
            success = step_methods[args.step]()
            sys.exit(0 if success else 1)
        else:
            print(f"Invalid step number: {args.step}")
            sys.exit(1)
    
    elif args.all:
        results = tester.run_all_steps()
        # Exit with error code if any step failed
        sys.exit(0 if all(results.values()) else 1)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
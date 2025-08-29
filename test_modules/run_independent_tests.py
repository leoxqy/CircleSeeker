#!/usr/bin/env python3
"""
Independent Module Test Runner
Runs each module's independent test in isolation.
"""

import subprocess
import sys
from pathlib import Path
import argparse
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_module_test(module_dir):
    """Run independent test for a specific module"""
    test_script = module_dir / "test_independent.py"
    
    if not test_script.exists():
        logger.warning(f"No independent test script found: {test_script}")
        return None
    
    logger.info(f"Running test for {module_dir.name}...")
    
    try:
        result = subprocess.run(
            [sys.executable, str(test_script)],
            cwd=str(module_dir),
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per module
        )
        
        if result.returncode == 0:
            logger.info(f"✅ {module_dir.name}: SUCCESS")
            return True
        else:
            logger.error(f"❌ {module_dir.name}: FAILED")
            logger.error(f"STDOUT:\n{result.stdout}")
            logger.error(f"STDERR:\n{result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"⏱️  {module_dir.name}: TIMEOUT")
        return False
    except Exception as e:
        logger.error(f"💥 {module_dir.name}: ERROR - {e}")
        return False

def main():
    """Main test runner"""
    parser = argparse.ArgumentParser(description="Run independent module tests")
    parser.add_argument("--module", help="Specific module to test (e.g., step02_carousel)")
    parser.add_argument("--verbose", action="store_true", help="Show detailed output")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    base_dir = Path(__file__).parent
    
    # Define modules with independent tests
    modules_to_test = [
        "step00_make_blastdb",
        "step01_tidehunter",
        "step02_carousel", 
        "step03_run_blast",
        "step04_gatekeeper",
        "step05_trapeze",
        "step06_menagerie",
        "step07_cd_hit", 
        "step08_harmonizer",
        "step09_sieve",
        "step10_minimap2",
        "step11_mosdepth",
        "step12_contortionist"
    ]
    
    if args.module:
        if args.module in modules_to_test:
            modules_to_test = [args.module]
        else:
            logger.error(f"Unknown module: {args.module}")
            logger.info(f"Available modules: {', '.join(modules_to_test)}")
            return 1
    
    logger.info("=" * 80)
    logger.info("RUNNING INDEPENDENT MODULE TESTS")
    logger.info("=" * 80)
    logger.info(f"Testing {len(modules_to_test)} modules...")
    
    results = {}
    
    # Run tests for each module
    for module_name in modules_to_test:
        module_dir = base_dir / module_name
        
        if not module_dir.exists():
            logger.warning(f"Module directory not found: {module_dir}")
            results[module_name] = None
            continue
        
        results[module_name] = run_module_test(module_dir)
    
    # Summary
    logger.info("=" * 80)
    logger.info("INDEPENDENT TESTS SUMMARY")
    logger.info("=" * 80)
    
    passed = 0
    failed = 0
    skipped = 0
    
    for module_name, result in results.items():
        if result is True:
            logger.info(f"✅ {module_name}")
            passed += 1
        elif result is False:
            logger.info(f"❌ {module_name}")
            failed += 1
        else:
            logger.info(f"⚠️  {module_name} (skipped)")
            skipped += 1
    
    logger.info("-" * 80)
    logger.info(f"Total: {len(results)} | Passed: {passed} | Failed: {failed} | Skipped: {skipped}")
    
    if failed > 0:
        logger.error("Some tests failed!")
        return 1
    
    logger.info("All available tests passed!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
"""TideHunter wrapper."""

from pathlib import Path
from typing import Optional
from circleseeker.external.base import ExternalTool


class TideHunter(ExternalTool):
    """TideHunter tandem repeat finder."""
    
    tool_name = "TideHunter"
    
    def run_analysis(
        self,
        input_file: Path,
        output_file: Path,
        k: int = 16,
        w: int = 1,
        p: int = 100,
        P: int = 2000000,
        e: float = 0.1,
        f: int = 2
    ) -> None:
        """Run TideHunter analysis."""
        cmd = [
            self.tool_name,
            "-f", str(f),
            "-t", str(self.threads),
            "-k", str(k),
            "-w", str(w),
            "-p", str(p),
            "-P", str(P),
            "-e", str(e),
            str(input_file)
        ]
        
        # Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Run TideHunter with output redirection
        with open(output_file, 'w') as out_handle:
            import subprocess
            try:
                result = subprocess.run(
                    cmd,
                    stdout=out_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
                
                self.logger.info(f"TideHunter completed successfully")
                self.logger.info(f"Output saved to: {output_file}")
                
                # Log output file size
                output_size = output_file.stat().st_size
                self.logger.debug(f"Output file size: {output_size} bytes")
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"TideHunter failed with exit code {e.returncode}")
                self.logger.error(f"Error message: {e.stderr}")
                from circleseeker.exceptions import ExternalToolError
                raise ExternalToolError(
                    f"TideHunter failed",
                    command=cmd,
                    returncode=e.returncode,
                    stderr=e.stderr
                )


class TideHunterRunner(TideHunter):
    """Backward compatibility wrapper for existing code."""
    
    def __init__(self, num_threads=8):
        super().__init__(threads=num_threads)
        self.num_threads = num_threads
    
    def run(self, input_fasta, output_path):
        """Run TideHunter with legacy interface."""
        return self.run_analysis(
            input_file=Path(input_fasta),
            output_file=Path(output_path)
        )

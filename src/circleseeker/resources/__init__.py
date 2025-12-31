"""Resource files and configuration templates."""


def get_default_config() -> str:
    """Return default configuration YAML content."""
    return """# CircleSeeker Configuration File

# Input files (can be overridden by CLI arguments)
input_file: ~
reference: ~
output_dir: "circleseeker_output"
prefix: "sample"

# Feature flags
enable_xecc: true

# Runtime settings
runtime:
  log_level: "WARNING"
  log_file: ~
  tmp_dir: ".tmp_work"
  keep_tmp: false
  checkpoint_policy: "continue"
  enable_progress: true

# Performance settings
performance:
  threads: 8

# External tool parameters
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2
  minimap2_align:
    preset: "sr"
    max_target_seqs: 200
    additional_args: ""
  minimap2:
    preset: "map-hifi"
    additional_args: ""
  samtools: {}
"""

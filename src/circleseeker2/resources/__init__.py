"""Resource files and configuration templates."""

def get_default_config() -> str:
    """Return default configuration YAML content."""
    return '''# CircleSeeker2 Configuration File

# Input files (can be overridden by CLI arguments)
input_file: ~
reference: ~
output_dir: "circleseeker2_output"
prefix: "sample"

# Feature flags
enable_xecc: false

# Runtime settings
runtime:
  log_level: "INFO"
  log_file: ~
  tmp_dir: ".tmp"
  checkpoint_interval: 5
  keep_tmp: false

# Performance settings
performance:
  threads: 8
  max_memory: "16G"
  chunk_size: 10000
  parallel_jobs: 4
  stream_buffer_size: 65536
  enable_profiling: false

# Quality control parameters
quality:
  min_quality_score: 0.99
  min_coverage: 10
  min_eccdna_size: 100
  max_eccdna_size: 1000000
  min_alignment_length: 100
  min_identity: 99.0

# External tool parameters
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2
  blast:
    word_size: 100
    evalue: "1e-50"
    perc_identity: 99.0
    dbtype: "nucl"
  minimap2:
    preset: "map-hifi"
    additional_args: ""
  samtools: {}
  mosdepth:
    min_mapq: 1
'''
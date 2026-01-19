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
    c: 2

  tandem_to_ring:
    min_ave_match: 99.0

  um_classify:
    gap_threshold: 10.0
    # Coverage model parameters (fractions 0-1)
    theta_full: 0.95
    theta_u: 0.95
    theta_m: 0.95
    # U requires the 2nd-best locus coverage to be low; set 1.0 to disable.
    theta_u2_max: 0.05
    # Optional: require minimap2 MAPQ >= this threshold for Uecc (0 disables).
    mapq_u_min: 0
    # Mecc ambiguity filters (0 = disabled)
    mapq_m_ambiguous_threshold: 0
    mecc_identity_gap_threshold: 0
    # Secondary mapping veto for U calls
    u_secondary_min_frac: 0.01
    u_secondary_min_bp: 50
    u_contig_gap_bp: 1000
    u_secondary_max_ratio: 0.05
    u_high_coverage_threshold: 0.98
    u_high_mapq_threshold: 50
    theta_locus: 0.95
    pos_tol_bp: 50
    # Ambiguity interception (fractions 0-1)
    delta_uc: 0.05
    epsilon_mc: 0.05
    # Legacy keys (accepted; prefer theta_full)
    min_full_length_coverage: 95.0
    max_identity_gap_for_mecc: 5.0

  cecc_build:
    overlap_threshold: 0.95
    min_segments: 2
    # Gap tolerance on query (bp)
    edge_tolerance: 20
    tau_gap: 20
    position_tolerance: 50
    # Buffer around query mid-point for doubled-sequence detection (bp)
    half_query_buffer: 50
    locus_overlap_threshold: 0.95
    # Chain coverage threshold (prefer theta_chain; fraction 0-1)
    theta_chain: 0.95
    # Legacy key (percent 0-100)
    min_match_degree: 95.0
    max_rotations: 20

  # General alignment configuration (applies to run_alignment)
  alignment:
    aligner: "minimap2"
    min_identity: 99.0
    min_alignment_length: 50
    db_prefix: ~

  # Minimap2 settings for run_alignment step
  minimap2_align:
    preset: "map-hifi"
    max_target_seqs: 5
    additional_args: ""
    min_identity: 99.0
    identity_decay_per_10kb: 0.5
    min_identity_floor: 97.0
    split_by_length: false
    split_length: 5000
    preset_short: "map-hifi"
    preset_long: "map-hifi"

  # Minimap2 settings for inference mapping (Cyrcular)
  minimap2:
    preset: "map-hifi"
    additional_args: ""

  samtools: {}
"""

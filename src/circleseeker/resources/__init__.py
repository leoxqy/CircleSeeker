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
  um_classify:
    gap_threshold: 10.0
    # Coverage model parameters (fractions 0-1)
    theta_full: 0.95
    # Optional split thresholds (fractions 0-1)
    # - If unset, theta_full is used for both U and M.
    theta_u: 0.95
    theta_m: 0.95
    # U requires the 2nd-best locus coverage to be low; set 1.0 to disable.
    theta_u2_max: 0.05
    # Optional: require minimap2 MAPQ >= this threshold for Uecc (0 disables).
    mapq_u_min: 0
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
    locus_overlap_threshold: 0.95
    # Chain coverage threshold (prefer theta_chain; fraction 0-1)
    theta_chain: 0.95
    # Legacy key (percent 0-100)
    min_match_degree: 95.0
    max_rotations: 20
  minimap2_align:
    preset: "sr"
    max_target_seqs: 200
    additional_args: ""
  minimap2:
    preset: "map-hifi"
    additional_args: ""
  samtools: {}
"""

# RAM-backed tmp_dir (tmpfs) for CircleSeeker

This note explains how to place CircleSeeker temporary files under a RAM-backed
filesystem (for example `/dev/shm`) to reduce I/O time on slow disks.

## When it helps

- Heavy I/O steps (LAST, minimap2 + samtools sort, CD-HIT, large FASTA/TSV).
- HDD or network storage bottlenecks.

## Method 1: Using --turbo mode (Recommended)

The simplest way is to use the `--turbo` command-line argument:

```bash
circleseeker -i reads.fasta -r ref.fa -o results/ --turbo -t 16
```

**Turbo mode features**:
- Automatically checks `/dev/shm` space (default: requires 10GB)
- Creates temporary directory in `/dev/shm`
- Creates `.tmp` symlink in output directory for easy access to intermediate files
- Automatically falls back to normal mode if space is insufficient
- Only effective on Linux systems (macOS has no `/dev/shm`)

**Cleanup behavior**:
- Default: Deletes `/dev/shm` contents and symlink after completion
- `--keep-tmp`: Moves files from `/dev/shm` to `output_dir/tmp/` for retention

**Configuration file method**:

```yaml
runtime:
  turbo_mode: true
  turbo_min_space_gb: 10.0  # Minimum space requirement (GB)
  keep_tmp: false
```

## Method 2: Manual configuration (Advanced)

1) Create a job-specific temp directory in `/dev/shm`.

```bash
export CS_TMPDIR="/dev/shm/cs_${SLURM_JOB_ID:-$$}"
mkdir -p "$CS_TMPDIR"
```

2) Point `runtime.tmp_dir` to that absolute path (absolute paths are allowed).

```yaml
runtime:
  tmp_dir: /dev/shm/cs_JOBID
  keep_tmp: true
```

Note: YAML values are not shell-expanded. Replace `cs_JOBID` with the real path
created in step 1.

3) Keep `output_dir` on disk so final results persist.

```bash
circleseeker -i reads.fasta -r ref.fa -o results/ -c config.yaml -t 16
```

## Cleanup

Temporary files in `/dev/shm` do not persist after a reboot or job end.
If a job fails, manual cleanup is recommended:

```bash
rm -rf /dev/shm/cs_${SLURM_JOB_ID}
```

## Limitations and notes

- `/dev/shm` is volatile. Do not place `output_dir` there unless you are
  prepared to lose final results.
- If `/dev/shm` fills up, the pipeline will fail when writing temp files.

## LAST index handling

The LAST database (used by CeccBuild) is now created **next to the reference
genome** (e.g., `ref.fa` → `ref.lastdb.*`) and reused across runs. This
one-time indexing significantly speeds up repeated analyses on the same
reference.

LAST alignment temporary files (query sequences, alignment results) are still
written to `runtime.tmp_dir` and benefit from RAM-backed storage.

## Quick check

```bash
df -hT /dev/shm
```

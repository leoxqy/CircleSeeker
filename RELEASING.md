# Releasing to Git

This repo contains a runnable pipeline plus large test data. To publish a clean, lightweight Git repository:

## 1) Prepare a clean tree

The provided `.gitignore` excludes caches, temporary workdirs, and very large test artifacts. Optionally, run the cleanup script to remove any generated files from your working tree:

```bash
bash scripts/cleanup_for_git.sh
```

## 2) Initialize and push

```bash
git init
git add -A
git commit -m "chore: initial public release"
git branch -M main
git remote add origin <your-remote-url>
git push -u origin main
```

## 3) (Optional) Tag a release

```bash
git tag -a v2.0.0 -m "CircleSeeker2 2.0.0"
git push --tags
```

## 4) Packaging (sdist/wheel)

```bash
python -m pip install --upgrade build
python -m build
# dist/ now contains .tar.gz and .whl
```

## Notes

- Large test data under `tests/` is ignored by default to keep the repo small. Provide sample datasets in `examples/` (or host separately) if desired.
- The pipeline writes intermediates to `<output_dir>/.tmp_work` and finalizes outputs into `<output_dir>` at the end of a successful run.
- External tools (TideHunter, Blast+, minimap2, samtools, mosdepth, cd-hit) are not bundled. Document installation in README.


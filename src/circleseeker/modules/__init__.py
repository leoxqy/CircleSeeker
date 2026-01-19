"""CircleSeeker analysis modules.

This package provides the core analysis modules for eccDNA detection:

- tandem_to_ring: Convert tandem repeats to circular DNA candidates
- um_classify: Classify eccDNA into U/M/C categories
- cecc_build: Build CeccDNA from chimeric alignments
- umc_process: Process U/M/C classified eccDNA
- ecc_dedup: Deduplicate eccDNA entries
- read_filter: Filter reads based on various criteria
- ecc_unify: Merge and unify eccDNA results
- ecc_summary: Generate summary statistics
- iecc_curator: Curate inferred eccDNA

Legacy naming (for reference):
- carousel         -> tandem_to_ring
- gatekeeper       -> um_classify
- trapeze          -> cecc_build
- menagerie        -> umc_process
- harmonizer       -> ecc_dedup
- sieve            -> read_filter
"""

from . import tandem_to_ring
from . import um_classify
from . import cecc_build
from . import umc_process
from . import ecc_dedup
from . import read_filter
from . import ecc_unify
from . import ecc_summary
from . import iecc_curator

__all__ = [
    "tandem_to_ring",
    "um_classify",
    "cecc_build",
    "umc_process",
    "ecc_dedup",
    "read_filter",
    "ecc_unify",
    "ecc_summary",
    "iecc_curator",
]

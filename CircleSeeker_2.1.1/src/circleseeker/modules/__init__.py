"""CircleSeeker analysis modules.

Renaming mapping (Python-safe underscores):

- carousel         -> tandem_to_ring
- gatekeeper       -> um_classify
- trapeze          -> cecc_build
- menagerie        -> umc_process
- harmonizer       -> ecc_dedup
- sieve            -> read_filter
- playbill         -> report_generator
- contortionist    -> (removed; superseded by cyrcular_calling step)
- juggler          -> (removed; validation folded into cyrcular_calling outputs)

Old imports remain valid for active modules; retired entries listed below are
no longer exported. New module names should be used going forward in internal
imports.
"""

# Re-export alias modules for discoverability
from . import tandem_to_ring as _tandem_to_ring  # noqa: F401
from . import um_classify as _um_classify  # noqa: F401
from . import cecc_build as _cecc_build  # noqa: F401
from . import umc_process as _umc_process  # noqa: F401
from . import ecc_dedup as _ecc_dedup  # noqa: F401
from . import read_filter as _read_filter  # noqa: F401
# from . import report_generator as _report_generator  # noqa: F401  # Module removed

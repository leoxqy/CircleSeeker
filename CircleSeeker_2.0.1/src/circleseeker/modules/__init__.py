"""CircleSeeker analysis modules.

Alias mapping (Python-safe underscores):

- tandem_to_ring -> carousel
- um_classify    -> gatekeeper
- cecc_build     -> trapeze
- umc_process    -> menagerie
- ecc_dedup      -> harmonizer
- read_filter    -> sieve
- split_refine   -> contortionist
- circle_validate -> juggler
- report_generator -> playbill
- file_organizer   -> propmaster

Old imports remain valid; the new names provide clearer semantics.
"""

# Re-export alias modules for discoverability
from . import carousel as _carousel  # noqa: F401
from . import gatekeeper as _gatekeeper  # noqa: F401
from . import trapeze as _trapeze  # noqa: F401
from . import menagerie as _menagerie  # noqa: F401
from . import harmonizer as _harmonizer  # noqa: F401
from . import sieve as _sieve  # noqa: F401
from . import contortionist as _contortionist  # noqa: F401
from . import juggler as _juggler  # noqa: F401
from . import playbill as _playbill  # noqa: F401
from . import propmaster as _propmaster  # noqa: F401

"""Enhanced progress helpers with better terminal compatibility."""

from __future__ import annotations

import sys
import os
from typing import Iterable, TypeVar, Iterator, Optional

T = TypeVar("T")


def iter_progress(
    iterable: Iterable[T],
    total: Optional[int] = None,
    desc: Optional[str] = None,
    enabled: bool = True,
    file=None,
) -> Iterator[T]:
    """Wrap iterable with tqdm if available and enabled, else return as-is.

    Enhanced version with better terminal detection and fallback.
    """
    # Check if progress should be disabled
    if not enabled:
        return iter(iterable)

    # Disable in non-TTY environments (pipes, redirects)
    if not sys.stdout.isatty():
        return iter(iterable)

    # Check environment override
    if os.environ.get('CIRCLESEEKER_NO_PROGRESS', '').lower() in ('1', 'true', 'yes'):
        return iter(iterable)

    try:
        from tqdm import tqdm

        # Detect terminal capabilities
        term = os.environ.get('TERM', '')
        ncols = None

        # Use auto-detection for most terminals
        if term and term != 'dumb':
            ncols = None  # Let tqdm auto-detect
        else:
            ncols = 80  # Fixed width for dumb terminals

        # Better formatted progress bar
        formatted_desc = f"{desc}" if desc else ""

        # Use file parameter for output control
        output_file = file if file is not None else sys.stderr

        return iter(tqdm(
            iterable,
            total=total,
            desc=formatted_desc,
            ncols=ncols,
            file=output_file,
            ascii=False if 'utf' in sys.stdout.encoding.lower() else True,
            leave=True,
            dynamic_ncols=True,
            smoothing=0.1,  # Smoother updates
            miniters=1,      # Update every iteration
        ))
    except ImportError:
        # Fallback: simple percentage printer
        if desc and total:
            print(f"Starting: {desc} (0/{total})", file=sys.stderr)

            class SimpleProgress:
                def __init__(self, iterable, total, desc):
                    self.iterable = iterable
                    self.total = total
                    self.desc = desc
                    self.count = 0

                def __iter__(self):
                    for item in self.iterable:
                        self.count += 1
                        if self.count % max(1, self.total // 20) == 0 or self.count == self.total:
                            pct = (self.count / self.total) * 100
                            print(f"\r{self.desc}: {pct:.0f}% ({self.count}/{self.total})",
                                  end='', file=sys.stderr)
                        yield item
                    print(file=sys.stderr)  # New line after completion

            return iter(SimpleProgress(iterable, total, desc))
        else:
            return iter(iterable)
    except Exception:
        # Any other error, just return the iterable
        return iter(iterable)
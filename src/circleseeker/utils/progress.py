"""Progress helpers (tqdm integration with graceful fallback)."""

from __future__ import annotations

from typing import Iterable, TypeVar, Iterator, Optional

T = TypeVar("T")


def iter_progress(
    iterable: Iterable[T],
    total: Optional[int] = None,
    desc: Optional[str] = None,
    enabled: bool = True,
) -> Iterator[T]:
    """Wrap iterable with tqdm if available and enabled, else return as-is."""
    if not enabled:
        return iter(iterable)
    try:
        from tqdm import tqdm

        # Format with bullet point and proper alignment
        # Use shorter bar to align with status line
        formatted_desc = f"· {desc:<12} " if desc else ""
        return iter(
            tqdm(
                iterable,
                total=total,
                desc=formatted_desc,
                bar_format="{desc}: {percentage:3.0f}%|{bar:30}| {n_fmt}/{total_fmt}",
                ncols=80,
            )
        )  # Fixed width for alignment
    except Exception:
        return iter(iterable)

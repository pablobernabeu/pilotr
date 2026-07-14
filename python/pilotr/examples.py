"""Locate the design specifications shipped with pilotr.

The JSON files under ``pilotr/examples/`` are packaged copies of the
repository's ``spec/examples/``, shipped inside the wheel so they are reachable
from an installed copy. They are the same files that drive the R twin and the
no-code app, so a design authored once runs unchanged across all three.
"""

from __future__ import annotations

from importlib import resources
from typing import Optional, Union


def _examples_dir():
    return resources.files("pilotr") / "examples"


def _available() -> list[str]:
    return sorted(
        entry.name[:-5]
        for entry in _examples_dir().iterdir()
        if entry.name.endswith(".json")
    )


def pilotr_example(name: Optional[str] = None) -> Union[list[str], str]:
    """List the bundled example specifications, or return the path to one.

    pilotr ships one ready-to-run specification per design family, as JSON. These
    are the same files that drive the R twin and the no-code app.

    Parameters
    ----------
    name : str, optional
        The base name of an example, with or without the ``.json`` extension, for
        example ``"between_2group_gaussian"``. When ``None`` (the default), the
        available example names are returned instead of a path.

    Returns
    -------
    list of str, or str
        When ``name`` is ``None``, the available example names. Otherwise, the
        path to that example's JSON file, ready to pass to
        :func:`pilotr.simulate.load_spec`.
    """
    available = _available()
    if name is None:
        return available
    base = name[:-5] if name.endswith(".json") else name
    if base not in available:
        raise ValueError(
            f"Unknown example {name!r}. Available: {', '.join(available)}."
        )
    return str(_examples_dir() / f"{base}.json")

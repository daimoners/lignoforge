"""On-disk helpers for LignoForge projects (JSON files in plain folders)."""

from __future__ import annotations

import json
import os
import re
import tempfile
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1

_id_lock = threading.Lock()


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_json(path: Path) -> Any:
    with open(path) as fh:
        return json.load(fh)


def write_json(path: Path, data: Any, indent: int = 2) -> None:
    """Atomic write: a crash never leaves a half-written record."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(data, fh, indent=indent)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def next_id(folder: Path, prefix: str) -> str:
    """Next sequential id ``<prefix>0001`` inside *folder* (thread-safe)."""
    with _id_lock:
        folder.mkdir(parents=True, exist_ok=True)
        used = [int(m.group(1)) for p in folder.iterdir()
                if (m := re.fullmatch(rf"{prefix}(\d+)", p.name))]
        new = f"{prefix}{(max(used) + 1 if used else 1):04d}"
        (folder / new).mkdir()
        return new


def lignoforge_version() -> str:
    try:
        from importlib.metadata import version
        return version("lignoforge")
    except Exception:                       # noqa: BLE001
        return "unknown"

"""Disk-bounded contracts for exact streaming-table validation."""

from __future__ import annotations

import sqlite3
from types import TracebackType


class DiskBackedUniqueIds:
    """Exact uniqueness set backed by an automatically deleted SQLite file."""

    def __init__(self) -> None:
        self._connection: sqlite3.Connection | None = None

    def __enter__(self) -> "DiskBackedUniqueIds":
        # An empty SQLite filename creates a temporary on-disk database that
        # SQLite removes automatically when the connection closes.
        connection = sqlite3.connect("")
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute("PRAGMA cache_size=-8192")
        connection.execute(
            "CREATE TABLE identifiers (value TEXT PRIMARY KEY) WITHOUT ROWID"
        )
        self._connection = connection
        return self

    def add(self, value: str) -> bool:
        if self._connection is None:
            raise RuntimeError("disk-backed ID set is not open")
        try:
            self._connection.execute(
                "INSERT INTO identifiers(value) VALUES (?)", (value,)
            )
        except sqlite3.IntegrityError:
            return False
        return True

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        if self._connection is not None:
            self._connection.close()
            self._connection = None

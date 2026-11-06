"""Output-contained run scratch and replayable diagnostic rows.

``ScratchWorkspace`` owns one hidden directory below a workflow's output
directory for the private files of a run or of one query batch; annotation
preparation, LD-score runs, gene-index builds and scoring, partitioned-h2 and
quantile-h2 all use it. ``FrameSpool`` appends bounded frames to one file in
such a directory and replays them.

Diagnostic rows that a writer publishes share the ``DiagnosticRows``
protocol, whether they are spooled (``FrameSpool``), staged in a compressed
TSV file (``TsvDiagnostics``) or held in memory (``MemoryDiagnostics``):
``frames()`` replays the rows in bounded frames and ``write_to(path)``
publishes them as one TSV file.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import pickle
import shutil
import tempfile
from typing import Iterator, Protocol

import pandas as pd

from .paths import output_directory_path


class ScratchWorkspace:
    """Own a unique scratch directory below an explicit workflow destination.

    The directory is ``.ldsc-annotation-<random>`` immediately below the
    destination, which is created when missing. Closing is idempotent and
    removes only this owner's private directory. Public artifacts and source
    files are never owned by this object. Borrowers must not close their
    caller's workspace.
    """

    def __init__(self, output_dir: str | Path):
        destination = output_directory_path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)
        self.path = Path(tempfile.mkdtemp(prefix=".ldsc-annotation-", dir=destination))
        self.closed = False

    def require_open(self) -> None:
        if self.closed:
            raise ValueError("Annotation workspace is closed.")

    def close(self) -> None:
        if not self.closed:
            shutil.rmtree(self.path)
            self.closed = True

    def __enter__(self):
        self.require_open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class DiagnosticRows(Protocol):
    """Complete diagnostic rows that a writer replays or publishes."""

    def frames(self) -> Iterator[pd.DataFrame]:
        """Yield the rows in order as bounded frames."""

    def write_to(self, path) -> None:
        """Publish every row as one TSV file at ``path``, gzip-compressed for a ``.gz`` suffix."""


@dataclass
class FrameSpool:
    """Append bounded private frames and replay them without retaining records.

    One append stream holds independently serialized frames, keeping file count
    constant as chunks accumulate. Pickles are locally generated scratch, never
    user inputs. Descriptors retain no records or open handles. ``columns`` is
    the header ``write_to`` writes for a spool that received no rows.
    """

    path: Path
    n_parts: int = 0
    n_rows: int = 0
    columns: tuple[str, ...] = ()

    def append(self, frame):
        if frame.empty:
            return
        self.path.mkdir(parents=True, exist_ok=True)
        with (self.path / "frames.pkl").open("ab" if self.n_parts else "wb") as stream:
            pickle.dump(frame, stream, protocol=pickle.HIGHEST_PROTOCOL)
        self.n_parts += 1
        self.n_rows += len(frame)

    def frames(self):
        if self.n_parts:
            with (self.path / "frames.pkl").open("rb") as stream:
                for _ in range(self.n_parts):
                    yield pickle.load(stream)

    def write_to(self, path):
        """Publish all records, retaining only one append chunk at a time."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt") as stream:
            header = True
            for frame in self.frames():
                frame.to_csv(stream, sep="\t", index=False, header=header, na_rep="")
                header = False
            if header:
                pd.DataFrame(columns=list(self.columns)).to_csv(stream, sep="\t", index=False)


@dataclass(frozen=True)
class TsvDiagnostics:
    """Complete compressed diagnostics with an explicit private or public owner.

    A private artifact borrows its scratch workspace and is readable until
    that owner closes. Publication replaces it with the persistent destination
    and no workspace dependency. Reads never cache complete diagnostic tables.
    """

    path: Path
    workspace: ScratchWorkspace | None = None

    def frames(self, *, chunk_rows=65536):
        """Replay bounded rows while checking the private owner's lifetime."""
        if self.workspace is not None:
            self.workspace.require_open()
        yield from pd.read_csv(self.path, sep="\t", chunksize=chunk_rows)

    def write_to(self, path):
        """Copy the complete artifact without loading or reserializing rows."""
        if self.workspace is not None:
            self.workspace.require_open()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.resolve() != self.path.resolve():
            shutil.copyfile(self.path, path)


@dataclass(frozen=True)
class MemoryDiagnostics:
    """Complete diagnostic rows owned by an in-memory result."""

    frame: pd.DataFrame

    def frames(self, *, chunk_rows=65536):
        for start in range(0, len(self.frame), chunk_rows):
            yield self.frame.iloc[start:start + chunk_rows]

    def write_to(self, path):
        """Serialize only when the caller explicitly chooses a writing workflow."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.frame.to_csv(path, sep="\t", index=False, na_rep="")

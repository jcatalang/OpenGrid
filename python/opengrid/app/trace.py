"""The session's Python trace: everything done in the window or its console, as code.

The session journal (``commands/base.py``) records every command, undo and redo as a call of
the public API. The trace writes it to a file as it grows, so the operations done by hand
can be read (and replayed) as a script. Unlike the journal, it lasts as long as the window:
File > New and Open continue it (``og.new()``, ``og.load(...)``), and the window adds what
changes no model (exports, saves: ``note``). The file is a temporary one
(``<temp>/opengrid-traces/``, files older than a week removed at the next start).
"""

import os
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

TRACE_DIR = Path(tempfile.gettempdir()) / "opengrid-traces"
KEEP_SECONDS = 7 * 24 * 3600


def new_trace_path(directory: Path | None = None) -> Path:
    """A fresh file name in the trace folder (old traces in it are removed)."""
    if directory is None:
        directory = TRACE_DIR
    directory.mkdir(parents=True, exist_ok=True)
    now = time.time()
    for old in directory.glob("session-*.py"):
        try:
            if now - old.stat().st_mtime > KEEP_SECONDS:
                old.unlink()
        except OSError:
            pass
    stamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(now))
    return directory / f"session-{stamp}-{os.getpid()}.py"


class SessionTrace:
    """Follows a session's journal into a file; ``on_append(text)`` gets each addition."""

    def __init__(self, session, path: Path):
        self.session = session
        self.path = Path(path)
        self.on_append: Callable[[str], None] | None = None
        self._journal = None
        self._written = 0
        started = time.strftime("%Y-%m-%d %H:%M:%S")
        header = (
            f"# OpenGrid session trace, started {started}: every command run in the window or\n"
            "# its console, as the Python that redoes it (run this file to replay it).\n"
            "import opengrid\n\nog = opengrid.Session()\n"
        )
        self.path.write_text(header, encoding="utf-8")
        session.commands.add_listener(self._on_change)
        self.sync()

    def text(self) -> str:
        return self.path.read_text(encoding="utf-8")

    def sync(self) -> None:
        """Write the journal lines not written yet (after File > New / Open: from the new
        journal, which starts with the load)."""
        journal = self.session.commands.journal
        if journal is not self._journal:
            replaced = self._journal is not None
            self._journal = journal
            self._written = 0
            if replaced and not (journal.lines and journal.lines[0].startswith("og.load(")):
                self._append(["og.new()"])
        lines = journal.lines[self._written :]
        if lines:
            self._written = len(journal.lines)
            self._append(lines)

    def note(self, source: str) -> None:
        """Add an operation that changes no model (an export, a save)."""
        self.sync()
        self._append([source])

    def _on_change(self, _change) -> None:
        self.sync()

    def _append(self, lines: list[str]) -> None:
        text = "\n".join(lines) + "\n"
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(text)
        if self.on_append is not None:
            self.on_append(text)

"""Command pattern core: Command, CommandStack (undo/redo) and the Python journal.

A command mutates the store inside a transaction. The resulting :class:`Change` (before/after
snapshots of the touched entities) is what undo and redo restore, so redo never re-runs the
algorithm and both are exact. After a successful run, the command renders itself as Python
source against the public API; the journal concatenates those lines into a script that
rebuilds the session from scratch.
"""

from abc import ABC, abstractmethod
from collections.abc import Callable

from opengrid.model.store import Change, Store

JOURNAL_HEADER = "import opengrid\n\nog = opengrid.Session()\n"


class Command(ABC):
    label = "command"

    @abstractmethod
    def execute(self, store: Store) -> object:
        """Apply the change (inside an open transaction); return a result for the caller."""

    @abstractmethod
    def script(self) -> str:
        """Python source reproducing this command (valid after ``execute``)."""


class Journal:
    """Replayable Python script of every command, undo and redo of a session."""

    def __init__(self):
        self.lines: list[str] = []

    def append(self, source: str) -> None:
        self.lines.append(source)

    def clear(self) -> None:
        self.lines.clear()

    def text(self) -> str:
        body = ""
        if self.lines:
            body = "\n".join(self.lines) + "\n"
        return JOURNAL_HEADER + body


class CommandStack:
    def __init__(self, store: Store):
        self.store = store
        self.journal = Journal()
        self._undo: list[tuple[Command, Change]] = []
        self._redo: list[tuple[Command, Change]] = []
        self._listeners: list[Callable[[Change], None]] = []

    def add_listener(self, listener: Callable[[Change], None]) -> None:
        """Call ``listener(change)`` after every run, undo and redo."""
        self._listeners.append(listener)

    def remove_listener(self, listener: Callable[[Change], None]) -> None:
        self._listeners.remove(listener)

    def run(self, command: Command) -> object:
        """Execute a command; on any error the store is rolled back and the error re-raised."""
        from opengrid.model.sketch import check_ownership

        self.store.begin()
        try:
            result = command.execute(self.store)
            check_ownership(self.store, self.store.touched())
        except BaseException:
            self.store.rollback()
            raise
        change = self.store.commit()
        self._undo.append((command, change))
        self._redo.clear()
        self.journal.append(command.script())
        self._notify(change)
        return result

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo_label(self) -> str | None:
        if self._undo:
            return self._undo[-1][0].label
        return None

    def redo_label(self) -> str | None:
        if self._redo:
            return self._redo[-1][0].label
        return None

    def undo(self) -> None:
        if not self._undo:
            raise IndexError("nothing to undo")
        command, change = self._undo.pop()
        self.store.apply_states(change.before)
        self._redo.append((command, change))
        self.journal.append("og.undo()")
        self._notify(change)

    def redo(self) -> None:
        if not self._redo:
            raise IndexError("nothing to redo")
        command, change = self._redo.pop()
        self.store.apply_states(change.after)
        self._undo.append((command, change))
        self.journal.append("og.redo()")
        self._notify(change)

    def clear(self) -> None:
        """Forget the undo/redo history (e.g. after loading a project)."""
        self._undo.clear()
        self._redo.clear()

    def _notify(self, change: Change) -> None:
        for listener in list(self._listeners):
            listener(change)

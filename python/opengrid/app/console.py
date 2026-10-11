"""Embedded Python console executing in the application's process and namespace."""

import code
import contextlib
import io

from PySide6.QtCore import Qt
from PySide6.QtGui import QFontDatabase, QKeyEvent
from PySide6.QtWidgets import QLineEdit, QPlainTextEdit, QVBoxLayout, QWidget

PROMPT = ">>> "
CONTINUATION_PROMPT = "... "


class _HistoryLineEdit(QLineEdit):
    """Single-line input with Up/Down command history."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.history: list[str] = []
        self._history_index = 0

    def remember(self, line: str) -> None:
        if line.strip():
            self.history.append(line)
        self._history_index = len(self.history)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802 (Qt override)
        if event.key() == Qt.Key.Key_Up:
            if self._history_index > 0:
                self._history_index -= 1
                self.setText(self.history[self._history_index])
            return
        if event.key() == Qt.Key.Key_Down:
            if self._history_index < len(self.history) - 1:
                self._history_index += 1
                self.setText(self.history[self._history_index])
            else:
                self._history_index = len(self.history)
                self.clear()
            return
        super().keyPressEvent(event)


class PythonConsole(QWidget):
    """Interactive Python console widget.

    Commands run synchronously on the GUI thread; long operations must be started through
    the worker-thread API (Phase 3+), not from here.
    """

    def __init__(self, namespace: dict | None = None, parent=None):
        super().__init__(parent)
        if namespace is None:
            namespace = {}
        self.namespace = namespace
        self._interpreter = code.InteractiveConsole(locals=self.namespace)
        self._more = False

        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        self.output = QPlainTextEdit(self)
        self.output.setReadOnly(True)
        self.output.setFont(font)
        self.input = _HistoryLineEdit(self)
        self.input.setFont(font)
        self.input.setPlaceholderText("Python command")
        self.input.returnPressed.connect(self._on_return)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.output)
        layout.addWidget(self.input)

        self.write("OpenGrid Python console. Namespace: " + ", ".join(sorted(self.namespace)))

    def write(self, text: str) -> None:
        self.output.appendPlainText(text.rstrip("\n"))

    def run_line(self, line: str) -> str:
        """Execute one line as if typed; return the captured stdout/stderr text."""
        if self._more:
            prompt = CONTINUATION_PROMPT
        else:
            prompt = PROMPT
        self.write(prompt + line)
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
            self._more = self._interpreter.push(line)
        captured = buffer.getvalue()
        if captured:
            self.write(captured)
        return captured

    def _on_return(self) -> None:
        line = self.input.text()
        self.input.remember(line)
        self.input.clear()
        self.run_line(line)

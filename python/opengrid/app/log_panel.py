"""The Log tab: what the program did and reported (Python logging, the "opengrid" logger at
INFO and above), one line each: the time and the message, warnings and errors in colour, a
message repeated at once counted on its line instead of written again."""

import html
import logging
import time

from PySide6.QtCore import QObject, Signal
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import QPlainTextEdit

LEVEL_STYLE = {
    logging.WARNING: ("#e8c547", "warning: "),
    logging.ERROR: ("#ff6b6b", "not done: "),
}


class _Emitter(QObject):
    message = Signal(int, str, float)


class QtLogHandler(logging.Handler):
    """logging.Handler forwarding records to a Qt signal (safe to call from worker threads)."""

    def __init__(self):
        super().__init__(level=logging.INFO)
        self.emitter = _Emitter()

    def emit(self, record: logging.LogRecord) -> None:
        try:
            text = record.getMessage()
        except (TypeError, ValueError):
            text = str(record.msg)
        self.emitter.message.emit(record.levelno, text, record.created)


class LogPanel(QPlainTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setMaximumBlockCount(10000)
        self._last: tuple[int, str] | None = None
        self._repeats = 0
        self.handler = QtLogHandler()
        self.handler.emitter.message.connect(self.add)

    def add(self, level: int, text: str, created: float | None = None) -> None:
        """Append a message (or count a repeat of the last one)."""
        if created is None:
            created = time.time()
        stamp = time.strftime("%H:%M:%S", time.localtime(created))
        if self._last == (level, text):
            self._repeats += 1
            self._replace_last(self._line(stamp, level, f"{text} ({self._repeats + 1} times)"))
            return
        self._last = (level, text)
        self._repeats = 0
        self.appendHtml(self._line(stamp, level, text))

    @staticmethod
    def _line(stamp: str, level: int, text: str) -> str:
        body = html.escape(text).replace("\n", "<br>&nbsp;&nbsp;&nbsp;&nbsp;")
        style = None
        if level >= logging.ERROR:
            style = LEVEL_STYLE[logging.ERROR]
        elif level >= logging.WARNING:
            style = LEVEL_STYLE[logging.WARNING]
        line = f'<span style="color:#8b949e;">{stamp}</span>&nbsp;&nbsp;'
        if style is None:
            return line + body
        colour, prefix = style
        return line + f'<span style="color:{colour};">{prefix}{body}</span>'

    def _replace_last(self, line: str) -> None:
        if self.document().blockCount() <= 1:
            self.clear()
            self.appendHtml(line)
            return
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.select(QTextCursor.SelectionType.BlockUnderCursor)
        cursor.removeSelectedText()
        self.appendHtml(line)

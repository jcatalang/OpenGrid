"""The message bar above the view: every warning and error the program reports (records of the
"opengrid" logger at WARNING and above, from any thread) is shown there, where it is seen, not
only in the Log tab. Messages that arrive together (one operation) are shown together; the bar
closes with its button, or by itself after a reading time (longer for longer text; not while
the pointer is over it). A record logged with ``extra={"on_screen": False}`` (shown by a dialog
already) is left to the Log tab."""

import logging

from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QToolButton, QVBoxLayout

MAX_LINES = 6  # more are counted, and listed in the Log tab
BURST_SECONDS = 1.0  # messages this close together are one burst (one operation)
COLOURS = {
    logging.WARNING: ("#4a3c12", "#e8c547", "Warning"),
    logging.ERROR: ("#4f1c21", "#ff6b6b", "Not done"),
}


class _Emitter(QObject):
    record = Signal(int, str)


class BarHandler(logging.Handler):
    """Forwards WARNING and ERROR records to the bar (safe from worker threads)."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.emitter = _Emitter()

    def emit(self, record: logging.LogRecord) -> None:
        if getattr(record, "on_screen", True) is False:
            return
        try:
            text = record.getMessage()
        except (TypeError, ValueError):
            text = str(record.msg)
        self.emitter.record.emit(record.levelno, text)


class MessageBar(QFrame):
    """The bar; ``log_requested`` asks the window to show the Log tab."""

    log_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("messageBar")
        self._lines: list[str] = []
        self._level = logging.WARNING
        self._burst = QTimer(self)  # while running, new messages join the shown ones
        self._burst.setSingleShot(True)
        self._burst.setInterval(int(BURST_SECONDS * 1000))
        self._close_timer = QTimer(self)
        self._close_timer.setSingleShot(True)
        self._close_timer.timeout.connect(self.hide)
        self._title = QLabel()
        self._title.setStyleSheet("font-weight: bold;")
        self._text = QLabel()
        self._text.setWordWrap(True)
        self._text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        log_button = QToolButton()
        log_button.setText("Log")
        log_button.setToolTip("Show every message in the Log tab")
        log_button.clicked.connect(self.log_requested)
        close_button = QToolButton()
        close_button.setText("✕")
        close_button.setToolTip("Close")
        close_button.clicked.connect(self.hide)
        text_column = QVBoxLayout()
        text_column.setSpacing(2)
        text_column.addWidget(self._title)
        text_column.addWidget(self._text)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(10, 6, 6, 6)
        layout.addLayout(text_column, 1)
        layout.addWidget(log_button, 0, Qt.AlignmentFlag.AlignTop)
        layout.addWidget(close_button, 0, Qt.AlignmentFlag.AlignTop)
        self.handler = BarHandler()
        self.handler.emitter.record.connect(self.post)
        self.hide()

    def lines(self) -> list[str]:
        """The messages shown (for tests)."""
        return list(self._lines)

    def post(self, level: int, text: str) -> None:
        """Show a message (a new burst replaces the messages shown before)."""
        if not self._burst.isActive() or not self.isVisible():
            self._lines = []
            self._level = level
        self._burst.start()
        for line in text.splitlines():
            line = line.strip()
            if line and line not in self._lines:
                self._lines.append(line)
        self._level = max(self._level, level)
        self._show()

    def _show(self) -> None:
        level = logging.WARNING
        if self._level >= logging.ERROR:
            level = logging.ERROR
        background, border, title = COLOURS[level]
        self.setStyleSheet(
            f"QFrame#messageBar {{ background: {background}; border-left: 4px solid {border}; }}"
            f" QFrame#messageBar QLabel {{ color: #f2f2f2; }}"
        )
        shown = self._lines[:MAX_LINES]
        text = "\n".join(shown)
        if len(self._lines) > MAX_LINES:
            text += f"\n… and {len(self._lines) - MAX_LINES} more (see the Log tab)"
        self._title.setText(title)
        self._title.setStyleSheet(f"font-weight: bold; color: {border};")
        self._text.setText(text)
        self.show()
        self._close_timer.start(self._reading_ms(text))

    @staticmethod
    def _reading_ms(text: str) -> int:
        return int(min(45_000, max(8_000, 4_000 + 60 * len(text))))

    def enterEvent(self, event) -> None:  # noqa: N802 (Qt)
        self._close_timer.stop()
        super().enterEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.isVisible():
            self._close_timer.start(4_000)
        super().leaveEvent(event)

"""Log / message panel connected to Python's logging module."""

import logging

from PySide6.QtCore import QObject, Signal
from PySide6.QtWidgets import QPlainTextEdit


class _Emitter(QObject):
    message = Signal(str)


class QtLogHandler(logging.Handler):
    """logging.Handler forwarding records to a Qt signal (safe to call from worker threads)."""

    def __init__(self):
        super().__init__()
        self.emitter = _Emitter()
        self.setFormatter(
            logging.Formatter(
                "%(asctime)s %(levelname)-7s %(name)s: %(message)s", datefmt="%H:%M:%S"
            )
        )

    def emit(self, record: logging.LogRecord) -> None:
        self.emitter.message.emit(self.format(record))


class LogPanel(QPlainTextEdit):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setMaximumBlockCount(10000)
        self.handler = QtLogHandler()
        self.handler.emitter.message.connect(self.appendPlainText)

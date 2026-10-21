"""The Script panel: the session's Python trace (``app/trace.py``) as it grows, beside the
console. Copy it, save it as a script, open the file in an editor."""

from pathlib import Path

from PySide6.QtCore import Qt, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFontDatabase, QGuiApplication, QTextCursor
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class ScriptPanel(QWidget):
    """Read-only view of the trace file; ``save_requested`` asks the window for a path."""

    save_requested = Signal()

    def __init__(self, trace, parent=None):
        super().__init__(parent)
        self.trace = trace
        self.view = QPlainTextEdit(self)
        self.view.setReadOnly(True)
        self.view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.view.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self._highlighter = _python_highlighter(self.view)
        self.path_label = QLabel(self)
        self.path_label.setWordWrap(True)
        self.path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.copy_button = QPushButton("Copy")
        self.copy_button.setToolTip("Copy the whole script")
        self.copy_button.clicked.connect(self.copy)
        self.save_button = QPushButton("Save as…")
        self.save_button.setToolTip("Save the script (a replayable Python file)")
        self.save_button.clicked.connect(self.save_requested)
        self.open_button = QPushButton("Open file")
        self.open_button.setToolTip("Open the trace file in the system's editor for .py files")
        self.open_button.clicked.connect(self.open_file)
        buttons = QHBoxLayout()
        buttons.addWidget(self.copy_button)
        buttons.addWidget(self.save_button)
        buttons.addWidget(self.open_button)
        buttons.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.path_label)
        layout.addWidget(self.view, 1)
        layout.addLayout(buttons)
        self.set_trace(trace)

    def set_trace(self, trace) -> None:
        self.trace = trace
        trace.on_append = self.append
        self.path_label.setText(
            f"Everything done in this window, as Python (also <code>trace_path</code> in the "
            f"console):<br><code>{trace.path}</code>"
        )
        self.view.setPlainText(trace.text())
        self._scroll_to_end()

    def text(self) -> str:
        return self.view.toPlainText()

    def append(self, text: str) -> None:
        cursor = self.view.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.insertText(text)
        self._scroll_to_end()

    def copy(self) -> None:
        QGuiApplication.clipboard().setText(self.trace.text())

    def open_file(self) -> bool:
        return QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.trace.path))))

    def _scroll_to_end(self) -> None:
        bar = self.view.verticalScrollBar()
        bar.setValue(bar.maximum())


def _python_highlighter(view: QPlainTextEdit):
    """Pygments colours through qtconsole's highlighter, when installed (else plain text)."""
    try:
        from pygments.lexers import PythonLexer
        from qtconsole.pygments_highlighter import PygmentsHighlighter
    except ImportError:
        return None
    highlighter = PygmentsHighlighter(view.document(), lexer=PythonLexer())
    if view.palette().window().color().lightness() < 128:
        highlighter.set_style("monokai")
    return highlighter

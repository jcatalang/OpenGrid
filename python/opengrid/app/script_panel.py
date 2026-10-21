"""The Script panel: the session's Python trace (``app/trace.py``) as it grows, beside the
console, and new code to run in the window. The lines already run are kept as they are
(read-only); lines written after them run with Ctrl+Enter or Run, in the window's model,
and then join the trace as written. Copy it, save it as a script, open the file in an
editor."""

from pathlib import Path

from PySide6.QtCore import QEvent, Qt, QUrl, Signal
from PySide6.QtGui import (
    QColor,
    QDesktopServices,
    QFontDatabase,
    QGuiApplication,
    QKeySequence,
    QTextCharFormat,
    QTextCursor,
    QTextFormat,
)
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMenu,
    QPlainTextEdit,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

PENDING_DARK = "#1f2a36"  # background of the lines not run yet
PENDING_LIGHT = "#eef4fb"


class ScriptEditor(QPlainTextEdit):
    """A script whose first part (what has run) cannot be edited; the rest is new code.
    ``run_requested`` on Ctrl+Enter / Shift+Enter."""

    run_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.committed_end = 0  # characters [0, committed_end) have run
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.setUndoRedoEnabled(False)  # its own undo would undo the trace's lines
        self.setAcceptDrops(False)
        self.textChanged.connect(self._mark_pending)

    # --------------------------------------------------------------- content
    def pending(self) -> str:
        return self.toPlainText()[self.committed_end :]

    def set_content(self, committed: str, pending: str = "") -> None:
        self.blockSignals(True)
        self.setPlainText(committed + pending)
        self.blockSignals(False)
        self.committed_end = len(committed)
        self._mark_pending()

    def append_committed(self, text: str) -> None:
        """Lines that ran (from the window or from here) go before the code not run yet."""
        cursor = QTextCursor(self.document())
        cursor.setPosition(self.committed_end)
        self.blockSignals(True)
        cursor.insertText(text)
        self.blockSignals(False)
        self.committed_end += len(text)
        self._mark_pending()

    def replace_pending(self, text: str) -> None:
        cursor = QTextCursor(self.document())
        cursor.setPosition(self.committed_end)
        cursor.movePosition(QTextCursor.MoveOperation.End, QTextCursor.MoveMode.KeepAnchor)
        cursor.insertText(text)
        self._mark_pending()

    # ------------------------------------------------- the read-only first part
    def _touches_committed(self) -> bool:
        cursor = self.textCursor()
        return min(cursor.selectionStart(), cursor.position()) < self.committed_end

    def _to_end(self) -> None:
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.setTextCursor(cursor)

    def keyPressEvent(self, event) -> None:  # noqa: N802 (Qt override)
        key = event.key()
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and event.modifiers() & (
            Qt.KeyboardModifier.ControlModifier
            | Qt.KeyboardModifier.MetaModifier
            | Qt.KeyboardModifier.ShiftModifier
        ):
            self.run_requested.emit()
            return
        if event.matches(QKeySequence.StandardKey.Copy) or event.matches(
            QKeySequence.StandardKey.SelectAll
        ):
            super().keyPressEvent(event)
            return
        edits = (
            bool(event.text())
            or key in (Qt.Key.Key_Backspace, Qt.Key.Key_Delete)
            or event.matches(QKeySequence.StandardKey.Paste)
            or event.matches(QKeySequence.StandardKey.Cut)
        )
        if edits and self._touches_committed():
            if key in (Qt.Key.Key_Backspace, Qt.Key.Key_Delete) or event.matches(
                QKeySequence.StandardKey.Cut
            ):
                return  # nothing that ran is deleted
            self._to_end()  # typing over what ran: written at the end instead
        cursor = self.textCursor()
        if (
            key == Qt.Key.Key_Backspace
            and not cursor.hasSelection()
            and cursor.position() <= self.committed_end
        ):
            return
        super().keyPressEvent(event)

    def insertFromMimeData(self, source) -> None:  # noqa: N802 (Qt override)
        if self._touches_committed():
            self._to_end()
        self.insertPlainText(source.text())

    def contextMenuEvent(self, event) -> None:  # noqa: N802 (Qt override)
        menu = QMenu(self)
        menu.addAction("Copy", self.copy).setEnabled(self.textCursor().hasSelection())
        menu.addAction("Paste", self.paste)
        menu.addAction("Select all", self.selectAll)
        menu.addSeparator()
        menu.addAction("Run the new lines (Ctrl+Enter)", self.run_requested.emit)
        menu.exec(event.globalPos())

    def event(self, event) -> bool:
        # Plain keys are typing here, not the window's single-key shortcuts.
        if event.type() == QEvent.Type.ShortcutOverride:
            plain = event.modifiers() & ~(
                Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.KeypadModifier
            )
            if plain == Qt.KeyboardModifier.NoModifier:
                event.accept()
                return True
        return super().event(event)

    def _mark_pending(self) -> None:
        """A tint behind the lines not run yet."""
        if len(self.toPlainText()) <= self.committed_end:
            self.setExtraSelections([])
            return
        tint = QTextCharFormat()
        tint.setBackground(QColor(PENDING_LIGHT))
        if self.palette().base().color().lightness() < 128:
            tint.setBackground(QColor(PENDING_DARK))
        tint.setProperty(QTextFormat.Property.FullWidthSelection, True)
        lines = []
        block = self.document().findBlock(self.committed_end)
        while block.isValid():
            line = QTextEdit.ExtraSelection()
            line.format = tint
            line.cursor = QTextCursor(block)
            lines.append(line)
            block = block.next()
        self.setExtraSelections(lines)


class ScriptPanel(QWidget):
    """The trace and the code to run; ``save_requested`` asks the window for a path and
    ``runner(code) -> (code not run, message)`` runs code in the window."""

    save_requested = Signal()

    def __init__(self, trace, parent=None):
        super().__init__(parent)
        self.trace = trace
        self.runner = None
        self.editor = ScriptEditor(self)
        self.editor.run_requested.connect(self.run)
        self._highlighter = _python_highlighter(self.editor)
        self.path_label = QLabel(self)
        self.path_label.setWordWrap(True)
        self.path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.message = QLabel(self)
        self.message.setWordWrap(True)
        self.run_button = QPushButton("Run")
        self.run_button.setToolTip(
            "Run the lines written after the trace in the window (Ctrl+Enter); they then join "
            "the trace"
        )
        self.run_button.clicked.connect(self.run)
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
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.copy_button)
        buttons.addWidget(self.save_button)
        buttons.addWidget(self.open_button)
        buttons.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.addWidget(self.path_label)
        layout.addWidget(self.editor, 1)
        layout.addWidget(self.message)
        layout.addLayout(buttons)
        self.set_trace(trace)

    @property
    def view(self) -> ScriptEditor:
        return self.editor

    def set_trace(self, trace) -> None:
        self.trace = trace
        trace.on_append = self.append
        self.path_label.setText(
            "Everything done in this window, as Python (<code>trace_path</code> in the "
            "console). Write new lines at the end and press Ctrl+Enter to run them in the "
            f"window.<br><code>{trace.path}</code>"
        )
        self.editor.set_content(trace.text())
        self._scroll_to_end()

    def text(self) -> str:
        return self.editor.toPlainText()

    def append(self, text: str) -> None:
        self.editor.append_committed(text)
        self._scroll_to_end()

    def run(self) -> bool:
        """Run the lines not run yet; what fails stays, with the error below."""
        code = self.editor.pending()
        if not code.strip() or self.runner is None:
            return False
        rest, message = self.runner(code)
        self.editor.replace_pending(rest)
        self.message.setText(message)
        self.message.setStyleSheet("")
        if rest.strip():
            self.message.setStyleSheet("color: #ff6b6b;")
        self._scroll_to_end()
        return not rest.strip()

    def copy(self) -> None:
        QGuiApplication.clipboard().setText(self.trace.text())

    def open_file(self) -> bool:
        return QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(self.trace.path))))

    def _scroll_to_end(self) -> None:
        bar = self.editor.verticalScrollBar()
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

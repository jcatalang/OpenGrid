"""Embedded consoles executing in the application's process and namespace.

``make_console`` gives an IPython console (Jupyter's Qt console on an in-process kernel:
completion, ``?`` help, magics, history) when ``qtconsole`` and ``ipykernel`` are installed,
else a plain Python console (``code.InteractiveConsole``). Both run the code on the GUI thread
and share the interface the window and the tests use: ``namespace``, ``run_line``,
``input_text`` / ``clear_input`` and ``focus_input``.
"""

import code
import contextlib
import io
import os
import re

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QFontDatabase, QKeyEvent
from PySide6.QtWidgets import (
    QApplication,
    QLineEdit,
    QPlainTextEdit,
    QVBoxLayout,
    QWidget,
)

PROMPT = ">>> "
CONTINUATION_PROMPT = "... "
_ANSI = re.compile(r"\x1b\[[0-9;]*m")
# One in-process kernel for the whole process (IPython's shell is a singleton anyway): the
# consoles of all windows are clients of it. It lives until the process ends.
_SHARED: dict = {}


def ipython_available() -> bool:
    """True if the IPython console can be used (qtconsole and ipykernel importable)."""
    try:
        import ipykernel.inprocess  # noqa: F401
        import qtconsole.inprocess  # noqa: F401
    except ImportError:
        return False
    return True


def make_console(namespace: dict, banner: str = "", parent=None) -> QWidget:
    """The IPython console if available, else the plain Python one."""
    if os.environ.get("OPENGRID_PLAIN_CONSOLE") is None and ipython_available():
        return IPythonConsole(namespace, banner, parent)
    return PythonConsole(namespace, banner, parent)


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
    """Plain interactive Python console (the fallback without IPython).

    Commands run synchronously on the GUI thread; long operations must be started through
    the worker-thread API, not from here.
    """

    kind = "python"

    def __init__(self, namespace: dict | None = None, banner: str = "", parent=None):
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

        if not banner:
            banner = "OpenGrid Python console. Namespace: " + ", ".join(sorted(self.namespace))
        self.write(banner)

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

    def input_text(self) -> str:
        return self.input.text()

    def clear_input(self) -> None:
        self.input.clear()

    def focus_input(self) -> QWidget:
        self.input.setFocus()
        return self.input

    def shutdown(self) -> None:
        """Nothing to stop (the IPython console stops its kernel)."""

    def _on_return(self) -> None:
        line = self.input.text()
        self.input.remember(line)
        self.input.clear()
        self.run_line(line)


class IPythonConsole(QWidget):
    """IPython in the application's process: Jupyter's Qt console on an in-process kernel.

    The kernel starts when a console is first shown or used (a window that is never shown,
    as in the tests, starts none); it is one per process, as IPython's shell is: each console
    binds its own namespace (``og``, ``window`` …) again when it gets the focus or runs code,
    so with several windows the names point at the window being typed in. Plain keys typed in
    it stay in it (the window's single-key shortcuts do not take them).
    """

    kind = "ipython"

    def __init__(self, namespace: dict | None = None, banner: str = "", parent=None):
        super().__init__(parent)
        if namespace is None:
            namespace = {}
        self.namespace = namespace
        self.banner = banner
        self.widget = None  # qtconsole's RichJupyterWidget, made at the first use
        self._manager = None
        self._client = None
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt override)
        self.start()
        super().showEvent(event)

    def start(self) -> None:
        """Start the kernel and the console widget (once)."""
        if self.widget is not None:
            return
        from qtconsole.rich_jupyter_widget import RichJupyterWidget

        manager = _kernel_manager()
        client = manager.client()
        client.start_channels()
        widget = RichJupyterWidget(self)
        widget.banner = self.banner.rstrip("\n") + "\n\n"
        widget.kernel_banner = ""
        if _palette_is_dark(self):
            widget.set_default_style(colors="linux")
        widget.kernel_manager = manager
        widget.kernel_client = client
        widget._control.installEventFilter(self)
        self._manager = manager
        self._client = client
        self.widget = widget
        self._layout.addWidget(widget)
        self.bind()

    def bind(self) -> None:
        """Point the shared IPython shell at this console's namespace."""
        if self._manager is None:
            return
        shell = self._manager.kernel.shell
        stale = [k for k, v in self.namespace.items() if shell.user_ns.get(k) is not v]
        if stale:
            shell.push(dict(self.namespace))

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt override)
        if event.type() == QEvent.Type.FocusIn:
            self.bind()
        if event.type() == QEvent.Type.ShortcutOverride:
            # Keys without Ctrl / Alt / Meta are typing (a line edit claims them the same way);
            # qtconsole handles its own Ctrl shortcuts.
            plain = event.modifiers() & ~(
                Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.KeypadModifier
            )
            if plain == Qt.KeyboardModifier.NoModifier:
                event.accept()
                return True
        return super().eventFilter(watched, event)

    def _wait_for_prompt(self) -> None:
        """The first prompt comes with the kernel's first reply (an event away)."""
        for _ in range(100):
            if self.widget._prompt:
                return
            QApplication.processEvents()

    def run_line(self, source: str) -> str:
        """Execute source as if typed (shown in the console); return its output as text
        (printed text, the value of an expression, a traceback without colours)."""
        self.start()
        self._wait_for_prompt()
        self.bind()
        captured: list[str] = []

        def on_message(message) -> None:
            kind = message["header"]["msg_type"]
            content = message["content"]
            if kind == "stream":
                captured.append(content["text"])
            elif kind == "execute_result":
                captured.append(content["data"].get("text/plain", "") + "\n")
            elif kind == "error":
                captured.append(_ANSI.sub("", "\n".join(content["traceback"])) + "\n")

        channel = self._client.iopub_channel
        channel.message_received.connect(on_message)
        try:
            self.widget.execute(source)
            QApplication.processEvents()
        finally:
            channel.message_received.disconnect(on_message)
        return "".join(captured)

    def input_text(self) -> str:
        if self.widget is None:
            return ""
        return self.widget.input_buffer

    def clear_input(self) -> None:
        if self.widget is not None:
            self.widget.input_buffer = ""

    def focus_input(self) -> QWidget:
        self.start()
        self._wait_for_prompt()
        self.widget._control.setFocus()
        self.bind()
        return self.widget._control

    def shutdown(self) -> None:
        """Detach from the shared kernel (the window is closing); the kernel stays."""
        if self._client is None:
            return
        self._client.stop_channels()
        frontends = self._manager.kernel.frontends
        if self._client in frontends:
            frontends.remove(self._client)
        self._client = None


def _kernel_manager():
    """The process's in-process kernel (started at the first call)."""
    manager = _SHARED.get("manager")
    if manager is None:
        # The kernel's debugger warns about frozen modules on every start otherwise.
        os.environ.setdefault("PYDEVD_DISABLE_FILE_VALIDATION", "1")
        from qtconsole.inprocess import QtInProcessKernelManager

        manager = QtInProcessKernelManager()
        manager.start_kernel(show_banner=False)
        manager.kernel.gui = "qt"
        manager.kernel.shell.banner1 = ""
        manager.kernel.shell.banner2 = ""
        manager.kernel.shell.enable_tip = False
        _SHARED["manager"] = manager
    return manager


def _palette_is_dark(widget: QWidget) -> bool:
    return widget.palette().window().color().lightness() < 128

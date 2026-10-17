"""Long operations off the GUI thread.

``run_blocking`` runs a function in a worker thread while the GUI thread keeps processing
events in a local loop; after a short grace period a window-modal busy dialog says what is
running (with a Cancel button when the operation can be cancelled). The C++ core and Netgen
release the GIL, so the window stays responsive (it repaints, the dialog animates) but takes
no input into the model meanwhile. The result, or the exception, comes back on the GUI thread.
The caller must not touch Qt from the function; the window suspends its model listener around
it (MainWindow._run_long), so commands run by the function (a CAD import) only record what
changed and the window refreshes once afterwards.
"""

import threading

from PySide6.QtCore import QEventLoop, QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import QProgressDialog

GRACE_MS = 250  # the dialog appears only for operations longer than this


class _Emitter(QObject):
    done = Signal()


def run_blocking(parent, label: str, function, cancel=None, grace_ms: int = GRACE_MS):
    """Run ``function()`` in a thread and return its result (or re-raise its exception),
    processing GUI events meanwhile; ``cancel()`` (optional) is called when the user presses
    the dialog's Cancel button."""
    box: dict = {}
    emitter = _Emitter()
    loop = QEventLoop(parent)
    emitter.done.connect(loop.quit, Qt.ConnectionType.QueuedConnection)

    def target():
        try:
            box["result"] = function()
        except BaseException as error:  # handed back to the GUI thread
            box["error"] = error
        finally:
            emitter.done.emit()

    thread = threading.Thread(target=target, name=f"opengrid: {label}", daemon=True)
    dialog: QProgressDialog | None = None

    def show_dialog():
        nonlocal dialog
        if not thread.is_alive():
            return
        button = ""
        if cancel is not None:
            button = "Cancel"
        dialog = QProgressDialog(label, button, 0, 0, parent)
        dialog.setWindowTitle("OpenGrid")
        dialog.setWindowModality(Qt.WindowModality.WindowModal)
        dialog.setMinimumDuration(0)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)
        if cancel is None:
            dialog.setCancelButton(None)
        else:
            dialog.canceled.connect(cancel)
        dialog.show()

    thread.start()
    QTimer.singleShot(grace_ms, show_dialog)
    loop.exec()
    thread.join()
    if dialog is not None:
        dialog.close()
        dialog.deleteLater()
    if "error" in box:
        raise box["error"]
    return box.get("result")

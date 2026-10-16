"""Entry point: ``python -m opengrid [project.ogp ...]`` or the ``opengrid`` console script.

Projects named on the command line are opened; on macOS the Finder hands documents over as
file-open events (double-clicked .ogp files, also while the application runs).
"""

import logging
import sys
from pathlib import Path


def _file_open_filter(window):
    from PySide6.QtCore import QEvent, QObject

    class FileOpenFilter(QObject):
        """Opens the documents of the application's file-open events in the window."""

        def eventFilter(self, watched, event):  # noqa: N802 - Qt's name
            if event.type() == QEvent.Type.FileOpen:
                path = event.file()
                if not path:
                    path = event.url().toLocalFile()
                if path:
                    window.open_document(path)
                return True
            return super().eventFilter(watched, event)

    return FileOpenFilter(window)


def main(argv=None) -> int:
    from PySide6.QtWidgets import QApplication

    from opengrid.app.main_window import MainWindow

    if argv is None:
        argv = sys.argv
    handler = logging.StreamHandler(sys.stderr)  # the launcher app keeps it in a log file
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    logging.getLogger("opengrid").addHandler(handler)
    app = QApplication.instance()
    if app is None:
        app = QApplication(argv)
    app.setApplicationName("OpenGrid")
    window = MainWindow()
    opener = _file_open_filter(window)
    app.installEventFilter(opener)
    window.show()
    for arg in argv[1:]:
        if arg.startswith("-psn"):  # older macOS: a process serial number from the Finder
            continue
        if Path(arg).suffix.lower() == ".ogp":
            window.open_document(arg)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())

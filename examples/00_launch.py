"""Phase 0 demo: launch the OpenGrid main window.

python examples/00_launch.py                  # interactive
python examples/00_launch.py --screenshot out.png   # save a screenshot and quit
"""

import argparse
import sys
from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication

from opengrid.app.main_window import MainWindow


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--screenshot", type=Path, default=None, help="save a PNG of the main window and exit"
    )
    args = parser.parse_args()

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    window.console.run_line("print('core version', opengrid.core_version())")

    if args.screenshot is not None:

        def save_and_quit():
            window.viewport.render()
            args.screenshot.parent.mkdir(parents=True, exist_ok=True)
            window.grab().save(str(args.screenshot))
            print(f"saved {args.screenshot}")
            window.close()
            app.quit()

        QTimer.singleShot(1500, save_and_quit)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())

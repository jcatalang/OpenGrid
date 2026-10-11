"""Elliptic smoother dock: residual history plot, run and cancel.

The smoother runs in a Python thread; the C++ solver releases the GIL, so the GUI stays
responsive. The inputs (initial grid and settings) are taken on the GUI thread and the
result is handed back through a queued Qt signal and stored in the session's cache on the
GUI thread, so the model is never touched from the worker.
"""

import threading

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QObject, QTimer, Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from opengrid import _core
from opengrid.model.domains import run_smoother


class _Emitter(QObject):
    done = Signal(object)


class SmootherDock(QWidget):
    finished = Signal(str, object)  # domain name, SmoothingResult (or exception)

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.domain_name: str | None = None
        self._progress = None
        self._thread: threading.Thread | None = None
        self._key = None
        self._emitter = _Emitter()
        self._emitter.done.connect(self._on_done)

        layout = QVBoxLayout(self)
        self.figure = Figure(figsize=(4, 2.2), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(240, 160)
        self.axes = self.figure.add_subplot(111)
        layout.addWidget(self.canvas, 1)
        row = QHBoxLayout()
        self.status = QLabel("Select a domain with elliptic smoothing")
        self.status.setWordWrap(True)
        self.run_button = QPushButton("Run smoother")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.run_button.clicked.connect(self.run_selected)
        self.cancel_button.clicked.connect(self.cancel)
        row.addWidget(self.status, 1)
        row.addWidget(self.run_button)
        row.addWidget(self.cancel_button)
        layout.addLayout(row)
        self._timer = QTimer(self)
        self._timer.setInterval(100)
        self._timer.timeout.connect(self._update_plot)

    @property
    def running(self) -> bool:
        return self._thread is not None

    def set_domain(self, name: str | None) -> None:
        if self.running:
            return
        self.domain_name = name
        self._plot([])
        if name is None:
            self.status.setText("Select a domain with elliptic smoothing")
            self.run_button.setEnabled(False)
            return
        entity = self.session.entity(name)
        if entity.kind != "domain":
            self.domain_name = None
            self.status.setText(f"{name}: elliptic smoothing applies to structured domains only")
            self.run_button.setEnabled(False)
            return
        generator = entity.params["generator"]
        self.run_button.setEnabled(generator["elliptic"] is not None)
        self.status.setText(f"{name}: {self.session.domain_status(name)}")
        cached = self.session.domains.smoothing_result(self.session.entity(name).id)
        if cached is not None:
            self._plot(cached.residuals)

    def run_selected(self) -> bool:
        if self.domain_name is None or self.running:
            return False
        return self.start(self.domain_name)

    def start(self, name: str) -> bool:
        entity = self.session.entity(name)
        elliptic = entity.params["generator"]["elliptic"]
        if elliptic is None:
            self.status.setText(f"{name} has no elliptic smoothing configured")
            return False
        self.domain_name = name
        self._key = self.session.domains.key(entity.id)
        initial = self.session.domains.initial_grid(entity.id)
        self._progress = _core.SmoothingProgress()
        progress = self._progress

        def work():
            try:
                result = run_smoother(initial, elliptic, progress)
            except Exception as error:  # reported on the GUI thread
                result = error
            self._emitter.done.emit(result)

        self._thread = threading.Thread(target=work, name="opengrid-smoother", daemon=True)
        self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.status.setText(f"{name}: smoothing…")
        self._timer.start()
        self._thread.start()
        return True

    def cancel(self) -> None:
        if self._progress is not None:
            self._progress.cancel()

    def wait(self, timeout: float = 60.0) -> None:
        """Block until the worker has finished (tests and scripts)."""
        if self._thread is not None:
            self._thread.join(timeout)

    def _on_done(self, result) -> None:
        self._timer.stop()
        thread = self._thread
        self._thread = None
        if thread is not None:
            thread.join()
        self.cancel_button.setEnabled(False)
        self.run_button.setEnabled(True)
        name = self.domain_name
        if isinstance(result, Exception):
            self.status.setText(f"{name}: smoothing failed: {result}")
            self.finished.emit(name, result)
            return
        self._plot(result.residuals)
        if result.cancelled:
            self.status.setText(f"{name}: cancelled after {result.iterations} iterations")
        else:
            self.session.domains.store_result(self._key, result)
            self.status.setText(f"{name}: {self.session.domain_status(name)}")
        self.finished.emit(name, result)

    def _update_plot(self) -> None:
        if self._progress is not None:
            self._plot(self._progress.residuals())

    def _plot(self, residuals) -> None:
        self.axes.clear()
        residuals = np.asarray(residuals)
        if residuals.size and residuals[0] > 0:
            self.axes.semilogy(residuals / residuals[0], color="#e4572e")
        self.axes.set_xlabel("iteration")
        self.axes.set_ylabel("residual / initial")
        self.canvas.draw_idle()

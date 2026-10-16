"""Elliptic smoother dock: residual history plot, run and cancel, for a structured domain or
a smoothing group (a domain in a group runs its group).

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
from opengrid.model.groups import run_group


class _Emitter(QObject):
    done = Signal(object)


class SmootherDock(QWidget):
    finished = Signal(str, object)  # domain or group name, Smoothing/GroupResult (or exception)

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.domain_name: str | None = None
        self._progress = None
        self._thread: threading.Thread | None = None
        self._key = None
        self._group = None  # (params, layout) while a group runs
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
        group = self._group_of(name)
        if group is not None:
            self.domain_name = group
            self.run_button.setEnabled(True)
            self.status.setText(f"{group}: {self.session.group_status(group)}")
            cached = self.session.groups.result(self.session.entity(group).id)
            if cached is not None:
                self._plot(cached.residuals)
            return
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

    def _group_of(self, name: str) -> str | None:
        """The smoothing group a name stands for: the group itself, or a member domain's."""
        entity = self.session.entity(name)
        if entity.kind == "smoothing_group":
            return name
        if entity.kind not in ("domain", "block"):
            return None
        try:
            if entity.kind == "block":
                group = self.session.groups.group_of_block(entity.id)
            else:
                group = self.session.groups.group_of(entity.id)
        except Exception:  # a domain in several groups: reported when it is evaluated
            return None
        if group is None:
            return None
        return self.session.store.get(group).name

    def run_selected(self) -> bool:
        if self.domain_name is None or self.running:
            return False
        return self.start(self.domain_name)

    def start(self, name: str) -> bool:
        group = self._group_of(name)
        if group is not None:
            return self._start_group(group)
        entity = self.session.entity(name)
        self._group = None
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

    def _start_group(self, name: str) -> bool:
        try:
            key, layout, params = self.session.groups.prepare(self.session.entity(name).id)
        except Exception as error:
            self.status.setText(f"{name}: {error}")
            return False
        self.domain_name = name
        self._key = key
        self._group = (params, layout)
        self._progress = _core.SmoothingProgress()
        progress = self._progress

        def work():
            try:
                result = run_group(layout, params, progress)
            except Exception as error:  # reported on the GUI thread
                result = error
            self._emitter.done.emit(result)

        self._thread = threading.Thread(target=work, name="opengrid-group", daemon=True)
        self.run_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        count = len(params["domains"]) + len(params["blocks"])
        self.status.setText(f"{name}: smoothing {count} members…")
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
        group = self._group
        self._group = None
        if group is not None:
            params, layout = group
            result = self.session.groups.finish(params, layout, result)
        self._plot(result.residuals)
        if result.cancelled:
            self.status.setText(f"{name}: cancelled after {result.iterations} iterations")
        elif group is not None:
            self.session.groups.store_result(self._key, result)
            self.status.setText(f"{name}: {self.session.group_status(name)}")
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

"""Quality dock: metric selection, histogram, statistics and worst-cell locator."""

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

# Metric -> (label, True if larger values are worse)
METRICS = {
    "scaled_jacobian": ("Scaled Jacobian", False),
    "skewness": ("Equiangle skewness", True),
    "aspect_ratio": ("Aspect ratio", True),
    "min_angle": ("Min included angle [deg]", False),
    "max_angle": ("Max included angle [deg]", True),
    "non_orthogonality": ("Non-orthogonality [deg]", True),
    "area_ratio": ("Neighbour area ratio", True),
}


def worst_cell(values: np.ndarray, larger_is_worse: bool) -> tuple[int, int]:
    """(j, i) of the worst cell (first one on ties)."""
    if larger_is_worse:
        flat = int(np.argmax(values))
    else:
        flat = int(np.argmin(values))
    j, i = np.unravel_index(flat, values.shape)
    return int(j), int(i)


class QualityDock(QWidget):
    metric_changed = Signal(str)
    locate_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        self.metric = QComboBox()
        for key, (label, _worse) in METRICS.items():
            self.metric.addItem(label, key)
        self.metric.currentIndexChanged.connect(
            lambda _i: self.metric_changed.emit(self.metric_key())
        )
        self.locate = QPushButton("Show worst cell")
        self.locate.clicked.connect(self.locate_requested.emit)
        row.addWidget(self.metric, 1)
        row.addWidget(self.locate)
        layout.addLayout(row)
        self.figure = Figure(figsize=(4, 2.2), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumSize(240, 160)
        self.axes = self.figure.add_subplot(111)
        layout.addWidget(self.canvas, 1)
        self.stats = QLabel("No domains")
        self.stats.setWordWrap(True)
        layout.addWidget(self.stats)

    def metric_key(self) -> str:
        return self.metric.currentData()

    def show_values(self, values: np.ndarray | None, scope: str) -> None:
        self.axes.clear()
        if values is None or values.size == 0:
            self.stats.setText("No domains")
            self.canvas.draw_idle()
            return
        label, larger_is_worse = METRICS[self.metric_key()]
        finite = values[np.isfinite(values)]
        low = float(np.min(finite))
        high = float(np.max(finite))
        if high - low <= 1e-12 * max(1.0, abs(high)):
            # All cells equal (e.g. a uniform grid): give the single bar a visible width.
            pad = 0.05 * max(1.0, abs(high))
            low -= pad
            high += pad
        self.axes.hist(finite, bins=40, range=(low, high), color="#4c8bf5")
        self.axes.set_xlabel(label)
        self.axes.set_ylabel("cells")
        self.canvas.draw_idle()
        worst = "min"
        worst_value = float(np.min(finite))
        if larger_is_worse:
            worst = "max"
            worst_value = float(np.max(finite))
        text = (
            f"{scope}: {values.size} cells — min {np.min(finite):.4g}, mean {np.mean(finite):.4g}, "
            f"max {np.max(finite):.4g} (worst = {worst} = {worst_value:.4g})"
        )
        self.stats.setText(text)

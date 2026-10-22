"""Quality dock: metric selection, the mesh coloured by it (or only its cells worse than a
threshold: in 3D also those inside blocks and volumes), histogram, statistics and worst-cell
locator."""

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

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
# Metric -> the threshold offered for "only cells worse than" (cells past it are suspect;
# inflation layers are meant to have large aspect ratios).
THRESHOLDS = {
    "scaled_jacobian": 0.3,
    "skewness": 0.8,
    "aspect_ratio": 1000.0,
    "min_angle": 20.0,
    "max_angle": 150.0,
    "non_orthogonality": 60.0,
    "area_ratio": 2.5,
}


def worse_than(values: np.ndarray, key: str, threshold: float) -> np.ndarray:
    """Mask of the cells worse than the threshold for the metric (beyond it, its bad side)."""
    if METRICS[key][1]:
        return values > threshold
    return values < threshold


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
    colouring_toggled = Signal(bool)  # colour the mesh by the metric
    threshold_changed = Signal(object)  # only cells worse than this (a float), or None

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        self.metric = QComboBox()
        for key, (label, _worse) in METRICS.items():
            self.metric.addItem(label, key)
        self.metric.currentIndexChanged.connect(lambda _i: self._on_metric())
        self.locate = QPushButton("Show worst cell")
        self.locate.clicked.connect(self.locate_requested.emit)
        row.addWidget(self.metric, 1)
        row.addWidget(self.locate)
        layout.addLayout(row)
        self.colour = QCheckBox("Colour the mesh")
        self.colour.setToolTip("Colour every cell of the domains, blocks and volumes by the metric")
        self.colour.toggled.connect(self.colouring_toggled.emit)
        self.only_worse = QCheckBox("Only cells worse than")
        self.only_worse.setToolTip(
            "Fill only the cells past this value (in 3D also those inside blocks and volumes), "
            "the rest as lines: where the distorted cells are"
        )
        self.threshold = QDoubleSpinBox()
        self.threshold.setDecimals(3)
        self.threshold.setRange(-1e9, 1e9)
        self.threshold.setKeyboardTracking(False)
        self.threshold.setValue(THRESHOLDS[self.metric_key()])
        self.only_worse.toggled.connect(lambda _on: self._emit_threshold())
        self.threshold.valueChanged.connect(lambda _v: self._emit_threshold())
        show = QHBoxLayout()
        show.addWidget(self.colour)
        show.addStretch(1)
        show.addWidget(self.only_worse)
        show.addWidget(self.threshold)
        layout.addLayout(show)
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

    def threshold_value(self) -> float | None:
        """The "only cells worse than" value, or None when off."""
        if not self.only_worse.isChecked():
            return None
        return float(self.threshold.value())

    def _on_metric(self) -> None:
        self.reset_threshold()
        self.metric_changed.emit(self.metric_key())

    def _emit_threshold(self) -> None:
        self.threshold_changed.emit(self.threshold_value())

    def reset_threshold(self) -> None:
        """The threshold offered for the current metric (after the metric changed)."""
        self.threshold.blockSignals(True)
        self.threshold.setValue(THRESHOLDS[self.metric_key()])
        self.threshold.blockSignals(False)

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
        threshold = self.threshold_value()
        if threshold is not None:
            bad = int(np.count_nonzero(worse_than(finite, self.metric_key(), threshold)))
            text += f"; {bad} worse than {threshold:g}"
            self.axes.axvline(threshold, color="#ff6b6b", linewidth=1.5)
            self.canvas.draw_idle()
        self.stats.setText(text)

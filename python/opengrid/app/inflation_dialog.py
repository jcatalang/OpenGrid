"""The "Inflation layers" dialog: boundary layers grown from the selected curves into the domains
they bound (Session.set_inflation)."""

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
)

from opengrid.app.forms import form_layout


class InflationDialog(QDialog):
    """Settings of the layers (prefilled from ``current``, an inflation dict, or defaults);
    ``values()`` after Apply, ``removed`` after "Remove layers"."""

    def __init__(self, curves: list[str], current=None, suggested_first=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Inflation layers")
        self.removed = False
        if current is None:
            current = {
                "first_height": None,
                "growth_rate": 1.15,
                "layers": 30,
                "total_thickness": None,
                "transition": True,
            }
        shown = ", ".join(curves[:6])
        if len(curves) > 6:
            shown += f" and {len(curves) - 6} more"
        intro = QLabel(f"Layers grow from {shown} into the domains they bound.")
        intro.setWordWrap(True)
        self.first_height = QLineEdit(_text(current.get("first_height")))
        hint = "required"
        if suggested_first is not None:
            hint = f"{suggested_first:.3g} (a twentieth of the curve spacing)"
        self.first_height.setPlaceholderText(hint)
        self.suggested_first = suggested_first
        self.growth = QLineEdit(_text(current.get("growth_rate", 1.15)))
        self.layers = QSpinBox()
        self.layers.setRange(1, 1000)
        self.layers.setValue(int(current.get("layers", 30)))
        self.thickness = QLineEdit(_text(current.get("total_thickness")))
        self.thickness.setPlaceholderText("none")
        self.transition = QCheckBox("stop at the core cell size")
        self.transition.setChecked(bool(current.get("transition", True)))
        self.slide = QCheckBox("run into the neighbouring curves (no steps)")
        self.slide.setChecked(True)
        self.slide.setToolTip(
            "The layers run along the boundary curves next to these up to them, instead of "
            "tapering into steps (those curves get points at the layer heights near the walls)"
        )
        form = form_layout()
        form.addRow("first height", self.first_height)
        form.addRow("growth rate", self.growth)
        form.addRow("layers", self.layers)
        form.addRow("total thickness", self.thickness)
        form.addRow(self.transition)
        form.addRow(self.slide)
        self.error = QLabel("")
        self.error.setStyleSheet("color: #ff6b6b;")
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Apply | QDialogButtonBox.StandardButton.Cancel
        )
        remove = QPushButton("Remove layers")
        remove.setToolTip("These curves stop growing layers")
        buttons.addButton(remove, QDialogButtonBox.ButtonRole.DestructiveRole)
        buttons.button(QDialogButtonBox.StandardButton.Apply).clicked.connect(self._apply)
        buttons.rejected.connect(self.reject)
        remove.clicked.connect(self._remove)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(form)
        layout.addWidget(self.error)
        layout.addWidget(buttons)

    def values(self) -> dict:
        """The settings, for Session.set_inflation; raises ValueError for bad input."""
        first = _number(self.first_height.text(), "first height", optional=True)
        if first is None:
            first = self.suggested_first
        if first is None:
            raise ValueError("first height: required")
        values = {
            "first_height": first,
            "growth_rate": _number(self.growth.text(), "growth rate"),
            "layers": self.layers.value(),
            "transition": self.transition.isChecked(),
            "slide": self.slide.isChecked(),
        }
        thickness = _number(self.thickness.text(), "total thickness", optional=True)
        if thickness is not None:
            values["total_thickness"] = thickness
        return values

    def _apply(self) -> None:
        try:
            self.values()
        except ValueError as error:
            self.error.setText(str(error))
            return
        self.accept()

    def _remove(self) -> None:
        self.removed = True
        self.accept()


def _text(value) -> str:
    if value is None:
        return ""
    return f"{value:g}"


def _number(text: str, label: str, optional: bool = False):
    text = text.strip()
    if not text:
        if optional:
            return None
        raise ValueError(f"{label}: required")
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{label}: not a number ({text!r})") from None

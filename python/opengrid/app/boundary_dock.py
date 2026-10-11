"""Boundary dock: patch/zone creation from the selection and the boundary validator report."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from opengrid.model.patches import PATCH_TYPES

# Categorical colours for patches (cycled); unassigned boundary connectors and interfaces have
# their own colours.
PATCH_PALETTE = (
    "#4c8bf5",
    "#2ee59d",
    "#ffb86b",
    "#c678dd",
    "#56b6c2",
    "#e5c07b",
    "#98c379",
    "#f78c6c",
)
UNASSIGNED_COLOR = "#ff3b30"
INTERFACE_COLOR = "#7f8c8d"


class BoundaryDock(QWidget):
    patch_requested = Signal(str, str)  # name, patch type
    add_requested = Signal(str)  # patch name
    remove_requested = Signal(str)
    zone_requested = Signal(str)
    validate_requested = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.name = QLineEdit()
        self.name.setPlaceholderText("patch or zone name")
        self.patch_type = QComboBox()
        self.patch_type.addItems(list(PATCH_TYPES))
        form.addRow("name", self.name)
        form.addRow("patch type", self.patch_type)
        layout.addLayout(form)
        row = QHBoxLayout()
        self.new_patch = QPushButton("New patch")
        self.add = QPushButton("Add to patch")
        self.remove = QPushButton("Remove from patch")
        self.new_zone = QPushButton("New zone")
        for button in (self.new_patch, self.add, self.remove, self.new_zone):
            row.addWidget(button)
        layout.addLayout(row)
        self.new_patch.clicked.connect(
            lambda: self.patch_requested.emit(
                self.name.text().strip(), self.patch_type.currentText()
            )
        )
        self.add.clicked.connect(lambda: self.add_requested.emit(self.name.text().strip()))
        self.remove.clicked.connect(lambda: self.remove_requested.emit(self.name.text().strip()))
        self.new_zone.clicked.connect(lambda: self.zone_requested.emit(self.name.text().strip()))
        validate_row = QHBoxLayout()
        self.validate = QPushButton("Validate boundary")
        self.validate.clicked.connect(self.validate_requested.emit)
        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        validate_row.addWidget(self.validate)
        validate_row.addWidget(self.summary, 1)
        layout.addLayout(validate_row)
        self.issues = QListWidget()
        layout.addWidget(self.issues, 1)

    def show_report(self, report) -> None:
        self.issues.clear()
        state = "OK - ready to export"
        if not report.ok:
            state = f"{len(report.issues)} problem(s)"
        self.summary.setText(
            f"{len(report.boundary)} boundary connectors, {len(report.interfaces)} interfaces, "
            f"{len(report.unassigned)} unassigned: {state}"
        )
        for issue in report.issues:
            self.issues.addItem(str(issue))
        for interface in report.interfaces:
            self.issues.addItem(f"interface: {interface}")

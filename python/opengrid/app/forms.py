"""Form layouts for the panels and dialogs."""

from PySide6.QtWidgets import QFormLayout


def form_layout(parent=None) -> QFormLayout:
    """A QFormLayout whose fields take the width available (AllNonFixedFieldsGrow, as Qt's
    Fusion style does on Linux). macOS's default keeps fields at their preferred width, so a
    word-wrapped label was given the height of fewer, longer lines than it gets, and its text
    was cut."""
    if parent is None:
        form = QFormLayout()
    else:
        form = QFormLayout(parent)
    form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
    return form

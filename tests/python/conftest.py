import os

# Allow GUI tests to run without a display (Linux CI under xvfb still works with this unset).
if os.environ.get("OG_GUI_ONSCREEN") is None and os.environ.get("DISPLAY") is None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

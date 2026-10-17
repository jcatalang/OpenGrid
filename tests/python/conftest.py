import os
import sys

# Allow GUI tests to run without a display (Linux CI under xvfb still works with this unset).
# On macOS Qt never uses DISPLAY (XQuartz exports one in every shell), so the tests run
# offscreen there unless asked otherwise: on-screen windows get resized by the window manager
# and the real-event tests' display coordinates no longer match.
_NO_DISPLAY = sys.platform == "darwin" or os.environ.get("DISPLAY") is None
if os.environ.get("OG_GUI_ONSCREEN") is None and _NO_DISPLAY:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

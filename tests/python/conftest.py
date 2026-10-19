import os
import sys

# GUI tests run on Qt's offscreen platform on macOS (Qt never uses DISPLAY there, XQuartz exports
# one in every shell; on-screen windows get resized by the window manager and the real-event
# tests' display coordinates no longer match) and where there is no display. On Linux with a
# display (xvfb-run) they use X11: the offscreen platform there can abort in pyvistaqt when a
# window closes (QOpenGLWidget is not supported offscreen; seen on openSUSE, 2026-10-19). The
# real-event harness activates its window, which a bare X server allows (tests/python/
# test_gui_events.py). OG_GUI_ONSCREEN keeps the platform as it is.
_NO_DISPLAY = sys.platform == "darwin" or os.environ.get("DISPLAY") is None
if os.environ.get("OG_GUI_ONSCREEN") is None and _NO_DISPLAY:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

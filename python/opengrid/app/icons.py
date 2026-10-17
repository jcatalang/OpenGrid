"""Toolbar icons drawn with QPainter: simple line glyphs in the palette's text colour, so they
need no image files and follow light and dark themes. ``icon(name)`` caches one QIcon per
name and colour, rendered at twice the nominal size for HiDPI screens."""

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPen, QPixmap

SIZE = 24  # nominal icon size in points
_cache: dict[tuple[str, str], QIcon] = {}


def icon(name: str, color: QColor | None = None) -> QIcon:
    """The icon of a glyph in GLYPHS (an empty icon for unknown names)."""
    if color is None:
        from PySide6.QtGui import QPalette
        from PySide6.QtWidgets import QApplication

        app = QApplication.instance()
        color = QColor("#d0d4d8")
        if app is not None:
            color = app.palette().color(QPalette.ColorRole.ButtonText)
    key = (name, color.name())
    if key not in _cache:
        draw = GLYPHS.get(name)
        result = QIcon()
        if draw is not None:
            for ratio in (1, 2):
                pixmap = QPixmap(SIZE * ratio, SIZE * ratio)
                pixmap.setDevicePixelRatio(ratio)
                pixmap.fill(Qt.GlobalColor.transparent)
                painter = QPainter(pixmap)
                painter.setRenderHint(QPainter.RenderHint.Antialiasing)
                painter.scale(SIZE, SIZE)  # glyphs are drawn in the unit square
                pen = QPen(color)
                pen.setWidthF(0.09)
                pen.setCapStyle(Qt.PenCapStyle.RoundCap)
                pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
                pen.setCosmetic(False)
                painter.setPen(pen)
                painter.setBrush(Qt.BrushStyle.NoBrush)
                draw(painter, color)
                painter.end()
                result.addPixmap(pixmap)
        _cache[key] = result
    return _cache[key]


def _poly(painter, points, close=False):
    path = QPainterPath(QPointF(*points[0]))
    for x, y in points[1:]:
        path.lineTo(x, y)
    if close:
        path.closeSubpath()
    painter.drawPath(path)


def _dot(painter, color, x, y, r=0.07):
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    painter.drawEllipse(QPointF(x, y), r, r)
    painter.restore()


def _thin(painter, factor=0.6):
    pen = painter.pen()
    pen.setWidthF(pen.widthF() * factor)
    painter.setPen(pen)


def _select(painter, color):
    painter.save()
    painter.setBrush(color)
    _poly(
        painter,
        [
            (0.25, 0.15),
            (0.25, 0.8),
            (0.42, 0.63),
            (0.55, 0.9),
            (0.66, 0.84),
            (0.53, 0.58),
            (0.78, 0.58),
        ],
        close=True,
    )
    painter.restore()


def _point(painter, color):
    _dot(painter, color, 0.5, 0.5, 0.14)


def _line(painter, color):
    _poly(painter, [(0.18, 0.82), (0.82, 0.18)])
    _dot(painter, color, 0.18, 0.82)
    _dot(painter, color, 0.82, 0.18)


def _arc(painter, color):
    path = QPainterPath()
    path.arcMoveTo(QRectF(0.15, 0.25, 0.7, 0.9), 160)
    path.arcTo(QRectF(0.15, 0.25, 0.7, 0.9), 160, -140)
    painter.drawPath(path)
    _dot(painter, color, 0.5, 0.25)


def _polyline(painter, color):
    points = [(0.15, 0.75), (0.35, 0.3), (0.55, 0.65), (0.85, 0.25)]
    _poly(painter, points)
    for x, y in points:
        _dot(painter, color, x, y)


def _spline(painter, color):
    path = QPainterPath(QPointF(0.12, 0.7))
    path.cubicTo(0.35, 0.0, 0.55, 1.0, 0.88, 0.3)
    painter.drawPath(path)
    for x, y in ((0.12, 0.7), (0.5, 0.5), (0.88, 0.3)):
        _dot(painter, color, x, y, 0.06)


def _nurbs(painter, color):
    painter.save()
    _thin(painter, 0.5)
    pen = painter.pen()
    pen.setStyle(Qt.PenStyle.DashLine)
    painter.setPen(pen)
    control = [(0.12, 0.8), (0.3, 0.15), (0.7, 0.15), (0.88, 0.8)]
    _poly(painter, control)
    painter.restore()
    path = QPainterPath(QPointF(0.12, 0.8))
    path.cubicTo(0.3, 0.15, 0.7, 0.15, 0.88, 0.8)
    painter.drawPath(path)
    painter.save()
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    for x, y in control:
        painter.drawRect(QRectF(x - 0.06, y - 0.06, 0.12, 0.12))
    painter.restore()


def _freehand(painter, color):
    path = QPainterPath(QPointF(0.12, 0.6))
    path.cubicTo(0.2, 0.1, 0.4, 0.9, 0.5, 0.45)
    path.cubicTo(0.6, 0.1, 0.7, 0.8, 0.88, 0.4)
    painter.drawPath(path)


def _sketch(painter, color):
    _poly(painter, [(0.18, 0.22), (0.82, 0.18), (0.84, 0.8), (0.2, 0.82)], close=True)
    painter.save()
    _thin(painter)
    _poly(painter, [(0.5, 0.2), (0.52, 0.81)])
    _poly(painter, [(0.19, 0.52), (0.83, 0.49)])
    painter.restore()
    for x, y in ((0.18, 0.22), (0.82, 0.18), (0.84, 0.8), (0.2, 0.82)):
        _dot(painter, color, x, y, 0.08)


def _finish(painter, color):
    painter.save()
    pen = painter.pen()
    pen.setWidthF(0.13)
    painter.setPen(pen)
    _poly(painter, [(0.18, 0.52), (0.42, 0.76), (0.84, 0.28)])
    painter.restore()


def _cancel(painter, color):
    painter.save()
    pen = painter.pen()
    pen.setWidthF(0.13)
    painter.setPen(pen)
    _poly(painter, [(0.24, 0.24), (0.76, 0.76)])
    _poly(painter, [(0.76, 0.24), (0.24, 0.76)])
    painter.restore()


def _connectors(painter, color):
    painter.save()
    _thin(painter)
    _poly(painter, [(0.12, 0.6), (0.88, 0.4)])
    painter.restore()
    for k in range(5):
        t = k / 4
        _dot(painter, color, 0.12 + 0.76 * t, 0.6 - 0.2 * t, 0.07)


def _chain(painter, color):
    painter.save()
    _thin(painter)
    _poly(painter, [(0.1, 0.72), (0.45, 0.55), (0.9, 0.3)])
    painter.restore()
    for t in (0.0, 0.25, 0.5):
        _dot(painter, color, 0.1 + 0.35 * t / 0.5, 0.72 - 0.17 * t / 0.5, 0.06)
    for t in (0.25, 0.5, 0.75, 1.0):
        _dot(painter, color, 0.45 + 0.45 * t, 0.55 - 0.25 * t, 0.06)


def _structured(painter, color):
    painter.drawRect(QRectF(0.15, 0.15, 0.7, 0.7))
    painter.save()
    _thin(painter)
    for k in (1, 2):
        v = 0.15 + 0.7 * k / 3
        _poly(painter, [(v, 0.15), (v, 0.85)])
        _poly(painter, [(0.15, v), (0.85, v)])
    painter.restore()


def _unstructured(painter, color):
    painter.save()
    _thin(painter, 0.8)
    points = [(0.15, 0.85), (0.5, 0.15), (0.85, 0.85), (0.5, 0.6), (0.3, 0.4), (0.7, 0.4)]
    _poly(painter, [points[0], points[1], points[2]], close=True)
    _poly(painter, [points[0], points[3]])
    _poly(painter, [points[2], points[3]])
    _poly(painter, [points[1], points[3]])
    _poly(painter, [points[4], points[3]])
    _poly(painter, [points[5], points[3]])
    painter.restore()


def _cube(painter, top=(0.5, 0.15)):
    # An isometric cube in the unit square.
    cx, cy = 0.5, 0.52
    dx, dy = 0.3, 0.17
    h = 0.36
    front_top = [
        (cx, cy - h),
        (cx + dx, cy - h + dy),
        (cx, cy - h + 2 * dy),
        (cx - dx, cy - h + dy),
    ]
    _poly(painter, front_top, close=True)
    _poly(painter, [(cx - dx, cy - h + dy), (cx - dx, cy + dy)])
    _poly(painter, [(cx + dx, cy - h + dy), (cx + dx, cy + dy)])
    _poly(painter, [(cx, cy - h + 2 * dy), (cx, cy + 2 * dy)])
    _poly(painter, [(cx - dx, cy + dy), (cx, cy + 2 * dy), (cx + dx, cy + dy)])


def _block(painter, color):
    _cube(painter)


def _extrude(painter, color):
    _poly(painter, [(0.2, 0.62), (0.5, 0.75), (0.8, 0.62), (0.5, 0.49)], close=True)
    painter.save()
    _thin(painter, 0.8)
    _poly(painter, [(0.2, 0.62), (0.2, 0.84), (0.5, 0.97), (0.8, 0.84), (0.8, 0.62)])
    _poly(painter, [(0.5, 0.75), (0.5, 0.97)])
    painter.restore()
    _poly(painter, [(0.5, 0.42), (0.5, 0.08)])
    _poly(painter, [(0.36, 0.22), (0.5, 0.08), (0.64, 0.22)])


def _volume(painter, color):
    a, b, c, d = (0.15, 0.8), (0.85, 0.72), (0.55, 0.2), (0.48, 0.55)
    _poly(painter, [a, b, c], close=True)
    for p in (a, b, c):
        _poly(painter, [p, d])


def _arrow_arc(painter, color, clockwise: bool):
    rect = QRectF(0.2, 0.25, 0.6, 0.6)
    path = QPainterPath()
    if clockwise:
        path.arcMoveTo(rect, 200)
        path.arcTo(rect, 200, -220)
        tip = (
            0.2 + 0.6 * (0.5 + 0.5 * math.cos(math.radians(200))),
            0.25 + 0.6 * (0.5 - 0.5 * math.sin(math.radians(200))),
        )
        head = [(tip[0] + 0.18, tip[1] - 0.02), tip, (tip[0] + 0.05, tip[1] - 0.2)]
    else:
        path.arcMoveTo(rect, -20)
        path.arcTo(rect, -20, 220)
        tip = (
            0.2 + 0.6 * (0.5 + 0.5 * math.cos(math.radians(-20))),
            0.25 + 0.6 * (0.5 - 0.5 * math.sin(math.radians(-20))),
        )
        head = [(tip[0] - 0.18, tip[1] - 0.02), tip, (tip[0] - 0.05, tip[1] - 0.2)]
    painter.drawPath(path)
    _poly(painter, head)


def _undo(painter, color):
    _arrow_arc(painter, color, clockwise=True)


def _redo(painter, color):
    _arrow_arc(painter, color, clockwise=False)


def _reset_camera(painter, color):
    painter.save()
    _thin(painter, 0.9)
    for sx in (1, -1):
        for sy in (1, -1):
            x = 0.5 + sx * 0.36
            y = 0.5 + sy * 0.36
            _poly(painter, [(x - sx * 0.18, y), (x, y), (x, y - sy * 0.18)])
    painter.restore()
    painter.drawEllipse(QPointF(0.5, 0.5), 0.13, 0.13)


def _zoom_selection(painter, color):
    painter.drawEllipse(QPointF(0.42, 0.42), 0.24, 0.24)
    painter.save()
    pen = painter.pen()
    pen.setWidthF(0.13)
    painter.setPen(pen)
    _poly(painter, [(0.6, 0.6), (0.86, 0.86)])
    painter.restore()
    _dot(painter, color, 0.42, 0.42, 0.06)


def _letters(text: str):
    """Bold letters filling the icon (drawn in device units: text does not scale well under
    the unit-square transform)."""

    def draw(painter, color):
        painter.save()
        painter.resetTransform()
        font = painter.font()
        font.setPixelSize(int(SIZE * 0.5))
        font.setBold(True)
        painter.setFont(font)
        painter.setPen(color)
        painter.drawText(QRectF(0, 0, SIZE, SIZE), Qt.AlignmentFlag.AlignCenter, text)
        painter.restore()

    return draw


def _cube_face(which: str):
    """The isometric cube with one face filled: "xy" its top, "xz" its front, "yz" its
    right side (the view presets look at those faces)."""
    cx, cy = 0.5, 0.52
    dx, dy = 0.3, 0.17
    h = 0.36
    top = [(cx, cy - h), (cx + dx, cy - h + dy), (cx, cy - h + 2 * dy), (cx - dx, cy - h + dy)]
    front = [(cx - dx, cy - h + dy), (cx, cy - h + 2 * dy), (cx, cy + 2 * dy), (cx - dx, cy + dy)]
    side = [(cx, cy - h + 2 * dy), (cx + dx, cy - h + dy), (cx + dx, cy + dy), (cx, cy + 2 * dy)]
    face = {"xy": top, "xz": front, "yz": side}[which]

    def draw(painter, color):
        _cube(painter)
        painter.save()
        fill = QColor(color)
        fill.setAlphaF(0.55)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(fill)
        path = QPainterPath(QPointF(*face[0]))
        for x, y in face[1:]:
            path.lineTo(x, y)
        path.closeSubpath()
        painter.drawPath(path)
        painter.restore()

    return draw


def _iso(painter, color):
    _cube(painter)
    _dot(painter, color, 0.5, 0.16, 0.07)


def _nav3d(painter, color):
    painter.drawEllipse(QPointF(0.5, 0.5), 0.22, 0.22)
    painter.save()
    painter.translate(0.5, 0.5)
    painter.rotate(-25)
    painter.drawEllipse(QPointF(0.0, 0.0), 0.4, 0.13)
    painter.restore()
    _dot(
        painter,
        color,
        0.5 + 0.4 * math.cos(math.radians(25)) * 0.95,
        0.5 - 0.4 * math.sin(math.radians(25)) * 0.95,
        0.07,
    )


def _perspective(painter, color):
    _poly(painter, [(0.3, 0.3), (0.7, 0.3), (0.86, 0.8), (0.14, 0.8)], close=True)
    painter.save()
    _thin(painter)
    _poly(painter, [(0.3, 0.3), (0.5, 0.12), (0.7, 0.3)])
    _poly(painter, [(0.5, 0.12), (0.5, 0.8)])
    painter.restore()


def _zoom_window(painter, color):
    painter.save()
    pen = painter.pen()
    pen.setStyle(Qt.PenStyle.DashLine)
    _thin(painter, 0.8)
    painter.setPen(pen)
    painter.drawRect(QRectF(0.18, 0.22, 0.52, 0.44))
    painter.restore()
    painter.drawEllipse(QPointF(0.62, 0.6), 0.16, 0.16)
    _poly(painter, [(0.74, 0.72), (0.88, 0.86)])


def _grid(painter, color):
    painter.save()
    _thin(painter, 0.5)
    for k in range(1, 5):
        v = 0.1 + 0.8 * k / 5
        _poly(painter, [(v, 0.1), (v, 0.9)])
        _poly(painter, [(0.1, v), (0.9, v)])
    painter.restore()
    _poly(painter, [(0.1, 0.58), (0.9, 0.58)])
    _poly(painter, [(0.42, 0.1), (0.42, 0.9)])


def _group(painter, color):
    painter.save()
    _thin(painter, 0.8)
    painter.drawRect(QRectF(0.14, 0.3, 0.36, 0.5))
    painter.drawRect(QRectF(0.5, 0.2, 0.36, 0.6))
    painter.restore()
    path = QPainterPath(QPointF(0.2, 0.55))
    path.cubicTo(0.4, 0.4, 0.6, 0.7, 0.8, 0.5)
    painter.drawPath(path)


def _network(painter, color):
    path = QPainterPath(QPointF(0.12, 0.35))
    path.cubicTo(0.4, 0.1, 0.6, 0.6, 0.88, 0.3)
    painter.drawPath(path)
    path = QPainterPath(QPointF(0.12, 0.7))
    path.cubicTo(0.4, 0.5, 0.6, 0.95, 0.88, 0.65)
    painter.drawPath(path)
    _poly(painter, [(0.3, 0.12), (0.32, 0.88)])
    _poly(painter, [(0.68, 0.12), (0.7, 0.88)])


GLYPHS = {
    "select": _select,
    "point": _point,
    "line": _line,
    "arc": _arc,
    "polyline": _polyline,
    "spline": _spline,
    "nurbs": _nurbs,
    "freehand": _freehand,
    "sketch": _sketch,
    "finish": _finish,
    "cancel": _cancel,
    "connectors": _connectors,
    "chain": _chain,
    "structured": _structured,
    "unstructured": _unstructured,
    "network": _network,
    "group": _group,
    "block": _block,
    "extrude": _extrude,
    "volume": _volume,
    "undo": _undo,
    "redo": _redo,
    "reset-camera": _reset_camera,
    "zoom-selection": _zoom_selection,
    "view-xy": _letters("XY"),
    "view-xz": _letters("XZ"),
    "view-yz": _letters("YZ"),
    "view-iso": _iso,
    "navigation-3d": _nav3d,
    "perspective": _perspective,
    "zoom-window": _zoom_window,
    "grid": _grid,
}

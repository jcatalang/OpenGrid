"""Plot3D structured multi-block grids (whole or multi-grid; 2D or 3D with nk = 1; ASCII or
Fortran-unformatted binary).

Binary files are Fortran sequential unformatted, little-endian, with 4-byte record markers,
int32 dimensions and float64 coordinates; records: [nblocks] (multi-grid only), [dims of all
blocks], then one record per block holding all x, then all y (then all z), with i varying
fastest. ASCII files contain the same numbers separated by white space.
"""

import struct
from pathlib import Path

import numpy as np


def write_plot3d(
    path,
    grids: list[np.ndarray],
    binary: bool = True,
    three_d: bool = False,
    multi_grid: bool = True,
) -> None:
    """Write grids of shape (nj, ni, 3). Two-dimensional files store x and y only."""
    if not grids:
        raise ValueError("nothing to export")
    if not multi_grid and len(grids) != 1:
        raise ValueError("a single-grid ('whole') Plot3D file holds exactly one block")
    dims = []
    blocks = []
    for grid in grids:
        grid = np.asarray(grid, dtype=float)
        nj, ni, _ = grid.shape
        components = [grid[:, :, 0].reshape(-1), grid[:, :, 1].reshape(-1)]
        if three_d:
            dims.extend([ni, nj, 1])
            components.append(grid[:, :, 2].reshape(-1))
        else:
            dims.extend([ni, nj])
        blocks.append(np.concatenate(components))
    path = Path(path)
    if binary:
        with path.open("wb") as f:
            if multi_grid:
                _record(f, np.array([len(grids)], dtype="<i4").tobytes())
            _record(f, np.array(dims, dtype="<i4").tobytes())
            for values in blocks:
                _record(f, values.astype("<f8").tobytes())
        return
    lines = []
    if multi_grid:
        lines.append(str(len(grids)))
    per_block = len(dims) // len(grids)
    for b in range(len(grids)):
        lines.append(" ".join(str(d) for d in dims[b * per_block : (b + 1) * per_block]))
    for values in blocks:
        for start in range(0, len(values), 4):
            lines.append(" ".join(repr(float(v)) for v in values[start : start + 4]))
    path.write_text("\n".join(lines) + "\n")


def _record(f, payload: bytes) -> None:
    marker = struct.pack("<i", len(payload))
    f.write(marker)
    f.write(payload)
    f.write(marker)


def read_plot3d(
    path, three_d: bool = False, multi_grid: bool = True, binary: bool | None = None
) -> list[np.ndarray]:
    """Read grids written by :func:`write_plot3d` (or compatible files); returns arrays of
    shape (nj, ni, 3) (z = 0 for 2D files). ``binary=None`` detects the format."""
    path = Path(path)
    raw = path.read_bytes()
    if binary is None:
        binary = _looks_binary(raw)
    per_block = 2
    if three_d:
        per_block = 3
    if binary:
        records = _records(raw)
        position = 0
        if multi_grid:
            n_blocks = int(np.frombuffer(records[0], dtype="<i4")[0])
            position = 1
        else:
            n_blocks = 1
        dims = np.frombuffer(records[position], dtype="<i4").reshape(n_blocks, per_block)
        data = [
            np.frombuffer(r, dtype="<f8") for r in records[position + 1 : position + 1 + n_blocks]
        ]
    else:
        tokens = raw.decode("ascii").split()
        cursor = 0
        if multi_grid:
            n_blocks = int(tokens[0])
            cursor = 1
        else:
            n_blocks = 1
        dims = np.array([int(t) for t in tokens[cursor : cursor + n_blocks * per_block]]).reshape(
            n_blocks, per_block
        )
        cursor += n_blocks * per_block
        data = []
        for d in dims:
            count = int(np.prod(d)) * per_block
            data.append(np.array([float(t) for t in tokens[cursor : cursor + count]]))
            cursor += count
    grids = []
    for d, values in zip(dims, data, strict=True):
        ni = int(d[0])
        nj = int(d[1])
        count = ni * nj
        grid = np.zeros((nj, ni, 3))
        grid[:, :, 0] = values[:count].reshape(nj, ni)
        grid[:, :, 1] = values[count : 2 * count].reshape(nj, ni)
        if three_d:
            grid[:, :, 2] = values[2 * count : 3 * count].reshape(nj, ni)
        grids.append(grid)
    return grids


def _looks_binary(raw: bytes) -> bool:
    if len(raw) < 4:
        return False
    marker = struct.unpack("<i", raw[:4])[0]
    return marker in (4, 8, 12) and len(raw) >= 8 + marker


def _records(raw: bytes) -> list[bytes]:
    records = []
    position = 0
    while position < len(raw):
        (size,) = struct.unpack("<i", raw[position : position + 4])
        payload = raw[position + 4 : position + 4 + size]
        (end,) = struct.unpack("<i", raw[position + 4 + size : position + 8 + size])
        if end != size:
            raise ValueError(f"corrupt Fortran record at byte {position}")
        records.append(payload)
        position += 8 + size
    return records

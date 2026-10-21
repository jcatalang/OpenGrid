"""Importers and exporters.

The VTKHDF writers need h5py, loaded on first use: OpenGrid imports and runs without it
(only the VTKHDF / HDF5 exports then say it is missing).
"""

from opengrid.io.airfoil import AirfoilData, read_airfoil_dat
from opengrid.io.csv_points import read_points_csv
from opengrid.io.plot3d import read_plot3d, write_plot3d


def __getattr__(name: str):
    if name in ("write_vtkhdf_collection", "write_vtkhdf_grid"):
        from opengrid.io import vtkhdf

        return getattr(vtkhdf, name)
    raise AttributeError(f"module 'opengrid.io' has no attribute {name!r}")


__all__ = [
    "AirfoilData",
    "read_airfoil_dat",
    "read_plot3d",
    "read_points_csv",
    "write_plot3d",
    "write_vtkhdf_collection",
    "write_vtkhdf_grid",
]

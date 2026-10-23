OpenGrid for Linux (x86_64), self-contained
============================================

This folder holds everything OpenGrid needs: Python, Qt, VTK, Netgen, OpenCASCADE, HDF5,
CGNS and OpenGrid itself. Nothing else has to be installed: no conda, no compiler, no root.

Requirements of the machine
  * Linux on x86_64 (Intel or AMD) with glibc 2.28 or newer (the file GLIBC_MIN says the
    exact minimum): openSUSE Leap 15.3 and later, SLES 15 SP3 and later, RHEL / Rocky / Alma 8
    and later, Ubuntu 20.04 and later. Check with: ldd --version
  * For the window: an X11 or Wayland desktop with OpenGL (Mesa), which desktop installs
    have. On a server without a display, scripts run as they are; the window needs X
    forwarding (ssh -X) or a virtual display (xvfb-run).

Start
  tar -xzf OpenGrid-<version>-linux-x86_64.tar.gz
  cd OpenGrid-<version>-linux-x86_64
  ./opengrid                         # the application (the first start takes a few seconds)
  ./opengrid examples/output/12a_far_field.ogp    # ... with a project (after running 12a)
  ./opengrid-python examples/12a_airfoil_cases.py # a script with OpenGrid's Python
  ./install.sh                       # optional: menu entry, .ogp files open in OpenGrid,
                                     # ~/.local/bin/opengrid and opengrid-python

Keep the folder where you first start it: the first start fixes the environment's paths to
that place. To move it, extract the archive again in the new place.

Scripting: opengrid-python is a normal Python 3.12 with numpy, h5py, meshio and pyvista; see
examples/README.md and the scripting guide (docs/SCRIPTING.md in the repository).

If the window does not open
  * "Could not load the Qt platform plugin xcb": the X11 helper libraries are missing:
      sudo zypper install libxkbcommon-x11-0 libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1
  * OpenGL errors (GLX, "failed to create context"): install Mesa:
      sudo zypper install Mesa-libGL1 Mesa-libEGL1 Mesa-dri
  * Over ssh: connect with ssh -X (or -Y); or render off screen with xvfb-run.
Build from source instead: README.md of the repository, "Building from source", Linux.

Licences
OpenGrid is MIT-licensed (LICENSE). The folder licenses/ holds the licence texts of every
third-party package in env/ (Python, Qt, VTK, Netgen, OpenCASCADE, HDF5, CGNS, ...), listed with
their versions in licenses/README.txt.

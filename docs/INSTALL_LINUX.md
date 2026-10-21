# OpenGrid on Linux

Three ways, from the simplest. All of them work for a normal user, without root, on x86_64
(Intel or AMD); they were checked on openSUSE Leap 15.6.

| Way | You need | Takes |
|---|---|---|
| 1. The self-contained bundle | the archive (about 0.8 GB), copied to the machine | a minute |
| 2. The installer script, from source | internet (direct or through a proxy), 3 GB of disk | 10–20 minutes |
| 3. By hand, from source | internet (direct or through a proxy), 3 GB of disk | 10–20 minutes |

**A machine without internet** (common on SLES servers): way 1 only. Ways 2 and 3 download
about 1.5 GB of conda-forge packages; without access to conda-forge they cannot work, whatever
the downloader.

## 1. The self-contained bundle (no compilation)

`OpenGrid-<version>-linux-x86_64.tar.gz` holds Python, Qt, VTK, Netgen, OpenCASCADE, HDF5,
CGNS and OpenGrid. Nothing else is installed: no conda, no compiler.

```bash
tar -xzf OpenGrid-0.0.1-linux-x86_64.tar.gz
cd OpenGrid-0.0.1-linux-x86_64
./opengrid                                   # the application
./opengrid-python examples/12a_airfoil_cases.py   # a script with OpenGrid's Python
./install.sh                                 # optional: menu entry, .ogp files, ~/.local/bin links
```

The first start takes a few seconds: it fixes the environment's paths to the folder's place,
so extract the archive where you will keep it (to move it, extract it again).

To get it onto a machine without internet, copy it from any computer that has the file:
`scp OpenGrid-0.0.1-linux-x86_64.tar.gz user@machine:` (or a USB disk). Nothing else is
needed: the bundle never downloads anything.

The machine needs glibc 2.28 or newer (`ldd --version`; openSUSE Leap 15.3 and later, SLES 15
SP3 and later) and, for the window, an X11 or Wayland desktop with OpenGL. The launcher checks
the glibc version first and says so when it is too old. Scripts run without a display.

**Getting the bundle.** On GitHub: Actions > **Linux bundle** > Run workflow. It builds the
bundle on GitHub's x86_64 Linux machines, tests it in clean openSUSE Leap 15.6, SLES 15 SP6
and SP3 and Rocky Linux 8 containers, and runs the C++ and Python test suites (about half an
hour). The run's artifact **OpenGrid-linux-x86_64** is a zip holding the `.tar.gz`.
On any Linux x86_64 machine with internet: `OG_NATIVE=1 packaging/linux/build_bundle.sh`
(no Docker; `packaging/linux/test_bundle.sh [archive] [image]` then checks it in a clean
container).

## 2. The installer script (compiles on the machine)

From a checkout of the repository:

```bash
git clone https://github.com/jcatalang/OpenGrid.git
cd OpenGrid
packaging/linux/install_from_source.sh          # add --test to run the test suites as well
opengrid                                        # ~/.local/bin/opengrid
```

It downloads with `curl`, `wget` or Python, whichever the machine has; behind a proxy, set
`export https_proxy=http://host:port` first (conda uses it too). A Miniforge installer copied
by hand is used with `MINIFORGE_INSTALLER=/path/to/Miniforge3-Linux-x86_64.sh`. When the
downloads fail it says so and points to the bundle.

The script installs Miniforge into a private folder (`~/.local/share/opengrid/miniforge`, or
`$OPENGRID_HOME`) unless it is there already, creates the `opengrid` environment from
`environment.yml` with conda-forge's GCC, builds OpenGrid into it (editable: after a
`git pull`, run the script again) and adds the launchers `opengrid` and `opengrid-python` and
a menu entry. It does not use or change any other conda installation.

## 3. By hand

These are the script's steps, for when something needs adapting.

**Miniforge** (skip if you have a conda; any recent one works):
```bash
curl -fsSLO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
source ~/miniforge3/etc/profile.d/conda.sh
```

**The environment** (Python 3.12, Qt, VTK, IPython for the console, Netgen, OpenCASCADE,
HDF5, CGNS, Eigen, CMake, Ninja, the test tools) **and a compiler:**
```bash
cd OpenGrid
conda env create -f environment.yml
conda install -n opengrid -c conda-forge cxx-compiler
conda activate opengrid
export CXX=$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++
```
A C++20 compiler is needed: GCC 11 or newer, or Clang 15 or newer. The distribution's own may
be too old (openSUSE Leap 15's default GCC is 7), so the environment gets conda-forge's GCC.
`CXX` points to its prefixed name, which finds its own archiver and linker; with the plain
`c++` CMake needs `ar` from the system's binutils and stops with a message saying so when
there is none. A system GCC works too: `sudo zypper install gcc13-c++` and `export CXX=g++-13`
(Leap 15.6 has GCC 11 to 15 as `gcc11-c++` … `gcc15-c++`).

**Build and install** (into the environment, editable):
```bash
pip install --no-build-isolation -e .
python -m opengrid                    # the application
```

**Tests** (optional; on a server without a display, under a virtual one):
```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
sudo zypper install xvfb-run xauth which   # once, for the GUI tests without a display
xvfb-run -a pytest -q
```

## System packages

The conda environment brings all libraries but the OpenGL driver and, for the window, a few X11
helpers that every desktop install has. On a minimal or server install:

```bash
sudo zypper install Mesa-libGL1 Mesa-libEGL1 Mesa-dri libxkbcommon-x11-0   # openSUSE / SLES
sudo apt install libgl1 libegl1 libxkbcommon-x11-0                         # Debian / Ubuntu
```

## Troubleshooting

- **`Could not load the Qt platform plugin "xcb"`:** X11 helper libraries are missing:
  `sudo zypper install libxkbcommon-x11-0 libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1`.
- **OpenGL errors (`GLX`, `failed to create context`):** install Mesa (above). Over ssh, connect
  with `ssh -X`; otherwise run under `xvfb-run`.
- **`GCC ... is too old (need >= 11)`:** CMake found the system's old compiler: set `CXX` as
  above.
- **`No 'ar' archiver found`:** set `CXX` to the prefixed compiler, or `sudo zypper install
  binutils`.
- **`version GLIBC_2.xx not found` (the bundle):** the machine's glibc is older than the
  bundle's minimum; build from source (way 2 or 3) instead.
- **A different conda is active and things get mixed:** start a new shell, or
  `conda deactivate` until no environment is active, before running the installer.

## Checking a distribution without installing anything

`packaging/linux/check_build.sh [image]` builds and tests OpenGrid in a clean Docker container
(openSUSE Leap 15.6 on x86_64 by default), the way way 3 does it, and reports PASS or FAIL. On a
Mac with Apple silicon the x86_64 container is emulated: it works, slowly.

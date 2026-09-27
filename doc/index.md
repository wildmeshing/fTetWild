# fTetWild Documentation

fTetWild converts an arbitrary triangle soup into a high-quality tetrahedral mesh.
It is robust to imperfect input (self-intersections, gaps, non-manifold geometry) and
operates entirely in floating-point arithmetic.

**Paper**: Yixin Hu, Teseo Schneider, Bolun Wang, Denis Zorin, Daniele Panozzo.
*Fast Tetrahedral Meshing in the Wild*. ACM Transactions on Graphics (SIGGRAPH 2020).
[`doc/paper.pdf`](paper.pdf)

---

## Contents

### Getting started
| | |
|---|---|
| [Dependencies](dependencies.md) | System requirements, FetchContent libraries, vendored code |
| [Library API](library_api.md) | CMake integration, C++ API, minimal example, parameter reference |

### Algorithm
| | |
|---|---|
| [Overview](overview.md) | The four pipeline phases, key properties, comparison with TetWild |
| [Phase 1 — Preprocessing](preprocessing.md) | Input simplification, parallel edge collapsing, envelope checking |
| [Phase 2 — Triangle Insertion](triangle_insertion.md) | Background mesh, snapping, tet-subdivision table, retry loop |
| [Phase 3 — Mesh Improvement](mesh_improvement.md) | AMIPS energy, edge split/collapse/swap, parallel vertex smoothing |
| [Phase 4 — Filtering & Booleans](filtering_and_booleans.md) | Winding number, flood-fill, CSG trees, manifold extraction |

### Reference
| | |
|---|---|
| [Data Structures](data_structures.md) | `Mesh`, `MeshVertex`, `MeshTet`, `Parameters`, `Statistics` |
| [Envelope & AABB](envelope_and_aabb.md) | `AABBWrapper`, sampling-based and exact envelope checks |

---

## Quick build reference

```bash
# No system packages required by default — all deps are fetched automatically.

# Configure and build
cd build
bash configure.sh          # Release + TBB, long double fallback (no GMP needed)
cmake --build . --parallel $(nproc)

# Optional: enable exact GMP rational arithmetic (requires libgmp-dev)
bash configure.sh --with-gmp

# Run tests
ctest --output-on-failure

# Mesh a surface
./FloatTetwild_bin -i input.obj -o output.msh
```

See `build/configure.sh --help` for all configuration options.

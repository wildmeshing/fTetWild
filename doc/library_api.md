# Using fTetWild as a C++ Library

fTetWild exposes a single-function API for integration into other C++ projects.

---

## CMake Integration

Add fTetWild as a subdirectory (or via FetchContent) and link against the `FloatTetwild` target:

```cmake
add_subdirectory(fTetWild)
target_link_libraries(my_target PRIVATE FloatTetwild)
```

All transitive dependencies (libigl, geogram, TBB, spdlog, GMP, etc.) are pulled in automatically through the CMake target.

The headers are installed under the `floattetwild/` include prefix (matching the `#include <floattetwild/...>` style used in the source).

---

## API

**Header**: `<floattetwild/FloatTetwild.h>`

```cpp
namespace floatTetWild {

int tetrahedralization(
    GEO::Mesh&       sf_mesh,       // [in]  input surface mesh in geogram format
    Parameters       params,        // [in]  algorithm parameters (copied by value)
    Eigen::MatrixXd& VO,            // [out] output vertices, shape (N, 3)
    Eigen::MatrixXi& TO,            // [out] output tetrahedra, shape (M, 4)
    int              boolean_op  = -1,     // [in]  -1 = none, 0 = union, 1 = intersection, 2 = diff
    bool             skip_simplify = false // [in]  skip preprocessing step
);

}
```

Returns `EXIT_SUCCESS` (0) on success, `EXIT_FAILURE` (1) if the input is empty or malformed.

---

## Minimal Example

```cpp
#include <floattetwild/FloatTetwild.h>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/Logger.hpp>
#include <floattetwild/Parameters.h>

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <Eigen/Dense>

int main() {
    // Initialize geogram (required once per process)
    GEO::initialize();

    // Initialize fTetWild's spdlog logger
    floatTetWild::Logger::init(/*use_cout=*/true, /*log_file=*/"");

    // Load input surface mesh into geogram
    GEO::Mesh sf_mesh;
    GEO::mesh_load("input.obj", sf_mesh);

    // Configure parameters
    floatTetWild::Parameters params;
    params.eps_rel             = 1e-3;   // envelope size: 0.1% of bbox diagonal
    params.ideal_edge_length_rel = 0.05; // target edge length: 5% of bbox diagonal
    params.stop_energy         = 10.0;   // stop when max AMIPS < 10
    params.max_its             = 80;     // max optimization iterations

    // Run
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    int result = floatTetWild::tetrahedralization(sf_mesh, params, V, T);

    if (result != EXIT_SUCCESS) {
        // V and T contain an empty or partial mesh
        return 1;
    }

    // V: (N×3) matrix of vertex positions
    // T: (M×4) matrix of tetrahedron vertex indices (0-based)
    printf("Vertices: %ld  Tetrahedra: %ld\n", V.rows(), T.rows());
    return 0;
}
```

---

## Loading Input with MeshIO

To avoid the geogram `mesh_load` dependency for formats other than what geogram supports natively, use fTetWild's own `MeshIO` which handles `.off`, `.obj`, `.stl`, `.ply`:

```cpp
#include <floattetwild/MeshIO.hpp>

std::vector<floatTetWild::Vector3>  input_vertices;
std::vector<floatTetWild::Vector3i> input_faces;
GEO::Mesh sf_mesh;
std::vector<int> input_tags;

bool ok = floatTetWild::MeshIO::load_mesh(
    "input.off", input_vertices, input_faces, sf_mesh, input_tags
);
```

The `sf_mesh` output is what you pass to `tetrahedralization()`.

---

## Writing Output

Use `MeshIO::write_mesh()` to write the internal `Mesh` object to `.msh` (binary by default) or `.mesh`:

```cpp
#include <floattetwild/MeshIO.hpp>

// After tetrahedralization, if you have the internal Mesh object:
floatTetWild::MeshIO::write_mesh("output.msh", mesh, /*is_surface=*/false);

// To write a surface-only mesh (for mesh repair use case):
Eigen::MatrixXd V_sf;
Eigen::MatrixXi F_sf;
floatTetWild::get_surface(mesh, V_sf, F_sf);
igl::write_triangle_mesh("surface.obj", V_sf, F_sf);
```

For reading/writing `.msh` files independently of the algorithm, use `PyMesh::MshLoader` and `PyMesh::MshSaver` from `src/external/`.

---

## Parameters Reference

See [`data_structures.md`](data_structures.md) for the full `Parameters` field reference, and the CLI flags in the `README.md` for the command-line equivalents of each parameter.

Key defaults:

| Parameter | Default | Effect |
|---|---|---|
| `eps_rel` | 1e-3 | Envelope size as fraction of bbox diagonal |
| `ideal_edge_length_rel` | 0.05 | Target edge length as fraction of bbox diagonal |
| `ideal_edge_length_abs` | 0 | Absolute target edge length (overrides relative) |
| `stop_energy` | 10 | Stop when max AMIPS < this |
| `max_its` | 80 | Max optimization iterations |
| `smooth_open_boundary` | false | Smooth open boundary faces |
| `manifold_surface` | false | Guarantee manifold output surface |
| `disable_filtering` | false | Keep all tets (no exterior removal) |
| `use_floodfill` | false | Use flood-fill instead of winding number |
| `coarsen` | false | Aggressively coarsen output |

---

## Geogram Initialization

`GEO::initialize()` must be called exactly once before using any geogram functionality (including fTetWild). On Linux/macOS, also set:

```cpp
#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif
GEO::initialize();
```

To suppress geogram's console output, redirect its logger (as done in `src/main.cpp`):

```cpp
GEO::Logger* geo_logger = GEO::Logger::instance();
geo_logger->unregister_all_clients();
// optionally re-register a custom client that forwards to spdlog
geo_logger->set_pretty(false);
```

---

## Thread Safety

The `Statistics` singleton (`src/Statistics.h`) is thread-safe via a mutex. The `Mesh` object and all algorithm functions are **not** thread-safe — do not call `tetrahedralization()` concurrently on the same `Mesh`. Concurrent calls on independent `Mesh` objects are safe.

TBB parallelism is used internally during preprocessing and vertex smoothing. The number of threads is controlled by `params.num_threads` (default: hardware concurrency).

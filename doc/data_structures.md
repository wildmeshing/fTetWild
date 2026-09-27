# Data Structures

All types live in `namespace floatTetWild`.

---

## Scalar and Linear Algebra Types (`src/Types.hpp`)

```cpp
typedef double Scalar;          // float if FLOAT_TETWILD_USE_FLOAT is defined

typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
typedef Eigen::Matrix<Scalar, 2, 1> Vector2;
typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;
typedef Eigen::Matrix<int, 4, 1>    Vector4i;
typedef Eigen::Matrix<int, 3, 1>    Vector3i;
typedef Eigen::Matrix<int, 2, 1>    Vector2i;
```

Numerical tolerances defined as preprocessor constants:
- `SCALAR_ZERO` = 1e-8 (doubles) / 1e-6 (floats)
- `SCALAR_ZERO_2` = 1e-16 (doubles)

---

## MeshVertex (`src/Mesh.hpp`)

Represents a vertex of the output tetrahedral mesh.

```cpp
class MeshVertex {
    Vector3 pos;                 // 3D position
    std::vector<int> conn_tets;  // indices of incident tetrahedra

    bool is_on_surface   = false;
    bool is_on_boundary  = false;  // open boundary of input
    bool is_on_cut       = false;
    bool is_on_bbox      = false;  // vertex is on bounding box
    bool is_outside      = false;

    bool is_removed  = false;  // lazy-deleted
    bool is_freezed  = false;  // not moved during optimization

    Scalar sizing_scalar = 1;  // local sizing field multiplier
    Scalar scalar        = 0;  // winding number / output value
};
```

**Lazy deletion**: vertices are never physically erased from `Mesh::tet_vertices`. Instead `is_removed = true` marks them as free slots. `Mesh::v_empty_start` is the cursor to the first free slot, reset by `reset_v_empty_start()`.

---

## MeshTet (`src/Mesh.hpp`)

Represents a tetrahedron. Each tet stores its 4 vertex indices and per-face metadata.

```cpp
class MeshTet {
    Vector4i indices;  // global vertex indices [v0, v1, v2, v3]

    // Per-face flags (4 faces, indexed j = 0..3; face j is opposite vertex j)
    std::array<char, 4> is_surface_fs;  // NOT_SURFACE / KNOWN_SURFACE / KNOWN_NOT_SURFACE
    std::array<char, 4> is_bbox_fs;     // NOT_BBOX (-1) or tag id (≥0)
    std::array<int,  4> opp_t_ids;      // opposite tet id (OPP_T_ID_UNKNOWN / OPP_T_ID_BOUNDARY / index)
    std::array<char, 4> surface_tags;   // input face tag per face

    Scalar quality   = 0;      // conformal AMIPS energy
    Scalar scalar    = 0;      // winding number at centroid
    bool is_removed  = false;  // lazy-deleted
    bool is_outside  = false;  // outside the input surface
};
```

Face indexing convention: **face j is the face opposite to vertex j**, i.e., the face formed by the three vertices that are _not_ `indices[j]`. Helper methods:
- `find(v)` — returns the local index (0–3) of vertex `v`
- `find_opp(v0, v1, v2)` — returns the local index of the vertex not in {v0,v1,v2}

Sentinel values:
```cpp
#define NOT_SURFACE       SCHAR_MAX    // face not yet classified as surface
#define KNOWN_NOT_SURFACE -SCHAR_MAX/2
#define KNOWN_SURFACE      SCHAR_MAX/2
#define NOT_BBOX          -1
#define OPP_T_ID_UNKNOWN  -2
#define OPP_T_ID_BOUNDARY -1
#define MAX_ENERGY         1e50
```

---

## Mesh (`src/Mesh.hpp`)

The top-level mesh container.

```cpp
class Mesh {
    std::vector<MeshVertex> tet_vertices;
    std::vector<MeshTet>    tets;
    Parameters params;

    int  t_empty_start = 0;  // first free tet slot
    int  v_empty_start = 0;  // first free vertex slot
    bool is_limit_length       = true;
    bool is_closed             = true;
    bool is_input_all_inserted = false;
    bool is_coarsening         = false;
};
```

Key helpers:
- `get_v_num()` / `get_t_num()` — count of non-removed vertices / tets
- `get_max_energy()` / `get_avg_energy()` — AMIPS statistics over live tets
- `reset_t_empty_start()` / `reset_v_empty_start()` — recompute empty-slot cursors
- `one_ring_vertex_sets(threshold, concurrent_sets, serial_set)` — partition vertices into independent sets for TBB parallel smoothing
- `partition(n_parts, tets_id)` — spatial partitioning for parallel operations

---

## Parameters (`src/Parameters.h`)

Holds all algorithm parameters. Computed values are populated by `Parameters::init(bbox_diag_length)` which must be called once after loading the input.

### User-set parameters

| Field | Default | Description |
|---|---|---|
| `eps_rel` | 1e-3 | ε = eps_rel × bbox diagonal |
| `ideal_edge_length_rel` | 1/20 | ℓ = ideal_edge_length_rel × bbox diagonal |
| `ideal_edge_length_abs` | 0 | Absolute ℓ (overrides relative if > 0) |
| `stop_energy` | 10 | Stop optimization when max AMIPS < this |
| `max_its` | 80 | Maximum optimization iterations |
| `smooth_open_boundary` | false | Laplacian-smooth open boundary faces |
| `manifold_surface` | false | Force manifold output surface |
| `disable_filtering` | false | Skip exterior tet removal |
| `use_floodfill` | false | Use flood-fill instead of winding number |
| `use_general_wn` | false | Generalized winding number |
| `coarsen` | false | Aggressively coarsen output |
| `num_threads` | max | TBB thread count |

### Derived parameters (set by `init()`)

| Field | Formula |
|---|---|
| `ideal_edge_length` | `bbox_diag_length × ideal_edge_length_rel` |
| `eps_input` | `bbox_diag_length × eps_rel` |
| `eps` | `eps_input` adjusted for stage |
| `split_threshold` | `ideal_edge_length × 4/3` |
| `collapse_threshold` | `ideal_edge_length × 4/5` |
| `eps_simplification` | `eps × 0.8` |
| `eps_coplanar` | `min(eps × 0.2, bbox_diag × 1e-6)` |
| `min_edge_length` | `bbox_diag × eps_rel` |

---

## Statistics (`src/Statistics.h`)

Thread-safe singleton collecting timing and quality data for each pipeline stage. Written to a CSV file at the end of execution. Stage IDs:

```
0 = init        1 = preprocessing   2 = tetrahedralization
3 = cutting     4 = optimization    5 = winding_number
6 = splitting   7 = collapsing      8 = swapping   9 = smoothing
```

CSV columns: `stage_id, time_s, v_count, t_count, max_energy, avg_energy, uninserted_faces, peak_mem_mb`

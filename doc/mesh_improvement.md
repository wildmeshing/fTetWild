# Phase 3: Mesh Improvement

**Source**: `src/MeshImprovement.h/.cpp`, `src/EdgeSplitting.h/.cpp`, `src/EdgeCollapsing.h/.cpp`, `src/EdgeSwapping.h/.cpp`, `src/VertexSmoothing.h/.cpp`, `src/LocalOperations.h/.cpp`  
**Entry point**: `floatTetWild::optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, ops)`

---

## Purpose

After triangle insertion, the mesh contains many poorly-shaped elements — especially around inserted faces where the tet-subdivision table cuts tets into small sub-tets. Phase 3 iteratively improves element quality by applying four types of local topology/geometry operations while maintaining the two validity invariants (no inverted tets, tracked surface within ε-envelope).

---

## Quality Metric: Conformal AMIPS Energy

Every operation uses the **conformal AMIPS 3D energy** to measure element quality:

```
AMIPS(T) = tr(Jᵀ J) / det(J)^(2/3)
```

where J is the Jacobian of the affine map from a regular reference tetrahedron to T.

- Minimum value: **3** (achieved by a regular tetrahedron)
- Range: [3, +∞)
- Scale-invariant: same quality for any scaling of T
- Differentiable: enables gradient-based vertex smoothing

**Numerical instability fix**: When AMIPS(T) > 10⁸, standard floating-point evaluation is non-deterministic (different vertex permutations give wildly different values — up to 4 orders of magnitude difference). The implementation detects this via `is_energy_unstable()` and falls back to higher-precision arithmetic for the intermediate cubed-energy computation, then takes the cubic root. This prevents unnecessary over-refinement.

Two backends are available for the fallback, selected at compile time via `FLOAT_TETWILD_USE_GMP`:
- **Default (OFF)**: `long double` (80-bit extended precision on x86/x86_64). No extra dependency.
- **GMP (ON)**: `triwild::Rational` backed by GMP arbitrary-precision rationals. Exact on all platforms, requires `libgmp-dev`.

---

## The Four Local Operations

### 1. Edge Splitting (`src/EdgeSplitting.cpp`)

Splits an edge [v₀, v₁] by inserting a new midpoint vertex. All incident tetrahedra are replaced by two new tets each.

- **Triggered when**: edge length > `split_threshold` = 4/3 × ℓ
- **Accepted when**: no new tet is inverted, no new tet exceeds a quality threshold, and if any split face was tracked surface, all new tracked faces remain within ε-envelope
- **Priority queue**: edges sorted by length (longest first)

### 2. Edge Collapsing (`src/EdgeCollapsing.cpp`)

Collapses an edge [v₀, v₁] by merging one vertex onto the other (or to a new midpoint). All tets sharing the collapsed edge are removed.

- **Triggered when**: edge length < `collapse_threshold` = 4/5 × ℓ
- **Rejection conditions**: collapse creates inverted tets, violates envelope, or removes a surface feature
- **Priority queue**: edges sorted by length (shortest first)
- **Surface vertices**: if v₀ is on the tracked surface, it may only move to positions that keep the surface within ε

### 3. Edge Swapping (`src/EdgeSwapping.cpp`)

Replaces an edge with a better-placed edge in the same "diamond" of surrounding tets (2→3 flip, 3→2 flip, or similar). This changes topology without moving vertex positions.

- **Triggered when**: swapping would reduce the maximum AMIPS energy among affected tets
- **Exact predicates** verify that no element is inverted after the swap
- The full set of swaps is performed round-robin until no further improvement is possible

### 4. Vertex Smoothing (`src/VertexSmoothing.cpp`)

Moves each vertex to a new position that minimizes the maximum AMIPS energy of its incident tets. This is a local unconstrained optimization using Newton's method on the AMIPS energy.

- **Parallelized** using a graph-coloring scheme: vertices are partitioned into independent sets (no two vertices in the same set share an incident tet) using `Mesh::one_ring_vertex_sets()`. Each set is processed in parallel with TBB, then the next set serially.
- **Surface vertices** are constrained: only moves that keep all tracked faces within ε are accepted.
- **Boundary vertices** (on open boundaries) are projected back to the boundary after smoothing.

---

## Optimization Loop

```
for iter = 1 to max_its:
    apply edge_splitting
    apply edge_collapsing
    apply edge_swapping
    apply vertex_smoothing

    if iter % 3 == 0:
        retry_failed_triangle_insertions()

    max_e = get_max_energy()
    if max_e < stop_energy:
        break

    update_scaling_field(mesh, max_e)  # if sizing field enabled
```

The `ops` parameter to `optimization()` is a 4-bit mask `{split, collapse, swap, smooth}` that enables/disables individual operations (used internally for debugging and for coarsening mode).

---

## Scaling Field

If a background sizing mesh is provided (`--bg-mesh`), the `sizing_scalar` of each vertex is interpolated from the background mesh's vertex scalar field. This scalar multiplicatively adjusts the split/collapse thresholds at each vertex, enabling adaptive mesh density.

If coarsening mode is enabled (`--coarsen`), the algorithm applies `apply_coarsening()` after the main optimization loop, aggressively collapsing edges to reduce the mesh to its coarsest valid form.

---

## Acceptance Criteria for All Operations

Every proposed local operation is rolled back if:
1. Any resulting tetrahedron is inverted (checked with exact orient3d predicates).
2. Any resulting tetrahedron has AMIPS energy worse than the element it replaces (for smoothing, the condition is the neighborhood's max energy).
3. Any tracked surface face is moved outside the ε-envelope after the operation.
4. Any bounding-box face would be disturbed.

---

## Empty Slot Management

The mesh uses lazy deletion: removed vertices and tets are flagged with `is_removed = true` rather than actually erased. Periodically (in `cleanup_empty_slots()`), when the fraction of removed slots exceeds a threshold (default 70%), the vectors are compacted to reclaim memory. The empty-slot cursors `t_empty_start` and `v_empty_start` are reset after each compaction via `reset_t_empty_start()` / `reset_v_empty_start()`.

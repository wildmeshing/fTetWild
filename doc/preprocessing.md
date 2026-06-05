# Phase 1: Input Preprocessing

**Source**: `src/Simplification.h/.cpp`  
**Entry point**: `floatTetWild::simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify)`

---

## Purpose

The raw input triangle soup often contains near-duplicate vertices, degenerate or nearly-degenerate triangles, and small topological noise. Preprocessing simplifies the input while keeping the result inside the **ε-prep envelope** (= 0.8 × ε by default), leaving room for the insertion snapping tolerance and for surface vertices to move during mesh improvement.

The preprocessing envelope is intentionally smaller than the final ε (80% of ε) for two reasons:
1. It preserves room for the snapping tolerance δ used during triangle insertion.
2. It prevents surface vertices from being too close to the envelope boundary, giving them freedom to move during mesh quality optimization.

---

## Steps

### 1. Remove Duplicate Vertices (`remove_duplicates`)

Merge input vertices whose distance is below `SCALAR_ZERO` (ε_zero = 1e-8). After merging, degenerate faces (with two identical vertex indices) are discarded.

### 2. Edge Collapsing (`collapsing`)

Iteratively collapse short edges of the input surface mesh, subject to two conditions:

1. **Manifold constraint**: the edge must be a manifold edge (at most two incident triangles), and all vertex-adjacent edges must also be manifold. This prevents destroying topological structure.
2. **Envelope constraint**: collapsing the edge must not move any triangle outside the `eps_simplification` envelope (verified using the AABB tree).

This step is **parallelized** using a 2-coloring strategy (see Parallel Coloring below).

### 3. Edge Swapping (`swapping`)

Swap edges to improve the triangulation of the input surface. Only swaps that keep all faces inside the envelope are accepted.

### 4. Flattening (`flattening`)

Merge nearly-coplanar vertices that are very close to each other but survived the collapsing step.

---

## Parallel Coloring for Edge Collapsing

The parallelization strategy marks a safe independent set of edges in each pass:

1. Color all input triangles **white**.
2. Mark an edge as **parallel-independent** if _all_ vertex-adjacent triangles are white. Color those triangles **black**.
3. All marked parallel-independent edges are collapsed simultaneously using TBB.
4. Iterate until fewer than 0.01% of original vertices are removed in a pass.

This guarantees no two simultaneously collapsed edges share a vertex or adjacent triangle, making the parallel execution safe.

---

## Envelope Checking

The `AABBWrapper` is used to sample each modified triangle and check whether the samples lie within the AABB tree's proximity threshold. The AABB tree is built from the original input surface `sf_mesh` and queried using `MeshFacetsAABBWithEps::is_in_envelope()`.

When `FLOAT_TETWILD_WITH_EXACT_ENVELOPE` is enabled (compile flag), the `FastEnvelope` library replaces the sampling-based check with an exact envelope containment test.

---

## Output

After preprocessing:
- `input_vertices` — reduced and deduplicated 3D vertex positions
- `input_faces` — reduced triangle index list
- `input_tags` — per-face integer tags (used for boolean/CSG operations)

These are passed directly to the Delaunay tetrahedralization and triangle insertion stages.

# Envelope and AABB Structures

**Source**: `src/AABBWrapper.h/.cpp`, `src/external/mesh_AABB.h/.cpp`

---

## The ε-Envelope

The ε-envelope is the central robustness mechanism in fTetWild. It is a tube of radius ε centered on the input surface: every face of the output mesh that corresponds to an input triangle is guaranteed to lie within this tube.

The envelope allows the algorithm to:
- Tolerate self-intersections, gaps, and degeneracies in the input (all imperfections smaller than ε are resolved automatically).
- Use floating-point arithmetic throughout without exact arithmetic, since small errors are absorbed by ε.
- Perform snapping during triangle insertion by moving vertices onto the insertion plane without violating correctness.

**Default ε**: `1e-3 × d` where d is the diagonal of the bounding box of the input.

---

## AABBWrapper (`src/AABBWrapper.h`)

`AABBWrapper` is the central proximity-query object. It wraps multiple AABB trees and exposes uniform query APIs for the two envelope implementations.

### Trees Maintained

```cpp
class AABBWrapper {
    const GEO::Mesh& sf_mesh;           // original input surface
    GEO::Mesh        b_mesh;            // background mesh (input vertices + grid points)
    GEO::Mesh        tmp_b_mesh;        // temporary background mesh

    MeshFacetsAABBWithEps sf_tree;      // AABB for sf_mesh
    shared_ptr<MeshFacetsAABBWithEps> b_tree;     // AABB for b_mesh
    shared_ptr<MeshFacetsAABBWithEps> tmp_b_tree; // AABB for tmp_b_mesh

#ifdef NEW_ENVELOPE
    // Exact envelope instances (when FLOAT_TETWILD_WITH_EXACT_ENVELOPE is ON)
    fastEnvelope::FastEnvelope sf_tree_exact;
    fastEnvelope::FastEnvelope sf_tree_exact_simplify;
    fastEnvelope::FastEnvelope b_tree_exact;
    fastEnvelope::FastEnvelope tmp_b_tree_exact;
#endif
};
```

### Key Query Operations

All geometry queries in the mesh improvement and insertion stages go through `AABBWrapper`:

- **`is_out_envelope(triangle, tree, params)`** — returns true if the triangle lies outside the ε-envelope (used to reject local operations that would violate the envelope).
- **`get_sf_diag()`** — returns the bounding box diagonal of `sf_mesh`, used to compute ε and ℓ from their relative values.

### Initialization Sequence

```cpp
// 1. Wrap the input surface
AABBWrapper tree(sf_mesh);

// 2. Initialize the exact envelope (only if NEW_ENVELOPE defined)
tree.init_sf_tree(input_vertices, input_faces, params.eps);

// 3. After preprocessing: build the background mesh tree
tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
```

---

## MeshFacetsAABBWithEps (`src/external/mesh_AABB.h`)

A custom AABB tree built on top of geogram's `GEO::MeshFacetsAABB`. Extended with an epsilon tolerance for proximity queries.

**Envelope containment check**: To test if a triangle T is inside the ε-envelope, the algorithm:
1. Samples T with a uniform point set at spacing proportional to ε.
2. For each sample point, queries the AABB for the nearest input face.
3. If all sample distances to the nearest face are ≤ ε, the triangle is inside the envelope.

Sampling introduces a conservative over-estimate: the sampling error is explicitly compensated so that the check is never falsely accepting (i.e., it may reject a triangle that is actually inside the envelope, but never accepts one that is outside).

---

## Exact Envelope (Optional)

When compiled with `-DFLOAT_TETWILD_WITH_EXACT_ENVELOPE=ON`, the preprocessor macro `NEW_ENVELOPE` is defined and the `FastEnvelope` library (github: wangbolun300/fast-envelope) is used instead of the sampling-based check.

`FastEnvelope` computes an exact answer to "is triangle T inside the ε-envelope?" using a BSP-tree over the thickened input faces. This eliminates false negatives (triangles rejected by the conservative sampling check) and can improve both output quality and running time for inputs with many near-envelope triangles.

The tradeoff is a significantly longer CMake configure/build time (the FastEnvelope library must be built) and a small runtime overhead for building the exact envelope structure.

**Two exact trees are built for the input surface**:
- `sf_tree_exact` — built with ε
- `sf_tree_exact_simplify` — built with 0.8ε (for the preprocessing stage, which uses a smaller envelope)

---

## Envelope in Local Operations

In `src/LocalOperations.h`:

```cpp
bool is_out_envelope(Mesh& mesh, int v_id, const Vector3& new_pos,
                     const AABBWrapper& tree);
bool is_out_boundary_envelope(const Mesh& mesh, int v_id, const Vector3& new_pos,
                              const AABBWrapper& tree);
```

Before any local operation (split, collapse, swap, smooth) moves a vertex, these functions verify that all tracked surface faces incident to the vertex remain within ε after the move. If either returns true, the operation is rolled back.

The key internal function is `sample_triangle_and_check_is_out()`, which samples a triangle and queries the AABB for the nearest face — all with a cached `prev_facet` hint to accelerate repeated queries on nearby triangles.

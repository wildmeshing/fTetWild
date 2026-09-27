# Phase 2: Incremental Triangle Insertion

**Sources**: `src/FloatTetDelaunay.h/.cpp`, `src/TriangleInsertion.h/.cpp`, `src/CutMesh.h/.cpp`, `src/auto_table.hpp/.cpp`  
**Entry points**:
- `FloatTetDelaunay::tetrahedralize(...)` — builds the initial background mesh
- `insert_triangles(...)` — performs the incremental insertion loop

---

## 2a. Background Mesh Generation

Before any triangle insertion, fTetWild creates an initial tetrahedral mesh that fills the input's expanded bounding box.

**Steps**:
1. Expand the bounding box by 2ε in all directions (so surface points have room to move).
2. Add 8 bounding-box corner vertices.
3. Add additional grid points at spacing ℓ/20 (where ℓ = `ideal_edge_length`), skipping any point within ε of an input face. These uniform points improve the initial element quality by filling the space evenly.
4. Run Geogram's Delaunay 3D tetrahedralization (`GEO::Delaunay_3d`) on the combined point set.
5. Populate the `Mesh` structure and initialize adjacency (`opp_t_ids`).

The result is a valid, inversion-free tetrahedral mesh that does not yet conform to any input triangles.

---

## 2b. Incremental Triangle Insertion

### Overview

Each input triangle T is inserted one at a time into the current `Mesh M`. A failed insertion (one that would create an inverted or degenerate element) is rolled back, and the triangle is deferred for later retry after the mesh quality has improved.

Insertion of a single triangle T has three main stages:

```
┌──────────────────────┐
│ Find cut tetrahedra  │  Identify TI = {tets that T cuts}
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Snapping             │  Move mesh vertices within δ of T's plane P
└──────────┬───────────┘  onto P (if doing so doesn't invert any tet)
           │
           ▼
┌──────────────────────┐
│ Table-based          │  Subdivide all tets in TI using the
│ tet-subdivision      │  precomputed tet-subdivision table
└──────────────────────┘
```

---

### Finding Cut Tetrahedra (TI)

A triangle T *cuts* a tetrahedron T if:
- T is completely contained inside T, **or**
- T cuts through at least one face of T (the intersection contains interior points of both T and the face).

Pure edge/vertex touches that contain no interior points are excluded. Exact predicates (Shewchuk) are used to determine containment and face-triangle intersection.

The set TI is grown iteratively: any tet adjacent to a vertex that was snapped and lies in T's cutting region is added.

---

### Snapping

Snapping reduces insertion failures by moving mesh vertices that are within a distance δ of T's plane P onto P itself, as long as this does not create inverted tetrahedra:

**δ values**:
- First insertion pass: δ = max(ε_zero, 10⁻³ × ε) — larger, more aggressive
- Subsequent passes: δ = ε_zero = 1e-8 — minimal

**Snapping algorithm (4 steps)**:
1. Find all vertices of TI within distance δ of P → set Vδ.
2. Move each vertex in Vδ to its closest point on P, if no element would be inverted. If a vertex cannot be moved to P, instead some vertices of the "cover face" F are snapped to the vertex.
3. For each vertex in Vδ, expand TI with any neighbor tets that (a) are cut by P and (b) have a face that projects over T.
4. Repeat until no new tets are added.

If moving a boundary vertex of the cover F would invalidate the cover (by changing F's boundary), TI is first expanded to include the 1-ring neighborhood until the affected vertex becomes an interior vertex of F.

---

### Table-Based Tet Subdivision (`src/auto_table.hpp/.cpp`)

After snapping, every tet T in TI is subdivided according to which of its 6 edges are cut by P. The precomputed **tet-subdivision table** maps an (primary_index, secondary_index) pair to a list of sub-tets.

**Primary index (I)**: 6-bit binary string indicating which edges are cut. 64 combinations exist, of which 23 are geometrically impossible (e.g. 5 or 6 cut edges). The 41 realizable configurations form 7 symmetry classes:

| Class | Description |
|---|---|
| (1) | No cut edges |
| (2) | 1 edge cut |
| (3) | 2 adjacent edges cut (sharing a vertex) |
| (4) | 2 opposite edges cut |
| (5) | 3 edges forming a path |
| (6) | 3 edges on one face (triangle of cuts) |
| (7) | 4 edges cut |

Classes (4) and (6) can only occur on neighboring tets of TI (not on tets that T directly cuts through).

**Secondary index (II)**: When two edges on a face are cut, two triangulations of that face exist. The secondary index selects one. The choice is deterministic and unique: for face [v₀, v₁, v₂] with intersection points p₁ and p₂, use the triangulation containing edge [p₂, v₁] if `label(v₁) > label(v₂)`, otherwise use [p₁, v₂]. This rule is based solely on global vertex ordering, ensuring consistency across adjacent tets without any look-ahead.

All generated sub-tetrahedra are checked for volume > ε³_zero. If any sub-tet is degenerate, the entire insertion is rolled back.

---

### Open-Boundary Edge Preservation

After a triangle is inserted, shared edges between non-coplanar adjacent triangles are naturally preserved (the plane of the second triangle will cut through the cover of the first). However, **open-boundary edges** (edges with only one incident non-coplanar triangle) require explicit handling.

For each open-boundary edge e of triangle T:
1. Project e and the cover F of T onto the best-fitting plane P' of F's vertices.
2. Compute the intersections of e's projection with the face projections in 2D.
3. Lift the 2D intersection points back to 3D on the faces of F.
4. Further subdivide the affected neighboring tetrahedra using the same table-based scheme.

If this fails numerically, the triangle insertion is rolled back and deferred.

---

### Retry Loop

Insertion is attempted for all input triangles sequentially. Any triangle that cannot be inserted (due to near-degenerate elements or numerical failures) is added to a "failed" list. During mesh improvement (Phase 3), every 3 optimization iterations the insertion of failed triangles is reattempted on the improved mesh. This continues until all triangles are inserted or the optimization terminates.

In practice on the Thingi10k dataset of 10,000 models, all triangles are always successfully inserted.

---

### Validity Invariant

At every point during triangle insertion (and throughout the entire algorithm), the mesh is maintained in a **valid** state:
1. Every tetrahedron has positive volume (checked with exact Shewchuk predicates).
2. Every successfully inserted triangle (the *tracked surface*) lies within ε of the input.

Any operation that would violate either condition is immediately rolled back.

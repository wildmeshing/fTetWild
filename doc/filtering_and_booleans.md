# Phase 4: Filtering and Boolean Operations

**Source**: `src/MeshImprovement.h/.cpp`, `src/CSGTreeParser.hpp/.cpp`  
**Key functions**: `filter_outside()`, `filter_outside_floodfill()`, `boolean_operation()`, `manifold_surface()`, `get_surface()`

---

## Purpose

After mesh improvement, the volumetric mesh fills the entire expanded bounding box. Phase 4 classifies each tetrahedron as *inside* or *outside* the input surface and removes the outside ones, producing the final mesh conforming to the input boundary.

---

## Interior Classification Methods

Three strategies are available, selected by command-line flags:

### 1. Winding Number (default)

For each non-removed tetrahedron, compute the **winding number** of its centroid with respect to the tracked surface. The winding number measures how many times the surface wraps around a point:
- Integer values (0 or ±1) indicate points clearly inside or outside.
- The threshold is 0.5: tetrahedra with |winding number| > 0.5 are kept as interior.

Two winding number variants exist:
- **Fast winding number** (`filter_outside(mesh)`) — uses the tracked surface, which has already been inserted into the mesh. This is the default.
- **General winding number** (`--use-general-wn`) — uses `filter_outside(mesh, input_vertices, input_faces)`, evaluating with respect to the original input surface rather than the tracked surface.

Winding numbers correctly classify interiors of non-manifold and non-orientable inputs, unlike signed distance.

### 2. Flood Fill (`--use-floodfill`)

`filter_outside_floodfill(mesh)` starts from tetrahedra adjacent to the bounding box (which are known to be outside) and floods inward, stopping at tracked surface faces. Any tet reachable from the bounding box without crossing a surface face is marked outside.

This is simpler and faster than winding numbers but requires the inserted surface to be watertight. It is useful when the input is known to be closed and the winding number is expensive.

### 3. Skip Filtering (`--disable-filtering`)

No tetrahedra are removed. The full bounding-box mesh is output. Useful for debugging insertion correctness.

---

## Open Boundary Handling (`--smooth-open-boundary`)

For inputs with open boundaries (non-closed surfaces), the filtering step must not remove tetrahedra on the open side. The `smooth_open_boundary()` function applies a few passes of Laplacian smoothing to the triangles on the open side of the surface before exterior tets are removed, producing a smoother cap over open regions.

After smoothing, exterior tets are removed in the same way as the default winding number approach.

---

## Boolean / CSG Operations

fTetWild supports approximate Boolean operations on triangle soups (union, intersection, difference) that do not need to be closed or manifold.

### Simple Two-Mesh Boolean (`--op 0/1/2`)

For two-mesh operations:
1. Input meshes are tagged (face tags 1 and 2 in `input_tags`).
2. Triangle insertion tracks the provenance of each face (which input mesh it came from).
3. After mesh improvement, `boolean_operation(mesh, op)` evaluates the winding number of each tet's centroid separately for each input mesh.
4. Tets are retained based on the Boolean logic:
   - **Union (0)**: inside mesh 1 OR inside mesh 2
   - **Intersection (1)**: inside mesh 1 AND inside mesh 2
   - **Difference (2)**: inside mesh 1 AND NOT inside mesh 2

### CSG Trees (`--csg tree.json`)

A JSON CSG tree allows arbitrary combinations of Boolean operations on multiple input meshes.

**JSON format**: A recursive tree where leaf nodes are mesh filenames and internal nodes are operations:
```json
{
  "operation": "union",
  "left": { "mesh": "object1.stl" },
  "right": {
    "operation": "difference",
    "left": { "mesh": "object2.obj" },
    "right": { "mesh": "object3.off" }
  }
}
```

`CSGTreeParser::get_meshes()` extracts the list of meshes and assigns integer IDs. After mesh improvement, `boolean_operation(mesh, csg_tree_with_ids, meshes)` computes a winding number for each mesh independently and evaluates the tree to decide which tets to retain.

---

## Surface Extraction

After filtering, the boundary faces of the remaining tets form the output surface.

### Standard (`get_surface`)

Collects all non-removed tet faces that are not shared by two non-removed tets (i.e., boundary faces). These faces form the output surface mesh (V_sf, F_sf).

### Manifold Output (`--manifold-surface`)

The raw output surface may be non-manifold at edges shared by more than two tets, or at vertices where the surface branches. `manifold_surface()` resolves this:
1. `manifold_edges()` — detects non-manifold edges and splits them.
2. `manifold_vertices()` — duplicates non-manifold vertices so the surface is globally manifold (using the algorithm of Attene et al. 2009). Note: duplicate vertices share identical 3D positions; only the topological structure differs.

The manifold extraction is needed when using the output as a repaired surface mesh (mesh repair application).

---

## Correct Tracked Surface Orientation

Before filtering, `correct_tracked_surface_orientation(mesh, tree)` re-examines each tracked surface face and corrects its normal orientation based on proximity to the original input surface. This ensures the winding number computation is consistent.

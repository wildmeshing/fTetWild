# fTetWild Algorithm Overview

fTetWild (Fast Tetrahedral Meshing in the Wild) converts an arbitrary triangle soup into a high-quality tetrahedral mesh. It requires no assumption about input quality: self-intersections, gaps, duplicate vertices, degenerate faces and non-manifold regions are all handled automatically.

**Reference**: Yixin Hu, Teseo Schneider, Bolun Wang, Denis Zorin, Daniele Panozzo. *Fast Tetrahedral Meshing in the Wild*. ACM Transactions on Graphics (SIGGRAPH 2020).

---

## Core Idea

The algorithm operates within an ε-envelope: a tube of radius ε around the input surface. The output mesh boundary is guaranteed to lie inside this envelope, so imperfections smaller than ε are resolved automatically. The default ε is 10⁻³ × d, where d is the diagonal of the bounding box.

Unlike TetWild (its predecessor), fTetWild works entirely in **floating-point arithmetic**. It maintains an inversion-free tetrahedral mesh at every algorithmic step, eliminating the expensive exact rational arithmetic that made TetWild slow.

---

## The Four Phases

```
Input surface mesh
       │
       ▼
┌─────────────────────┐
│ 1. Preprocessing    │  Simplify input, merge near-duplicate vertices,
│    (Simplification) │  collapse short edges while staying in ε-envelope.
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│ 2. Background mesh  │  Delaunay tetrahedralization on preprocessed
│  + Triangle Insert. │  points + grid points. Then insert input
│    (incremental)    │  triangles one at a time using a table-based
│                     │  tet-subdivision scheme (snapping + cut).
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│ 3. Mesh Improvement │  Iterative local operations (split, collapse,
│    (optimization)   │  swap, smooth) minimizing conformal AMIPS
│                     │  energy. Retries failed insertions every 3
│                     │  iterations.
└────────┬────────────┘
         │
         ▼
┌─────────────────────┐
│ 4. Filtering        │  Remove tetrahedra outside the surface using
│    (winding number  │  winding number or flood-fill. Boolean CSG
│     or flood-fill)  │  operations are applied here.
└────────┬────────────┘
         │
         ▼
Output tetrahedral mesh (.msh / .mesh)
```

---

## Key Properties

| Property | Value |
|---|---|
| Always produces valid output | Yes (inversion-free, floating-point) |
| Tolerates imperfect input | Yes (self-intersections, gaps, non-manifold) |
| Envelope guarantee | Boundary within ε of input |
| Default ε | 10⁻³ × bbox diagonal |
| Default target edge length ℓ | bbox diagonal / 20 |
| Mesh quality metric | Conformal AMIPS energy ∈ [3, +∞), optimal = 3 |
| Default stop energy | 10 |
| Max optimization iterations | 80 |
| Parallelism | TBB (preprocessing + vertex smoothing) |

---

## Comparison With TetWild

| | TetWild | fTetWild |
|---|---|---|
| Arithmetic | Rational (exact) | Floating-point |
| Triangle insertion | BSP (all at once) | Incremental (one at a time) |
| Valid FP output guaranteed | No | Yes |
| Speed (avg on Thingi10k) | 360 s | 49.8 s (with parallelism) |
| Success rate on 10k models | 99.89% | 99.97% |

---

## Source Layout

| File(s) | Phase |
|---|---|
| `src/Simplification.h/.cpp` | Phase 1: preprocessing |
| `src/FloatTetDelaunay.h/.cpp` | Phase 2: background mesh generation |
| `src/TriangleInsertion.h/.cpp`, `src/CutMesh.h/.cpp` | Phase 2: triangle insertion |
| `src/MeshImprovement.h/.cpp` | Phases 3 & 4: optimization + filtering |
| `src/EdgeSplitting`, `EdgeCollapsing`, `EdgeSwapping`, `VertexSmoothing` | Phase 3: local operations |
| `src/LocalOperations.h/.cpp` | Shared geometry helpers (energy, predicates) |
| `src/AABBWrapper.h/.cpp` | Envelope queries (wraps geogram AABB) |
| `src/FloatTetwild.h/.cpp` | Library entry point |
| `src/main.cpp` | CLI entry point |
| `src/Mesh.hpp`, `src/Parameters.h`, `src/Types.hpp` | Data structures |

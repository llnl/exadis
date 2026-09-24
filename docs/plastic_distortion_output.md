# Plastic distortion tensor output

Adds three accumulated plastic-deformation output fields to the driver:
`Eptot`, `Wptot`, and `Pdist`. They are summed from the per-step increments
(`system->dEp`, `system->dWp`) in `ExaDiSApp::update_mechanics`, persisted in
restart files, and emitted as output columns when requested.

## The fields

| Field   | Meaning                          | Symmetry      | Columns |
|---------|----------------------------------|---------------|---------|
| `Eptot` | Accumulated plastic strain       | symmetric     | 9 (full 3×3) |
| `Wptot` | Accumulated plastic spin/rotation| antisymmetric | 9 (full 3×3) |
| `Pdist` | Full plastic distortion `βp`     | general       | 9 (full 3×3) |

Relationship:

```
Pdist = Eptot + Wptot          (βp = εp + ωp)
Eptot = ½ (βp + βpᵀ)           symmetric part  → plastic strain
Wptot = ½ (βp − βpᵀ)           antisymmetric part → plastic rotation
```

## Column layout

Each field is written as a full 3×3 tensor in **row-major** order:

```
xx xy xz  yx yy yz  zx zy zz
```

So the header for `Eptot` is `Epxx Epxy Epxz Epyx Epyy Epyz Epzx Epzy Epzz`
(`Wp…` for `Wptot`, `Bp…` for `Pdist`).

## Enabling the output

Add the property name to the driver's output property list (case-insensitive):

| Property name        | Field   |
|----------------------|---------|
| `Eptot` / `eptot`    | `Eptot` |
| `Wptot` / `wptot`    | `Wptot` |
| `Pdist` / `pdist`    | `Pdist` |

## Interpreting values

- **`Eptot` (plastic strain).** Symmetric. Diagonal = normal plastic strains
  along x/y/z; off-diagonal = engineering shear (½γ) components. For a single
  active slip system you'll see the resolved shear dominate one off-diagonal
  pair. Trace ≈ 0 — dislocation glide is volume-preserving (plastic
  incompressibility).
- **`Wptot` (plastic spin).** Antisymmetric (zero diagonal, `ij = −ji`). It is
  the lattice rotation accumulated by plastic flow — relevant for texture
  evolution and for the stress counter-rotation applied when
  `ctrl.rotation` is on.
- **`Pdist` (`βp`).** The full displacement-gradient contribution from
  plasticity. Use it when you need strain and rotation together (e.g. comparing
  against a continuum `βp`), rather than reconstructing it from the two halves.

## Restart

`Eptot` and `Wptot` are written to and read back from the restart file
(`write_restart` / `read_restart`), so accumulation continues seamlessly across
restarts. `Pdist` is derived on the fly (`Eptot + Wptot`) and is not stored
separately.

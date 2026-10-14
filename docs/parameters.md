# Parameters

The notebooks keep every choice in one configuration cell. The web form uses the same choices, with the labels below. Downstream CCP4, Coot, and PyMOL code reads those values and does not need to be edited for a new site or a new contour.

Shared crystallographic meaning is on the [background](concepts.md) page. Each workflow page says which of these choices change that map.

## Names

| Notebook | Web form | Meaning |
| --- | --- | --- |
| `DARK_ID` | Dark name, or omit-model name | Label used in output file names |
| `LIGHT_ID` | Light name, or observed-data name | Label used in output file names |
| `PDB_ID` | Structure name | Label for the single-structure map |
| `CCP4_ROOT` | CCP4 directory | Suite prefix. Default `/Applications/ccp4-9` |
| `RESOLUTION_SOURCE` | Resolution | `"sf"` reads the MTZ range. `"explicit"` uses the two cutoffs |
| `LOWRES` | Low-resolution cutoff (Å) | Larger number. Ignored when the source is `"sf"` |
| `DMIN` | High-resolution cutoff (Å) | Smaller number. Ignored when the source is `"sf"` |
| `MAP_BORDER_A` | Map border (Å) | `extends` padding around the model. Default 5 |
| `PHASE_SOURCE` | Phase source | Fo−Fo only. Default `light_validation` |
| `DELTA_B` | Sharpening ΔB (Å²) | 2mFo−DFc only. Default −40 |
| `REFMAC_NCYC` | REFMAC cycles | Omit map only. `0` does not move atoms |
| `SITE_DISTANCE_RESTRAINTS` | Distance restraints | Omit map. One pair per line |
| `SITE_DISTANCE_SIGMA` | Distance sigma (Å) | Default 0.50 |
| `SITE_ANGLE_RESTRAINTS` | Angle restraints | Omit map. One angle per line. Sigma in degrees |
| `SITE_ANGLE_SIGMA` | Angle sigma (degrees) | Default 25 |
| `USE_LIBCHECK` | Widen distances in a libcheck library | Omit map. Default off |
| `LIB_MONOMERS` | Monomer names for libcheck | Empty means the site residue plus other hetero residues |
| `SITE_CHAIN`, `SITE_RESI` | Chain, residue | Centre of the clip and of the figure |
| `SITE_LABEL` | Label | File-name tag such as `OEC` or `heme` |
| `SITE_COPIES` | Additional copies | Further `(chain, residue)` pairs. The form always includes the main site |
| `EXTRA_RESIDUES` | Neighbouring residues | Sticks drawn beside the site |
| `EXPORT_RADIUS_A` | Clip radius (Å) | Coot sphere. Default 10 |
| `CARVE_A` | Carve radius (Å) | PyMOL skin around displayed atoms. Default 3 |
| `SIGMA_LEVELS` | Sigma levels | Contour list |
| `SIGMA_POSITIVE`, `SIGMA_NEGATIVE` | Positive levels only, negative levels only | Difference maps. Blank means the shared list |
| `CONTOUR_MODE` | Contour | `"sigma"` or `"absolute"` |
| `surface_representation` | Surface | `"mesh"` or `"solid"` |
| `SHOW_OEC_BONDS` | Draw Mn/Ca–O bonds | Skipped when the atom names are absent |
| `DRAW_DEFAULT_OEC_DISTANCES` | Draw the default distances | Mn1–O5 and Glu OE2 distances when those atoms exist |
| `PNG_WIDTH`, `PNG_HEIGHT` | Figure width, figure height | Pixels, from 200 to 4000 |
| `RAY_TRACE` | Ray-trace the figure | Falls back to an untraced render if the display cannot ray-trace |
| `RECOMPUTE_MAP` | Recompute the map | Off reuses an extended CCP4 map with the same name |
| `RECOMPUTE_FRAGMENT` | Recompute the Coot fragment | Off reuses the `.map` clip |
| `PYMOL_VIEW` | — | Notebook camera matrix. The web application uses the browser viewer |
| `INSPECT_IN_PYMOL` | — | Notebook only. Opens the PyMOL GUI after the PNG |

## Resolution

`"sf"` is the default in every notebook.

- Fo−Fo uses the common shell of the two structure-factor MTZs: the poorer low-resolution limit and the poorer high-resolution limit.
- The omit map uses the observed structure-factor MTZ, after `cif2mtz` and before `ctruncate`.
- The 2mFo−DFc map uses the structure-factor MTZ as the cutoff. The coefficients that are transformed are `FWT` / `PHWT`, and the log warns if `DMIN` is finer than the `FWT` file.

`"explicit"` requires `LOWRES > DMIN > 0`. The explicit numbers stored in the notebooks, used only in that mode, are 20.0–2.10 Å (Fo−Fo), 30.0–1.89 Å (omit map), and 30.0–2.09 Å (2mFo−DFc).

## Site

`SITE_CHAIN` is a single character, as stored in the PDB file. `SITE_RESI` is the residue sequence number, which changes between entries that contain the same chemical group.

Neighbouring residues are one per line:

```text
A 189 GLU189
```

The third field is a label used in the scene. It must start with a letter and contain only letters, digits, and underscores.

Additional copies are one `chain residue` pair per line. Coot writes one fragment per copy. The figure uses the first, which is the main site.

The residue scan lists non-water `HETATM` records and hides common lipids and pigments (`LMT`, `LHG`, `SQD`, `DGD`, `BCR`, `CLA`, `CHL`, `PHO`, `STE`, `LMG`, `PL9`, `FME`, `BCT`). A cluster stored as `ATOM` records will not appear. Set the chain and residue from the file in that case.

## Contours

Difference-map levels are positive numbers. The negative surface is drawn at the same magnitudes unless the negative list is set. The Fo−Fo notebook uses `3, 4, 5`. The omit-map notebook uses `3.5`.

The 2mFo−DFc list may contain a negative level, which is drawn red. The notebook uses `2.0, 3.5` and a solid surface. A classical unsharpened envelope is nearer `1` to `1.5`.

`"sigma"` normalises the map so the root-mean-square is 1. `"absolute"` uses raw map units. `"mesh"` is `isomesh`. `"solid"` is a transparent `isosurface`.

Colours, from the lowest magnitude to the highest:

| Map | Positive or density | Negative |
| --- | --- | --- |
| Fo−Fo and mFo−DFc | forest, cyan, blue, then density | red, salmon, magenta, then hotpink |
| 2mFo−DFc | marine, blue, deepblue, then density | red, and only if a level is negative |

## Omit-map restraints

Distance lines accept two, three, or four fields:

```text
MN1 O5
MN3 O5 2.4
MN4 O5 2.4 0.5
```

Angle lines accept three, four, or five fields. The angle and the sigma are degrees:

```text
MN3 O5 MN4
MN3 O5 CA1 90
O5 CA1 O1 90 25
```

A missing distance or angle keeps the value in the starting model. A missing sigma uses the form-wide sigma. External distances replace the dictionary bond for that atom pair in `SITE_COPIES` only. A `libcheck` library changes every copy of the residue type. Angles are external in both modes.

The configuration-cell defaults are the O5 distances to Mn1, Mn3, Mn4, and Ca, and these angles:

```text
MN3 O5 MN4
MN3 O5 CA1
MN4 O5 CA1
O5 CA1 O1
O2 MN3 O5
O3 MN3 O5
O4 MN3 O5
O4 MN4 O5
```

## Re-runs

Building the MTZ files, running REFMAC, and transforming the map are the slow steps. Clipping and rendering are cheaper.

- Leave map recomputation off to reuse an extended CCP4 file of the same name.
- Leave fragment recomputation off to reuse the Coot clip while you change colours, contours, or the camera.
- PyMOL runs on every iteration.
- Changing $\Delta B$, the resolution shell, the cycle count, or the restraints changes the map. Turn recomputation on for that run. The 2mFo−DFc file name already contains $\Delta B$ and the resolution tag, so an old sharpening is not silently reused under a new name.

## What the output directory contains

Every job writes logs beside the maps. The names follow the labels you set.

| Workflow | Extended map | Figure |
| --- | --- | --- |
| Fo−Fo | `{light}_minus_{dark}_{low}-{high}_extended.ccp4` | `{light}_minus_{dark}_{site}_{chain}{resi}_diffmap.png` |
| mFo−DFc | `{light}_DELFWT_{dark}_{low}-{high}_extended.ccp4` | `{light}_DELFWT_{dark}_{site}_{chain}{resi}_diffmap.png` |
| 2mFo−DFc | `{id}_2fofc_{btag}_{low}-{high}_extended.ccp4` | `{id}_2fofc_{btag}_{site}_{chain}{resi}.png` |

The resolution tag replaces the decimal point with `p`. The omit-map job also writes `{dark}_refined.pdb`, the external-restraint text, and, when requested, `{dark}_lib.cif`. The web application adds `webapp_manifest.json` listing the figure, the extended map, the fragment, and the PyMOL script.

# Parameters

Every choice is a field on the workflow form. CCP4, Coot, and PyMOL read those values, so a new site or a new contour does not require an edit to the programs.

Shared crystallographic meaning is on the [background](concepts.md) page. Each workflow page says which of these choices change that map.

## Names

| Form field | Meaning |
| --- | --- |
| Dark name, or omit-model name | Label used in output file names |
| Light name, or observed-data name | Label used in output file names |
| Structure name | Label for the single-structure map |
| CCP4 directory | Suite prefix. Default `/Applications/ccp4-9` |
| Resolution | Structure-factor file reads the MTZ range. Explicit uses the two cutoffs |
| Low-resolution cutoff (Å) | Larger number. Ignored when the source is the structure-factor file |
| High-resolution cutoff (Å) | Smaller number. Ignored when the source is the structure-factor file |
| FFT grid resolution (Å) | Omit map only. Blank samples at the high-resolution limit. A smaller number extends the `fft` grid with zero amplitude |
| Map border (Å) | `extends` padding around the model. Default 5 |
| Phase source | Fo−Fo only. Default is the light validation file |
| Sharpening ΔB (Å²) | 2mFo−DFc only. Default −40 |
| REFMAC cycles | Omit map only. `0` does not move atoms |
| Distance restraints | Omit map. One pair per line |
| Distance sigma (Å) | Default 0.50 |
| Angle restraints | Omit map. One angle per line. Sigma in degrees |
| Angle sigma (degrees) | Default 25 |
| Widen distances in a libcheck library | Omit map. Default off |
| Monomer names for libcheck | Empty means the site residue plus other hetero residues |
| Chain, residue | Centre of the clip and of the figure |
| Label | File-name tag such as `OEC` or `heme` |
| Additional copies | Further chain and residue pairs. The form always includes the main site |
| Neighbouring residues | Sticks drawn beside the site |
| Clip radius (Å) | Coot sphere. Default 10 |
| Carve radius (Å) | PyMOL skin around displayed atoms. Default 3 |
| PNG sigma levels | Contours drawn in the PyMOL PNG |
| Positive levels only, negative levels only | Difference maps. Blank means the PNG sigma levels |
| Viewer sigma minimum, viewer sigma maximum | Inclusive range of the one isosurface in the browser. The minimum must be lower than the maximum |
| Contour | Sigma, after the map is normalised, or absolute map units |
| Surface | Mesh or solid, for the PyMOL PNG. The viewer has its own solid/mesh switch |
| Draw Mn/Ca–O bonds | Skipped when the atom names are absent |
| Draw the default distances | Mn1–O5 and Glu OE2 distances when those atoms exist |
| Figure width, figure height | Pixels, from 200 to 4000 |
| Ray-trace the figure | Falls back to an untraced render if the display cannot ray-trace |
| Recompute the map | Off reuses an extended CCP4 map with the same name |
| Recompute the Coot fragment | Off reuses the `.map` clip |

## Resolution

The structure-factor file is the default resolution source on every form.

- Fo−Fo uses the common shell of the two structure-factor MTZs: the poorer low-resolution limit and the poorer high-resolution limit.
- The omit map uses the observed structure-factor MTZ, after `cif2mtz` and before `ctruncate`. Its FFT grid resolution, when set to a smaller number of ångströms, samples that shell on the grid of the finer limit. Reciprocal-space points beyond the coefficients contribute zero amplitude. A number coarser than the coefficient limit is rejected. The extended-map name gains `_grid` and that limit, so the finer map is not reused in place of the ordinary one.
- The 2mFo−DFc map uses the structure-factor MTZ as the cutoff. The coefficients that are transformed are `FWT` / `PHWT`, and the log warns if the high-resolution cutoff is finer than the `FWT` file.

Explicit resolution requires the low-resolution cutoff to be larger than the high-resolution cutoff, and both to be positive. The explicit numbers stored on the forms, used only in that mode, are 20.0–2.10 Å (Fo−Fo), 30.0–1.89 Å (omit map), and 30.0–2.09 Å (2mFo−DFc).

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

Difference-map PNG levels are positive numbers. The negative surface is drawn at the same magnitudes unless the negative list is set. The Fo−Fo form uses `3, 4, 5`. The omit-map form uses `3.5`.

The 2mFo−DFc list may contain a negative level, which is drawn red. The form uses `2.0, 3.5` and a solid surface. A classical unsharpened envelope is nearer `1` to `1.5`.

Sigma mode normalises the map so the root-mean-square is 1. Absolute mode uses raw map units. Mesh is `isomesh`. Solid is a transparent `isosurface`.

The viewer sigma minimum and maximum bound a different contour: the single isosurface in the browser. Difference maps default to \(-5\) and \(5\). The sharpened map defaults to \(1\) and \(5\). The starting viewer level is the PNG level closest to zero that lies inside that range. Blue/red and temperature colouring, and the viewer’s own solid/mesh switch, are chosen on the result panel after the job finishes. They are not form fields.

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

The form defaults are the O5 distances to Mn1, Mn3, Mn4, and Ca, and these angles:

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
| mFo−DFc | `{light}_DELFWT_{dark}_{low}-{high}[_grid{limit}]_extended.ccp4` | `{light}_DELFWT_{dark}_{site}_{chain}{resi}_diffmap.png` |
| 2mFo−DFc | `{id}_2fofc_{btag}_{low}-{high}_extended.ccp4` | `{id}_2fofc_{btag}_{site}_{chain}{resi}.png` |

The resolution tag replaces the decimal point with `p`. The omit-map job also writes `{dark}_refined.pdb`, the external-restraint text, and, when requested, `{dark}_lib.cif`. The web application adds `webapp_manifest.json` listing the figure, the extended map, the fragment, and the PyMOL script.

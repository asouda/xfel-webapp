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
| Resolution | Subframe on every workflow. The source is the structure-factor file or explicit cutoffs |
| Low-resolution cutoff (Å) | Shown only for explicit cutoffs. Larger number |
| High-resolution cutoff (Å) | Shown only for explicit cutoffs. Smaller number |
| Set a custom FFT grid resolution for the final map (Å) | Switch on every workflow. Off samples at the high-resolution limit. On, a smaller number extends the `fft` grid with zero amplitude |
| Extended map border (Å) | In the Resolution subframe. Padding around the model, in ångströms. Default 5 |
| Phase source | Fo−Fo only. Default is the light validation file |
| Sharpening ΔB (Å²) | 2mFo−DFc only. Default −40 |
| REFMAC refinement | Omit map and reconstruction. The section starts off, which sets zero cycles. On shows cycles, B-factor refinement, refinement type, restraints, and the custom dictionary |
| B-factor refinement | `REFI BREF`. ISOTropic (default), ANISotropic, MIXEd, or OVERall |
| Refinement type | `REFI TYPE`. RESTrained (default) or UNREstrained |
| REFMAC cycles | Shown while REFMAC refinement is on. `0` does not move atoms. The section off forces `0` |
| Number of atoms to reconstruct | Reconstruction only. How many heavy atoms to insert. Default 1 |
| Map check interval | Reconstruction only. `0` shows the viewer at the end. A positive integer pauses after that many search iterations |
| Neighborhood size | Reconstruction only, in Density peak search. Odd integer, the side of the local-maximum cube. Default 3 |
| Spline maximum | Reconstruction only, in Density peak search. Off by default. On places the atom at the spline maximum |
| Spline order | Reconstruction only. Integer from 0 through 5. Default 3. Shown only while the spline maximum is on |
| External harmonic restraints | Reconstruction only. Off by default. The sigma field is shown only while the switch is on. Holds atoms already added at the coordinates in the current omit model |
| Harmonic sigma (Å) | Reconstruction only. Default 0.1. Shown only while the harmonic switch is on, and used only when the cycle count is greater than zero |
| Distance restraints | Omit map and reconstruction. The lines are shown only while the switch is on. One pair per line. A distance and sigma on the line replace the fields for that pair. The omit-map switch starts on. The reconstruction switch starts off |
| Distance sigma (Å) | Default 0.50. Used when a distance line does not give its own sigma |
| Angle restraints | Omit map and reconstruction. The lines are shown only while the switch is on. One angle per line. An angle and sigma on the line, in degrees, replace the fields for that angle. The omit-map switch starts on. The reconstruction switch starts off |
| Angle sigma (degrees) | Default 25. Used when an angle line does not give its own sigma |
| Dihedral restraints | Omit map and reconstruction. The lines are shown only while the switch is on. One torsion per line, four atoms. An angle and sigma on the line, in degrees, replace the fields for that torsion. The switch starts off |
| Dihedral sigma (degrees) | Default 10. Used when a dihedral line does not give its own sigma |
| Custom dictionary | Omit map and reconstruction, inside REFMAC refinement. Off by default. The file fields are shown only while that switch is on, and only while the REFMAC section is on. Read only when the cycle count is greater than zero, and read before the external restraints |
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

The structure-factor file is the default resolution source on every form. The two cutoff fields are hidden until explicit cutoffs are chosen. The FFT-grid number is hidden until its switch is on. Extended map border stays visible in the same section.

- Fo−Fo uses the common shell of the two structure-factor MTZs: the poorer low-resolution limit and the poorer high-resolution limit.
- The omit map and the reconstruction use the observed reflection MTZ. An mmCIF is converted with `cif2mtz` first; an MTZ supplied on the form is that file. Resolution is read before `ctruncate`. The reconstruction repeats that shell on every iteration.
- The 2mFo−DFc map uses the structure-factor MTZ as the cutoff. The coefficients that are transformed are `FWT` / `PHWT`, and the log warns if the high-resolution cutoff is finer than the `FWT` file.
- Every workflow can set a custom FFT grid for the final map. The switch starts off, and the grid then matches the high-resolution limit of the coefficients. A smaller number samples that shell on the grid of the finer limit. Reciprocal-space points beyond the coefficients contribute zero amplitude. A number coarser than the coefficient limit is rejected. The extended-map name gains `_grid` and that limit, so the finer map is not reused in place of the ordinary one.

Explicit resolution requires the low-resolution cutoff to be larger than the high-resolution cutoff, and both to be positive. The explicit numbers stored on the forms, used only in that mode, are 20.0–2.10 Å (Fo−Fo), 30.0–1.89 Å (omit map and reconstruction), and 30.0–2.09 Å (2mFo−DFc).

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

Difference-map PNG levels are positive numbers. The negative surface is drawn at the same magnitudes unless the negative list is set. The Fo−Fo form uses `3, 4, 5`. The omit-map and reconstruction forms use `3.5`.

The 2mFo−DFc list may contain a negative level, which is drawn red. The form uses `2.0, 3.5` and a solid surface. A classical unsharpened envelope is nearer `1` to `1.5`.

Sigma mode normalises the map so the root-mean-square is 1. Absolute mode uses raw map units. Mesh is `isomesh`. Solid is a transparent `isosurface`.

The viewer sigma minimum and maximum bound a different contour: the single isosurface in the browser. Difference maps default to \(-5\) and \(5\). The sharpened map defaults to \(1\) and \(5\). The starting viewer level is the PNG level closest to zero that lies inside that range. Blue/red and temperature colouring, and the viewer’s own solid/mesh switch, are chosen on the result panel after the job finishes. They are not form fields.

Colours, from the lowest magnitude to the highest:

| Map | Positive or density | Negative |
| --- | --- | --- |
| Fo−Fo and mFo−DFc | forest, cyan, blue, then density | red, salmon, magenta, then hotpink |
| 2mFo−DFc | marine, blue, deepblue, then density | red, and only if a level is negative |

## Reconstruction search

The neighborhood size must be a positive odd integer. The spline maximum is off by default, and the spline order field is hidden until that switch is on. The order must be from 0 through 5. Leaving the switch off places each atom at the voxel centre. The map-check interval must be zero or positive. The number of atoms must be at least 1.

Harmonic restraints sit inside REFMAC refinement. The section starts off, which hides that switch and forces zero cycles. The harmonic switch itself is off by default, and its sigma field is hidden until the switch is on. A sigma of 0.1 Å is the default when the field is left blank. The restraint is still skipped when the cycle count is 0.

The 1 Å warning and the 3 Å stop are not parameters. Hydrogens, deuterium, and atoms with occupancy 0 or less are never candidates.

The reconstruction form does not show the switches for Mn/Ca–O bonds, ray tracing, map recomputation, or fragment recomputation. Every iteration recomputes the map and the clip, ray-traces the figure, and draws bonds from covalent radii.

## Omit-map restraints

These fields are inside REFMAC refinement, on the omit map and the reconstruction. The section starts off. While it is off the fields are hidden, nothing from the section is submitted, and the cycle count is 0. Turn the section on to set cycles, `REFI BREF`, `REFI TYPE`, the restraints, and a custom dictionary. `REFI BREF` is ISOTropic, ANISotropic, MIXEd, or OVERall, and the default is ISOTropic. `REFI TYPE` is RESTrained or UNREstrained, and the default is RESTrained. The Fo−Fo zero-cycle job is not this section. It stays `REFI TYPE REST` and `BREF ANIS`.

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

Dihedral lines accept four, five, or six fields. The angle and the sigma are degrees, from −180 to 180:

```text
MN1 O5 MN4 CA1
MN1 O5 MN4 CA1 180
MN1 O5 MN4 CA1 180 10
```

The torsion is the angle between the plane of the first three atoms and the plane of the last three. A missing distance, angle, or dihedral keeps the value in the starting model. A missing sigma uses the form-wide sigma for that restraint. A distance and sigma on the line, such as `MN1 O5 2.40 0.50`, replace the form-wide values for that pair. The same applies to an angle line and a dihedral line. Each restraint group is applied only while its switch is on. External restraints replace the CCP4 monomer library for those atoms in `SITE_COPIES` only.

A custom dictionary is a `LIBIN` line. The file fields are shown only while its switch is on, and turning the switch off uses no dictionary. The file is passed only when the switch is on and the cycle count is greater than zero. Zero cycles keeps the CCP4 monomer library. The external restraints are read after the dictionary, so they replace matching distances, angles, and dihedrals in that file.

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
- PyMOL runs for Fo−Fo, the omit map, and the sharpened map. In a reconstruction it writes a script every iteration and runs that script only when the iteration pauses for a map check.
- Changing $\Delta B$, the resolution shell, the FFT grid, the cycle count, or the restraints changes the map. Turn recomputation on for that run. The 2mFo−DFc file name contains $\Delta B$. A finer FFT grid adds `_grid` and that limit on every workflow, so the finer map is a different file.
- The reconstruction recomputes the map and the fragment on every iteration. Each iteration writes a new `r-iteration-NN` directory, so an earlier iteration is not overwritten.

## What the output directory contains

Every job writes logs beside the maps. The names follow the labels you set.

| Workflow | Extended map | Figure |
| --- | --- | --- |
| Fo−Fo | `{light}_minus_{dark}[_grid{limit}]_extended.ccp4` | `{light}_minus_{dark}_{site}_{chain}{resi}_diffmap.png` |
| mFo−DFc | `{light}_DELFWT_{dark}[_grid{limit}]_extended.ccp4` | `{light}_DELFWT_{dark}_{site}_{chain}{resi}_diffmap.png` |
| Reconstruction | The omit-map name, inside `r-iteration-01`, `r-iteration-02`, and so on | The omit-map figure name, in the same iteration directory |
| 2mFo−DFc | `{id}_2fofc_{btag}[_grid{limit}]_extended.ccp4` | `{id}_2fofc_{btag}_{site}_{chain}{resi}.png` |

`_grid{limit}` is present only when the FFT grid is finer than the coefficients. The limit replaces the decimal point with `p`, as in `_grid1p20`. The omit-map job also writes `{dark}_refined.pdb` and the external-restraint text. A reconstruction iteration writes those same omit-map files, plus `reconstruction.log` and, when a peak was found, `density_sections_stepNN.txt`. When an atom is added it also writes `{dark}_refined_added_{atom}_stepNN.pdb`, the refined coordinates with that atom at the peak, and the next iteration starts from that file. PDB, CCP4, map, PNG, and that density-section file in the iteration directory have `_stepNN` before the extension. The output directory also receives `{dark}_working.pdb` and, at the end, `{dark}_reconstructed_final.pdb`, a copy of the last iteration’s refined coordinates. The last iteration inserts no atom. The web application adds `webapp_manifest.json` listing the figure, the extended map, the fragment, and the PyMOL script.

# Sharpened 2mFo−DFc map

This workflow builds the weighted $2F_o-F_c$ electron-density map of one structure, sharpens it by a $B$-factor shift $\Delta B$, clips a fragment around a site, and draws $\sigma$ isosurfaces in PyMOL with parts of the model superimposed.

The notebook is `notebooks/XFEL_2fo-fc_map.ipynb`. The web application runs the same steps from the 2Fo−Fc page.

The map is the electron density of that structure. Mesh that covers an atom is the expected signal. Difference peaks belong to the [Fo−Fo](fo-fo.md) and [omit-map](omit-map.md) workflows.

## Worked example

The configuration cell, and the web form, use PDB `9z35`: the room-temperature XFEL structure of the two-flash, $S_3$-rich state of photosystem II, determined at 1.93 Å from anomalous data at 9.5 keV.

| File | Name in the configuration |
| --- | --- |
| Coordinates | `9z35.pdb`, from `data/9z35/` |
| Structure factors | `9z35-sf.cif` |
| Validation $2F_o-F_c$ coefficients | `9z35_validation_2fo-fc_map_coef.cif` |

The site is chain `A` residue **601**, ligand `OEY`, the major oxygen-evolving complex in this entry, drawn with D1-Glu189 (`A 189`). The notebook introduction also discusses the Ibrahim *et al.* 2F structure `6w1v`, where the corresponding cluster is `OEY` A416 and the validation coefficients stop near 2.09 Å. Residue numbers do not transfer between those entries. Use the residue list for the file you actually load.

The figure is sharpened with $\Delta B = -40$ Å² and drawn as solid surfaces at $2.0\sigma$ and $3.5\sigma$ of that sharpened map.

## Background

The wwPDB validation coefficient file stores the EDS weighted map

$$
\rho_{2mF_o-DF_c}(\mathbf{r}) \;\propto\; \sum_{hkl} F_{\mathrm{WT}}\, e^{i\varphi_{\mathrm{WT}}}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}},
$$

where $F_{\mathrm{WT}}$ (`FWT`) is already $2m|F_o|-D|F_c|$ and $\varphi_{\mathrm{WT}}$ (`PHWT`) is the model phase. The companion `FOM` column is the figure of merit that went into that weighting. It is left unused here, so the weighting is applied once.

### Sharpening

Before the transform, CCP4 `fft` rescales each amplitude:

$$
F_{\mathrm{used}} = F_{\mathrm{WT}}\,\exp\!\big(-\Delta B\,(\sin\theta/\lambda)^2\big).
$$

$\Delta B$ is the temperature factor in `SCALE F1`. The default $\Delta B = -40$ Å² makes the exponent positive at higher resolution, so those terms are boosted and the map looks sharper. The same factor is $\exp\!\big(+B_{\mathrm{sharpen}}(\sin\theta/\lambda)^2\big)$ with $B_{\mathrm{sharpen}} = -\Delta B$, which is 40 Å² for the default. $\Delta B = 0$ leaves `FWT` unchanged. A more negative $\Delta B$ sharpens more.

Because $(\sin\theta/\lambda)^2 = 1/(4d^2)$, the multiplier at the high-resolution limit $d_{\min}$ is $\exp(-\Delta B / 4d_{\min}^2)$. The map step prints that number. σ levels are measured on the sharpened map, so a $1.5\sigma$ surface at $-40$ Å² is a different surface from $1.5\sigma$ with no sharpening.

### How to read the figure

- Mesh sitting on atoms is accounted-for electron density.
- A gap over an atom, or a blob with no atom, is the usual $2mF_o-DF_c$ hint that the model and the data disagree. The difference map that would show that disagreement directly is the [omit-model workflow](omit-map.md).
- $2F_o-F_c$ maps are normally viewed near $1\sigma$. Levels of $3$–$5\sigma$ are the usual range for difference maps and will often leave an unsharpened $2mF_o-DF_c$ mesh empty. The worked example contours the sharpened map at 2.0 and 3.5σ, which is a stronger cutoff of a map whose high-resolution terms have been boosted.

## Files

Three files. Paths left unset are `{pdbid}.pdb`, `{pdbid}-sf.cif`, and `{pdbid}_validation_2fo-fc_map_coef.cif` inside the input directory. For `9z35` that directory is `data/9z35/`.

| Role | Notebook variable | Web form |
| --- | --- | --- |
| Coordinates | `PDB_XYZ` | Coordinates. Drawn on the map, and used as the `extends` mask. |
| Structure factors | `PDB_SF` | Structure factors. With the resolution source set to the structure-factor file, `LOWRES` and `DMIN` are the range of the MTZ made from this file. |
| Validation coefficients | `PDB_VALIDATION` | Validation 2Fo−Fc coefficients. `FWT` and `PHWT` are the map. |

`PDB_ID` is a short label in output names. It does not have to be a PDB code. The structure-factor MTZ is not scaled into the map. Its resolution range is what the `"sf"` resolution source uses as the Fourier cutoff.

## Steps

The programs run in this order: `cif2mtz`, `fft` with `SCALE F1` set to $\Delta B$, `extends`, Coot (as a regular map), and PyMOL.

The output map name includes both $\Delta B$ and the resolution tag, so a new sharpening or a new shell is a new file. `RECOMPUTE_MAP` must still be on the first time you ask for that file. An older extended map is reused only when the name matches and recomputation is off.

### 1. Check the programs and the files

The programs are `cif2mtz`, `fft`, `extends`, and `mtzdump`. Coot comes from the CCP4 tree. PyMOL must be on `PATH`. The three input files must already exist.

### 2. Inspect the coordinates and choose the site

The residue list comes from the coordinate file. For photosystem II, look for `OEX` or `OEY`. The workflow checks the chosen residue, prints its atom names, and prints the centroid used for the clip. Put neighbouring protein residues in `EXTRA_RESIDUES`; the default is D1-Glu189.

On the web form, **List hetero residues** reads these coordinates.

### 3. Build the sharpened map

**3a. Convert both reflection files.** `cif2mtz` writes `{id}-sf.mtz` and, when the map is rebuilt, `{id}_validation_2fo-fc_map_coef.mtz`. The structure-factor conversion runs even when an extended map is reused, so the resolution range can still be reported.

**3b. Read the shell and the coefficient labels.** `mtzdump` supplies the resolution range. With `RESOLUTION_SOURCE = "sf"` that range becomes `LOWRES` and `DMIN`. The notebook’s explicit fallback is 30.0–2.09 Å. The validation MTZ is checked for `FWT` and `PHWT`, and its own resolution limits are printed.

If the high-resolution cutoff is finer than the `FWT` limit, the extra shells are not in the coefficient file. The log warns about that. `DMIN` should be no finer than the `FWT` limit. For `6w1v` the notebook records that limit at about 2.09 Å. For another entry, use the limit the map step prints.

**3c. Sharpen and transform.** The amplitude multiplier at `DMIN` is printed, then `fft` writes `{id}_2fofc_{btag}_{tag}.ccp4`:

```text
TITLE {id} weighted 2mFo-DFc sharpened by dB=<ΔB>
RESOLUTION <low> <high>
LABIN F1=FWT PHI=PHWT
SCALE F1 1.0 <ΔB>
```

`SCALE F1 1.0 ΔB` is the factor in the sharpening equation: scale 1, temperature factor $\Delta B$.

**3d. Cover the model.** `extends` pads the map by `MAP_BORDER_A` (default 5 Å) around the coordinates and writes the extended CCP4 file.

The log states that the blue mesh is sharpened 2mFo−DFc density for this entry, at the $\Delta B$ you set.

### 4. Export the site fragment

Coot loads the extended map as a regular map. The difference-map flag is 0. It centres on each residue in `SITE_COPIES` and writes a clip of radius `EXPORT_RADIUS_A` (default 10 Å). The figure is drawn for the first copy. The Coot contour does not change the exported voxels.

### 5. Contour and render

PyMOL loads the coordinates and the fragment and builds one isosurface per entry in `SIGMA_LEVELS`. The configuration cell uses `[2.0, 3.5]` and `surface_representation = "solid"`. The site residue and every extra residue are shown as sticks and, for the site, spheres. The rest of the model stays hidden.

| Lowest $\sigma$ | Middle | Highest $\sigma$ |
| --- | --- | --- |
| marine | blue | deepblue |

A negative number in the level list is drawn red. The default list is positive density only. A fourth positive level uses the colour `density`.

`CARVE_A` (default 3 Å) keeps each surface near the displayed atoms. Mn/Ca–oxygen bonds, without an O6 pair in this workflow’s bond list, are added when those atom names exist.

Outputs:

| File | Contents |
| --- | --- |
| `{id}_2fofc_{btag}_{site}_{chain}{resi}.png` | The figure |
| `pymol_{id}_2fofc_{btag}_{site}.pml` | The scene |
| `{id}_2fofc_{btag}_{tag}_extended.ccp4` | The map |

`btag` encodes $\Delta B$, so `-40` Å² and `0` Å² never overwrite each other.

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| `DELTA_B` | Sharpening. `-40` boosts high resolution. `0` leaves `FWT` as deposited. More negative sharpens more |
| Resolution source, `LOWRES`, `DMIN` | Fourier shell. The high-resolution cutoff should respect the `FWT` limit printed in the log |
| `SIGMA_LEVELS` | Contours of the sharpened map. `[1.2]` is a single, more classical envelope. The worked example uses `2` and `3.5` |
| Surface | `"mesh"` is a wireframe. The worked example uses `"solid"` |
| `RECOMPUTE_MAP` | Turn on after changing $\Delta B$ or the resolution cutoffs, even though the file name also changes |

## Reading the default figure

For `9z35` the mesh is the validation $2mF_o-DF_c$ map sharpened by $\Delta B = -40$ Å², clipped around `OEY` A601, with D1-Glu189 from the extra-residue list. Density that wraps those atoms is the $2F_o-F_c$ envelope of this structure.

If the mesh is empty, lower the contour (try `1.0`). If the cage swallows the model, decrease `CARVE_A` or raise the cutoff.

## Adapting it

1. Provide `{pdbid}.pdb`, `{pdbid}-sf.cif`, and `{pdbid}_validation_2fo-fc_map_coef.cif`, or point the three path variables at your files.
2. Set `PDB_ID` and `DELTA_B`. Leave the resolution source on `"sf"`, or set `LOWRES` and `DMIN` no finer than the `FWT` limit.
3. Run the residue list, then set the chain, residue, and neighbours.
4. Clear any saved camera unless you intend to reuse it.
5. Recompute the map the first time, and again whenever $\Delta B$ or the resolution cutoffs change. Leave recomputation off while you tune σ and the view.

The extended map can be opened with the coordinates:

```text
coot {extended_map} {pdb}
```

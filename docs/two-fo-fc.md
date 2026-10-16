# Sharpened 2mFo−DFc map

This workflow builds the weighted \(2mF_o-DF_c\) electron-density map of one structure from its deposited validation coefficients, sharpens it by a \(B\)-factor shift \(\Delta B\), clips a fragment around a site, and draws \(\sigma\) isosurfaces with parts of the model superimposed.

The map is the electron density of that structure. Mesh that covers an atom is the expected signal. Difference peaks belong to the [Fo−Fo](fo-fo.md) and [omit-map](omit-map.md) workflows. The [theoretical background](theory.md) explains why `FWT` is already \(2m|F_o|-D|F_c|\) and what a validation coefficient file stores. The steps below are the programs that turn those coefficients into the worked-example figure.

## Worked example

The form is set up on PDB `9z35`, the room-temperature XFEL structure of the two-flash, \(S_3\)-rich state of photosystem II, determined at 1.93 Å from anomalous data at 9.5 keV. The repository holds the coordinate file, the structure-factor mmCIF, and both validation coefficient files under `data/9z35/`.

| Role | File | Download |
| --- | --- | --- |
| Coordinates | `9z35.pdb`. Drawn on the map, and used as the mask for `extends`. | [9z35.pdb](https://github.com/asouda/xfel-webapp/raw/main/data/9z35/9z35.pdb) |
| Structure factors | `9z35-sf.cif`. Anomalous intensities. Used here for the resolution range. | [9z35-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/9z35/9z35-sf.cif.gz) |
| Validation \(2F_o-F_c\) coefficients | `9z35_validation_2fo-fc_map_coef.cif`. `FWT` and `PHWT` are the map. | [9z35 2Fo−Fc coefficients](https://github.com/asouda/xfel-webapp/raw/main/data/9z35/9z35_validation_2fo-fc_map_coef.cif.gz) |
| Validation \(F_o-F_c\) coefficients | Present in the same directory. This workflow does not read it. | [9z35 Fo−Fc coefficients](https://github.com/asouda/xfel-webapp/raw/main/data/9z35/9z35_validation_fo-fc_map_coef.cif.gz) |

The mmCIF files are gzip-compressed. Point the form at the `.gz` path; a name that ends in `.gz` is decompressed with `gunzip` before `cif2mtz`.

The site is chain `A` residue **601**, ligand `OEY`, the major oxygen-evolving complex in this entry, drawn with D1-Glu189 (`A 189`). Residue numbers do not transfer from `6w1v`, where the corresponding cluster is `OEY` A416. Use **List hetero residues** on the file you actually load.

The figure is sharpened with \(\Delta B = -40\) Å² and drawn as solid surfaces at \(2.0\sigma\) and \(3.5\sigma\) of that sharpened map. The viewer slider for this workflow defaults to the range \(1\) to \(5\sigma\).

### How the 9z35 map is calculated

The programs run in this order: `cif2mtz`, `mtzdump`, `fft` with `SCALE F1` set to \(\Delta B\), `extends`, Coot, and PyMOL. There is no REFMAC job and no scaling of a second dataset.

**1. Check the programs and the files.** The programs are `cif2mtz`, `fft`, `extends`, and `mtzdump`. Coot comes from the CCP4 tree. PyMOL must be on `PATH`. The three input files must already exist.

**2. Inspect the coordinates and choose the site.** **List hetero residues** reads `9z35.pdb`. Confirm `OEY` A601, and keep D1-Glu189 in the neighbouring-residue list. The job checks the residue, prints its atom names, and prints the centroid used for the clip.

**3. Convert both reflection files.** `cif2mtz` writes `9z35-sf.mtz` from the structure-factor mmCIF and `9z35_validation_2fo-fc_map_coef.mtz` from the coefficient mmCIF. The structure-factor conversion is repeated even when an extended map is reused, so the resolution range can still be reported. The coefficient MTZ contains `FWT`, `PHWT`, and `FOM`. `FWT` is already the weighted amplitude \(2m|F_o|-D|F_c|\). `PHWT` is the model phase. `FOM` is the figure of merit that was used to build `FWT`; the transform leaves it unread.

**4. Read the shell and the coefficient labels.** `mtzdump` supplies the resolution range of `9z35-sf.mtz`. With the resolution source set to that file, the range becomes the low- and high-resolution cutoffs. The explicit pair stored for this workflow is 30.0–2.09 Å. The validation MTZ is checked for `FWT` and `PHWT`, and its own resolution limits are printed.

If the high-resolution cutoff is finer than the `FWT` limit, the extra shells are not in the coefficient file and the log says so. The cutoff should be no finer than the `FWT` limit the map step prints. For `9z35` the structure was determined at 1.93 Å; follow the limit printed for the coefficient file you loaded.

**5. Sharpen and transform.** Before the sum, `fft` rescales each amplitude:

\[
F_{\mathrm{used}} = F_{\mathrm{WT}}\, \exp\!\big(-\Delta B\, (\sin\theta/\lambda)^2\big).
\]

\(\Delta B\) is the temperature factor in `SCALE F1`. The default \(\Delta B = -40\) Å² makes the exponent positive at higher resolution, so those terms are boosted and the map looks sharper. The same factor is \(\exp\!\big(+B_{\mathrm{sharpen}}(\sin\theta/\lambda)^2\big)\) with \(B_{\mathrm{sharpen}} = -\Delta B\), which is 40 Å² for the default. \(\Delta B = 0\) leaves `FWT` unchanged. A more negative \(\Delta B\) sharpens more.

Because \((\sin\theta/\lambda)^2 = 1/(4d^2)\), the multiplier at the high-resolution limit \(d_{\min}\) is \(\exp(-\Delta B / 4 d_{\min}^2)\). The log prints that number. Sigma levels are measured on the sharpened map, so a \(1.5\sigma\) surface at \(-40\) Å² is a different surface from \(1.5\sigma\) with no sharpening.

`fft` then evaluates

\[
\rho(\mathbf{r}) \propto \sum_{hkl} F_{\mathrm{used}}(hkl)\, e^{i\, \mathrm{PHWT}(hkl)}\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}}
\]

and writes `9z35_2fofc_{btag}_{tag}.ccp4`:

```text
TITLE 9z35 weighted 2mFo-DFc sharpened by dB=<ΔB>
RESOLUTION <low> <high>
LABIN F1=FWT PHI=PHWT
SCALE F1 1.0 <ΔB>
```

`SCALE F1 1.0 ΔB` is scale 1 and temperature factor \(\Delta B\).

**6. Cover the model.** `extends` pads the map by the map border (default 5 Å) around `9z35.pdb` and writes the extended CCP4 file.

**7. Export the site fragment.** Coot loads the extended map as a regular map. The difference-map flag is 0. It centres on `OEY` A601 and writes a clip of the chosen radius (default 10 Å). The Coot contour does not change the exported voxels.

**8. Contour and render.** PyMOL loads `9z35.pdb` and the fragment and builds one isosurface per PNG sigma level. The worked example uses \(2.0\) and \(3.5\), as solid surfaces. The site residue is shown as sticks and spheres. Each neighbouring residue is shown as sticks. The rest of the model stays hidden.

| Lowest \(\sigma\) | Middle | Highest \(\sigma\) |
| --- | --- | --- |
| marine | blue | deepblue |

A negative number in the level list is drawn red. The default list is positive density only. A fourth positive level uses the colour `density`. The carve radius (default 3 Å) keeps each surface near the displayed atoms.

The browser viewer draws one isosurface of the same fragment, between the viewer sigma minimum and maximum. For this workflow those defaults are \(1\) and \(5\). Blue/red colouring paints positive density blue. Temperature colouring runs from cold at the low end of the slider to hot at the high end. The PNG remains the PyMOL figure at \(2.0\sigma\) and \(3.5\sigma\).

Outputs:

| File | Contents |
| --- | --- |
| `9z35_2fofc_{btag}_{site}_{chain}{resi}.png` | The PyMOL figure |
| `pymol_9z35_2fofc_{btag}_{site}.pml` | The scene |
| `9z35_2fofc_{btag}_{tag}_extended.ccp4` | The map |

`btag` encodes \(\Delta B\), so \(-40\) Å² and \(0\) Å² never overwrite each other. The output name also includes the resolution tag. `RECOMPUTE_MAP` must still be on the first time you ask for that file.

## Background

The wwPDB validation coefficient file stores the EDS weighted map

\[
\rho_{2mF_o-DF_c}(\mathbf{r}) \propto \sum_{hkl} F_{\mathrm{WT}}\, e^{i\varphi_{\mathrm{WT}}}\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}},
\]

where \(F_{\mathrm{WT}}\) (`FWT`) is already \(2m|F_o|-D|F_c|\) and \(\varphi_{\mathrm{WT}}\) (`PHWT`) is the model phase. The columns in the `9z35` file are `_refln.pdbx_FWT`, `_refln.pdbx_PHWT`, and `_refln.fom`.

### How to read the figure

- Mesh sitting on atoms is accounted-for electron density.
- A gap over an atom, or a blob with no atom, is the usual \(2mF_o-DF_c\) hint that the model and the data disagree. The difference map that shows that disagreement directly is the [omit-model workflow](omit-map.md).
- \(2F_o-F_c\) maps are normally viewed near \(1\sigma\). Levels of \(3\)–\(5\sigma\) are the usual range for difference maps and will often leave an unsharpened \(2mF_o-DF_c\) mesh empty. The worked example contours the sharpened map at 2.0 and 3.5σ, a stronger cutoff of a map whose high-resolution terms have been boosted.

## Files

Three files. The structure-factor MTZ is not scaled into the map. Its resolution range is what the structure-factor resolution source uses as the Fourier cutoff.

| Role | Form field |
| --- | --- |
| Coordinates | Coordinates |
| Structure factors | Structure factors |
| Validation coefficients | Validation 2Fo−Fc coefficients |

The structure name is a short label in output file names. It does not have to be a PDB code.

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| Sharpening \(\Delta B\) | \(-40\) boosts high resolution. \(0\) leaves `FWT` as deposited. More negative sharpens more |
| Resolution source and the two cutoffs | Fourier shell. The high-resolution cutoff should respect the `FWT` limit printed in the log |
| PNG sigma levels | Contours of the sharpened map. `1.2` is a single envelope nearer the classical \(1\sigma\) view. The worked example uses \(2\) and \(3.5\) |
| Viewer sigma minimum and maximum | Range of the one browser isosurface. Defaults are \(1\) and \(5\) for this workflow |
| Surface | `mesh` is a wireframe. The worked example uses `solid` |
| Recompute the map | Turn on after changing \(\Delta B\) or the resolution cutoffs |

## Reading the figure

For `9z35` the mesh is the validation \(2mF_o-DF_c\) map sharpened by \(\Delta B = -40\) Å², clipped around `OEY` A601, with D1-Glu189 from the neighbouring-residue list. Density that wraps those atoms is the \(2F_o-F_c\) envelope of this structure.

If the mesh is empty, lower the contour (try `1.0`). If the cage swallows the model, decrease the carve radius or raise the cutoff.

## Adapting it

1. Provide the coordinate file, the structure-factor mmCIF, and the validation \(2F_o-F_c\) coefficient file.
2. Set the structure name and \(\Delta B\). Leave the resolution source on the structure-factor file, or set the cutoffs no finer than the `FWT` limit.
3. List hetero residues, then set the chain, residue, and neighbours.
4. Recompute the map the first time, and again whenever \(\Delta B\) or the resolution cutoffs change. Leave recomputation off while you tune σ and the view.

The extended map can be opened with the coordinates:

```text
coot {extended_map} {coordinates}
```

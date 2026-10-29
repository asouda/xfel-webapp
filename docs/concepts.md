# Shared background

The three workflows differ in which amplitudes they transform. They share the crystallographic objects those amplitudes come from, the files you supply, and the clipping and drawing steps at the end.

How a structure factor is measured, how an amplitude and a model phase are obtained, and what a validation coefficient file stores are on the [theoretical background](theory.md) page.

## A map is a Fourier sum

An electron-density map is the Fourier transform of complex structure factors. Each reflection $hkl$ contributes an amplitude and a phase:

$$
\rho(\mathbf{r}) \;\propto\; \sum_{hkl} |F(hkl)|\, e^{i\varphi(hkl)}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}}.
$$

The amplitude $|F|$ is measured, or calculated from a model. The phase $\varphi$ is not measured in a typical MX experiment. Every map on this site borrows its phases from a model, either directly from a zero-cycle or refined REFMAC job, or from the phase column already stored in a wwPDB validation coefficient file.

That borrowing is why the maps are interpretations as well as transforms. A difference peak is only as meaningful as the phases and the scale that went into it.

## Three coefficient choices

**Fo−Fo.** Two measured amplitude sets are subtracted and given one reference phase $\varphi_{\mathrm{ref}}$:

$$
\Delta\rho(\mathbf{r}) \;\propto\; \sum_{hkl} \big(|F_{\mathrm{LIGHT}}|-|F_{\mathrm{DARK}}|\big)\, e^{i\varphi_{\mathrm{ref}}}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}}.
$$

This is the difference-Fourier approximation. The true phase of the difference vector $F_{\mathrm{LIGHT}}-F_{\mathrm{DARK}}$ is unknown. The approximation holds when the crystals are isomorphous and the structural change is local, so the reference phases still point roughly the right way.

**mFo−DFc.** One observed amplitude $|F_o|$ is compared with the amplitude $|F_c|$ calculated from a model, using the model phase $\varphi_c$. REFMAC weights the difference by the figure of merit $m$ and by $D$, the factor that damps the calculated amplitude according to coordinate error:

$$
\Delta\rho(\mathbf{r}) \;\propto\; \sum_{hkl} \big(m|F_o|-D|F_c|\big)\, e^{i\varphi_c}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}}.
$$

REFMAC writes that complex coefficient as `DELFWT` / `PHDELWT`. The workflow Fourier-transforms those columns. It does not recompute $m|F_o|-D|F_c|$ by hand.

**2mFo−DFc.** The same weighting, arranged so that density covers the model when the model and the data agree:

$$
\rho(\mathbf{r}) \;\propto\; \sum_{hkl} \big(2m|F_o|-D|F_c|\big)\, e^{i\varphi_c}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}}.
$$

The wwPDB validation file already stores this coefficient as `FWT` / `PHWT`. The companion `FOM` column is the figure of merit that went into the weighting. The workflow leaves `FOM` unused, so the figure of merit is applied once.

## Isomorphism, occupancy, and noise

A pump–probe or time-resolved XFEL experiment often produces two datasets from closely related crystals: a dark or resting state and a light or activated state. They can be subtracted only when they are isomorphous — the same space group and a very similar cell. Scale, overall $B$ factor, and completeness still differ between collections. The Fo−Fo workflow removes that global difference with `scaleit` before the Fourier transform. An overall SCALEIT $R$ factor around 0.10–0.15 is typical of an isomorphous XFEL pair that still yields a usable map. A much larger $R$ usually means the crystals are not isomorphous, or that one dataset is poor.

Activated-state occupancy in these experiments is often only 10–30 %. The difference signal is a fraction of a fully occupied atom, and the maps are noisy. Contouring at $3\sigma$ shows weaker features and more noise. Contouring at $4$–$5\sigma$ keeps the strongest peaks. Report the level you used and look at more than one level. There is no single correct contour.

## What a σ level means here

`CONTOUR_MODE = "sigma"` asks PyMOL to normalise the CCP4 map so that its root-mean-square value is 1. A contour of `3` is then $3\sigma$ of that map. `"absolute"` turns the normalisation off and treats the numbers as raw map units.

The σ you contour is the σ of the map that was written. For the sharpened 2mFo−DFc workflow, that is the sharpened map. Changing $\Delta B$ changes the features and the root-mean-square, so $1.5\sigma$ after sharpening is a different surface from $1.5\sigma$ with $\Delta B = 0$.

The Coot σ slider is a display setting inside Coot. The voxels written to the fragment are unchanged. Contouring happens later, in PyMOL.

## Files from the PDB

Each workflow starts from files you put on the form. The worked examples use files distributed by the [RCSB PDB](https://www.rcsb.org), and the application does not download them for you. Copies of the example files are in this repository; the reflection mmCIF files there are gzip-compressed. A path or an upload whose name ends in `.gz` is decompressed with `gunzip` before the CCP4 programs read it. PDB and mmCIF are used as stored.

| File | Typical name | Role |
| --- | --- | --- |
| Coordinates | `{pdbid}.pdb` or mmCIF | The model. REFMAC, `extends`, Coot, and PyMOL all read it. PDB and mmCIF coordinates are both accepted. |
| Structure factors | `{pdbid}-sf.cif` | Deposited reflections. Fo−Fo and the omit map use the amplitudes (after conversion). The 2mFo−DFc workflow uses this file for the resolution range. |
| Validation $2F_o-F_c$ coefficients | `{pdbid}_validation_2fo-fc_map_coef.cif` | `FWT` and `PHWT`. This file is the 2mFo−DFc map. In the Fo−Fo workflow it is only a phase source. |

Validation coefficients are also published at

```text
https://files.rcsb.org/pub/pdb/validation_reports/xx/{pdbid}/{pdbid}_validation_2fo-fc_map_coef.cif.gz
```

Deposited structure factors for the `6w1o` and `6w1v` examples are intensities, `_refln.intensity_meas` and `_refln.intensity_sigma`, with the free set marked by `_refln.status` `f`. `cif2mtz` stores those as `I`, `SIGI`, and `FREE`. The omit-map workflow converts intensities to amplitudes with `ctruncate` and copies `FREE` onto the amplitude MTZ. The Fo−Fo workflow lets a zero-cycle REFMAC job place the amplitudes in a REFMAC-style MTZ. The `9z35` structure-factor file stores anomalous intensities and `_refln.pdbx_r_free_flag`. The omit-map workflow converts that Friedel pair to mean amplitudes. The sharpened-map workflow uses the same file for the resolution range.

Column names are read from each MTZ with `mtzdump` and matched against a short list of candidates (`FP`, `F`, `FOBS`, `SIGFP`, `PHWT`, `PHIC`, `DELFWT`, `FWT`, and so on), so a slightly different mmCIF layout still runs.

## Resolution limits

`RESOLUTION_SOURCE = "sf"` reads the resolution range from the structure-factor MTZ after `cif2mtz`. `mtzdump` reports that range. The mmCIF often has no `pdbx_d_res_low` or `pdbx_d_res_high`, so the MTZ is the source used here.

For two datasets, the common shell is the poorer low-resolution limit and the poorer high-resolution limit.

`RESOLUTION_SOURCE = "explicit"` ignores the MTZ range and uses `LOWRES` and `DMIN`. `LOWRES` is the larger number in ångströms. `DMIN` must be coarser than or equal to the high-resolution limit of the coefficients that are actually transformed. Very low-resolution terms are sensitive to bulk solvent and scaling. An explicit `LOWRES` is how you cut them when the MTZ range is wider than you want.

## The shared tail of every workflow

After the map exists, the three workflows do the same three things, and then open the fragment in the browser.

### Extend the map

`fft` writes an asymmetric-unit map. `extends` expands it so that it covers the coordinate file plus a border, `MAP_BORDER_A` (default 5 Å). The coordinate file is the light model, the refined omit model, or the single structure, depending on the workflow.

### Clip a site fragment

Headless Coot (`--no-graphics --script`) loads that extended map, centres on the centroid of each residue in `SITE_COPIES`, and writes a sphere of radius `EXPORT_RADIUS_A` (default 10 Å). Difference workflows load the map as a difference map. The 2mFo−DFc workflow loads it as a regular map.

The generated Coot script is Python 2, because the CCP4 `coot` wrapper still ships a Python 2 interpreter.

Leaving fragment recomputation off reuses an existing clip. The fragment is what PyMOL and the browser contour, so a change of radius requires a new clip.

### Draw the figure

PyMOL loads the coordinate file and the fragment.

- Surface `mesh` draws a wireframe with `isomesh`. Surface `solid` draws a transparent surface with `isosurface`.
- The carve radius (default 3 Å) hides the surface farther than this distance from the displayed atoms, so the cage stays on the site.
- The site residue is drawn as sticks and spheres. Each neighbouring residue is drawn as sticks. The rest of the model stays hidden.
- Mn/Ca–oxygen bonds are added when those atom names exist. A few default cluster distances are drawn the same way. A structure without those names simply skips the bonds.
- The scene is written as a `.pml` file next to the PNG.

PyMOL is always re-run. Contour, colour, and camera changes do not require a new CCP4 map. Leaving map recomputation off reuses the extended CCP4 map and skips the slow steps.

### The browser viewer

The page then loads the same model and the same clipped map with NGL. The viewer draws one isosurface. The PNG sigma levels are the contours in the PyMOL figure; the slider does not change that PNG.

The slider runs from the viewer sigma minimum to the viewer sigma maximum on the form. The number box sets any sigma inside that range. The starting level is the PNG level closest to zero that still sits inside the range. Difference maps default to \(-5\) through \(5\). The sharpened \(2mF_o-DF_c\) map defaults to \(1\) through \(5\).

Two colourings are available. Blue/red paints positive sigma blue and negative sigma red, darker farther from zero, and paints zero grey. Temperature paints the low end of the slider cold (deep blue) and the high end hot (deep red), through cyan, yellow, and orange. Solid and mesh switch the isosurface style, and the transparency slider runs from opaque to clear.

The edit buttons act on the model in the viewer. Add bond and Add distance each take two atom clicks. Delete removes a stick, or the bond and distance between two atoms. Add label shows an atom name, and a second click removes it. Distance labels turn the ångström lengths on the dashed lines on or off. Save PNG writes the orientation currently on screen. Reset view returns the camera to the site.

## Programs

| Program | Where it is used |
| --- | --- |
| `cif2mtz` | Reflection mmCIF to MTZ, in every workflow |
| `mtzdump` | Column names and resolution limits |
| `refmac5` | Zero-cycle amplitude preparation (Fo−Fo), or the omit-model job that writes `DELFWT` |
| `cad` | Combine amplitude and phase columns (Fo−Fo), or copy `FREE` (omit map) |
| `scaleit` | Scale light amplitudes onto the dark reference (Fo−Fo only) |
| `ctruncate` | Intensities to amplitudes (omit map only) |
| `libcheck` | Optional monomer library when omit-map distances are widened in the dictionary |
| `fft` | The Fourier transform, including sharpening for 2mFo−DFc |
| `extends` | Pad the map around the model |
| Coot | Clip the site |
| PyMOL | Contour and render |

# Omit-model mFo−DFc map

This workflow calculates the weighted difference map \(mF_o-DF_c\) for an omit or reference model against one set of observed amplitudes, clips a fragment around a site, and draws the \(\pm n\sigma\) surfaces on the coordinates written by REFMAC.

One REFMAC job does the comparison. It reads the omit coordinates and the observed amplitudes, runs for the number of cycles on the form, and writes those coordinates together with \(F_o\) (`FP`), \(F_c\) (`FC` or `FC_ALL`), and the weighted difference `DELFWT` / `PHDELWT`. The map is the Fourier transform of `DELFWT` / `PHDELWT`. No second dataset is scaled in, so `scaleit` is not used. `cad` copies a free-R flag from the intensity MTZ onto the amplitudes when that column is present.

The [theoretical background](theory.md) defines \(m\), \(D\), and a refinement cycle. The steps below are the programs that compute them for the worked example.

## Worked example

Fig. 2a and Fig. 2b of Wang, Armstrong, and Batista, [Proc. Natl. Acad. Sci. U.S.A. 118, e2023982118 (2021)](https://doi.org/10.1073/pnas.2023982118), compare an omit model of the oxygen-evolving complex with the observed amplitudes of the Ibrahim *et al.* two-flash dataset `6w1v`. The inputs are the `6w1v` coordinates, the `6w1v` structure-factor mmCIF, and the single-conformation omit model `6w1v_omitO6_single_confA.pdb`.

| Role | What the form needs | Download |
| --- | --- | --- |
| Omit coordinates | `6w1v_omitO6_single_confA.pdb`. O6 is absent. This is the model REFMAC reads, and the figure is drawn on the coordinates REFMAC writes. | [6w1v_omitO6_single_confA.pdb](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v_omitO6_single_confA.pdb) |
| Observed structure factors | `6w1v-sf.cif`. Intensities of the 2F, 200 ms dataset. \(F_o\) comes from this file. | [6w1v-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-sf.cif.gz) |
| Light coordinates | The deposited `6w1v` model, used as the dataset’s coordinate file on the form. | [bundle 1](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle1.pdb), [bundle 2](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle2.pdb) |

The structure-factor file in the repository is gzip-compressed. Point the form at the `.gz` file; a name ending in `.gz` is decompressed with `gunzip`. The deposited `6w1v` coordinates are stored as two PDB bundles because the entry exceeds one classic PDB file. Bundle 1 contains the cluster.

In `6w1v_omitO6_single_confA.pdb` the cluster that corresponds to the major site is a single conformation, ligand `OEX`, chain `A`, residue **415**, with occupancy 1 and with O6 removed. The deposited `6w1v` entry is different: residue **415** is the minor conformation `OEX` and residue **416** is the major conformation `OEY`, which includes O6. Set the site from the omit file, chain `A` residue `415`, and list D1-Glu189 (`A 189`) as a neighbour. Confirm both with **List hetero residues** before trusting a residue number copied from `6w1v`. The default PNG contour for this workflow is a single mesh at \(3.5\sigma\).

Positive density in the figure is density in the `6w1v` observations that this omit model does not include. The model drawn is the REFMAC coordinate file. With zero cycles, the default, that file is the omit model rewritten with hydrogens and with the atoms unmoved.

### How the Fig. 2 map is calculated

The programs run in this order: `cif2mtz`, `mtzdump`, `ctruncate`, `cad`, `refmac5`, `fft`, `extends`, Coot, and PyMOL. `libcheck` is added only when distance restraints are widened in a monomer library.

**1. Check the programs and the files.** The programs are `cif2mtz`, `ctruncate`, `cad`, `refmac5`, `fft`, `extends`, and `mtzdump`, plus `libcheck` when that option is on. Coot comes from the CCP4 tree. PyMOL must be on `PATH`. The log prints the cycle count before REFMAC starts, so a zero-cycle job is visible in advance.

**2. Inspect the omit model and choose the site.** **List hetero residues** reads `6w1v_omitO6_single_confA.pdb`. Select `OEX` A415. The job checks that chain and residue, prints the atom names, and prints the centroid of the clip.

**3. Convert the observed mmCIF.** `cif2mtz` writes `6w1v-sf.mtz`. The columns at this stage are still the deposited intensities. For this dataset those are `_refln.intensity_meas` and `_refln.intensity_sigma`, stored in the MTZ as `I` and `SIGI`.

**4. Choose the resolution shell.** With the resolution source set to the structure-factor file, `mtzdump` sets the low- and high-resolution cutoffs from `6w1v-sf.mtz`. With an explicit shell, the two numbers on the form are used. The explicit pair stored for this workflow is 30.0–1.89 Å, which matches the high-resolution limit recorded for the `6w1v` intensities. The high-resolution cutoff should be no finer than that limit.

**5. Turn intensities into amplitudes.** `ctruncate` reads `I` and `SIGI` and writes `6w1v-sf-truncate.mtz` with amplitudes `F` and `SIGF`:

```text
ctruncate -hklin 6w1v-sf.mtz -hklout 6w1v-sf-truncate.mtz -colin '/*/*/[I,SIGI]'
```

The program applies the French–Wilson estimate: for each reflection it returns the expected \(|F|\) given the measured intensity and its uncertainty, including the case where the measured intensity is negative. Those amplitudes are \(|F_o|\).

**6. Copy the free-R flag.** `cad` takes every column from the amplitude MTZ and the `FREE` column from the intensity MTZ, and writes `6w1v-sf-truncate-free.mtz`:

```text
LABIN FILE 1 ALL
LABIN FILE 2 E1=FREE
```

`FREE` marks reflections that REFMAC must leave out of the refinement target. In `6w1v-sf.cif` that set is the status value `f` on `_refln.status`. `cif2mtz` writes it as `FREE`, and `cad` copies it onto the amplitudes so \(R_{\mathrm{free}}\) is still computed on the held-out reflections.

**7. Calculate \(F_c\), or refine, with REFMAC.** The job reads the amplitude MTZ and `6w1v_omitO6_single_confA.pdb` and writes `6w1v_omitO6_single_confA_vs_6w1v.mtz` and `6w1v_omitO6_single_confA_refined.pdb`.

Each cycle computes

\[
F_c(hkl) = \sum_j \mathrm{occ}_j\, f_j(|\mathbf{s}|)\, e^{-B_j |\mathbf{s}|^{2}/4}\, e^{2\pi i\, \mathbf{h}\cdot\mathbf{r}_j} + F_{\mathrm{solvent}}
\]

from the current coordinates, scales \(|F_o|\) to \(|F_c|\), and forms the weighted difference

\[
m|F_o| - D|F_c|
\]

with phase \(\varphi_c = \arg(F_c)\). REFMAC stores that complex number as `DELFWT` and `PHDELWT`. \(m\) is the figure of merit. \(D\) down-weights \(F_c\) as the coordinate error grows. Because \(F_o\) and \(F_c\) are scaled inside this job, a second scaling program is unnecessary.

`NCYC 0` stops after that calculation. Atoms do not move. A larger cycle count repeats the shift of coordinates and \(B\) factors against the likelihood target plus the geometry restraints, and the figure is drawn on the shifted model. The keywords always include

```text
NCYC <cycles>
MAKE HYDR Y
MAKE HOUT Y
```

`MAKE HYDR` / `MAKE HOUT` generate riding hydrogens and write them to the coordinate file. Zero cycles still produce a coordinate file.

External restraints, when the lists on the form are non-empty, are written to a text file and placed in front of those keywords. The default distance pairs are Mn1–O5, Mn3–O5, Mn4–O5, and Ca1–O5, with a sigma of 0.50 Å. The default angles are the angles at O5 among Mn3, Mn4, and Ca, the O5–Ca–O1 angle, and the angles at Mn3 and Mn4 that include O5, with a sigma of 25°. With the library option off, the distances are external restraints and replace the dictionary bond for those atom pairs in the listed copies only. With the library option on, `libcheck` widens those distances in a monomer library that applies to every copy of that residue type, and the angles stay external. An empty distance or angle list adds no restraint of that kind.

The MTZ is checked for `FP` (or `F`, `FOBS`), `FC_ALL` (or `FC`), `DELFWT`, and `PHDELWT`. The transform uses the weighted coefficient.

**8. Fourier-transform the difference.** `fft` writes `6w1v_DELFWT_6w1v_omitO6_single_confA_{tag}.ccp4`:

```text
RESOLUTION <low> <high>
LABIN F1=DELFWT PHI=PHDELWT
```

This is the sum

\[
\Delta\rho(\mathbf{r}) \propto \sum_{hkl} \mathrm{DELFWT}(hkl)\, e^{i\, \mathrm{PHDELWT}(hkl)}\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}}.
\]

Positive density is present in the `6w1v` observations and absent from the omit model. Negative density is placed by the model where the observations are weaker.

**9. Cover the REFMAC model.** `extends` pads the map by the map border (default 5 Å) around the refined coordinate file and writes the extended CCP4 map.

**10. Export the site fragment.** Coot loads the extended map as a difference map and the REFMAC coordinates, centres on `OEX` A415, and writes a clip of the chosen radius (default 10 Å). The Coot contour does not change the exported voxels.

**11. Contour and render.** PyMOL loads the REFMAC coordinates and the fragment. Colours follow the Fo−Fo figure.

| Sign | Lowest contour | Middle | Highest contour |
| --- | --- | --- | --- |
| Positive, in the data and missing from the model | forest | cyan | blue |
| Negative, in the model and missing from the data | red | salmon | magenta |

Cluster bonds use the Mn/Ca–oxygen pairs without an O6 pair, matching an omit model from which O6 was removed. Bonds are added only when the atom names exist.

The browser viewer draws one isosurface of the same fragment. Its slider runs from the viewer sigma minimum to the maximum (default \(-5\) to \(5\)). Blue/red colouring uses blue for positive density and red for negative density. The PNG stays at the PyMOL levels, \(3.5\sigma\) unless you change that list.

Outputs:

| File | Contents |
| --- | --- |
| `6w1v_DELFWT_{omit}_{site}_{chain}{resi}_diffmap.png` | The PyMOL figure |
| `pymol_6w1v_DELFWT_{omit}_{site}.pml` | The scene |
| `{omit}_refined.pdb` | Coordinates after the requested cycles |
| `*_extended.ccp4` | The map |

A map left by an older `cad` / `scaleit` pipeline has a different file name and is not reused. Turn map recomputation on to overwrite the extended map after any change of cycles, restraints, or resolution.

## Background

The sign convention is data minus model.

- Positive density (forest, cyan, blue at increasing \(\sigma\)) means electrons in the observations that the omit model does not account for.
- Negative density (red, salmon, magenta) means electrons the omit model places where the observations are weaker.
- An atom present in the data and left out of the model is a positive peak.

XFEL difference maps are noisy, and the activated-state occupancy is often only 10–30%. A contour of \(3\sigma\) shows weaker features and more noise. Contours of \(4\)–\(5\sigma\) keep the strongest peaks. The worked example draws a single mesh at \(3.5\sigma\). Report the level you used.

The point of an omit map is to keep the feature of interest out of \(F_c\), so that it can appear in the difference. Refinement will try to fill unexplained density by moving nearby atoms. The distance and angle restraints are how selected cluster geometry is kept from collapsing while those atoms are still allowed to move. A normal dictionary bond sigma is about 0.02 Å, so the default 0.50 Å is loose. `refmac5` reads the angle sigma in degrees, so the number on the form is written unchanged.

## Files

| Role | Form field | Required for the map |
| --- | --- | --- |
| Omit coordinates | Omit-model coordinates | Yes. This is the model REFMAC reads. |
| Observed structure factors | Observed structure factors | Yes. |
| Light coordinates | Light coordinates | Yes, as the coordinate file of the dataset being compared. The figure is drawn on the REFMAC output. |

The omit model’s own structure-factor file is not an input. \(F_o\) comes from the observed dataset. The two name fields are labels in the output files. They do not have to be PDB codes.

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| REFMAC cycles | `0` leaves atoms fixed. A larger integer refines against the observed amplitudes |
| Distance restraints and distance sigma | Which site distances are loose, and how loose |
| Angle restraints and angle sigma | External angles, sigma in degrees |
| Widen distances in a libcheck library | External distances on the listed copies, or a widened monomer library for every copy of that residue type |
| Resolution source and the two cutoffs | Shell of the `DELFWT` transform |
| PNG sigma levels | The worked example uses a single \(3.5\sigma\) mesh |
| Viewer sigma minimum and maximum | Range of the one browser isosurface |

## Reading the figure

The model in the figure is the REFMAC coordinate set. With zero cycles it is the omit model rewritten by REFMAC, hydrogens included, atoms unmoved. A positive peak where `6w1v` has an atom and the omit model does not is the expected \(mF_o-DF_c\) signal. Whether the 2F data support an inserted oxo or water at the O6 position is the question discussed by Wang, Armstrong, and Batista; this workflow reproduces the map they compare.

If the mesh is empty, lower the contour (try `3`). If the cage swallows the model, decrease the carve radius or raise the contour.

## Adapting it

1. Provide the omit coordinates, the observed structure factors, and the light coordinates.
2. Set the two labels. Leave the resolution source on the structure-factor file, or set the two cutoffs yourself.
3. List hetero residues on the omit model, then set the chain, residue, and neighbours.
4. Set the cycle count and the restraint lists. Leave the library option off to pass distances straight to REFMAC. Turn it on to widen them in a `libcheck` library.
5. Recompute the map the first time, and again after any change of cycles, restraints, or resolution. Leave recomputation off while you tune σ and the view.

The extended map can be opened with the REFMAC coordinates:

```text
coot {extended_map} {omit}_refined.pdb
```

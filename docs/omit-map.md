# Omit-model mFo−DFc map

This workflow calculates the weighted difference map \(mF_o-DF_c\) for an omit or reference model against one set of observed amplitudes, clips a fragment around a site, and draws the \(\pm n\sigma\) surfaces on the coordinates written by REFMAC.

One REFMAC job does the comparison. It reads the omit coordinates and the observed amplitudes, runs for the number of cycles on the form, and writes those coordinates together with \(F_o\) (`FP`), \(F_c\) (`FC` or `FC_ALL`), and the weighted difference `DELFWT` / `PHDELWT`. The map is the Fourier transform of `DELFWT` / `PHDELWT`. No second dataset is scaled in, so `scaleit` is not used. `cad` copies a free-R flag from the intensity MTZ onto the amplitudes when that column is present.

The [theoretical background](theory.md) defines \(m\), \(D\), and a refinement cycle. The steps below are the programs that compute them for the worked example.

## Worked example

Fig. 2a and Fig. 2b of Wang, Armstrong, and Batista, [Proc. Natl. Acad. Sci. U.S.A. 118, e2023982118 (2021)](https://doi.org/10.1073/pnas.2023982118), compare an omit model of the oxygen-evolving complex with the observed amplitudes of the Ibrahim *et al.* two-flash dataset `6w1v`. The inputs are the `6w1v` coordinates, the `6w1v` structure-factor mmCIF, and the single-conformation omit model `6w1v_omitO6_single_confA.pdb`.

| Role | What the form needs | Download |
| --- | --- | --- |
| Omit coordinates | `6w1v_omitO6_single_confA.pdb`. O6 is absent. This is the model REFMAC reads, and the figure is drawn on the coordinates REFMAC writes. | [6w1v_omitO6_single_confA.pdb](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v_omitO6_single_confA.pdb) |
| Observed structure factors | `6w1v-sf.cif`, or an MTZ of the same reflections. Intensities of the 2F, 200 ms dataset. \(F_o\) comes from this file. | [6w1v-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-sf.cif.gz) |
| Light coordinates | The deposited `6w1v` model, used as the dataset’s coordinate file on the form. | [bundle 1](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle1.pdb), [bundle 2](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle2.pdb) |

The structure-factor file in the repository is gzip-compressed. Point the form at the `.gz` file; a name ending in `.gz` is decompressed with `gunzip`. The deposited `6w1v` coordinates are stored as two PDB bundles because the entry exceeds one classic PDB file. Bundle 1 contains the cluster.

In `6w1v_omitO6_single_confA.pdb` the cluster that corresponds to the major site is a single conformation, ligand `OEX`, chain `A`, residue **415**, with occupancy 1 and with O6 removed. The deposited `6w1v` entry is different: residue **415** is the minor conformation `OEX` and residue **416** is the major conformation `OEY`, which includes O6. Set the site from the omit file, chain `A` residue `415`, and list D1-Glu189 (`A 189`) as a neighbour. Confirm both with **List hetero residues** before trusting a residue number copied from `6w1v`. The default PNG contour for this workflow is a single mesh at \(3.5\sigma\).

Positive density in the figure is density in the `6w1v` observations that this omit model does not include. The model drawn is the REFMAC coordinate file. With zero cycles, the default, that file is the omit model rewritten with hydrogens and with the atoms unmoved.

### How the Fig. 2 map is calculated

The programs run in this order: `cif2mtz`, `mtzdump`, `ctruncate`, `cad`, `refmac5`, `fft`, Coot, and PyMOL. The map is then padded around the model in the app.

**1. Check the programs and the files.** The programs are `cif2mtz`, `ctruncate`, `cad`, `refmac5`, `fft`, and `mtzdump`. Coot comes from the CCP4 tree. PyMOL must be on `PATH`. The log prints the cycle count before REFMAC starts, so a zero-cycle job is visible in advance.

**2. Inspect the omit model and choose the site.** **List hetero residues** reads `6w1v_omitO6_single_confA.pdb`. Select `OEX` A415. The job checks that chain and residue, prints the atom names, and prints the centroid of the clip.

**3. Convert the observed mmCIF, or use an MTZ.** `cif2mtz` writes `6w1v-sf.mtz`. The columns at this stage are still the deposited intensities. For this dataset those are `_refln.intensity_meas` and `_refln.intensity_sigma`, stored in the MTZ as `I` and `SIGI`. An MTZ supplied on the form is copied to that file and `cif2mtz` is not run. Intensities in it still follow the next steps. Amplitudes already named `F`/`SIGF`, `FP`/`SIGFP`, `FOBS`/`SIGFOBS`, or `FMEAN`/`SIGFMEAN` skip `ctruncate`. A `FREE` column is the free set; without one, every reflection is used.

**4. Choose the resolution shell.** The Resolution section starts on the structure-factor file. `mtzdump` then sets the low- and high-resolution cutoffs from `6w1v-sf.mtz`. Choose explicit cutoffs and the two fields appear. The explicit pair stored for this workflow is 30.0–1.89 Å, which matches the high-resolution limit recorded for the `6w1v` intensities. That shell is the range of `DELFWT` coefficients. The FFT-grid switch in the same section starts off, which samples the map on the grid for this shell. Extended map border, default 5 Å, is the padding used in step 9.

**5. Turn intensities into amplitudes.** `ctruncate` reads `I` and `SIGI` and writes `6w1v-sf-truncate.mtz` with amplitudes `F` and `SIGF`:

```text
ctruncate -hklin 6w1v-sf.mtz -hklout 6w1v-sf-truncate.mtz -colin '/*/*/[I,SIGI]'
```

The program applies the French–Wilson estimate: for each reflection it returns the expected \(|F|\) given the measured intensity and its uncertainty, including the case where the measured intensity is negative. Those amplitudes are \(|F_o|\).

A deposition that stores Friedel intensities instead of a mean, such as `9z35`, has `I(+)`, `SIGI(+)`, `I(-)`, and `SIGI(-)` after `cif2mtz`. The same step then runs:

```text
ctruncate -hklin 9z35-sf.mtz -hklout 9z35-sf-truncate.mtz -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'
```

That writes the mean amplitudes as `FMEAN` and `SIGFMEAN`, together with the Friedel amplitudes. Mean intensities are used when both layouts are present.

**6. Copy the free-R flag.** For a mean-intensity file, `cad` takes every column from the amplitude MTZ and the `FREE` column from the intensity MTZ, and writes `6w1v-sf-truncate-free.mtz`:

```text
LABIN FILE 1 ALL
LABIN FILE 2 E1=FREE
```

`FREE` marks reflections that REFMAC must leave out of the refinement target. In `6w1v-sf.cif` that set is the status value `f` on `_refln.status`. `cif2mtz` writes it as `FREE`, and `cad` copies it onto the amplitudes so \(R_{\mathrm{free}}\) is still computed on the held-out reflections.

For a Friedel-intensity file, `cad` copies `FMEAN` and `SIGFMEAN` only, renames them to `F` and `SIGF`, and adds `FREE`. REFMAC then sees the same mean-amplitude columns as in the `6w1v` run.

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
MAKE -
    HYDROGEN NO -
    HOUT NO -
    PEPTIDE NO -
    CISPEPTIDE YES -
    SSBRIDGE YES -
    SYMMETRY YES -
    SUGAR YES -
    CONNECTIVITY NO -
    LINK NO
REFI -
    TYPE REST -
    RESI MLKF -
    METH CGMAT -
    BREF ISOT
SCAL -
    TYPE SIMP -
    LSSC -
    ANISO -
    EXPE
SOLVENT YES
WEIGHT -
    AUTO
NCYC <cycles>
```

Riding hydrogens are not built and are not written to the coordinate file. The REFMAC refinement section starts off, which forces `NCYC 0`. Turn it on to set the cycle count. `REFI TYPE` is `REST` (restrained, the default) or `UNRE` (unrestrained). `REFI BREF` is `ISOT` (isotropic, the default), `ANIS`, `MIXE`, or `OVER`. The block above is that default. Scaling is simple, anisotropic, and experimental, bulk solvent is included, and the X-ray weight is automatic. A line that ends in `-` continues the same keyword. Zero cycles still produce a coordinate file. The residual is maximum likelihood (`RESI MLKF`) with the conjugate-gradient method (`METH CGMAT`).

While REFMAC refinement is on, external restraints are written to a text file and placed after a custom dictionary and in front of the fixed keywords. Each group is used only while its switch is on. The distance and angle switches start on. The dihedral switch starts off. The default distance pairs are Mn1–O5, Mn3–O5, Mn4–O5, and Ca1–O5, with a sigma of 0.50 Å. The default angles are the angles at O5 among Mn3, Mn4, and Ca, the O5–Ca–O1 angle, and the angles at Mn3 and Mn4 that include O5, with a sigma of 25°. A dihedral line names four atoms, and its default sigma is 10°. A line such as `MN1 O5 2.40 0.50` uses that distance and that sigma. A line such as `MN3 O5` keeps the distance in the starting model and uses the Distance sigma field. Angles and dihedrals work the same way, in degrees. These external restraints replace the CCP4 monomer library for the same atoms, in the listed copies only. A switch that is off adds no restraint of that kind. The protein harmonic switch also starts off. Turn it on to write one `exte harm` line for every protein atom in every chain, except atoms in the site residue, the additional copies, and the neighbouring residues. Waters, ligands, and other non-protein residues are left out. The sigma field is shown only while that switch is on, and the default is 0.1 Å. Those lines are written only when the cycle count is greater than zero. REFMAC holds each listed atom at the coordinates in the input model. A custom dictionary is a `LIBIN` line. Its file fields are shown only while REFMAC refinement is on and its own switch is on, and it is read only when that switch is on and the cycle count is greater than zero. The external restraints follow that line, so they replace matching distances, angles, and dihedrals in the dictionary.

The MTZ is checked for `FP` (or `F`, `FOBS`), `FC_ALL` (or `FC`), `DELFWT`, and `PHDELWT`. The transform uses the weighted coefficient.

**8. Fourier-transform the difference.** `fft` writes `6w1v_DELFWT_6w1v_omitO6_single_confA_{tag}.ccp4`:

```text
RESOLUTION <low> <grid>
LABIN F1=DELFWT PHI=PHDELWT
```

`<grid>` is the high-resolution cutoff when the FFT-grid switch is off. Turn the switch on and a smaller number replaces `<grid>` only. `fft` takes the maximum Miller indices from `<grid>` and the cell, then samples the asymmetric unit at about one third of that spacing. Indices inside the extended box that are not in the MTZ are omitted from the sum, so their amplitude is zero. The coefficients that were measured are unchanged. When the grid is finer, the map name gains `_grid` and that limit, for example `_grid1p20`.

This is the sum

\[
\Delta\rho(\mathbf{r}) \propto \sum_{hkl} \mathrm{DELFWT}(hkl)\, e^{i\, \mathrm{PHDELWT}(hkl)}\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}}.
\]

Positive density is present in the `6w1v` observations and absent from the omit model. Negative density is placed by the model where the observations are weaker.

When the whole-cell grid is too large for `fft` (the three sampling intervals no longer fit in its 32-bit allocator, so the program stops with `CCPALC: length <=0` or a segmentation fault), the app evaluates this sum only on a cube centered on the site. The cube is wide enough for the clip radius. Reflections beyond the measured coefficients still contribute nothing. The file is a P1 map of that cube: the cell edges are scaled with the grid so the crystallographic grid indices still fall on the model coordinates. The sigma stored on that cube is the full-cell root-mean-square of the coefficients.

**9. Cover the REFMAC model.** The app pads the map by the extended map border (default 5 Å) around the refined coordinate file and writes the extended CCP4 map. The sigma on that file stays the value from `fft`. A grid that `fft` cannot allocate skips this pad: the site cube is already the map Coot reads.

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

A map left by an older `cad` / `scaleit` pipeline has a different file name and is not reused. Turn map recomputation on to overwrite the extended map after any change of cycles, restraints, resolution, or FFT grid resolution for a final map. A finer grid writes a different file name, so it does not replace the map sampled at the coefficient limit.

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
| Observed structure factors | Observed structure factors, mmCIF or MTZ | Yes. |
| Light coordinates | Light coordinates | Yes, as the coordinate file of the dataset being compared. The figure is drawn on the REFMAC output. |

The omit model’s own structure-factor file is not an input. \(F_o\) comes from the observed dataset. The two name fields are labels in the output files. They do not have to be PDB codes.

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| REFMAC refinement | Starts off, which forces zero cycles. On shows cycles, B-factor refinement, refinement type, the restraints, the protein harmonic restraint, and the custom dictionary |
| B-factor refinement | `REFI BREF`. ISOTropic (default), ANISotropic, MIXEd, or OVERall |
| Refinement type | `REFI TYPE`. RESTrained (default) or UNREstrained |
| REFMAC cycles | Shown while the section is on. `0` leaves atoms fixed. A larger integer refines against the observed amplitudes |
| Distance restraints and distance sigma | Inside the section. Switch starts on. External distances. A sigma on a line replaces the field for that pair |
| Angle restraints and angle sigma | Inside the section. Switch starts on. External angles, sigma in degrees. A sigma on a line replaces the field for that angle |
| Dihedral restraints and dihedral sigma | Inside the section. Switch starts off. External torsions, four atoms, sigma in degrees. The default sigma is 10 |
| Protein harmonic restraints and sigma | Inside the section. Switch starts off. Holds every protein atom, in every chain, except the site residue, additional copies, and neighbouring residues. Default sigma 0.1 Å. Applied only when the cycle count is greater than zero |
| Custom dictionary | Inside the section. Switch starts off. Extra monomer library, used only when the switch is on and the cycle count is greater than zero. External restraints are read after it |
| Resolution source | Shell of the `DELFWT` coefficients. Explicit cutoffs show the two fields |
| Custom FFT grid for the final map | Off uses that shell. On, a smaller number samples the finer grid, with zero amplitude beyond the coefficients |
| Extended map border | Padding around the refined model. Default 5 Å |
| PNG sigma levels | The worked example uses a single \(3.5\sigma\) mesh |
| Viewer sigma minimum and maximum | Range of the one browser isosurface |

## Reading the figure

The model in the figure is the REFMAC coordinate set. With zero cycles it is the omit model rewritten by REFMAC, with the atoms unmoved and no riding hydrogens added. A positive peak where `6w1v` has an atom and the omit model does not is the expected \(mF_o-DF_c\) signal. Whether the 2F data support an inserted oxo or water at the O6 position is the question discussed by Wang, Armstrong, and Batista; this workflow reproduces the map they compare. The [reconstruction](reconstruction.md) repeats this calculation and places the missing reference atom closest to each peak.

If the mesh is empty, lower the contour (try `3`). If the cage swallows the model, decrease the carve radius or raise the contour.

## Adapting it

1. Provide the omit coordinates, the observed structure factors, and the light coordinates.
2. Set the two labels. Leave the resolution source on the structure-factor file, or switch to explicit cutoffs and set the two fields. Leave the FFT-grid switch off unless the final map should be sampled more finely than the coefficients.
3. List hetero residues on the omit model, then set the chain, residue, and neighbours.
4. Leave REFMAC refinement off for zero cycles. Turn it on to set the cycle count, `REFI BREF`, `REFI TYPE`, and the restraint lists. A line with its own distance and sigma uses those values. A line with only the atom names uses the sigma fields. The protein harmonic switch holds every protein atom outside the site, the additional copies, and the neighbouring residues.
5. Recompute the map the first time, and again after any change of cycles, restraints, the protein harmonic restraint, B-factor model, refinement type, resolution, or FFT grid. Leave recomputation off while you tune σ and the view.

The extended map can be opened with the REFMAC coordinates:

```text
coot {extended_map} {omit}_refined.pdb
```

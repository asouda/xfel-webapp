# Omit-model mFo−DFc map

This workflow calculates the weighted difference map $mF_o-DF_c$ for an omit or reference model against one set of observed amplitudes, clips a fragment around a site, and draws the $\pm n\sigma$ surfaces on the refined coordinates.

The notebook is `notebooks/XFEL_fo-fc_diffmap_refinement.ipynb`. The web application runs the same steps from the omit-map page.

One REFMAC job does the whole comparison. It reads the omit coordinates and the light amplitudes, runs for `REFMAC_NCYC` cycles, and writes the refined coordinates together with $F_o$ (`FP`), $F_c$ (`FC` or `FC_ALL`), and the weighted difference `DELFWT` / `PHDELWT`. The map is the Fourier transform of `DELFWT` / `PHDELWT`. No second dataset is scaled in, so `scaleit` is not used. `cad` only copies the `FREE` flag from the intensity MTZ onto the amplitudes.

## Worked example

The configuration cell uses a single-conformation omit model of the Ibrahim *et al.* 2F photosystem II dataset.

| Role | Name | What it supplies |
| --- | --- | --- |
| Omit model | `6w1v_omitO6_single_confB` | Coordinates. O6 is absent from this single conformation. |
| Observed data | `6w1v` | 2F, 200 ms structure factors. $F_o$ always comes from this dataset. |

The notebook introduction places the comparison with Fig. 2B of Wang, Armstrong, and Batista, PNAS, 2021, which uses that 2F dataset. The closing notes of the notebook also describe the same kind of figure for the deposited single-conformation model `11yo` against the `6w1v` amplitudes. In either case the model drawn in the figure is the REFMAC output, and positive density is density in the observations that the omit model does not include.

The site in the configuration cell is chain `A` residue **416**, the `OEY` cluster with O6 removed. Neighbouring sticks are D1-Glu189. Confirm the residue from the hetero-residue list before trusting a residue number copied from another PDB entry.

## Background

The map is the REFMAC weighted difference

$$
\Delta\rho(\mathbf{r}) \;\propto\; \sum_{hkl} \big(m|F_o|-D|F_c|\big)\, e^{i\varphi_c}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}},
$$

with $|F_o|$ the observed light amplitude and $|F_c|$, $\varphi_c$ calculated from the omit model. $m$ is the figure of merit. $D$ down-weights $F_c$ as the model’s coordinate error grows. REFMAC stores the weighted coefficient as `DELFWT` with phase `PHDELWT`. Because $F_o$ and $F_c$ were scaled inside that one job, a second program does not have to put them on a common scale.

### How to read the sign

- Positive density (forest, cyan, blue at increasing $\sigma$) means electrons in the observations that the omit model does not account for.
- Negative density (red, salmon, magenta) means electrons the omit model places where the observations are weaker.
- An atom present in the data and left out of the model — for example an oxygen removed from a single-conformation cluster — is a positive peak.

This is a different sign convention from a casual reading of “light minus dark.” Here the subtraction is data minus model. The colours match the Fo−Fo figure, and the caption is the one to trust: positive density is in the light observations and absent from the omit model.

XFEL difference maps are noisy, and the activated-state occupancy is often only 10–30 %. A contour of $3\sigma$ shows weaker features and more noise. Contours of $4$–$5\sigma$ keep the strongest peaks. The configuration cell draws a single mesh at $3.5\sigma$. Report the level you used.

### Zero cycles or a short refinement

`REFMAC_NCYC = 0`, which is what the configuration cell and the web form use, calculates $F_c$ and the difference map without moving atoms. A positive cycle count refines coordinates and $B$ factors of the omit model against the observed amplitudes, and the figure is drawn on that refined model.

The point of an omit map is to keep the feature of interest out of $F_c$, so that it can appear in the difference. Refinement will try to fill unexplained density by moving nearby atoms. The restraints below are how the notebook keeps selected cluster distances and angles from collapsing while still allowing them to move.

### Distance and angle restraints

`SITE_DISTANCE_RESTRAINTS` lists atom pairs in the site residue. Each entry is one of

- `(atom_1, atom_2)`, which keeps the distance currently in the model
- `(atom_1, atom_2, distance_Å)`
- `(atom_1, atom_2, distance_Å, sigma_Å)`

The default pairs are Mn1–O5, Mn3–O5, Mn4–O5, and Ca1–O5. An empty list adds no distance restraints.

`SITE_DISTANCE_SIGMA` (default 0.50 Å) is the standard deviation used when an entry does not set its own. A normal dictionary bond is about 0.02 Å, so 0.50 Å is loose and those atoms can move.

`SITE_ANGLE_RESTRAINTS` lists angles in the same residue. Each entry is `(atom_1, vertex, atom_3)`, optionally followed by the angle and a sigma. The angle is at the middle atom. A missing angle keeps the angle in the starting model. The defaults are the angles at O5 among Mn3, Mn4, and Ca, the O5–Ca–O1 angle, and the angles at Mn3 and Mn4 that include O5. An empty list adds no angle restraints.

`SITE_ANGLE_SIGMA` (default 25) is in degrees. `refmac5` reads that sigma in degrees, so the number is written unchanged.

`USE_LIBCHECK` chooses where the distance restraints go.

- `False` (the default) passes them to `refmac5` as external distances. REFMAC then uses that restraint instead of the dictionary bond for the same atom pair, and only in the residues listed in `SITE_COPIES`.
- `True` runs `libcheck` and widens those distances in a monomer library. The library applies to every copy of that residue type. The written file has the same layout as a CCP4 monomer dictionary. If a residue’s coordinates are not connected well enough for `libcheck` (this happens for the lipid LMG), that residue is copied from the CCP4 monomer library instead.

Angle restraints stay external either way. `LIB_MONOMERS = None` sends the site residue type and every other non-water hetero residue through `libcheck`. A list of residue names limits that set.

## Files

| Role | Notebook variable | Web form | Required for the map |
| --- | --- | --- | --- |
| Omit coordinates | `DARK_PDB` | Omit-model coordinates | Yes. This is the model REFMAC reads. |
| Observed structure factors | `LIGHT_SF` | Observed structure factors | Yes. Intensities `I`, `SIGI`, as deposited for `6w1v`. |
| Light coordinates | `LIGHT_PDB` | Light coordinates | Used for the residue list. The figure is drawn on the refined omit model. |

The omit model’s own structure-factor file is not an input. $F_o$ comes from the light dataset. `DARK_ID` and `LIGHT_ID` are labels in output names (`6w1v_omitO6_single_confB`, `6w1v`, `11yo`, …). They do not have to be PDB codes.

## Steps

The programs run in this order: `cif2mtz`, `ctruncate`, `cad` (to copy `FREE`), `refmac5`, `fft` of `DELFWT`, `extends`, Coot, and PyMOL.

A map left by an older `cad` / `scaleit` pipeline has a different file name and is not reused. Set `RECOMPUTE_MAP` to overwrite the extended map for this omit model, dataset, and resolution tag.

### 1. Check the programs and the files

The programs are `cif2mtz`, `ctruncate`, `cad`, `refmac5`, `fft`, `extends`, and `mtzdump`. `libcheck` is required only when `USE_LIBCHECK` is on. Coot comes from the CCP4 tree. PyMOL must be on `PATH`.

### 2. Inspect the omit model and choose the site

The residue list is built from the omit coordinates. For photosystem II, look for `OEX` or `OEY`. The workflow checks that `SITE_CHAIN` / `SITE_RESI` exists in the omit model, prints the atom names, and prints the centroid. The log also prints the cycle count, so a zero-cycle job is visible before REFMAC starts.

On the web form, **List hetero residues** reads the omit-model coordinates.

### 3. Build the mFo−DFc map

**3a. Convert the observed mmCIF.** `cif2mtz` writes `{light}-sf.mtz`. The columns at this stage are the deposited intensities.

**3b. Choose the resolution shell.** With `RESOLUTION_SOURCE = "sf"`, `LOWRES` and `DMIN` are the resolution range of that MTZ. With `"explicit"`, the numbers in the configuration are used. The notebook’s explicit fallback is 30.0–1.89 Å. `DMIN` should be no finer than the high-resolution limit of the observed data.

**3c. Intensities to amplitudes.** `ctruncate` reads `I` and `SIGI` and writes `{light}-sf-truncate.mtz` with amplitudes `F` and `SIGF`:

```text
ctruncate -hklin {light}-sf.mtz -hklout {light}-sf-truncate.mtz -colin '/*/*/[I,SIGI]'
```

**3d. Copy the free-R flag.** `cad` takes all columns from the amplitude MTZ and the `FREE` column from the intensity MTZ, and writes `{light}-sf-truncate-free.mtz`:

```text
LABIN FILE 1 ALL
LABIN FILE 2 E1=FREE
```

**3e. Refine, or just calculate, with REFMAC.** The job reads the amplitude MTZ and the omit coordinates and writes `{dark}_vs_{light}.mtz` and `{dark}_refined.pdb`.

External restraints, when the lists are non-empty, are written to `{dark}_external_restraints.txt` and placed in front of the REFMAC keywords. With `USE_LIBCHECK`, the distance terms are omitted from that file and a `{dark}_lib.cif` library is passed as `libin` instead. The keywords always include

```text
NCYC <cycles>
MAKE HYDR Y
MAKE HOUT Y
```

`MAKE HYDR` / `MAKE HOUT` tell REFMAC to generate hydrogens and to write them to the coordinate file. Zero cycles still produce a coordinate file; the log calls it the refined model even when nothing moved.

The MTZ is checked for `FP` (or `F`, `FOBS`), `FC_ALL` (or `FC`), `DELFWT`, and `PHDELWT`. The transform uses the weighted coefficient, not a raw $F_o-F_c$.

**3f. Fourier transform.** `fft` writes `{light}_DELFWT_{dark}_{tag}.ccp4`:

```text
RESOLUTION <low> <high>
LABIN F1=DELFWT PHI=PHDELWT
```

**3g. Cover the refined model.** `extends` pads the map by `MAP_BORDER_A` (default 5 Å) around `{dark}_refined.pdb` and writes the extended CCP4 map.

The log states the interpretation: positive density is present in the light observations and absent from the omit model.

### 4. Export the site fragment

Coot loads the extended map as a difference map and the refined coordinates, centres on each residue in `SITE_COPIES`, and writes a clip of radius `EXPORT_RADIUS_A` (default 10 Å). As in the other workflows, the Coot contour does not change the exported voxels.

### 5. Contour and render

PyMOL loads the refined coordinates and the fragment. Positive and negative levels default to the same list; the configuration cell sets that list to `3.5`. Colours follow the Fo−Fo figure.

| Sign | Lowest contour | Middle | Highest contour |
| --- | --- | --- | --- |
| Positive, in the data and missing from the model | forest | cyan | blue |
| Negative, in the model and missing from the data | red | salmon | magenta |

Cluster bonds in this workflow use the Mn/Ca–oxygen pairs without O6, matching an omit model from which O6 was removed. Bonds are added only when the atom names exist.

Outputs:

| File | Contents |
| --- | --- |
| `{light}_DELFWT_{dark}_{site}_{chain}{resi}_diffmap.png` | The figure |
| `pymol_{light}_DELFWT_{dark}_{site}.pml` | The scene |
| `{dark}_refined.pdb` | Coordinates after `REFMAC_NCYC` cycles |
| `*_extended.ccp4` | The map |

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| `REFMAC_NCYC` | `0` leaves atoms fixed. A larger integer refines against the observed amplitudes |
| Distance restraints and `SITE_DISTANCE_SIGMA` | Which site distances are loose, and how loose. 0.50 Å is far wider than a dictionary bond |
| Angle restraints and `SITE_ANGLE_SIGMA` | External angles, sigma in degrees |
| `USE_LIBCHECK` | External distances on the listed copies, or a widened monomer library for every copy of that residue type |
| Resolution source, `LOWRES`, `DMIN` | Shell of the `DELFWT` transform |
| `SIGMA_LEVELS` | The configuration cell uses a single $3.5\sigma$ mesh |

## Reading the default figure

The model in the figure is the REFMAC coordinate set. With zero cycles it is the omit model rewritten by REFMAC, hydrogens included, atoms unmoved. Positive density is electron density in the 2F observations that this single-conformation model does not include. A positive peak on an atom present in `6w1v` and absent from the omit model is the expected $mF_o-DF_c$ signal.

A strong positive peak sitting exactly on O6 is not required. Whether the 2F data support an inserted oxo or water there is debated. The workflow reproduces the map.

If the mesh is empty, lower the contour (try `3`). If the cage swallows the model, decrease `CARVE_A` or raise the contour.

## Adapting it

1. Provide the omit coordinates, the observed structure factors, and the light coordinates used for the residue list.
2. Set the two labels. Leave the resolution source on `"sf"`, or set `LOWRES` and `DMIN` yourself.
3. Run the residue list, then set the chain, residue, and neighbours.
4. Set the cycle count and the restraint lists. Leave `USE_LIBCHECK` off to pass distances straight to REFMAC. Turn it on to widen them in a `libcheck` library. An empty distance or angle list adds no restraints of that kind.
5. Recompute the map the first time, and again after any change of cycles, restraints, or resolution. Leave recomputation off while you tune σ and the view.

The extended map can be opened with the refined coordinates:

```text
coot {extended_map} {dark}_refined.pdb
```

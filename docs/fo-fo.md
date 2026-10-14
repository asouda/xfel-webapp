# Isomorphous Fo−Fo difference map

This workflow turns a pair of isomorphous crystal structures — a dark or reference state and a light or activated state — into a signed $|F_{\mathrm{LIGHT}}|-|F_{\mathrm{DARK}}|$ difference map, clips a fragment around a site, and draws the $\pm n\sigma$ surfaces in PyMOL.

The notebook is `notebooks/XFEL_fo-fo_diffmap.ipynb`. The web application runs the same steps from the Fo−Fo page.

## Worked example

The default pair is photosystem II from Ibrahim *et al.*, Science, 2020.

| Role | PDB | State | What it supplies |
| --- | --- | --- | --- |
| Dark | `6w1o` | 0F, dark-stable $S_1$ | Reference amplitudes and coordinates |
| Light | `6w1v` | 2F, 200 ms | Activated amplitudes. The figure is drawn on this model. Two oxygen-evolving-complex conformations are present. |

The same steps work for any other isomorphous pair: 1F against 0F, 2F against 1F, ligand-bound against apo. Change the identifiers, the files, the site, and the contour levels.

On `6w1v`, chain `A` residue **416** is ligand `OEY`, the major (75 %) 2F cluster conformation, which includes O6. Residue **415** is `OEX`, the minor conformation without O6. On `6w1o` the unique cluster is `OEX` at residue **414**. Those residue numbers belong to different entries. The residue list in step 2 is how you confirm which one you centred.

## Background

A time-resolved measurement gives two amplitude sets. The chemical signal is the difference electron density

$$
\Delta\rho(\mathbf{r}) \;\propto\; \sum_{hkl} \big(|F_{\mathrm{LIGHT}}|-|F_{\mathrm{DARK}}|\big)\, e^{i\varphi_{\mathrm{ref}}}\, e^{-2\pi i \mathbf{h}\cdot\mathbf{r}}.
$$

The phases $\varphi_{\mathrm{ref}}$ come from a complete model, or from a $2F_o-F_c$ map, rather than from the unknown phase of the difference vector. The sum is a useful picture of a local change when the two crystals share a space group and a very similar cell. A large domain motion or a different cell breaks the approximation, because the reference phases then describe a different structure.

### How to read the sign

- Positive density, coloured forest, then cyan, then blue as $|\sigma|$ increases, means more electrons in the light dataset. An atom appeared there, or moved there.
- Negative density, coloured red, then salmon, then magenta, means more electrons in the dark dataset. An atom left that position.
- A real shift is a negative peak on the old position paired with a positive peak on the new position.

### Why the amplitudes are scaled first

Subtracting raw amplitudes mixes the structural signal with differences of scale, overall $B$ factor, and completeness. `scaleit` applies an isotropic scale and $B$ so that the light amplitudes sit on the dark reference inside the chosen resolution shell. The log reports an overall $R$ factor. Values around 0.10–0.15 are typical for a usable isomorphous XFEL pair. A much larger $R$ usually means non-isomorphism or a data-quality problem. It does not, by itself, mean a larger chemical change.

### Which phases

A textbook difference-Fourier map uses phases from the reference (dark) model. This tutorial follows the published photosystem II comparison and, by default, takes `PHWT` from the **light** validation $2F_o-F_c$ coefficients, so the map uses the same phase convention as that figure. `PHASE_SOURCE` can switch to the dark validation file, or to `PHIC` / `PHWT` from either zero-cycle REFMAC job.

The validation coefficient file is a phase source. It is converted to MTZ so `cad` can copy one phase column. It is not the difference map.

## Files

Six files, three per structure. Paths left unset are built from the identifier inside the input directory: `{id}.pdb`, `{id}-sf.cif`, and `{id}_validation_2fo-fc_map_coef.cif`.

| Role | Notebook variable | Web form |
| --- | --- | --- |
| Dark coordinates | `DARK_PDB` | Dark coordinates |
| Light coordinates | `LIGHT_PDB` | Light coordinates |
| Dark structure factors | `DARK_SF` | Dark structure factors |
| Light structure factors | `LIGHT_SF` | Light structure factors |
| Dark validation coefficients | `DARK_VALIDATION` | Dark 2Fo−Fc coefficients |
| Light validation coefficients | `LIGHT_VALIDATION` | Light 2Fo−Fc coefficients |

`DARK_ID` and `LIGHT_ID` are short labels in output file names. They do not have to be PDB codes. Both `.pdb` and mmCIF coordinates are accepted.

## Steps

The programs run in this order: `cif2mtz`, then `refmac5` for zero cycles, `cad`, `scaleit`, `fft`, `extends`, Coot, and PyMOL.

Each CCP4 program writes a log in the output directory. Set `RECOMPUTE_MAP` to rebuild the extended map. While that file is kept, later runs reuse it and only redraw the figure.

### 1. Check the programs and the files

CCP4 is expected at `/Applications/ccp4-9` (`CCP4_ROOT` in the notebook, the CCP4 directory on the form). The workflow sources `bin/ccp4.setup-sh`, then looks for `cif2mtz`, `refmac5`, `cad`, `scaleit`, `fft`, `extends`, and `mtzdump`. Coot comes from that CCP4 tree. PyMOL must be on `PATH`.

The six input files must already exist. An empty or missing file stops the job before any map is written.

### 2. Inspect the light model and choose the site

The notebook lists non-water hetero residues in the light coordinates, skipping common lipids and pigments, and prints chain, residue number, residue name, atom count, and occupancy. For photosystem II, look for `OEX` or `OEY`. For another system, look for the chromophore, the ligand, or the metal site.

If the table is empty, the group may be stored as `ATOM` records. Open the file in a text editor or in Coot and set `SITE_CHAIN` and `SITE_RESI` from that residue.

The workflow then checks that the chosen residue exists, lists its atom names, and prints the centroid that Coot will centre on. The web form has a **List hetero residues** button that fills the chain and residue from this scan.

### 3. Build the difference map

**3a. Convert reflections to MTZ.** `cif2mtz` writes `{id}-sf.mtz` for both structure-factor files and `{id}_validation_2fo-fc_map_coef.mtz` for both coefficient files.

**3b. Calculate model terms without moving atoms.** `refmac5` runs on each structure with

```text
NCYC 0
```

and writes `{id}_calc.mtz` plus `{id}_checked.pdb`. Zero cycles put `FP` and `SIGFP` into a consistent REFMAC-style MTZ. The deposited coordinates are not refined.

**3c. Read the resolution shell.** With `RESOLUTION_SOURCE = "sf"`, `mtzdump` supplies the range of each structure-factor MTZ, and the job keeps the common shell. With `"explicit"`, `LOWRES` and `DMIN` are used as entered. The default explicit pair in the notebook is 20.0–2.10 Å, used only when the source is explicit.

**3d. Combine columns.** `cad` copies the light amplitude, the dark amplitude, and one phase into `cad_{light}_{dark}.mtz`. Amplitudes are labelled `F_LIGHT`, `SIGF_LIGHT`, `F_DARK`, and `SIGF_DARK`. The phase is labelled `PHI_REF`. The phase file and column depend on `PHASE_SOURCE`:

| `PHASE_SOURCE` | File | Column candidates |
| --- | --- | --- |
| `light_validation` (default) | Light validation MTZ | `PHWT`, `PHIC`, `PHI` |
| `dark_validation` | Dark validation MTZ | `PHWT`, `PHIC`, `PHI` |
| `light_refmac` | Light zero-cycle MTZ | `PHWT`, `PHIC`, `PHI` |
| `dark_refmac` | Dark zero-cycle MTZ | `PHWT`, `PHIC`, `PHI` |

**3e. Scale light onto dark.** `scaleit` reads the combined MTZ and writes `scaleit_{light}_{dark}.mtz`:

```text
RESOLUTION <low> <high>
REFINE ISOTROPIC
LABIN FP=F_DARK SIGFP=SIGF_DARK FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
```

Dark is the reference (`FP`). Light is treated as the derivative (`FPH1`), which is the SCALEIT convention for a second amplitude set. The log line to read is the overall $R$ factor.

**3f. Transform the difference.** `fft` writes `{light}_minus_{dark}_{tag}.ccp4`:

```text
RESOLUTION <low> <high>
LABIN F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_REF
```

`F1` minus `F2`, with `PHI`, is the signed difference Fourier.

**3g. Cover the light model.** `extends` pads that map by `MAP_BORDER_A` (default 5 Å) around the light coordinates and writes the `*_extended.ccp4` file. Later steps use the extended map.

The log states the sign convention explicitly: positive density is stronger in the light identifier, negative density is stronger in the dark identifier.

### 4. Export the site fragment

Coot loads the extended map as a difference map and the light coordinates, centres on the centroid of every `(chain, residue)` pair in `SITE_COPIES`, and writes a clip of radius `EXPORT_RADIUS_A` (default 10 Å). One fragment is written per copy. The figure uses the first.

The σ value passed into the Coot script only affects the on-screen contour inside Coot. The voxel values in the `.map` file are the clipped map.

### 5. Contour and render

PyMOL loads the light coordinates and the fragment. For each positive level it builds a surface at $+n\sigma$, and for each negative level a surface at $-n\sigma$. The default levels are 3, 4, and 5.

| Sign | Lowest contour | Middle | Highest contour |
| --- | --- | --- | --- |
| Positive, light stronger than dark | forest | cyan | blue |
| Negative, dark stronger than light | red | salmon | magenta |

A fourth level, if you add one, uses `density` and `hotpink`. `SIGMA_POSITIVE` and `SIGMA_NEGATIVE` override the shared list when you want, for example, $+5\sigma$ and $-3\sigma$.

`CARVE_A` (default 3 Å) keeps each surface near the displayed atoms. The site residue and `EXTRA_RESIDUES` are shown; the default extra residue is D1-Glu189 (`A 189`). Oxygen-evolving-complex bonds, including O6 when that atom exists, are drawn only if the atom names are present.

Outputs:

| File | Contents |
| --- | --- |
| `{light}_minus_{dark}_{site}_{chain}{resi}_diffmap.png` | The figure |
| `pymol_{light}_minus_{dark}_{site}.pml` | The scene |
| `*_extended.ccp4` | The map you can reopen in Coot |

In the notebook, `INSPECT_IN_PYMOL` opens the scene. Quit PyMOL when the camera is right; the notebook prints a view matrix you can paste into `PYMOL_VIEW`. The web application shows the same fragment in the browser and can save a PNG of that view.

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| `PHASE_SOURCE` | Switches between light validation phases (the published-figure default), dark validation phases (the classical reference), and either zero-cycle REFMAC phase |
| Resolution source, `LOWRES`, `DMIN` | Shell used by `scaleit` and `fft`. The high-resolution cutoff should be no finer than the worse dataset |
| `SIGMA_LEVELS` | Nested contours. `[4.0]` is a cleaner figure than `[3, 4, 5]` |
| `MAP_BORDER_A`, clip radius, carve | How much map is stored, how much is clipped, and how tightly the mesh hugs the atoms |
| `RECOMPUTE_MAP` | `False` reuses the extended map so contour and camera edits stay cheap |

The full list, including figure size and ray tracing, is on the [parameters](parameters.md) page. Shared contour and clipping behaviour is on the [background](concepts.md) page.

## Reading the default photosystem II figure

Around `OEY` A416 in `6w1v` minus `6w1o` you should see the Mn$_4$Ca cluster and the dangling Mn4, with O6 present in the major light conformation, and difference peaks on the cluster and on D1-Glu189. Glu189 is modelled with two alternate conformations in the 2F structure.

A strong positive peak sitting exactly on O6 is not required by this workflow. Whether the 2F data support an inserted oxo or water at that position is debated in the literature. The notebook reproduces the map.

If the mesh is empty, lower the contour (try `3` alone). If the cage swallows the model, decrease `CARVE_A` or raise the contour.

## Adapting it

1. Put the six files in the input directory, or point the six path variables at them.
2. Set `DARK_ID`, `LIGHT_ID`, and `PHASE_SOURCE`. Leave the resolution source on `"sf"` for the common shell, or switch to `"explicit"` and set `LOWRES` and `DMIN`.
3. Run the residue list, then set the chain, residue, and neighbouring residues.
4. Clear any saved camera unless you intend to reuse it.
5. Recompute the map the first time. Leave recomputation off while you tune σ and the view.

You can also open the extended map directly:

```text
coot {extended_map} {light_pdb}
```

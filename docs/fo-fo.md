# Isomorphous Fo−Fo difference map

This workflow turns a pair of isomorphous crystal structures — a dark or reference state and a light or activated state — into a signed \(|F_{\mathrm{LIGHT}}|-|F_{\mathrm{DARK}}|\) difference map, clips a fragment around a site, and draws the \(\pm n\sigma\) surfaces.

The calculation is \(|F|\) minus \(|F|\) with one borrowed phase. The [theoretical background](theory.md) explains where those amplitudes and phases come from. The steps below are the CCP4 programs that evaluate the sum for the worked example.

## Worked example

Fig. 1a of Wang, Armstrong, and Batista, [Proc. Natl. Acad. Sci. U.S.A. 118, e2023982118 (2021)](https://doi.org/10.1073/pnas.2023982118), is the isomorphous difference \(F_o(6\mathrm{w1v}/2\mathrm{F}) - F_o(6\mathrm{w1o}/0\mathrm{F})\), drawn on the two-flash model. The structures are the photosystem II series of Ibrahim *et al.*, [Proc. Natl. Acad. Sci. U.S.A. 117, 12624 (2020)](https://doi.org/10.1073/pnas.2000529117): `6w1o` is the dark-stable \(S_1\) state (0F) and `6w1v` is the 200 ms two-flash state (2F).

| Role | What the form needs | Download |
| --- | --- | --- |
| Dark coordinates | `6w1o` model. The oxygen-evolving complex in this entry is `OEX` A414. | [bundle 1](https://github.com/asouda/xfel-webapp/raw/main/data/6w1o/6w1o-pdb-bundle1.pdb), [bundle 2](https://github.com/asouda/xfel-webapp/raw/main/data/6w1o/6w1o-pdb-bundle2.pdb) |
| Light coordinates | `6w1v` model. The figure is drawn on this structure. | [bundle 1](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle1.pdb), [bundle 2](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle2.pdb) |
| Dark structure factors | `6w1o-sf.cif` | [6w1o-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/6w1o/6w1o-sf.cif.gz) |
| Light structure factors | `6w1v-sf.cif` | [6w1v-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-sf.cif.gz) |
| Dark validation coefficients | Phase source only if you switch away from the light file | [6w1o validation 2Fo−Fc](https://github.com/asouda/xfel-webapp/raw/main/data/6w1o/6w1o_validation_2fo-fc_map_coef.cif.gz) |
| Light validation coefficients | Default phase source, column `PHWT` | [6w1v validation 2Fo−Fc](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v_validation_2fo-fc_map_coef.cif.gz) |

`6w1o` and `6w1v` are larger than one classic PDB file, so the repository stores each coordinate set as the two RCSB bundles. Bundle 1 contains the oxygen-evolving complex used for this figure. The form takes one coordinate file per structure: point it at bundle 1, or at a single merged `6w1o.pdb` / `6w1v.pdb` of the same entry if you have one. The mmCIF files are gzip-compressed. Give the form the `.gz` path; a name that ends in `.gz` is decompressed with `gunzip` before `cif2mtz` reads it.

On `6w1v`, chain `A` residue **416** is ligand `OEY`, the major (75%) 2F cluster, which includes O6. Residue **415** is `OEX`, the minor conformation without O6. On `6w1o` the cluster is `OEX` at residue **414**. Set the site to chain `A`, residue `416`, and list D1-Glu189 (`A 189`) as a neighbouring residue. The published figure is contoured at \(+5\), \(+4\), and \(+3\sigma\) (blue, cyan, green) and at \(-4\) and \(-3\sigma\) (red, salmon). The form’s default PNG levels are \(3\), \(4\), and \(5\) for both signs.

The same programs work for any other isomorphous pair. Change the six files, the site, and the contour levels.

### How the Fig. 1a map is calculated

The programs run in this order: `cif2mtz`, `refmac5` for zero cycles, `mtzdump`, `cad`, `scaleit`, `fft`, `extends`, Coot, and PyMOL.

**1. Check the programs and the files.** CCP4 is expected at `/Applications/ccp4-9` unless the form names another suite prefix. The job sources `bin/ccp4.setup-sh` and looks for `cif2mtz`, `refmac5`, `cad`, `scaleit`, `fft`, `extends`, and `mtzdump`. Coot comes from that CCP4 tree. PyMOL must be on `PATH`. An empty or missing file stops the job before a map is written.

**2. Inspect the light model and choose the site.** The **List hetero residues** button reads the light coordinates, skips water and common lipids and pigments, and lists chain, residue number, residue name, atom count, and occupancy. For this example, confirm `OEY` A416, then keep that chain and residue. The job checks that the residue exists, lists its atom names, and prints the centroid Coot will centre on.

**3. Convert reflections to MTZ.** `cif2mtz` reads each structure-factor mmCIF and each validation mmCIF and writes `{id}-sf.mtz` and `{id}_validation_2fo-fc_map_coef.mtz`. For `6w1o` and `6w1v` the structure-factor loop holds intensities (`intensity_meas`, `intensity_sigma`). `cif2mtz` stores them as `I` and `SIGI`. The validation loop holds `FWT`, `PHWT`, and `FOM`. Nothing is scaled or transformed yet; the MTZ is the same list of reflections in the binary layout the later programs read.

**4. Place the amplitudes in a REFMAC MTZ, without moving atoms.** `refmac5` runs once on `6w1o` and once on `6w1v` with

```text
NCYC 0
```

and writes `{id}_calc.mtz` plus `{id}_checked.pdb`. Zero cycles compute \(F_c\) from the deposited coordinates and write the observed amplitudes as `FP` and `SIGFP`. The coordinates used for the figure stay at the deposited positions. This step is bookkeeping of the amplitude columns. The difference phase comes from the next step, not from this \(F_c\), unless the phase source is switched to a REFMAC phase.

**5. Read the resolution shell.** With the resolution source set to the structure-factor files, `mtzdump` reports the range of each `{id}-sf.mtz`, and the job keeps the common shell: the poorer low-resolution limit and the poorer high-resolution limit. With an explicit shell, the two cutoffs on the form are used instead. The explicit pair stored for this workflow is 20.0–2.10 Å.

**6. Combine the two amplitudes and one phase.** `cad` copies four amplitude columns and one phase into `cad_6w1v_6w1o.mtz`:

| Output label | Taken from |
| --- | --- |
| `F_LIGHT`, `SIGF_LIGHT` | Light `FP`, `SIGFP` |
| `F_DARK`, `SIGF_DARK` | Dark `FP`, `SIGFP` |
| `PHI_REF` | `PHWT` in the light validation MTZ, by default |

The phase file follows the phase-source field:

| Phase source | File | Column candidates |
| --- | --- | --- |
| Light validation (default) | `6w1v` validation MTZ | `PHWT`, `PHIC`, `PHI` |
| Dark validation | `6w1o` validation MTZ | `PHWT`, `PHIC`, `PHI` |
| Light REFMAC | `6w1v_calc.mtz` | `PHWT`, `PHIC`, `PHI` |
| Dark REFMAC | `6w1o_calc.mtz` | `PHWT`, `PHIC`, `PHI` |

`cad` does not subtract anything. It builds one reflection list so the next two programs see both datasets and the phase on the same \(hkl\).

**7. Scale the light amplitudes onto the dark amplitudes.** `scaleit` refines an isotropic scale and an isotropic \(B\) factor that put \(|F_{\mathrm{light}}|\) onto the \(|F_{\mathrm{dark}}|\) reference inside the resolution shell:

```text
RESOLUTION <low> <high>
REFINE ISOTROPIC
LABIN FP=F_DARK SIGFP=SIGF_DARK FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT
```

Dark is `FP`, the reference. Light is `FPH1`, the second amplitude set in the SCALEIT convention. The output is `scaleit_6w1v_6w1o.mtz`, with the light amplitudes replaced by their scaled values. The log line to read is the overall \(R\) factor. Values around 0.10–0.15 are typical for a usable isomorphous XFEL pair. A much larger \(R\) usually means the crystals differ by more than a local change, or that one dataset is poor.

**8. Fourier-transform the difference.** `fft` sums

\[
\Delta\rho(\mathbf{r}) \propto \sum_{hkl} \big(|F_{\mathrm{LIGHT}}| - |F_{\mathrm{DARK}}|\big)\, e^{i\varphi_{\mathrm{ref}}}\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}}
\]

over the same resolution shell and writes `6w1v_minus_6w1o_{tag}.ccp4`:

```text
RESOLUTION <low> <high>
LABIN F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_REF
```

`F1` minus `F2`, with phase `PHI`, is the signed difference Fourier. Positive density is stronger in `6w1v`. Negative density is stronger in `6w1o`. The map at this stage covers the asymmetric unit.

**9. Cover the light model.** `extends` repeats the asymmetric unit by the space-group symmetry and pads the grid by the map border (default 5 Å) around the light coordinates, writing `*_extended.ccp4`. Later steps use that extended map.

**10. Export the site fragment.** Headless Coot loads the extended map as a difference map and the light coordinates, centres on the centroid of `OEY` A416 (and on any extra copies listed on the form), and writes a sphere of radius equal to the clip radius (default 10 Å). The σ value passed into the Coot script affects only the contour drawn inside that Coot session. The voxels in the `.map` file are the clipped map.

**11. Contour and render.** PyMOL loads the light coordinates and the fragment. For each positive PNG level it builds a surface at \(+n\sigma\), and for each negative level a surface at \(-n\sigma\).

| Sign | Lowest contour | Middle | Highest contour |
| --- | --- | --- | --- |
| Positive, light stronger than dark | forest | cyan | blue |
| Negative, dark stronger than light | red | salmon | magenta |

A fourth level uses the colours `density` and `hotpink`. Separate positive and negative lists override the shared PNG list, so the figure can be \(+5\sigma\) and \(-3\sigma\), which is close to the published Fig. 1a colouring. The carve radius (default 3 Å) keeps each surface near the displayed atoms. `OEY` A416 and D1-Glu189 are shown; the rest of the model stays hidden. Oxygen-evolving-complex bonds, including O6 when that atom exists, are drawn when the atom names are present.

The browser viewer then loads the same fragment and draws one isosurface. The slider and the number box move that single contour between the viewer sigma minimum and maximum (default \(-5\) to \(5\)). Blue/red colouring paints positive density blue and negative density red, darker farther from zero. Temperature colouring paints the low end of the range cold and the high end hot. The saved PNG remains the PyMOL rendering at the PNG sigma levels.

Outputs:

| File | Contents |
| --- | --- |
| `6w1v_minus_6w1o_{site}_{chain}{resi}_diffmap.png` | The PyMOL figure |
| `pymol_6w1v_minus_6w1o_{site}.pml` | The scene |
| `*_extended.ccp4` | The map that can be reopened in Coot |

While an extended map of the same name is kept, a later run can reuse it and only redraw the figure. Turn recomputation on after changing the files, the phase source, or the resolution shell.

## Background

A time-resolved measurement gives two amplitude sets. The chemical signal is the difference electron density in step 8 above. The phases \(\varphi_{\mathrm{ref}}\) come from a refined model, through the validation coefficient file or through a zero-cycle REFMAC job, rather than from the unknown phase of the difference vector \(F_{\mathrm{light}}-F_{\mathrm{dark}}\). The sum is a useful picture of a local change when the two crystals share a space group and a very similar cell. A large domain motion or a different cell breaks the approximation, because the reference phases then describe a different structure.

Subtracting raw amplitudes mixes the structural signal with differences of scale, overall \(B\) factor, and completeness. That is what `scaleit` removes before `fft`.

### How to read the sign

- Positive density, coloured forest, then cyan, then blue as \(|\sigma|\) increases, means more electrons in the light dataset. An atom appeared there, or moved there.
- Negative density, coloured red, then salmon, then magenta, means more electrons in the dark dataset. An atom left that position.
- A real shift is a negative peak on the old position paired with a positive peak on the new position.

The validation coefficient file supplies `PHI_REF`. The difference map is the `fft` output, not that file.

## Files

Six files, three per structure. The labels on the form are short names used in the output files. They do not have to be PDB codes. Both `.pdb` and mmCIF coordinates are accepted.

| Role | Form field |
| --- | --- |
| Dark coordinates | Dark coordinates |
| Light coordinates | Light coordinates |
| Dark structure factors | Dark structure factors |
| Light structure factors | Light structure factors |
| Dark validation coefficients | Dark 2Fo−Fc coefficients |
| Light validation coefficients | Light 2Fo−Fc coefficients |

## Parameters that change this map

| Choice | Effect |
| --- | --- |
| Phase source | Light validation phases (the Fig. 1a default), dark validation phases (the classical reference), or either zero-cycle REFMAC phase |
| Resolution source and the two cutoffs | Shell used by `scaleit` and `fft`. The high-resolution cutoff should be no finer than the worse dataset |
| PNG sigma levels | Nested contours in the PyMOL figure. A single level is a cleaner figure than \(3\), \(4\), and \(5\) |
| Viewer sigma minimum and maximum | Range of the one isosurface in the browser. The PNG levels are independent |
| Map border, clip radius, carve | How much map is stored, how much is clipped, and how tightly the mesh hugs the atoms |
| Recompute the map | Off reuses the extended map so contour and camera edits stay cheap |

The full list is on the [parameters](parameters.md) page. Shared contour and clipping behaviour is on the [background](concepts.md) page.

## Reading the photosystem II figure

Around `OEY` A416 in `6w1v` minus `6w1o` you should see the Mn\(_4\)Ca cluster and the dangling Mn4, with O6 present in the major light conformation, and difference peaks on the cluster and on D1-Glu189. Glu189 is modelled with two alternate conformations in the 2F structure. Wang, Armstrong, and Batista contour this map at \(+3\) to \(+5\sigma\) and \(-3\) to \(-4\sigma\) and discuss the peak nearest the proposed oxo position.

If the mesh is empty, lower the contour (try `3` alone). If the cage swallows the model, decrease the carve radius or raise the contour.

## Adapting it

1. Point the six file fields at the new pair, or upload the files. A `.gz` name is decompressed automatically.
2. Set the two labels and the phase source. Leave the resolution source on the structure-factor files for the common shell, or switch to explicit cutoffs.
3. List hetero residues, then set the chain, residue, and neighbouring residues.
4. Recompute the map the first time. Leave recomputation off while you tune σ and the view.

The extended map can also be opened directly:

```text
coot {extended_map} {light_coordinates}
```

# Electron density maps from XFEL structure factors

This site explains the map workflows in the `xfel-webapp` web application. Each workflow collects coordinates and reflection files, runs CCP4, Coot, and PyMOL, and draws an isodensity map on the structure. The page keeps an interactive 3D view of that map and writes a PNG of the PyMOL figure around a site you choose. The reconstruction repeats the omit map and places one missing atom on each peak.

The [theoretical background](theory.md) covers structure factors, amplitudes, model phases, the Fourier transform, refinement cycles, the validation coefficient files, and the peak search that rebuilds an omit model. The pages linked below are the formula for each map, the files it needs, and the program sequence that builds it.

Three of the calculations are written out as Jupyter notebooks in this repository. The reconstruction has no separate notebook: it repeats the omit-map calculation.

- [Fo−Fo difference map](https://github.com/asouda/xfel-webapp/raw/main/notebooks/XFEL_fo-fo_diffmap.ipynb) (`notebooks/XFEL_fo-fo_diffmap.ipynb`)
- [Omit-model mFo−DFc map](https://github.com/asouda/xfel-webapp/raw/main/notebooks/XFEL_fo-fc_diffmap_refinement.ipynb) (`notebooks/XFEL_fo-fc_diffmap_refinement.ipynb`)
- [Sharpened 2mFo−DFc map](https://github.com/asouda/xfel-webapp/raw/main/notebooks/XFEL_2fo-fc_map.ipynb) (`notebooks/XFEL_2fo-fc_map.ipynb`)

## Choose a map

<div class="grid cards" markdown>

-   **Isomorphous Fo−Fo**

    ---

    Signed difference density between a dark reference dataset and a light, activated dataset. Positive density is stronger in the light data. Negative density is stronger in the dark data.

    [:octicons-arrow-right-24: Theory and steps](fo-fo.md)

-   **Omit-model mFo−DFc**

    ---

    Weighted difference between one set of observed amplitudes and an omit model. Positive density is in the observations and missing from the model. REFMAC can leave the coordinates fixed or refine them.

    [:octicons-arrow-right-24: Theory and steps](omit-map.md)

-   **Reconstruct from omit maps**

    ---

    The omit-model \(mF_o-DF_c\) map, repeated. Each iteration refines the omit model, places the missing reference atom closest to the highest peak into those refined coordinates, and the next map starts from that model.

    [:octicons-arrow-right-24: Theory and steps](reconstruction.md)

-   **Sharpened 2mFo−DFc**

    ---

    Electron density of a single structure, from the deposited validation coefficients, sharpened by a \(B\)-factor shift and drawn on the model.

    [:octicons-arrow-right-24: Theory and steps](two-fo-fc.md)

</div>

## What each map answers

| Question | Map | Observations | Model |
| --- | --- | --- | --- |
| Where did atoms move between two isomorphous states? | Fo−Fo | Light and dark amplitudes, scaled together | Phases from one refined model; atoms are not moved |
| What is in the data that this model leaves out? | mFo−DFc | One amplitude set | Omit coordinates, optionally refined against those amplitudes |
| Which missing reference atom sits on the highest omit-map peak? | mFo−DFc, one atom at a time | The same amplitude set, recalculated after each insertion | The omit model, with each new atom placed at the peak |
| Where is the electron density of this structure? | 2mFo−DFc | Validation coefficients of that structure | The same structure, superimposed on the map |

A chemically meaningful Fo−Fo feature is a negative peak on the old position next to a positive peak on the new one. An omit-map feature for an atom left out of the model is a positive peak. The reconstruction treats the strongest of those peaks as the next atom and repeats the map. A 2mFo−DFc mesh that wraps an atom is accounted-for density. Difference maps are usually inspected at 3–5σ, while a 2mFo−DFc map is usually inspected nearer 1–2σ of its own root-mean-square.

Phases, weighting, σ contours, the RCSB files, and the Coot, PyMOL, and browser-viewer tail shared by every workflow are on the [background](concepts.md) page. Form fields are on the [parameters](parameters.md) page.

## Worked examples

The examples are photosystem II, arranged to follow Wang, Armstrong, and Batista, [Proc. Natl. Acad. Sci. U.S.A. 118, e2023982118 (2021)](https://doi.org/10.1073/pnas.2023982118). Download links for the input files are on each workflow page. Reflection files in the repository are gzip-compressed; the application decompresses a file whose name ends in `.gz`.

| Workflow | Figure | Inputs | Site |
| --- | --- | --- | --- |
| Fo−Fo | Fig. 1a, \(F_o(6\mathrm{w1v})-F_o(6\mathrm{w1o})\) | `6w1v` and `6w1o` coordinates, both structure-factor mmCIF files, and both validation coefficient files | `OEY` A416 on the 2F model, plus D1-Glu189 |
| mFo−DFc | Fig. 2a and Fig. 2b | `6w1v` coordinates, `6w1v-sf.cif`, and the omit model `6w1v_omitO6_single_confA.pdb` | `OEX` A415 in the omit model, the single conformation with O6 absent |
| Reconstruction | The same omit comparison, one atom at a time | `6w1v` coordinates as the reference, `6w1v-sf.cif`, and `6w1v_omitO6_single_confB.pdb` | `OEY` A416 in the omit model, the single major conformation with O6 absent |
| 2mFo−DFc | Room-temperature two-flash map | `9z35` coordinates, structure factors, and validation \(2F_o-F_c\) coefficients | `OEY` A601, plus D1-Glu189 |

The `6w1o` and `6w1v` structures are the dark-stable \(S_1\) (0F) and 200 ms two-flash (2F) photosystem II models of Ibrahim *et al.*, [Proc. Natl. Acad. Sci. U.S.A. 117, 12624 (2020)](https://doi.org/10.1073/pnas.2000529117). Residue numbers belong to one PDB entry. `OEY` A416 in `6w1v` is a different residue record from `OEX` A414 in `6w1o`, from `OEX` A415 in the conformation-A omit model, from `OEY` A416 in the conformation-B omit model, and from `OEY` A601 in `9z35`. Every workflow can list hetero residues before the map is clipped.

## Where the calculation runs

Each workflow page runs the same five stages. The reconstruction repeats stages 3–5 once per atom, and can pause after a map so the current peak can be inspected before the next atom is added.

1. Check CCP4, Coot, PyMOL, and the input files. A path or an upload whose name ends in `.gz` is decompressed with `gunzip`. PDB and mmCIF are read as stored.
2. Inspect the coordinates and the chosen site.
3. Build the map with CCP4.
4. Export a radius-limited fragment with headless Coot.
5. Contour that fragment in PyMOL and write a PNG and a `.pml` script.

The page then opens the model and the clipped map in the browser. That viewer draws one isosurface. Its slider runs from the viewer sigma minimum to the viewer sigma maximum on the form, and the number box sets an exact level inside that range. Positive sigma is blue and negative sigma is red in the blue/red colouring, darker farther from zero. Temperature colouring runs from cold at the low end of the range to hot at the high end. The PNG sigma levels are a separate list: they are the contours in the PyMOL figure, and the slider does not change that PNG.

CCP4 is expected at `/Applications/ccp4-9` unless you set another suite prefix. Coot is the binary shipped with CCP4. PyMOL must be on `PATH`. Installation of those programs is described in the project `README.md` and on the application home page.

## Read this site

```bash
uv sync --group docs
uv run mkdocs serve
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000). `uv run mkdocs build` writes the static pages to `site/`.

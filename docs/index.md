# Electron density maps from XFEL structures

This site explains the three map workflows in the tutorial notebooks and in the local web application. Each workflow collects coordinates and reflection files, runs the same CCP4, Coot, and PyMOL programs, and draws an isodensity figure around a site you choose.

The notebooks are the detailed record of the calculation. The web application is the form that runs those steps. The pages here follow the notebook text: the formula for the map, the files it needs, and the program sequence that builds it.

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

-   **Sharpened 2mFo−DFc**

    ---

    Electron density of a single structure, from the deposited validation coefficients, sharpened by a $B$-factor shift and drawn on the model.

    [:octicons-arrow-right-24: Theory and steps](two-fo-fc.md)

</div>

## What each map answers

| Question | Map | Observations | Model |
| --- | --- | --- | --- |
| Where did atoms move between two isomorphous states? | Fo−Fo | Light and dark amplitudes, scaled together | Phases from one refined model; atoms are not moved |
| What is in the data that this model leaves out? | mFo−DFc | One amplitude set | Omit coordinates, optionally refined against those amplitudes |
| Where is the electron density of this structure? | 2mFo−DFc | Validation coefficients of that structure | The same structure, superimposed on the map |

A chemically meaningful Fo−Fo feature is a negative peak on the old position next to a positive peak on the new one. An omit-map feature for an atom left out of the model is a positive peak. A 2mFo−DFc mesh that wraps an atom is accounted-for density. The three figures use different contour ranges: difference maps are usually inspected at 3–5σ, while a 2mFo−DFc map is usually inspected nearer 1–2σ of its own root-mean-square.

The shared ideas — phases, weighting, σ contours, the RCSB files, and the Coot and PyMOL tail that every workflow shares — are on the [background](concepts.md) page. Names of the configuration variables and the matching web-form fields are on the [parameters](parameters.md) page.

## Worked examples

The notebooks are set up on photosystem II.

| Workflow | Default inputs | Site |
| --- | --- | --- |
| Fo−Fo | Dark `6w1o` (0F, dark-stable $S_1$) and light `6w1v` (2F, 200 ms), Ibrahim *et al.*, Science, 2020 | `OEY` A416, the major 2F oxygen-evolving complex, plus D1-Glu189 |
| mFo−DFc | Omit model `6w1v_omitO6_single_confB` against the `6w1v` amplitudes. The notebook introduction places this comparison with Fig. 2B of Wang, Armstrong, and Batista, PNAS, 2021 | `OEY` A416 in the omit model, with O6 absent |
| 2mFo−DFc | `9z35`, the room-temperature XFEL structure of the two-flash, $S_3$-rich state | `OEY` A601, plus D1-Glu189 |

Residue numbers are properties of one PDB entry. `OEY` A416 in `6w1v` is a different residue record from `OEX` A414 in `6w1o` and from `OEY` A601 in `9z35`. Every workflow lists hetero residues so the chain and residue can be checked before the map is clipped.

## Where the calculation runs

Both interfaces perform the same five stages.

1. Check CCP4, Coot, PyMOL, and the input files.
2. Inspect the coordinates and the chosen site.
3. Build the map with CCP4.
4. Export a radius-limited fragment with headless Coot.
5. Contour that fragment in PyMOL and write a PNG and a `.pml` script.

The web application also opens the model and the clipped map in the browser. Camera inspection in the notebook (`INSPECT_IN_PYMOL`) is the PyMOL window; in the application it is the on-page viewer.

CCP4 is expected at `/Applications/ccp4-9` unless you set another suite prefix. Coot is the binary shipped with CCP4. PyMOL must be on `PATH`. Installation of those programs is described in the project `README.md` and on the application home page.

## Read this site

```bash
uv sync --group docs
uv run mkdocs serve
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000). `uv run mkdocs build` writes the static pages to `site/`.

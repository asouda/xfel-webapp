# Construction of XFEL electron density maps

Local site web application for the three workflows:

- isomorphous Fo−Fo difference map
- omit-model mFo−DFc difference map, with optional refinement
- sharpened 2mFo−DFc electron density map

CCP4, Coot, and PyMOL do the calculation. This application collects the files and parameters, runs the same programs, and shows the isodensity figure.

The documentation site in `docs/` explains the theory and the step-by-step calculation for each workflow. Build and read it with:

```bash
uv sync --group docs
uv run mkdocs serve
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000).

## Environment

Python dependencies are declared in `pyproject.toml` and locked in `uv.lock`. The project needs Python 3.9 or newer (`.python-version` pins 3.9). [uv](https://docs.astral.sh/uv/) creates the virtual environment and installs Flask.

From this directory:

```bash
uv sync
```

That writes `.venv` and installs the locked dependencies. You do not need to activate the environment before `uv run`.

The map programs are separate from the Python environment:

1. **CCP4.** Install the [CCP4 suite](https://www.ccp4.ac.uk/download/). The default directory in the app is `/Applications/ccp4-9`. Point the CCP4 field at the suite prefix if yours is elsewhere; the app sources `bin/ccp4.setup-sh` and then calls `cif2mtz`, `refmac5`, `cad`, `scaleit`, `fft`, `extends`, `mtzdump`, `ctruncate`, and `libcheck`.
2. **Coot.** Use the binary shipped with CCP4, usually `$CCP4/bin/coot`. It runs without a window and writes a map fragment around the residue you choose.
3. **PyMOL.** The `pymol` executable must be on `PATH`. The conda-forge package is enough:
  ```bash
   conda install -c conda-forge pymol
  ```

The home page lists whether each program was found. Change the CCP4 directory there and choose **Recheck** after installing or moving the suite.

## Run the application

```bash
cd xfel-webapp
uv sync
uv run python xfel-app.py
```

Open [http://127.0.0.1:5050](http://127.0.0.1:5050). The server listens on `127.0.0.1`. Set `PORT` to use another port:

```bash
PORT=8080 uv run python xfel-app.py
```

Installation notes for CCP4, Coot, and PyMOL are also on that page. Uploaded files and intermediates (MTZ, logs, the extended map, the site fragment, the PNG, and the PyMOL script) are written to the output directory entered on the workflow form. The default is `runs/`.

## Directory structure

```
xfel-webapp/
├── xfel-app.py          Flask application
├── pyproject.toml       Project metadata, dependencies, and license
├── uv.lock              Locked dependency versions
├── .python-version      Python version used by uv (3.9)
├── xfel_web/            Calculation package
├── templates/           HTML pages
├── static/              CSS, JavaScript, and the NGL viewer
├── notebooks/           Tutorial notebooks the workflows follow
├── data/                Example coordinate and reflection files
├── tests/               Unit tests
├── docs/                Documentation site (MkDocs)
├── runtime/             Per-job working files (created at startup, not in git)
└── runs/                Default output directory (not in git)
```



### `xfel_web/`

Python package imported by `xfel-app.py`.

- `__init__.py` — package marker.
- `ccp4.py` — locates CCP4, Coot, and PyMOL, runs the suite after sourcing `ccp4.setup-sh`, and builds the PyMOL scene.
- `forms.py` — reads the workflow form into the dictionaries the pipelines expect.
- `pipelines.py` — the three calculations: Fo−Fo, omit-model mFo−DFc, and 2mFo−DFc.
- `restraints.py` — REFMAC external restraints and optional `libcheck` monomer libraries.



### `templates/`

Jinja templates rendered by Flask.

- `base.html` — shared page shell, navigation, and footer.
- `index.html` — workflow cards and the software check.
- `workflow.html` — file and parameter form, job log, figure, and interactive viewer for one workflow.



### `static/`

Browser assets.

- `css/app.css` — layout and styling.
- `js/home.js` — rechecks CCP4, Coot, and PyMOL from the home page.
- `js/workflow.js` — submits a job, polls its log, and drives the figure and viewer.
- `js/particles.js` — background animation on the page canvas.
- `vendor/ngl-2.5.0.js` — [NGL](https://github.com/nglviewer/ngl) build used to show the model and clipped map in the browser.



### `notebooks/`

The three tutorial Jupyter notebooks:

- `XFEL_fo-fo_diffmap.ipynb` — isomorphous Fo−Fo difference map
- `XFEL_fo-fc_diffmap_refinement.ipynb` — omit-model mFo−DFc map, with optional refinement
- `XFEL_2fo-fc_map.ipynb` — sharpened 2mFo−DFc map

The web workflows follow the same CCP4, Coot, and PyMOL steps.

### `data/`

Example inputs grouped by PDB id: `6w1o`, `6w1p`, `6w1v`, `8iri`, `9z35`, and `11yo`. A directory typically holds coordinates (PDB or mmCIF), structure-factor mmCIF (`*-sf.cif`), and validation map-coefficient files. Some entries also include a validation PDF, a validation mmCIF, or a chain-id mapping. The form defaults look for related files under the parent tutorial `data/` directory when those paths exist.

### `tests/`

`test_forms.py` checks form parsing, resolution selection, and restraint text. Run it with:

```bash
uv run python -m unittest tests.test_forms
```



### `docs/`

MkDocs site for the three workflows. The pages follow the tutorial notebooks: the map each workflow calculates, the files it needs, and the CCP4, Coot, and PyMOL steps.

```bash
uv sync --group docs
uv run mkdocs serve
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000). `uv run mkdocs build` writes the static site to `site/` (gitignored).

### `runtime/`

`runtime/jobs/` holds one working directory per submitted job. The application creates it on startup. It is listed in `.gitignore`.

### `runs/`

Default destination for maps, logs, figures, and PyMOL scripts when the form’s output directory is left at its preset. It is listed in `.gitignore` and stays empty until a workflow writes there.

## Authors

- Alexander Soudackov — [asoudakov@icloud.com](mailto:asoudakov@icloud.com)



## License

MIT License

Copyright (c) 2026 Alexander Soudackov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
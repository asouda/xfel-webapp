"""Local site for the three XFEL electron-density workflows."""

from __future__ import annotations

import json
import os
import threading
import time
import traceback
import uuid
from pathlib import Path

from flask import Flask, abort, jsonify, render_template, request, send_file
from werkzeug.utils import secure_filename

from xfel_web.ccp4 import LIPID_PIGMENT, probe_software, pymol_scene_links, scan_hetatm
from xfel_web.forms import FormError, decompress_gzip, parse_job
from xfel_web.pipelines import RUNNERS

WEB_ROOT = Path(__file__).resolve().parent
PROJECT = WEB_ROOT.parent
DATA = PROJECT / "data"
RUNTIME = WEB_ROOT / "runtime" / "jobs"
RUNTIME.mkdir(parents=True, exist_ok=True)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 800 * 1024 * 1024

_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()
_PROBE_CACHE: dict[str, tuple[float, dict]] = {}


def _existing(relative: str) -> str:
    path = (PROJECT / relative).resolve()
    return str(path) if path.is_file() else ""


def _output_default(name: str) -> str:
    return str((WEB_ROOT / "runs" / name).resolve())


WORKFLOWS = {
    "fofo": {
        "slug": "fofo",
        "nav": "Fo−Fo map",
        "title": "Isomorphous Fo−Fo difference map",
        "kicker": "XFEL_fo-fo_map_fig1.ipynb",
        "card": "Signed difference density between a dark reference dataset and a light, activated dataset.",
        "how": "cif2mtz and a zero-cycle REFMAC prepare both amplitude sets. cad and scaleit put them on one phase set, fft builds |F_light| − |F_dark|, and Coot plus PyMOL draw the ±nσ surfaces.",
        "scan_key": "light_pdb",
        "scan_label": "light coordinates",
        "files": [
            ("dark_pdb", "Dark coordinates", "Reference model, PDB or mmCIF."),
            ("light_pdb", "Light coordinates", "Activated model. The figure is drawn on this structure."),
            ("dark_sf", "Dark structure factors", "mmCIF reflection file, usually {id}-sf.cif."),
            ("light_sf", "Light structure factors", "mmCIF reflection file for the activated dataset."),
            ("dark_validation", "Dark 2Fo−Fc coefficients", "Validation map coefficients. Used only if phases come from the dark model."),
            ("light_validation", "Light 2Fo−Fc coefficients", "Validation FWT/PHWT file. The default phase source."),
        ],
        "paths": {
            "dark_pdb": _existing("data/6w1o.pdb"),
            "light_pdb": _existing("data/6w1v.pdb"),
            "dark_sf": _existing("data/6w1o-sf.cif"),
            "light_sf": _existing("data/6w1v-sf.cif"),
            "dark_validation": _existing("data/6w1o_validation_2fo-fc_map_coef.cif"),
            "light_validation": _existing("data/6w1v_validation_2fo-fc_map_coef.cif"),
        },
        "defaults": {
            "dark_id": "6w1o",
            "light_id": "6w1v",
            "phase_source": "light_validation",
            "site_chain": "A",
            "site_resi": "416",
            "site_label": "OEC",
            "extra_residues": "A 189 GLU189",
            "sigma_levels": "3, 4, 5",
            "surface": "mesh",
            "lowres": "20",
            "dmin": "2.10",
            "png_width": "1600",
            "png_height": "800",
        },
    },
    "refinement": {
        "slug": "refinement",
        "nav": "Omit mFo−DFc",
        "title": "Omit-model mFo−DFc map",
        "kicker": "XFEL_refinement.ipynb",
        "card": "Weighted difference map of an omit model against one set of observed amplitudes, with optional REFMAC.",
        "how": "ctruncate converts deposited intensities to amplitudes. REFMAC writes the weighted DELFWT coefficient for the number of cycles you set, then fft, Coot, and PyMOL produce the figure.",
        "scan_key": "dark_pdb",
        "scan_label": "omit-model coordinates",
        "files": [
            ("dark_pdb", "Omit-model coordinates", "Single-conformation model the map is calculated from."),
            ("light_pdb", "Light coordinates", "Coordinates of the dataset being compared. The figure itself uses the refined omit model."),
            ("light_sf", "Observed structure factors", "mmCIF with mean intensities I and SIGI, or Friedel intensities I(+) and I(-)."),
        ],
        "paths": {
            "dark_pdb": _existing("data/6w1v_omitO6_single_confB.pdb"),
            "light_pdb": _existing("data/6w1v.pdb"),
            "light_sf": _existing("data/6w1v-sf.cif"),
        },
        "defaults": {
            "dark_id": "6w1v_omitO6_single_confB",
            "light_id": "6w1v",
            "site_chain": "A",
            "site_resi": "416",
            "site_label": "OEC",
            "extra_residues": "A 189 GLU189",
            "sigma_levels": "3.5",
            "surface": "mesh",
            "refmac_ncyc": "0",
            "distance_sigma": "0.5",
            "angle_sigma": "25",
            "distance_restraints": "MN1 O5\nMN3 O5\nMN4 O5\nCA1 O5",
            "angle_restraints": "MN3 O5 MN4\nMN3 O5 CA1\nMN4 O5 CA1\nO5 CA1 O1\nO2 MN3 O5\nO3 MN3 O5\nO4 MN3 O5\nO4 MN4 O5",
            "lowres": "30",
            "dmin": "1.89",
            "png_width": "1600",
            "png_height": "800",
        },
    },
    "twofofc": {
        "slug": "twofofc",
        "nav": "2Fo−Fc map",
        "title": "Sharpened 2mFo−DFc map",
        "kicker": "XFEL_2fo-fc_map.ipynb",
        "card": "Weighted electron density from deposited validation coefficients, sharpened and drawn on the model.",
        "how": "fft sharpens the deposited FWT amplitudes by a B-factor shift and transforms them with PHWT. Coot clips the site, and PyMOL draws one isosurface at each σ level.",
        "scan_key": "pdb",
        "scan_label": "coordinates",
        "files": [
            ("pdb", "Coordinates", "Model drawn on the map and used to extend it."),
            ("sf", "Structure factors", "mmCIF reflections. Their resolution range is the default FFT shell."),
            ("validation", "Validation 2Fo−Fc coefficients", "FWT and PHWT. This file is the map, not a difference map."),
        ],
        "paths": {
            "pdb": _existing("data/9z35/9z35.pdb"),
            "sf": _existing("data/9z35/9z35-sf.cif"),
            "validation": _existing("data/9z35/9z35_validation_2fo-fc_map_coef.cif"),
        },
        "defaults": {
            "pdb_id": "9z35",
            "delta_b": "-40",
            "site_chain": "A",
            "site_resi": "601",
            "site_label": "OEC",
            "extra_residues": "A 189 GLU189",
            "sigma_levels": "2, 3.5",
            "viewer_sigma_min": "1",
            "viewer_sigma_max": "5",
            "surface": "solid",
            "lowres": "30",
            "dmin": "2.09",
            "png_width": "1600",
            "png_height": "800",
        },
    },
}

SHARED_DEFAULTS = {
    "ccp4_root": "/Applications/ccp4-9",
    "resolution_source": "sf",
    "map_border": "5",
    "export_radius": "10",
    "carve": "3",
    "contour_mode": "sigma",
    "viewer_sigma_min": "-5",
    "viewer_sigma_max": "5",
    "extra_copies": "",
}


def _job_path(job_id: str) -> Path:
    return RUNTIME / f"{job_id}.json"


def _save_job(job: dict) -> None:
    payload = {
        "id": job["id"],
        "workflow": job["workflow"],
        "status": job["status"],
        "error": job.get("error"),
        "created": job["created"],
        "log_path": job["log_path"],
        "result": job.get("result"),
    }
    _job_path(job["id"]).write_text(json.dumps(payload, indent=2) + "\n")


def _load_jobs() -> None:
    for path in RUNTIME.glob("*.json"):
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if data.get("status") == "running":
            data["status"] = "failed"
            data["error"] = "The server stopped while this calculation was running. Files already written were left in the output directory."
            path.write_text(json.dumps(data, indent=2) + "\n")
        _JOBS[data["id"]] = data


_load_jobs()


class _Log:
    def __init__(self, path: Path):
        self.path = path
        self.lock = threading.Lock()

    def __call__(self, message: str = "") -> None:
        text = str(message).rstrip("\n")
        with self.lock:
            with self.path.open("a", encoding="utf-8") as handle:
                handle.write(text + "\n")
        print(text, flush=True)


def _execute(job_id: str) -> None:
    with _JOBS_LOCK:
        job = _JOBS[job_id]
        job["status"] = "running"
        _save_job(job)
        config = job["config"]
        workflow = job["workflow"]
    emit = _Log(Path(job["log_path"]))
    try:
        result = RUNNERS[workflow](config, emit)
        with _JOBS_LOCK:
            job["result"] = result
            job["status"] = "succeeded"
            job["error"] = None
            _save_job(job)
        emit("")
        emit("Finished.")
        emit(result["caption"])
    except Exception as exc:
        emit("")
        emit(f"Failed: {exc}")
        emit(traceback.format_exc())
        with _JOBS_LOCK:
            job["status"] = "failed"
            job["error"] = str(exc)
            _save_job(job)


def _software(ccp4_root: str) -> dict:
    now = time.time()
    cached = _PROBE_CACHE.get(ccp4_root)
    if cached and now - cached[0] < 60:
        return cached[1]
    report = probe_software(ccp4_root)
    _PROBE_CACHE[ccp4_root] = (now, report)
    return report


@app.get("/")
def index():
    root = request.args.get("ccp4_root", "/Applications/ccp4-9")
    return render_template(
        "index.html",
        workflows=WORKFLOWS,
        software=_software(root),
        ccp4_root=root,
    )


@app.get("/workflow/<slug>")
def workflow(slug: str):
    spec = WORKFLOWS.get(slug)
    if spec is None:
        abort(404)
    defaults = {**SHARED_DEFAULTS, **spec["defaults"], "output_dir": _output_default(slug)}
    return render_template(
        "workflow.html",
        spec=spec,
        defaults=defaults,
        paths=spec["paths"],
        workflows=WORKFLOWS,
    )


@app.get("/api/software")
def software():
    root = request.args.get("ccp4_root", "/Applications/ccp4-9")
    return jsonify(_software(root))


@app.post("/api/scan")
def scan():
    upload = request.files.get("file")
    raw_path = (request.form.get("path") or "").strip()
    folder = WEB_ROOT / "runtime" / "scans"
    if upload and upload.filename:
        folder.mkdir(parents=True, exist_ok=True)
        filename = secure_filename(upload.filename) or "coordinates.pdb"
        dest = folder / f"{uuid.uuid4().hex}_{filename}"
        upload.save(dest)
        pdb_path = dest
    elif raw_path:
        pdb_path = Path(raw_path).expanduser()
        if not pdb_path.is_absolute():
            pdb_path = (PROJECT / pdb_path).resolve()
    else:
        return jsonify({"error": "Choose a coordinate file or enter its path."}), 400
    if not pdb_path.is_file():
        return jsonify({"error": f"Coordinate file not found: {pdb_path}"}), 400
    try:
        if pdb_path.name.lower().endswith(".gz"):
            folder.mkdir(parents=True, exist_ok=True)
            pdb_path = decompress_gzip(pdb_path, folder, "Coordinate file", dest_name=pdb_path.name)
        rows = scan_hetatm(pdb_path)
    except FormError as exc:
        return jsonify({"error": str(exc)}), 400
    except OSError as exc:
        return jsonify({"error": str(exc)}), 400
    visible = []
    for rec in rows:
        if rec["resname"] in LIPID_PIGMENT:
            continue
        visible.append(
            {
                "chain": rec["chain"],
                "resi": rec["resi"],
                "resname": rec["resname"],
                "natoms": len(rec["atoms"]),
                "occ": f"{rec['occ'][0]:.2f}–{rec['occ'][1]:.2f}",
                "atoms": " ".join(sorted(rec["atoms"])[:12]),
            }
        )
    return jsonify(
        {
            "file": pdb_path.name,
            "shown": len(visible),
            "total": len(rows),
            "rows": visible[:250],
        }
    )


@app.post("/api/jobs/<slug>")
def start_job(slug: str):
    if slug not in RUNNERS:
        abort(404)
    try:
        config = parse_job(slug, request.form, request.files, PROJECT)
    except FormError as exc:
        return jsonify({"error": str(exc)}), 400
    job_id = uuid.uuid4().hex
    log_path = RUNTIME / f"{job_id}.log"
    log_path.write_text("")
    job = {
        "id": job_id,
        "workflow": slug,
        "status": "queued",
        "error": None,
        "created": time.time(),
        "log_path": str(log_path),
        "config": config,
        "result": None,
    }
    with _JOBS_LOCK:
        _JOBS[job_id] = job
        _save_job(job)
    threading.Thread(target=_execute, args=(job_id,), daemon=True).start()
    return jsonify({"id": job_id, "output_dir": str(config["output_dir"])})


def _get_job(job_id: str) -> dict:
    if not job_id or len(job_id) != 32:
        abort(404)
    with _JOBS_LOCK:
        job = _JOBS.get(job_id)
    if job is None:
        abort(404)
    return job


def _script_links(result: dict) -> dict:
    """Bonds, distances, and map path taken from the PyMOL script for this job."""
    empty = {"map": None, "bonds": [], "distances": []}
    script = result.get("pymol_script")
    if not script:
        return empty
    path = Path(script)
    if not path.is_file():
        return empty
    return pymol_scene_links(path)


def _viewer_map_file(result: dict, links: dict) -> Path:
    """Serve the map the PyMOL script loads, which is the Coot site fragment."""
    fallback = Path((result.get("viewer") or {}).get("map", ""))
    loaded = links.get("map")
    if not loaded or not result.get("output_dir"):
        return fallback
    candidate = Path(loaded)
    try:
        candidate.resolve().relative_to(Path(result["output_dir"]).resolve())
    except ValueError:
        return fallback
    if candidate.is_file():
        return candidate
    return fallback


@app.get("/api/jobs/<job_id>")
def job_status(job_id: str):
    job = _get_job(job_id)
    try:
        offset = max(0, int(request.args.get("offset", "0")))
    except ValueError:
        offset = 0
    log_path = Path(job["log_path"])
    chunk = ""
    next_offset = offset
    if log_path.is_file():
        with log_path.open("r", encoding="utf-8", errors="replace") as handle:
            handle.seek(offset)
            chunk = handle.read(120_000)
            next_offset = handle.tell()
    result = job.get("result")
    public_result = None
    if result:
        links = _script_links(result)
        viewer = None
        if result.get("viewer"):
            viewer = {
                **result["viewer"],
                "model_url": f"/api/jobs/{job_id}/viewer/model",
                "map_url": f"/api/jobs/{job_id}/viewer/map",
                "map_file": _viewer_map_file(result, links).name,
                "bonds": links["bonds"],
                "distances": links["distances"],
            }
            # The sharpened 2mFo−DFc figure keeps density within this radius of
            # the model. The viewer uses the same radius on the Coot fragment.
            if job.get("workflow") == "twofofc" and links.get("carve"):
                viewer["carve"] = links["carve"]
        public_result = {
            "caption": result["caption"],
            "output_dir": result["output_dir"],
            "final_map": result["final_map"],
            "fragment": result["fragment"],
            "pymol_script": result["pymol_script"],
            "figure_url": f"/api/jobs/{job_id}/figure",
            "viewer": viewer,
            "files": result["files"],
        }
    return jsonify(
        {
            "id": job_id,
            "status": job["status"],
            "error": job.get("error"),
            "log": chunk,
            "offset": next_offset,
            "result": public_result,
        }
    )


@app.post("/api/jobs/<job_id>/snapshot")
def job_snapshot(job_id: str):
    """Store a PNG of the interactive view, after the user has moved the camera."""
    job = _get_job(job_id)
    if job.get("status") != "succeeded":
        return jsonify({"error": "The scene is not ready to save."}), 400
    result = job.get("result") or {}
    output_dir = Path(result.get("output_dir") or "")
    if not output_dir.is_dir():
        return jsonify({"error": "The output directory for this calculation is not available."}), 400
    upload = request.files.get("image")
    if upload is None:
        return jsonify({"error": "No image was sent."}), 400
    data = upload.read(40 * 1024 * 1024 + 1)
    if len(data) < 8 or data[:8] != b"\x89PNG\r\n\x1a\n":
        return jsonify({"error": "The image is not a PNG."}), 400
    if len(data) > 40 * 1024 * 1024:
        return jsonify({"error": "The image is too large to save."}), 400
    stamp = time.strftime("%Y%m%d-%H%M%S")
    millis = int((time.time() % 1) * 1000)
    dest = output_dir / f"viewer_{stamp}-{millis:03d}.png"
    while dest.exists():
        millis = (millis + 1) % 1000
        dest = output_dir / f"viewer_{stamp}-{millis:03d}.png"
    dest.write_bytes(data)
    with _JOBS_LOCK:
        files = result.setdefault("files", [])
        files.append({"role": "Viewer snapshot", "path": str(dest)})
        _save_job(job)
    return jsonify({"path": str(dest), "name": dest.name})


@app.get("/api/jobs/<job_id>/figure")
def job_figure(job_id: str):
    job = _get_job(job_id)
    result = job.get("result") or {}
    figure = result.get("figure")
    if not figure or not Path(figure).is_file():
        abort(404)
    return send_file(figure, mimetype="image/png")


@app.get("/api/jobs/<job_id>/viewer/<asset>")
def job_viewer_asset(job_id: str, asset: str):
    """Serve only the model and clipped map that belong to this finished job."""
    job = _get_job(job_id)
    scene = (job.get("result") or {}).get("viewer")
    if not scene or asset not in {"model", "map"}:
        abort(404)
    if asset == "map":
        path = _viewer_map_file(job.get("result") or {}, _script_links(job.get("result") or {}))
    else:
        path = Path(scene[asset])
    if not path.is_file():
        abort(404)
    mimetype = "chemical/x-pdb" if asset == "model" else "application/octet-stream"
    return send_file(path, mimetype=mimetype, conditional=True)


@app.errorhandler(404)
def not_found(_exc):
    if request.path.startswith("/api/"):
        return jsonify({"error": "Not found."}), 404
    return render_template("index.html", workflows=WORKFLOWS, software=_software("/Applications/ccp4-9"), ccp4_root="/Applications/ccp4-9"), 404


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5050"))
    app.run(host="127.0.0.1", port=port, threaded=True, debug=False)

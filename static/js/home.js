(() => {
  const form = document.getElementById("probe-form");
  const list = document.getElementById("status-list");
  if (!form || !list) return;

  function item(ok, name, path) {
    const li = document.createElement("li");
    li.className = ok ? "ok" : "bad";
    const span = document.createElement("span");
    span.textContent = name;
    const code = document.createElement("code");
    code.textContent = path || "not found";
    li.append(span, code);
    return li;
  }

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    const root = new FormData(form).get("ccp4_root");
    const response = await fetch(`/api/software?ccp4_root=${encodeURIComponent(root)}`);
    const data = await response.json();
    list.replaceChildren(
      item(data.ccp4_ok, "CCP4", data.ccp4_ok ? data.ccp4_setup : "setup script not found"),
      item(data.coot_ok, "Coot", data.coot_path),
      item(data.pymol_ok, "PyMOL", data.pymol_path || "not on PATH"),
      ...data.programs.map((program) => item(program.ok, program.name, program.path)),
    );
  });
})();

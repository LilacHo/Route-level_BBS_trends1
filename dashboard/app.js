// app.js -- wires up MANIFEST (data/manifest.js) and TABLES (data/tables.js),
// both loaded as plain globals via <script src> so this page works from a
// double-clicked index.html with no server / fetch() needed.

(function () {
  "use strict";

  const MODELS = ["base", "anthro"];
  const MODEL_LABEL = { base: "Base", anthro: "Anthro" };
  const SCENARIOS = ["rcp45", "rcp85"];
  const SCENARIO_LABEL = { rcp45: "RCP 4.5", rcp85: "RCP 8.5" };
  const GROUP_LEVELS = ["Contraction", "Stable", "Expansion"];

  const speciesByCode = {};
  MANIFEST.species.forEach((s) => { speciesByCode[s.code] = s; });

  const state = {
    group: "",       // "" = all groups
    code: MANIFEST.species[0].code,
    model: "base",
    scenario: "rcp45",
  };

  const el = {
    groupSelect: document.getElementById("group-select"),
    speciesFilter: document.getElementById("species-filter"),
    speciesSelect: document.getElementById("species-select"),
    modelToggle: document.getElementById("model-toggle"),
    scenarioSelect: document.getElementById("scenario-select"),
    heading: document.getElementById("species-heading"),
    rcpImg: document.getElementById("rcp-violin-img"),
    rcpWrap: document.getElementById("rcp-violin-wrap"),
    groupImg: document.getElementById("group-violin-img"),
    groupWrap: document.getElementById("group-violin-wrap"),
    mapImg: document.getElementById("map-img"),
    mapWrap: document.getElementById("map-wrap"),
    tableBody: document.getElementById("stats-tbody"),
    tableNote: document.getElementById("table-note"),
  };

  function groupLabel(g) {
    return g.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
  }

  function populateGroupSelect() {
    el.groupSelect.innerHTML = "";
    const optAll = document.createElement("option");
    optAll.value = "";
    optAll.textContent = "All groups";
    el.groupSelect.appendChild(optAll);
    MANIFEST.groups.forEach((g) => {
      const opt = document.createElement("option");
      opt.value = g;
      opt.textContent = groupLabel(g);
      el.groupSelect.appendChild(opt);
    });
  }

  function filteredSpecies() {
    const q = el.speciesFilter.value.trim().toLowerCase();
    return MANIFEST.species.filter((s) => {
      if (state.group && s.group !== state.group) return false;
      if (q && !(s.name.toLowerCase().includes(q) || s.code.toLowerCase().includes(q))) return false;
      return true;
    });
  }

  function populateSpeciesSelect(preserveSelection) {
    const list = filteredSpecies();
    const prev = state.code;
    el.speciesSelect.innerHTML = "";
    list.forEach((s) => {
      const opt = document.createElement("option");
      opt.value = s.code;
      opt.textContent = `${s.name} (${s.code})`;
      el.speciesSelect.appendChild(opt);
    });
    if (list.length === 0) {
      const opt = document.createElement("option");
      opt.textContent = "No species match";
      opt.disabled = true;
      el.speciesSelect.appendChild(opt);
      return;
    }
    const stillPresent = preserveSelection && list.some((s) => s.code === prev);
    state.code = stillPresent ? prev : list[0].code;
    el.speciesSelect.value = state.code;
  }

  function fmt(x) {
    if (x === null || x === undefined) return "—";
    return Number(x).toFixed(2);
  }

  function renderPanels() {
    const sp = speciesByCode[state.code];
    el.heading.innerHTML = `${sp.name} <span class="code">(${sp.code}, ${groupLabel(sp.group)})</span>`;

    setImg(el.rcpWrap, el.rcpImg, `img/violin_rcp/${sp.code}_${state.model}.webp`,
      "By-RCP violin not available.");

    if (sp.has_group_violin) {
      setImg(el.groupWrap, el.groupImg, `img/violin_group/${sp.code}_${state.model}.webp`,
        "By-group violin not available.");
    } else {
      setPlaceholder(el.groupWrap,
        "No routes fall in Contraction/Stable/Expansion for either scenario (all routes are in the raster's “never suitable” category) — this plot was never generated for this species.");
    }

    setImg(el.mapWrap, el.mapImg, `img/maps/${sp.code}_${state.model}_${state.scenario}.webp`,
      "Range-shift map not available.");
  }

  function setImg(wrap, imgEl, src, missingMsg) {
    imgEl.onerror = () => setPlaceholder(wrap, missingMsg);
    imgEl.onload = () => {
      wrap.querySelectorAll(".placeholder").forEach((p) => p.remove());
      imgEl.style.display = "";
    };
    imgEl.style.display = "";
    imgEl.src = src;
  }

  function setPlaceholder(wrap, msg) {
    wrap.querySelectorAll("img").forEach((i) => { i.style.display = "none"; });
    let ph = wrap.querySelector(".placeholder");
    if (!ph) {
      ph = document.createElement("div");
      ph.className = "placeholder";
      wrap.appendChild(ph);
    }
    ph.textContent = msg;
  }

  function renderTable() {
    const rowsByCategory = TABLES[state.code];
    el.tableBody.innerHTML = "";
    let anyData = false;

    GROUP_LEVELS.forEach((cat) => {
      const tr = document.createElement("tr");
      const th = document.createElement("th");
      th.scope = "row";
      th.innerHTML = `<span class="cat-dot cat-${cat}"></span>${cat}`;
      tr.appendChild(th);

      MODELS.forEach((model) => {
        const scenarioData = (rowsByCategory && rowsByCategory[model] && rowsByCategory[model][state.scenario]) || {};
        const row = scenarioData[cat];
        if (row) anyData = true;
        ["n", "mean", "sd"].forEach((field) => {
          const td = document.createElement("td");
          if (row) {
            td.textContent = field === "n" ? row.n : fmt(row[field]);
          } else {
            td.textContent = "—";
            td.className = "na-cell";
          }
          tr.appendChild(td);
        });
      });
      el.tableBody.appendChild(tr);
    });

    el.tableNote.textContent = anyData
      ? "n = routes, mean/sd = annual population trend (% per year), by SDM range-change category. — = no routes in that category for this model/scenario."
      : "No Contraction/Stable/Expansion route data for this species/scenario (all routes fall in the raster's “never suitable” category).";
  }

  function render() {
    renderPanels();
    renderTable();
  }

  // Events ---------------------------------------------------------------
  el.groupSelect.addEventListener("change", () => {
    state.group = el.groupSelect.value;
    populateSpeciesSelect(false);
    state.code = el.speciesSelect.value;
    render();
  });

  el.speciesFilter.addEventListener("input", () => {
    populateSpeciesSelect(true);
    state.code = el.speciesSelect.value;
    render();
  });

  el.speciesSelect.addEventListener("change", () => {
    state.code = el.speciesSelect.value;
    render();
  });

  el.modelToggle.querySelectorAll("button").forEach((btn) => {
    btn.addEventListener("click", () => {
      state.model = btn.dataset.model;
      el.modelToggle.querySelectorAll("button").forEach((b) => b.classList.toggle("active", b === btn));
      renderPanels();
    });
  });

  el.scenarioSelect.addEventListener("change", () => {
    state.scenario = el.scenarioSelect.value;
    render();
  });

  // Init -------------------------------------------------------------------
  populateGroupSelect();
  SCENARIOS.forEach((sc) => {
    const opt = document.createElement("option");
    opt.value = sc;
    opt.textContent = SCENARIO_LABEL[sc];
    el.scenarioSelect.appendChild(opt);
  });
  populateSpeciesSelect(false);
  render();
})();

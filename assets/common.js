/* Shared rendering helpers for the trunnion results pages. Data: window.TRUNNION (data.js). */
(function () {
  "use strict";
  const T = () => window.TRUNNION;

  const fmt = {
    n(v, d = 1) { return Number(v).toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d }); },
    um(v) { return fmt.n(v, 0) + " \u00b5m"; },
    arcmin(v) { return fmt.n(v, 1) + "\u2032"; },
  };

  function el(tag, attrs, children) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === "class") e.className = v;
      else if (k === "text") e.textContent = v;
      else e.setAttribute(k, v);
    }
    for (const c of children || []) e.append(c);
    return e;
  }

  /* Specification rows shown for one design: [label, value, unit, note] */
  function specRows(key) {
    const d = T().designs[key];
    return [
      ["Printed volume", fmt.n(d.volume_cm3, 0), "cm\u00b3", ""],
      ["Mass (PA-CF, 1.25 g/cm\u00b3)", fmt.n(d.mass_g, 0), "g", "100 % infill"],
      ["Worst-case compliance, load envelope", fmt.n(d.worst_compliance, 2), "N\u00b7mm", "lower is stiffer"],
      ["Seat deflection, 80 N vertical", fmt.n(d.defl_vertical, 0), "\u00b5m", ""],
      ["Seat deflection, 80 N horizontal", fmt.n(d.defl_horizontal, 0), "\u00b5m", ""],
      ["Seat deflection, 30 N axial", fmt.n(d.defl_axial, 0), "\u00b5m", ""],
      ["Axle tilt, 100 N belt pull", fmt.n(d.belt_arcmin, 1), "arcmin", "target \u2264 5"],
      ["Peak von Mises stress", fmt.n(d.vm_max, 1), "MPa", "all load cases"],
      ["Unsupported overhang (45\u00b0)", fmt.n(d.overhang_pct, 1), "%", "of down-facing area"],
      ["Separate solid bodies", String(d.bodies), "", ""],
    ];
  }

  function renderSpecTable(container, key) {
    const tbody = el("tbody");
    for (const [label, value, unit, note] of specRows(key)) {
      tbody.append(el("tr", {}, [
        el("th", { scope: "row", text: label }),
        el("td", { class: "num", text: value }),
        el("td", { class: "unit", text: unit }),
        el("td", { class: "note", text: note }),
      ]));
    }
    container.replaceChildren(el("table", { class: "spec" }, [
      el("thead", {}, [el("tr", {}, [
        el("th", { scope: "col", text: "Property" }),
        el("th", { scope: "col", class: "num", text: "Value" }),
        el("th", { scope: "col", text: "Unit" }),
        el("th", { scope: "col", text: "Note" }),
      ])]),
      tbody,
    ]));
  }

  const COMPARE_COLS = [
    ["Volume", (d) => fmt.n(d.volume_cm3, 0), "cm\u00b3"],
    ["Worst-case compliance", (d) => fmt.n(d.worst_compliance, 2), "N\u00b7mm"],
    ["Axial seat deflection", (d) => fmt.n(d.defl_axial, 0), "\u00b5m"],
    ["Axle tilt under belt", (d) => fmt.n(d.belt_arcmin, 1), "arcmin"],
    ["Peak stress", (d) => fmt.n(d.vm_max, 1), "MPa"],
    ["Overhang", (d) => fmt.n(d.overhang_pct, 1), "%"],
  ];

  function renderCompareTable(container, onSelect, selectedKey) {
    const data = T();
    const head = el("tr", {}, [el("th", { scope: "col", text: "Design" })]);
    for (const [label, , unit] of COMPARE_COLS) {
      head.append(el("th", { scope: "col", class: "num" }, [label, el("span", { class: "u", text: unit })]));
    }
    const tbody = el("tbody");
    for (const key of data.order) {
      const d = data.designs[key];
      const name = el("th", { scope: "row" }, [
        el("span", { class: "dname", text: d.label }),
        el("span", { class: "dsub", text: d.short }),
      ]);
      const row = el("tr", { "data-key": key }, [name]);
      if (key === data.recommended) row.classList.add("recommended");
      if (key === selectedKey) row.classList.add("selected");
      for (const [, get] of COMPARE_COLS) row.append(el("td", { class: "num", text: get(d) }));
      if (onSelect) {
        row.tabIndex = 0;
        row.addEventListener("click", () => onSelect(key));
        row.addEventListener("keydown", (e) => {
          if (e.key === "Enter" || e.key === " ") { e.preventDefault(); onSelect(key); }
        });
      }
      tbody.append(row);
    }
    container.replaceChildren(el("table", { class: "compare" }, [el("thead", {}, [head]), tbody]));
  }

  function renderGates(container, key) {
    const d = T().designs[key];
    const list = el("ul", { class: "gates" });
    for (const g of d.gates) {
      list.append(el("li", { class: g.pass ? "pass" : "fail" }, [
        el("span", { class: "gstate", text: g.pass ? "Pass" : "Fail" }),
        el("span", { class: "glabel", text: g.label }),
        el("span", { class: "gvalue", text: g.value }),
      ]));
    }
    container.replaceChildren(list);
  }

  /* Tilt dial: seat deflection (radius) against A-axis tilt (angle, -90..+90 deg). */
  function renderTiltDial(container, opts) {
    const o = Object.assign({ size: 320, accent: "#B3460E", baseline: "#6E7C8A", ink: "#2A2F35", line: "#D7DBDF" }, opts || {});
    const data = T();
    const ns = "http://www.w3.org/2000/svg";
    const W = o.size, H = Math.round(o.size * 0.62), cx = W / 2, cy = H - 26, R = Math.min(cx - 20, cy - 18);
    let maxD = 0;
    for (const k of data.order) for (const [, dv] of data.designs[k].tilt) maxD = Math.max(maxD, dv);
    const ringStep = maxD > 60 ? 20 : 10;
    const rMax = Math.ceil(maxD / ringStep) * ringStep;
    const svg = document.createElementNS(ns, "svg");
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.setAttribute("role", "img");
    svg.setAttribute("aria-label", "Bearing seat deflection versus A-axis tilt for each design");
    function add(tag, attrs, text) {
      const e = document.createElementNS(ns, tag);
      for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
      if (text) e.textContent = text;
      svg.append(e);
      return e;
    }
    const toXY = (deg, r) => {
      const a = (deg * Math.PI) / 180;
      return [cx + r * Math.sin(a), cy - r * Math.cos(a)];
    };
    for (let r = ringStep; r <= rMax; r += ringStep) {
      const rr = (r / rMax) * R;
      const [x1, y1] = toXY(-90, rr), [x2, y2] = toXY(90, rr);
      add("path", { d: `M${x1},${y1} A${rr},${rr} 0 0 1 ${x2},${y2}`, fill: "none", stroke: o.line, "stroke-width": 1 });
      add("text", { x: cx + 3, y: cy - rr + 11, "font-size": 10, fill: o.ink, "font-family": "inherit" }, `${r}`);
    }
    for (const deg of [-90, -45, 0, 45, 90]) {
      const [x, y] = toXY(deg, R + 2);
      add("line", { x1: cx, y1: cy, x2: x, y2: y, stroke: o.line, "stroke-width": 1 });
      const [tx, ty] = toXY(deg, R + 12);
      add("text", { x: tx, y: ty + 3, "font-size": 10, fill: o.ink, "text-anchor": "middle", "font-family": "inherit" }, `${deg}\u00b0`);
    }
    add("text", { x: cx, y: H - 6, "font-size": 10, fill: o.ink, "text-anchor": "middle", "font-family": "inherit" },
      "A-axis tilt; radius = seat deflection under 80 N (\u00b5m)");
    const keys = data.order.filter((k) => k !== data.recommended).concat([data.recommended]);
    for (const k of keys) {
      const pts = data.designs[k].tilt.map(([deg, dv]) => toXY(deg, (dv / rMax) * R).join(",")).join(" ");
      const isRec = k === data.recommended;
      add("polyline", {
        points: pts, fill: "none",
        stroke: isRec ? o.accent : o.baseline,
        "stroke-width": isRec ? 2.5 : 1.2,
        "stroke-dasharray": isRec ? "" : (k.startsWith("R0") ? "2 3" : k.startsWith("v0") ? "6 3" : ""),
        opacity: isRec ? 1 : 0.85,
      });
    }
    container.replaceChildren(svg);
    const legend = el("ul", { class: "legend" });
    for (const k of data.order) {
      const isRec = k === data.recommended;
      legend.append(el("li", { class: isRec ? "rec" : "" }, [el("span", { class: "swatch" }), data.designs[k].label]));
    }
    container.append(legend);
  }

  function renderFiles(container) {
    const list = el("ul", { class: "files" });
    for (const f of T().files) {
      const a = el("a", { href: f.href, download: "" }, [f.name]);
      list.append(el("li", {}, [a, el("span", { class: "fmeta", text: f.meta })]));
    }
    container.replaceChildren(list);
  }

  function renderAssumptions(container) {
    const list = el("ol", { class: "assumptions" });
    for (const a of T().assumptions) list.append(el("li", { text: a }));
    container.replaceChildren(list);
  }

  window.TrunnionUI = { fmt, renderSpecTable, renderCompareTable, renderGates, renderTiltDial, renderFiles, renderAssumptions };
})();

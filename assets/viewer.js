/* Minimal offline WebGL mesh viewer (no dependencies).
   createViewer(canvas, options) -> { show(meshKey), reset() }
   Meshes come from window.TRUNNION.meshes: base64 uint16 positions (quantised) + uint32 indices. */
(function () {
  "use strict";

  function decodeBase64(b64) {
    const bin = atob(b64);
    const bytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    return bytes.buffer;
  }

  function buildGeometry(mesh) {
    const q = new Uint16Array(decodeBase64(mesh.positions));
    const idx = new Uint32Array(decodeBase64(mesh.indices));
    const n = q.length / 3;
    const pos = new Float32Array(q.length);
    for (let i = 0; i < n; i++) {
      for (let a = 0; a < 3; a++) {
        pos[3 * i + a] = mesh.lo[a] + (q[3 * i + a] / 65535) * mesh.span;
      }
    }
    const nrm = new Float32Array(q.length);
    for (let t = 0; t < idx.length; t += 3) {
      const a = idx[t] * 3, b = idx[t + 1] * 3, c = idx[t + 2] * 3;
      const ux = pos[b] - pos[a], uy = pos[b + 1] - pos[a + 1], uz = pos[b + 2] - pos[a + 2];
      const vx = pos[c] - pos[a], vy = pos[c + 1] - pos[a + 1], vz = pos[c + 2] - pos[a + 2];
      const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      for (const k of [a, b, c]) {
        nrm[k] += nx; nrm[k + 1] += ny; nrm[k + 2] += nz;
      }
    }
    for (let i = 0; i < nrm.length; i += 3) {
      const l = Math.hypot(nrm[i], nrm[i + 1], nrm[i + 2]) || 1;
      nrm[i] /= l; nrm[i + 1] /= l; nrm[i + 2] /= l;
    }
    return { pos, nrm, idx };
  }

  const VS = `
    attribute vec3 p; attribute vec3 n;
    uniform mat4 mvp; uniform mat3 nm;
    varying vec3 vn;
    void main() { vn = nm * n; gl_Position = mvp * vec4(p, 1.0); }`;
  const FS = `
    precision mediump float;
    varying vec3 vn; uniform vec3 color;
    void main() {
      vec3 n = normalize(vn);
      vec3 l1 = normalize(vec3(0.4, 0.5, 0.8));
      vec3 l2 = normalize(vec3(-0.6, -0.2, 0.3));
      float d = 0.30 + 0.62 * max(dot(n, l1), 0.0) + 0.18 * max(dot(n, l2), 0.0);
      gl_FragColor = vec4(color * d, 1.0);
    }`;

  function mat4mul(a, b) {
    const o = new Float32Array(16);
    for (let r = 0; r < 4; r++)
      for (let c = 0; c < 4; c++) {
        let s = 0;
        for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
        o[c * 4 + r] = s;
      }
    return o;
  }

  function hexToRgb(hex) {
    const v = parseInt(hex.replace("#", ""), 16);
    return [((v >> 16) & 255) / 255, ((v >> 8) & 255) / 255, (v & 255) / 255];
  }

  window.createViewer = function (canvas, options) {
    const opts = Object.assign({ color: "#B3460E", background: "#FFFFFF" }, options || {});
    const gl = canvas.getContext("webgl", { antialias: true }) || canvas.getContext("experimental-webgl");
    if (!gl) {
      canvas.insertAdjacentHTML("afterend", '<p class="viewer-fallback">3D view needs WebGL. Open the STL in any slicer instead.</p>');
      return { show() {}, reset() {}, setColor() {} };
    }
    gl.getExtension("OES_element_index_uint");
    function compile(type, src) {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      return s;
    }
    let loc = null;
    let cache = {};
    function setupProgram() {
      const prog = gl.createProgram();
      gl.attachShader(prog, compile(gl.VERTEX_SHADER, VS));
      gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, FS));
      gl.linkProgram(prog);
      gl.useProgram(prog);
      loc = {
        p: gl.getAttribLocation(prog, "p"),
        n: gl.getAttribLocation(prog, "n"),
        mvp: gl.getUniformLocation(prog, "mvp"),
        nm: gl.getUniformLocation(prog, "nm"),
        color: gl.getUniformLocation(prog, "color"),
      };
      cache = {};
    }
    setupProgram();
    // Recover if the browser drops the GPU context (driver reset, too many tabs).
    canvas.addEventListener("webglcontextlost", (e) => e.preventDefault());
    canvas.addEventListener("webglcontextrestored", () => {
      gl.getExtension("OES_element_index_uint");
      setupProgram();
      draw();
    });
    let current = null;
    let yaw = -0.75, pitch = 0.32, dist = 330;
    let color = hexToRgb(opts.color);
    const bg = hexToRgb(opts.background);

    function upload(key) {
      if (cache[key]) return cache[key];
      const g = buildGeometry(window.TRUNNION.meshes[key]);
      const pb = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, pb);
      gl.bufferData(gl.ARRAY_BUFFER, g.pos, gl.STATIC_DRAW);
      const nb = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, nb);
      gl.bufferData(gl.ARRAY_BUFFER, g.nrm, gl.STATIC_DRAW);
      const ib = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, ib);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, g.idx, gl.STATIC_DRAW);
      cache[key] = { pb, nb, ib, count: g.idx.length };
      return cache[key];
    }

    function draw() {
      if (gl.isContextLost()) return;
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.clientWidth, h = canvas.clientHeight;
      if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
        canvas.width = Math.round(w * dpr);
        canvas.height = Math.round(h * dpr);
      }
      gl.viewport(0, 0, canvas.width, canvas.height);
      gl.clearColor(bg[0], bg[1], bg[2], 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      gl.enable(gl.DEPTH_TEST);
      if (!current) return;
      const m = upload(current);
      // model: centre the part (Z up in CAD -> Y up on screen)
      const cx = 60, cy = 20, cz = 71;
      const cy_ = Math.cos(yaw), sy_ = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
      // C maps CAD (Z up) to view space (Y up); then yaw about the vertical, then pitch.
      const C = new Float32Array([1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1]);
      const RY = new Float32Array([cy_, 0, -sy_, 0, 0, 1, 0, 0, sy_, 0, cy_, 0, 0, 0, 0, 1]);
      const RX = new Float32Array([1, 0, 0, 0, 0, cp, sp, 0, 0, -sp, cp, 0, 0, 0, 0, 1]);
      const R = mat4mul(RX, mat4mul(RY, C));
      const T = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -cx, -cy, -cz, 1]);
      const V = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -dist, 1]);
      const aspect = canvas.width / Math.max(canvas.height, 1);
      const f = 1 / Math.tan((30 * Math.PI) / 180 / 2), near = 10, far = 2000;
      const P = new Float32Array([
        f / aspect, 0, 0, 0,
        0, f, 0, 0,
        0, 0, (far + near) / (near - far), -1,
        0, 0, (2 * far * near) / (near - far), 0,
      ]);
      const MV = mat4mul(V, mat4mul(R, T));
      const MVP = mat4mul(P, MV);
      const NM = new Float32Array([R[0], R[1], R[2], R[4], R[5], R[6], R[8], R[9], R[10]]);
      gl.uniformMatrix4fv(loc.mvp, false, MVP);
      gl.uniformMatrix3fv(loc.nm, false, NM);
      gl.uniform3fv(loc.color, color);
      gl.bindBuffer(gl.ARRAY_BUFFER, m.pb);
      gl.enableVertexAttribArray(loc.p);
      gl.vertexAttribPointer(loc.p, 3, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ARRAY_BUFFER, m.nb);
      gl.enableVertexAttribArray(loc.n);
      gl.vertexAttribPointer(loc.n, 3, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, m.ib);
      gl.drawElements(gl.TRIANGLES, m.count, gl.UNSIGNED_INT, 0);
    }

    let dragging = false, lx = 0, ly = 0;
    canvas.addEventListener("pointerdown", (e) => {
      dragging = true; lx = e.clientX; ly = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      yaw += (e.clientX - lx) * 0.008;
      pitch = Math.max(-1.4, Math.min(1.4, pitch + (e.clientY - ly) * 0.008));
      lx = e.clientX; ly = e.clientY;
      draw();
    });
    canvas.addEventListener("pointerup", () => { dragging = false; });
    canvas.addEventListener("wheel", (e) => {
      e.preventDefault();
      dist = Math.max(160, Math.min(700, dist * (1 + Math.sign(e.deltaY) * 0.08)));
      draw();
    }, { passive: false });
    canvas.tabIndex = 0;
    canvas.addEventListener("keydown", (e) => {
      const step = 0.12;
      if (e.key === "ArrowLeft") yaw -= step;
      else if (e.key === "ArrowRight") yaw += step;
      else if (e.key === "ArrowUp") pitch = Math.min(1.4, pitch + step);
      else if (e.key === "ArrowDown") pitch = Math.max(-1.4, pitch - step);
      else return;
      e.preventDefault();
      draw();
    });
    window.addEventListener("resize", draw);
    if (window.ResizeObserver) {
      new ResizeObserver(() => draw()).observe(canvas);
    }

    return {
      show(key) { current = key; draw(); },
      reset() { yaw = -0.75; pitch = 0.32; dist = 330; draw(); },
      setColor(hex) { color = hexToRgb(hex); draw(); },
      view(name) {
        const views = { iso: [-0.75, 0.32], front: [0, 0], side: [-Math.PI / 2, 0], top: [0, 1.39] };
        [yaw, pitch] = views[name] || views.iso;
        draw();
      },
    };
  };
})();

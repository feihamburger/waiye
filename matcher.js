/* 无序边长 → 常见房屋模板识别
 * 支持：矩形(4)、L 形(6)、U 形(8)、T 形(8)
 * 不要求测量顺序；用边长关系做最小二乘匹配。
 */
(function (global) {
  "use strict";

  const TEMPLATES = {
    rect: {
      name: "矩形", n: 4,
      params: ["a", "b"],
      edges: [[1, 0], [0, 1], [1, 0], [0, 1]],
      build: p => [[0, 0], [p.a, 0], [p.a, p.b], [0, p.b]],
      valid: p => p.a > 0.2 && p.b > 0.2
    },
    L: {
      name: "L 形", n: 6,
      params: ["W", "H", "nw", "nh"],
      edges: [
        [1, 0, 0, 0],
        [0, 1, 0, -1],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, -1, 0],
        [0, 1, 0, 0]
      ],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H - p.nh], [p.W - p.nw, p.H - p.nh],
                   [p.W - p.nw, p.H], [0, p.H]],
      valid: p => p.W > 0.4 && p.H > 0.4 && p.nw > 0.1 && p.nh > 0.1
                   && p.W - p.nw > 0.1 && p.H - p.nh > 0.1
    },
    U: {
      name: "U 形", n: 8,
      params: ["W", "H", "rw", "hb", "lw"],
      edges: [
        [1, 0, 0, 0, 0],     // W
        [0, 1, 0, 0, 0],     // H
        [0, 0, 1, 0, 0],     // rw
        [0, 1, 0, -1, 0],    // H-hb
        [1, 0, -1, 0, -1],   // W-rw-lw
        [0, 1, 0, -1, 0],    // H-hb
        [0, 0, 0, 0, 1],     // lw
        [0, 1, 0, 0, 0]      // H
      ],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H], [p.W - p.rw, p.H],
                   [p.W - p.rw, p.hb], [p.lw, p.hb], [p.lw, p.H], [0, p.H]],
      valid: p => p.W > 1 && p.H > p.hb + 0.2 && p.rw > 0.2 && p.lw > 0.2 && p.hb > 0.1 && p.W - p.rw - p.lw > 0.2
    },
    stair: {
      name: "阶梯形", n: 8,
      params: ["W", "H", "w1", "w2", "h1", "h2"],
      edges: [
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [1, 0, -1, -1, 0, 0],
        [0, 1, 0, 0, 0, -1],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, 1],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0]
      ],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H], [p.W - p.w1 - p.w2, p.H],
                   [p.W - p.w1 - p.w2, p.h2], [p.w1, p.h2], [p.w1, p.h1], [0, p.h1]],
      valid: p => p.W > 1 && p.H > p.h2 + 0.2 && p.h2 > p.h1 + 0.1 && p.h1 > 0.1
                   && p.w1 > 0.1 && p.w2 > 0.1 && p.W - p.w1 - p.w2 > 0.1
    },
    T: {
      name: "T 形", n: 8,
      params: ["W", "H", "th", "sw"],
      edges: [
        [1, 0, 0, 0],
        [0, 1, -1, 0],
        [0.5, 0, 0, -0.5],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
        [0.5, 0, 0, -0.5],
        [0, 1, -1, 0]
      ],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H - p.th], [(p.W + p.sw) / 2, p.H - p.th],
                   [(p.W + p.sw) / 2, p.H], [(p.W - p.sw) / 2, p.H], [(p.W - p.sw) / 2, p.H - p.th], [0, p.H - p.th]],
      valid: p => p.W > 1 && p.H > p.th + 0.2 && p.th > 0.2 && p.sw > 0.2 && p.W - p.sw > 0.2
    }
  };

  function solveLS(A, y, k) {
    const n = A.length;
    const AtA = Array.from({ length: k }, () => new Array(k).fill(0));
    const Aty = new Array(k).fill(0);
    for (let i = 0; i < n; i++) {
      for (let r = 0; r < k; r++) {
        Aty[r] += A[i][r] * y[i];
        for (let c = 0; c < k; c++) AtA[r][c] += A[i][r] * A[i][c];
      }
    }
    for (let r = 0; r < k; r++) AtA[r][r] += 1e-9;
    // 高斯消元
    const M = AtA.map((row, i) => row.concat([Aty[i]]));
    for (let col = 0; col < k; col++) {
      let piv = col;
      for (let r = col + 1; r < k; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
      if (Math.abs(M[piv][col]) < 1e-12) return null;
      [M[col], M[piv]] = [M[piv], M[col]];
      const pv = M[col][col];
      for (let c = col; c <= k; c++) M[col][c] /= pv;
      for (let r = 0; r < k; r++) {
        if (r === col) continue;
        const f = M[r][col];
        if (!f) continue;
        for (let c = col; c <= k; c++) M[r][c] -= f * M[col][c];
      }
    }
    return M.map(row => row[k]);
  }

  function permutations(arr, cb) {
    const a = arr.slice();
    const n = a.length;
    (function rec(k) {
      if (k === n) { cb(a); return; }
      for (let i = k; i < n; i++) {
        [a[k], a[i]] = [a[i], a[k]];
        rec(k + 1);
        [a[k], a[i]] = [a[i], a[k]];
      }
    })(0);
  }

  function shoelace(v) {
    let s = 0;
    for (let i = 0; i < v.length; i++) {
      const p = v[i], q = v[(i + 1) % v.length];
      s += p[0] * q[1] - q[0] * p[1];
    }
    return Math.abs(s) / 2;
  }

  function segInter(a, b, c, d) {
    const o = (p, q, r) => (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    const o1 = o(a, b, c), o2 = o(a, b, d), o3 = o(c, d, a), o4 = o(c, d, b);
    return o1 * o2 < -1e-9 && o3 * o4 < -1e-9;
  }

  function simple(v) {
    const n = v.length;
    for (let i = 0; i < n; i++) {
      for (let j = i + 2; j < n; j++) {
        if (i === 0 && j === n - 1) continue;
        if (segInter(v[i], v[(i + 1) % n], v[j], v[(j + 1) % n])) return false;
      }
    }
    for (let i = 0; i < n; i++) {
      const a = v[i], b = v[(i + 1) % n], c = v[(i + 2) % n];
      const cr = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
      if (Math.abs(cr) < 1e-9) return false;
    }
    return shoelace(v) > 0.3;
  }

  function normKey(v) {
    const mx = Math.min(...v.map(p => p[0])), my = Math.min(...v.map(p => p[1]));
    return v.map(p => `${(p[0] - mx).toFixed(2)},${(p[1] - my).toFixed(2)}`).join(";");
  }

  function canonical(v) {
    let best = null;
    const T = [
      p => [p[0], p[1]], p => [p[1], -p[0]], p => [-p[0], -p[1]], p => [-p[1], p[0]],
      p => [-p[0], p[1]], p => [p[0], -p[1]], p => [p[1], p[0]], p => [-p[1], -p[0]]
    ];
    for (const t of T) {
      const a = v.map(t);
      for (const arr of [a, a.slice().reverse()]) {
        const k = normKey(arr);
        if (best === null || k < best) best = k;
      }
    }
    return best;
  }

  function dPointSeg(p, a, b) {
    const dx = b[0] - a[0], dy = b[1] - a[1];
    const L2 = dx * dx + dy * dy;
    if (L2 === 0) return Math.hypot(p[0] - a[0], p[1] - a[1]);
    let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(p[0] - (a[0] + t * dx), p[1] - (a[1] + t * dy));
  }
  function dPointPoly(p, v) {
    let m = Infinity;
    for (let i = 0; i < v.length; i++) m = Math.min(m, dPointSeg(p, v[i], v[(i + 1) % v.length]));
    return m;
  }
  function hausdorff(a, b) {
    let m = 0;
    for (const p of a) m = Math.max(m, dPointPoly(p, b));
    for (const p of b) m = Math.max(m, dPointPoly(p, a));
    return m;
  }
  function centroid(v) {
    return [v.reduce((s, p) => s + p[0], 0) / v.length, v.reduce((s, p) => s + p[1], 0) / v.length];
  }
  function shapeDist(a, b) {
    const ca = centroid(a);
    const T = [
      p => [p[0], p[1]], p => [p[1], -p[0]], p => [-p[0], -p[1]], p => [-p[1], p[0]],
      p => [-p[0], p[1]], p => [p[0], -p[1]], p => [p[1], p[0]], p => [-p[1], -p[0]]
    ];
    let best = Infinity;
    for (const t of T) {
      const bt = b.map(t);
      const cb = centroid(bt);
      const moved = bt.map(p => [p[0] + ca[0] - cb[0], p[1] + ca[1] - cb[1]]);
      best = Math.min(best, hausdorff(a, moved));
    }
    return best;
  }

  function match(lengths, opts) {
    const o = Object.assign({ tolMm: 12, maxCandidates: 6 }, opts || {});
    const vals = lengths.map(Number).filter(v => isFinite(v) && v > 0);
    const n = vals.length;
    if (![4, 6, 8].includes(n)) {
      return { error: "目前只支持 4 / 6 / 8 条边（矩形、L、U、T）", candidates: [] };
    }
    const found = new Map();

    for (const key of Object.keys(TEMPLATES)) {
      const t = TEMPLATES[key];
      if (t.n !== n) continue;
      if (key === "stair" && o.template !== "stair") continue;
      const idx = [...Array(n).keys()];
      permutations(idx, perm => {
        const L = perm.map(i => vals[i]);
        const x = solveLS(t.edges, L, t.params.length);
        if (!x) return;
        const p = {};
        t.params.forEach((name, j) => p[name] = x[j]);
        if (!t.valid(p)) return;
        const fit = t.edges.map(row => row.reduce((s, c, j) => s + c * x[j], 0));
        let s2 = 0;
        for (let i = 0; i < n; i++) s2 += (fit[i] - L[i]) ** 2;
        const res = Math.sqrt(s2 / n) * 1000;
        if (res > o.tolMm) return;
        const verts = t.build(p);
        if (!simple(verts)) return;
        // 按模板 + 参数（厘米级舍入）去重，避免毫米极差造成大量重复候选
        const key2 = key + "|" + t.params.map(name => p[name].toFixed(2)).join(",");
        const cand = {
          template: key, name: t.name, params: p,
          vertices: verts, residualMm: res, area: shoelace(verts),
          assignment: L.map((v, i) => ({ edge: i + 1, measured: v, fitted: fit[i] })),
          tolerances: o.tolMm
        };
        const prev = found.get(key2);
        if (!prev || res < prev.residualMm) found.set(key2, cand);
      });
    }

    // 按形状聚类（8 cm 内视为同一种），每种留残差最小的
    const all = [...found.values()].sort((a, b) => a.residualMm - b.residualMm);
    const clusters = [];
    for (const cand of all) {
      if (clusters.some(cl => shapeDist(cand.vertices, cl.vertices) < 0.08)) continue;
      clusters.push(cand);
    }
    // 先按残差（毫米级），残差相同则面积大的优先（真实房间的“缺口”通常更小）
    clusters.sort((a, b) => (Math.round(a.residualMm) - Math.round(b.residualMm)) || (b.area - a.area));
    const candidates = clusters.slice(0, o.maxCandidates);

    return { input: vals, candidates };
  }


  // ---------- 最少边（按对称）输入 ----------
  const MIN_TEMPLATES = {
    rect: {
      name: "矩形", n: 2, params: ["a", "b"],
      build: p => [[0, 0], [p.a, 0], [p.a, p.b], [0, p.b]],
      valid: p => p.a > 0.2 && p.b > 0.2
    },
    L: {
      name: "L 形", n: 4, params: ["W", "Hnh", "nw", "nh"],
      build: p => [[0, 0], [p.W, 0], [p.W, p.Hnh], [p.W - p.nw, p.Hnh],
                   [p.W - p.nw, p.Hnh + p.nh], [0, p.Hnh + p.nh]],
      valid: p => p.W > p.nw + 0.1 && p.nw > 0.1 && p.nh > 0.1 && p.Hnh > 0.1
    },
    U: {
      name: "U 形", n: 4, params: ["W", "H", "rw", "Hin"],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H], [p.W - p.rw, p.H],
                   [p.W - p.rw, p.H - p.Hin], [p.rw, p.H - p.Hin],
                   [p.rw, p.H], [0, p.H]],
      valid: p => p.W > 2 * p.rw + 0.1 && p.rw > 0.1 && p.Hin > 0.1 && p.H > p.Hin + 0.1
    },
    T: {
      name: "T 形", n: 4, params: ["W", "th", "Hsth", "sw"],
      build: p => [[0, 0], [p.W, 0], [p.W, p.Hsth], [(p.W + p.sw) / 2, p.Hsth],
                   [(p.W + p.sw) / 2, p.Hsth + p.th], [(p.W - p.sw) / 2, p.Hsth + p.th],
                   [(p.W - p.sw) / 2, p.Hsth], [0, p.Hsth]],
      valid: p => p.W > p.sw + 0.1 && p.sw > 0.1 && p.th > 0.1 && p.Hsth > 0.1
    },
    stair: {
      name: "阶梯形", n: 6, params: ["W", "H", "Wtop", "Hright", "w2", "h2h1"],
      build: p => [[0, 0], [p.W, 0], [p.W, p.H], [p.W - p.Wtop, p.H],
                   [p.W - p.Wtop, p.H - p.Hright],
                   [p.w2, p.H - p.Hright],
                   [p.w2, p.H - p.Hright - p.h2h1],
                   [0, p.H - p.Hright - p.h2h1]],
      valid: p => p.W > p.Wtop + p.w2 + 0.1 && p.Wtop > 0.1 && p.w2 > 0.1
                   && p.Hright > 0.1 && p.h2h1 > 0.1 && p.H > p.Hright + p.h2h1 + 0.1
    }
  };

  function matchMinimal(values, opts) {
    const o = Object.assign({ template: "auto", maxCandidates: 8 }, opts || {});
    const vals = values.map(Number).filter(v => isFinite(v) && v > 0);
    const n = vals.length;
    if (![2, 4, 6].includes(n)) {
      return { error: "最少边模式需要 2（矩形）、4（L/U/T）或 6（阶梯形）个值", input: vals, candidates: [] };
    }
    const out = [];
    for (const key of Object.keys(MIN_TEMPLATES)) {
      const t = MIN_TEMPLATES[key];
      if (t.n !== n) continue;
      if (o.template && o.template !== "auto" && o.template !== key) continue;
      const k = t.params.length;
      permutations([...Array(n).keys()], perm => {
        const chosen = perm.slice(0, k).map(i => vals[i]);
        const p = {};
        t.params.forEach((name, j) => p[name] = chosen[j]);
        if (!t.valid(p)) return;
        const verts = t.build(p);
        if (!simple(verts)) return;
        const area = shoelace(verts);
        if (area < 0.5) return;
        const lengths = verts.map((v, i) => {
          const w = verts[(i + 1) % verts.length];
          return Math.hypot(w[0] - v[0], w[1] - v[1]);
        });
        out.push({
          template: key, name: t.name, params: p, vertices: verts,
          residualMm: 0, area,
          assignment: lengths.map((L, i) => ({ edge: i + 1, fitted: L }))
        });
      });
    }
    out.sort((a, b) => b.area - a.area);
    const clusters = [];
    for (const cand of out) {
      if (clusters.some(cl => shapeDist(cand.vertices, cl.vertices) < 0.05)) continue;
      clusters.push(cand);
    }
    return { input: vals, candidates: clusters.slice(0, o.maxCandidates) };
  }

  const api = { match, matchMinimal, TEMPLATES, MIN_TEMPLATES };
  global.TemplateMatch = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);

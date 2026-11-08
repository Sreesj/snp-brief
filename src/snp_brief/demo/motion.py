"""The animated DNA helix: drawn in code (no stock image), always moving, driven by scroll.

`HELIX_JS` defines `vaHelix(canvas, options)`. It renders a white, faceted, double helix with
depth shading on a 2D canvas, rotating continuously. `options.scrollTarget` (a window or an
element) makes scrolling turn the helix and slide it back and forth. Scroll-reveal for
elements with the class `va-reveal` is included.

`injector(mode)` returns a tiny HTML component that adds two canvases to the Streamlit page
itself (the component iframe is same-origin): a large background helix and a small vertical
one in the sidebar. It is idempotent across reruns. `mode` is "landing" (full strength) or
"app" (faint, so results stay readable).
"""

from __future__ import annotations

import json

HELIX_JS = r"""
function vaHelix(canvas, opts) {
  opts = Object.assign({vertical: false, opacity: 0.6, radius: 0.16, turns: 3.2, speed: 0.35,
                        scrollTarget: null, density: 1.0}, opts || {});
  const ctx = canvas.getContext('2d');
  const win = canvas.ownerDocument.defaultView;
  let w = 0, h = 0, dpr = 1, phase = 0, scrollY = 0, drift = 0, last = 0;
  const N = Math.round(150 * opts.density);
  const noise = Array.from({length: N * 2 + 8}, (_, i) => 0.7 + 0.6 * Math.abs(Math.sin(i * 12.9898) * 43758.5453 % 1));
  function resize() {
    dpr = Math.min(win.devicePixelRatio || 1, 2);
    const r = canvas.getBoundingClientRect();
    w = Math.max(1, r.width); h = Math.max(1, r.height);
    canvas.width = w * dpr; canvas.height = h * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  function readScroll() {
    const t = opts.scrollTarget;
    if (!t) return 0;
    return (t.scrollY !== undefined ? t.scrollY : t.scrollTop) || 0;
  }
  function frame(now) {
    const dt = Math.min(0.05, (now - last) / 1000 || 0); last = now;
    phase += dt * opts.speed;
    const target = readScroll();
    scrollY += (target - scrollY) * 0.08;
    drift += ((Math.sin(scrollY / 420) * 0.12) - drift) * 0.06;
    draw();
    win.requestAnimationFrame(frame);
  }
  function draw() {
    ctx.clearRect(0, 0, w, h);
    const L = opts.vertical ? h * 1.15 : Math.hypot(w, h) * 1.05;
    const R = (opts.vertical ? w * 0.42 : Math.min(w, h)) * opts.radius;
    // Axis: diagonal from lower-left to upper-right (horizontal mode) or top to bottom (vertical).
    const ang = opts.vertical ? Math.PI / 2 : -0.32;
    const ax = Math.cos(ang), ay = Math.sin(ang), nx = -ay, ny = ax;
    const cx = w * (opts.vertical ? 0.5 : 0.58 + drift), cy = h * (opts.vertical ? 0.5 : 0.52);
    const spin = phase + scrollY / 260;
    const items = [];
    for (let i = 0; i < N; i++) {
      const s = i / (N - 1) - 0.5;
      const px = cx + ax * s * L, py = cy + ay * s * L;
      for (let k = 0; k < 2; k++) {
        const th = s * opts.turns * Math.PI * 2 + spin + k * Math.PI;
        const z = Math.sin(th), off = Math.cos(th) * R;
        const scale = 1 / (1 - z * 0.28);
        items.push({kind: 0, x: px + nx * off * scale, y: py + ny * off * scale, z, r: R * 0.13 * scale * noise[i * 2 + k], i, k});
      }
      if (i % 5 === 0) {
        const th = s * opts.turns * Math.PI * 2 + spin;
        const z1 = Math.sin(th), z2 = -z1, o1 = Math.cos(th) * R, o2 = -o1;
        const s1 = 1 / (1 - z1 * 0.28), s2 = 1 / (1 - z2 * 0.28);
        items.push({kind: 1, x1: px + nx * o1 * s1, y1: py + ny * o1 * s1, x2: px + nx * o2 * s2, y2: py + ny * o2 * s2,
                    z: Math.min(z1, z2) - 0.05, wdt: R * 0.07});
      }
    }
    items.sort((a, b) => a.z - b.z);
    for (const it of items) {
      const light = 0.35 + 0.65 * (it.z * 0.5 + 0.5);
      const g = Math.round(150 + 105 * light);
      ctx.globalAlpha = opts.opacity * (0.35 + 0.65 * light);
      if (it.kind === 1) {
        ctx.strokeStyle = `rgb(${g},${g},${g})`;
        ctx.lineWidth = Math.max(1, it.wdt);
        ctx.beginPath(); ctx.moveTo(it.x1, it.y1); ctx.lineTo(it.x2, it.y2); ctx.stroke();
      } else {
        // a small faceted blob: an irregular polygon, so the strand reads as sculpted, not smooth
        const sides = 6, rot = (it.i * 0.7 + it.k) % 6.283;
        ctx.fillStyle = `rgb(${g},${g},${g})`;
        ctx.beginPath();
        for (let j = 0; j < sides; j++) {
          const a = rot + j * Math.PI * 2 / sides;
          const rr = it.r * (0.75 + 0.35 * noise[(it.i + j * 7) % noise.length] / 1.3);
          const X = it.x + Math.cos(a) * rr, Y = it.y + Math.sin(a) * rr;
          j ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y);
        }
        ctx.closePath(); ctx.fill();
        ctx.globalAlpha *= 0.35;
        ctx.strokeStyle = 'rgb(20,20,22)'; ctx.lineWidth = 0.6; ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
  }
  resize();
  win.addEventListener('resize', resize);
  win.requestAnimationFrame(frame);
  return {resize, setOpacity: (o) => { opts.opacity = o; }};
}

function vaReveal(doc) {
  const win = doc.defaultView;
  if (!win.IntersectionObserver) { doc.querySelectorAll('.va-reveal').forEach(e => e.classList.add('va-in')); return; }
  const io = new win.IntersectionObserver((entries) => {
    entries.forEach(e => { if (e.isIntersecting) { e.target.classList.add('va-in'); io.unobserve(e.target); } });
  }, {threshold: 0.15});
  doc.querySelectorAll('.va-reveal:not(.va-in)').forEach(e => io.observe(e));
}
"""


INSTALL_JS = r"""
(function () {
  const doc = document, pwin = window;
  const main = doc.querySelector('[data-testid="stMain"]') || doc.querySelector('section.main') || pwin;
  function install() {
    if (!doc.getElementById('va-helix-bg')) {
      const shade = doc.createElement('div');
      shade.id = 'va-helix-shade';
      shade.style.cssText = 'position:fixed;inset:0;z-index:-1;pointer-events:none;' +
        'background:radial-gradient(70% 60% at 50% 38%, rgba(8,9,10,.62) 0%, rgba(8,9,10,.45) 45%, rgba(8,9,10,.85) 100%);';
      const c = doc.createElement('canvas');
      c.id = 'va-helix-bg';
      c.style.cssText = 'position:fixed;inset:0;width:100vw;height:100vh;z-index:-1;pointer-events:none;transition:opacity .8s ease;';
      doc.body.prepend(shade);
      doc.body.prepend(c);
      vaHelix(c, {scrollTarget: main, opacity: 0.6});
    }
    const side = doc.querySelector('[data-testid="stSidebar"]');
    if (side && !doc.getElementById('va-helix-side')) {
      const s = doc.createElement('canvas');
      s.id = 'va-helix-side';
      s.style.cssText = 'position:absolute;right:-24px;bottom:0;width:150px;height:62%;z-index:0;pointer-events:none;opacity:.2;';
      side.appendChild(s);
      vaHelix(s, {vertical: true, opacity: 0.9, radius: 0.5, turns: 2.4, speed: 0.25, density: 0.6,
                  scrollTarget: side.querySelector('[data-testid="stSidebarContent"]')});
    }
  }
  window.__vaSetMode = function (mode) {
    install();
    const bg = doc.getElementById('va-helix-bg');
    if (bg) bg.style.opacity = mode === 'landing' ? '1' : '0.3';
    vaReveal(doc);
    setTimeout(() => vaReveal(doc), 500);
    setTimeout(() => vaReveal(doc), 1500);
  };
})();
"""


def injector(mode: str = "landing") -> str:
    """A zero-height component that installs the helix into the parent Streamlit page.

    The drawing code is added to the parent page itself (once), so the animation keeps running
    when Streamlit rebuilds or removes this little iframe on a rerun.
    """
    code = json.dumps(HELIX_JS + INSTALL_JS)
    return (
        "<script>(function(){"
        "let w; try { w = window.parent; w.document; } catch (e) { return; }"
        "if (!w.__vaSetMode) { const s = w.document.createElement('script'); s.textContent = "
        + code
        + "; w.document.head.appendChild(s); }"
        f"w.__vaSetMode({json.dumps(mode)});"
        "})();</script>"
    )

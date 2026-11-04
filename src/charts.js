import {
  exact, compact, niceTicks, scaleLinear, scaleLog, logTicks, labelStep, waffle as waffleCells, stack, easeOut, spread, squarify, heat,
} from './lib.js';

const NS = 'http://www.w3.org/2000/svg';
const reduced = () => matchMedia('(prefers-reduced-motion: reduce)').matches;

export const el = (tag, attrs = {}, ...kids) => {
  const node = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v !== null && v !== undefined && v !== false) node.setAttribute(k, v);
  }
  for (const kid of kids.flat()) {
    if (kid !== null && kid !== undefined && kid !== false) node.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return node;
};

export const html = (tag, props = {}, ...kids) => {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') node.className = v;
    else node.setAttribute(k, v);
  }
  for (const kid of kids.flat()) {
    if (kid !== null && kid !== undefined && kid !== false) node.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return node;
};

const tipNode = () => document.getElementById('tip');

const clearMarks = () => {
  document.querySelectorAll('.chart .sel').forEach((n) => n.classList.remove('sel'));
  document.querySelectorAll('.chart .cross').forEach((c) => c.setAttribute('visibility', 'hidden'));
};

const place = (node, pageX, pageY) => {
  const w = node.offsetWidth;
  const h = node.offsetHeight;
  const vw = document.documentElement.clientWidth;
  const minX = window.scrollX + 8;
  const maxX = window.scrollX + vw - w - 8;
  const maxY = window.scrollY + window.innerHeight - h - 8;
  let left = pageX + 14;
  if (left > maxX) left = pageX - w - 14;
  let top = pageY + 14;
  if (top > maxY) top = pageY - h - 14;
  node.style.left = `${Math.max(minX, Math.min(left, maxX))}px`;
  node.style.top = `${Math.max(window.scrollY + 8, top)}px`;
};

export const tip = {
  show(lines, pageX, pageY, pinned = false) {
    const node = tipNode();
    const close = pinned ? html('button', { type: 'button', class: 'tip-close', 'aria-label': 'Close' }, '×') : null;
    close?.addEventListener('click', () => tip.hide());
    node.replaceChildren(...lines, ...(close ? [close] : []));
    node.hidden = false;
    node.dataset.pinned = String(pinned);
    place(node, pageX, pageY);
  },
  hide() {
    const node = tipNode();
    node.hidden = true;
    node.dataset.pinned = 'false';
    clearMarks();
  },
  get pinned() {
    const node = tipNode();
    return !node.hidden && node.dataset.pinned === 'true';
  },
};

const installGlobal = () => {
  if (document.body.dataset.tipReady) return;
  document.body.dataset.tipReady = '1';
  document.addEventListener('click', (event) => {
    if (tip.pinned && !event.target.closest('#tip')) tip.hide();
  });
  document.addEventListener('keydown', (event) => {
    if (event.key === 'Escape') tip.hide();
  });
};

const tipLink = (href, text) => html('a', { href, target: '_blank', rel: 'noopener' }, `${text ?? 'Open the source'} ↗`);

export const tipLines = (value, label, ...extra) => [
  html('b', {}, value),
  html('span', {}, label),
  ...extra.filter(Boolean).map((e) => {
    if (e.nodeType) return e;
    if (e.href) return tipLink(e.href, e.text);
    return html('span', {}, e);
  }),
];

export const nearestMark = (root, cx, cy, axis) => {
  let best = null;
  let bestD = Infinity;
  for (const m of root.querySelectorAll('[data-k]')) {
    if (m.classList.contains('hit')) continue;
    const r = m.getBoundingClientRect();
    if (!r.width && !r.height) continue;
    const d = axis === 'x'
      ? Math.abs(r.left + r.width / 2 - cx)
      : Math.hypot(Math.max(r.left - cx, 0, cx - r.right), Math.max(r.top - cy, 0, cy - r.bottom));
    if (d < bestD) {
      bestD = d;
      best = m;
    }
  }
  return bestD <= (axis === 'x' ? 40 : 36) ? best : null;
};

const hover = (root, lookup, { axis, steps } = {}) => {
  installGlobal();
  const resolve = (event) => {
    const direct = event.target.closest?.('[data-k]');
    if (direct && root.contains(direct)) return direct;
    return Number.isFinite(event.clientX) ? nearestMark(root, event.clientX, event.clientY, axis) : null;
  };
  const select = (mark) => {
    clearMarks();
    mark.classList.add('sel');
  };
  root.addEventListener('pointermove', (event) => {
    if (event.pointerType !== 'mouse' || tip.pinned) return;
    const mark = resolve(event);
    if (!mark) {
      tip.hide();
      return;
    }
    tip.show(lookup(mark.dataset.k, event, mark), event.pageX, event.pageY);
  });
  root.addEventListener('pointerleave', (event) => {
    if (event.pointerType === 'mouse' && !tip.pinned) tip.hide();
  });
  root.addEventListener('click', (event) => {
    const mark = resolve(event);
    if (!mark) return;
    event.stopPropagation();
    event.preventDefault();
    select(mark);
    tip.show(lookup(mark.dataset.k, event, mark), event.pageX, event.pageY, true);
  });
  root.setAttribute('tabindex', '0');
  root.setAttribute('aria-keyshortcuts', 'ArrowLeft ArrowRight ArrowUp ArrowDown Home End Enter');
  root.querySelectorAll('[tabindex]').forEach((n) => n !== root && n.removeAttribute('tabindex'));
  const marks = () => [...root.querySelectorAll('[data-k]')]
    .filter((m) => !m.classList.contains('hit') && !m.closest('[visibility="hidden"]') && m.getAttribute('visibility') !== 'hidden')
    .map((m, i) => ({ m, rank: Number(m.dataset.rank ?? i) }))
    .sort((a, b) => a.rank - b.rank)
    .map((x) => x.m);
  let cursor = -1;
  const showAt = (i, pinned) => {
    if (steps) {
      const hit = root.querySelector('.hit');
      const out = lookup(hit.dataset.k, { index: i }, hit);
      const cross = root.querySelector('.cross');
      const r = (cross ?? root).getBoundingClientRect();
      tip.show(out, r.left + window.scrollX, r.top + window.scrollY + r.height / 3, pinned);
      return;
    }
    const m = marks()[i];
    if (!m) return;
    select(m);
    const r = m.getBoundingClientRect();
    tip.show(lookup(m.dataset.k, { index: i }, m), r.right + window.scrollX, r.top + window.scrollY, pinned);
  };
  const MOVES = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1 };
  root.addEventListener('keydown', (event) => {
    const count = steps ?? marks().length;
    if (!count) return;
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      showAt(Math.max(cursor, 0), true);
      return;
    }
    if (event.key in MOVES) cursor = Math.max(0, Math.min(count - 1, cursor + MOVES[event.key]));
    else if (event.key === 'Home') cursor = 0;
    else if (event.key === 'End') cursor = count - 1;
    else return;
    event.preventDefault();
    showAt(cursor, tip.pinned);
  });
  root.addEventListener('focus', () => {
    cursor = Math.max(cursor, 0);
    if (!tip.pinned) showAt(cursor, false);
  });
  root.addEventListener('blur', () => {
    if (!tip.pinned) tip.hide();
  });
};

const widthOf = (mount) => Math.round(Math.min(Math.max(mount.clientWidth || 800, 280), 1120));

const measureLines = (svg) => {
  for (const path of svg.querySelectorAll('.line')) {
    const len = Math.ceil(path.getTotalLength ? path.getTotalLength() : 2000);
    path.style.setProperty('--len', `${len}`);
  }
};

const mountSvg = (mount, svg, ...before) => {
  mount.replaceChildren(...before, svg);
  measureLines(svg);
  return svg;
};

const legendOf = (items) => html('ul', { class: 'legend' }, items.map((s) => html('li', {}, html('i', { style: `background:${s.color}` }), s.label)));

const noteAt = (cx, text, W) => {
  const tw = text.length * 6.8;
  if (cx + 6 + tw <= W) return { x: cx + 6, 'text-anchor': 'start' };
  if (cx - 6 - tw >= 0) return { x: cx - 6, 'text-anchor': 'end' };
  return { x: Math.max(0, W - tw), 'text-anchor': 'start' };
};

const barPath = (x, y, w, h) => {
  const r = Math.min(3, w / 2, h);
  const b = y + h;
  return `M${x},${b}V${y + r}Q${x},${y} ${x + r},${y}H${x + w - r}Q${x + w},${y} ${x + w},${y + r}V${b}Z`;
};

const fitsIn = (text, px) => String(text).length <= Math.max(4, Math.floor(px / 6.6));

const rowLayout = (rows, labelW, rowH) => {
  const fits = rows.every((r) => fitsIn(r.label, labelW));
  return fits
    ? { L: labelW + 10, rowH, labelY: 4, markY: 0, cls: 'rows' }
    : { L: 14, rowH: rowH + 16, labelY: -8, markY: 5, cls: 'rows stacked' };
};

const axisTicks = (lo, hi) => {
  const span = Math.max(hi, -lo);
  const step = niceTicks(span, 4)[1] ?? 1;
  const out = [];
  for (let t = Math.floor(lo / step) * step; t <= Math.ceil(hi / step) * step + step / 2; t += step) out.push(Number(t.toPrecision(12)));
  return out;
};

export const columns = ({
  mount, rows, key, value, overlay, attr = 'year', label = (r) => key(r), highlight = () => false,
  annotations = [], aria, describe, height = 250, yFormat = compact,
}) => {
  const W = widthOf(mount);
  const T = annotations.length ? 22 + annotations.length * 18 : 12;
  const L = 46;
  const B = 24;
  const H = height + T;
  const ticks = niceTicks(Math.max(...rows.map(value)), 4);
  const top = ticks.at(-1) || 1;
  const band = (W - L) / rows.length;
  const bw = Math.max(1, Math.min(band - 2, 26));
  const y = scaleLinear([0, top], [H - B, T]);
  const xAt = (i) => L + i * band + (band - bw) / 2;
  const widest = Math.max(...rows.map((r) => String(label(r)).length));
  const step = labelStep(rows.length, W - L, Math.max(34, widest * 6.6 + 12));
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: t === 0 ? 'base' : 'grid', x1: L, x2: W, y1: y(t), y2: y(t) }),
      el('text', { class: 'tick', x: L - 8, y: y(t) + 4, 'text-anchor': 'end' }, yFormat(t)))),
    rows.map((r, i) => el('path', {
      class: `bar${highlight(r) ? ' hi' : ''}${r.partial ? ' part' : ''}`,
      d: barPath(xAt(i), y(value(r)), bw, Math.max(0, H - B - y(value(r)))),
      [`data-${attr}`]: key(r),
      'data-k': i,
      style: `--i:${i}`,
    })),
    overlay ? rows.map((r, i) => el('path', {
      class: 'bar2',
      d: barPath(xAt(i) + bw * 0.22, y(overlay(r)), bw * 0.56, Math.max(0, H - B - y(overlay(r)))),
    })) : null,
    rows.map((r, i) => (i % step === 0 ? el('text', {
      class: 'tick', x: xAt(i) + bw / 2, y: H - 6, 'text-anchor': 'middle',
    }, label(r)) : null)),
    annotations.map((a, n) => {
      const i = rows.findIndex((r) => key(r) === a.key);
      if (i < 0) return null;
      const cx = xAt(i) + bw / 2;
      const ty = 14 + n * 18;
      return el('g', {},
        el('path', { class: 'note-line', d: `M${cx},${ty + 4}V${y(value(rows[i])) - 5}` }),
        el('text', { class: 'note', y: ty + 4, ...noteAt(cx, a.text, W) }, a.text));
    }));
  hover(svg, (k) => describe(rows[Number(k)], rows[Number(k) - 1]), { axis: 'x' });
  return mountSvg(mount, svg);
};

export const dotplot = ({
  mount, rows, attr = 'row', scale = 'log', format = exact, describe, rowH = 24, href, linkText, labelShare = 0.34, aria,
}) => {
  const W = widthOf(mount);
  const R = 70;
  const T = 22;
  const lay = rowLayout(rows, Math.min(240, Math.max(96, W * labelShare)), rowH);
  const { L } = lay;
  const H = T + rows.length * lay.rowH + 6;
  const values = rows.map((r) => r.value).filter((v) => v > 0);
  const min = Math.min(...values);
  const max = Math.max(...values);
  const lo = scale === 'log' ? 10 ** Math.floor(Math.log10(min)) : 0;
  const hi = scale === 'log' ? 10 ** Math.ceil(Math.log10(max)) : niceTicks(max, 4).at(-1);
  const x = scale === 'log' ? scaleLog([lo, hi], [L, W - R]) : scaleLinear([0, hi], [L, W - R]);
  const spaced = (list) => list.reduce((kept, t) => {
    if (!kept.length || x(t) - x(kept.at(-1)) >= 40) kept.push(t);
    return kept;
  }, []);
  const ticks = spaced(scale === 'log' ? logTicks(lo, hi) : niceTicks(max, 4));
  const at = (v) => (v > 0 ? x(v) : L);
  const rowTop = (rank) => T + rank * lay.rowH + lay.rowH / 2;
  const rowsG = el('g', { class: lay.cls }, rows.map((r, i) => el('g', {
    [`data-${attr}`]: r.key,
    'data-row': attr === 'row' ? r.key : null,
    'data-k': i,
    'data-rank': i,
    transform: `translate(0 ${rowTop(i)})`,
  },
  el('rect', { class: 'row-hit', x: 0, y: -lay.rowH / 2, width: W, height: lay.rowH }),
  el('text', { class: 'lab', x: 0, y: lay.labelY }, r.label),
  el('line', { class: 'stem', x1: L, x2: at(r.value), y1: lay.markY, y2: lay.markY, style: `--i:${i}` }),
  el('circle', { class: `dot${r.hi ? ' hi' : ''}`, cx: at(r.value), cy: lay.markY, r: 4.5, style: `--i:${i}` }),
  el('text', { class: 'val', x: at(r.value) + 9, y: lay.markY + 4 }, format(r.value)))));
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: 'grid', x1: x(t), x2: x(t), y1: T - 6, y2: H }),
      el('text', { class: 'tick', x: x(t), y: 10, 'text-anchor': 'middle' }, compact(t)))),
    rowsG);
  hover(svg, (k) => {
    const r = rows[Number(k)];
    const out = describe(r);
    return href ? [...out, tipLink(href(r), linkText?.(r))] : out;
  });
  mountSvg(mount, svg);
  const groups = [...rowsG.children];
  return {
    order(keys) {
      keys.forEach((key, rank) => {
        const g = groups.find((n) => n.getAttribute(`data-${attr}`) === String(key));
        g.setAttribute('transform', `translate(0 ${rowTop(rank)})`);
        g.style.transform = `translate(0px, ${rowTop(rank)}px)`;
        g.dataset.rank = rank;
      });
    },
  };
};

export const hbars = ({
  mount, rows, segments, aria, describe, format = exact, tickFormat = compact, labelShare = 0.16, rowH = 22, showSign = false,
}) => {
  const W = widthOf(mount);
  const R = 78;
  const T = 22;
  const lay = rowLayout(rows, Math.min(150, Math.max(44, W * labelShare)), rowH);
  const { L } = lay;
  const H = T + rows.length * lay.rowH + 6;
  const totals = rows.map((r) => segments.reduce((n, s) => n + (r.values[s.key] ?? 0), 0));
  const ticks = axisTicks(Math.min(0, ...totals), Math.max(0, ...totals));
  const x = scaleLinear([ticks[0], ticks.at(-1)], [L, W - R]);
  const bh = Math.min(14, lay.rowH - 8);
  const text = (v) => (showSign && v > 0 ? `+${format(v)}` : v < 0 ? `−${format(-v)}` : format(v));
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: t === 0 ? 'base' : 'grid', x1: x(t), x2: x(t), y1: T - 6, y2: H }),
      el('text', { class: 'tick', x: x(t), y: 10, 'text-anchor': 'middle' }, t < 0 ? `−${tickFormat(-t)}` : tickFormat(t)))),
    el('g', { class: lay.cls }, rows.map((r, i) => {
      const mid = T + i * lay.rowH + lay.rowH / 2;
      const y = mid + lay.markY;
      let acc = 0;
      const segs = segments.map((s) => {
        const v = r.values[s.key] ?? 0;
        if (!v) return null;
        const from = v > 0 ? acc : v;
        const to = v > 0 ? acc + v : 0;
        if (v > 0) acc += v;
        return el('rect', {
          class: `seg${v < 0 ? ' neg' : ''}`, x: x(from), y: y - bh / 2, width: Math.max(0.5, x(to) - x(from)), height: bh, style: `fill:${s.color};--i:${i}`,
        });
      });
      const total = totals[i];
      return el('g', { 'data-row': r.key, 'data-k': i, 'data-rank': i },
        el('rect', { class: 'row-hit', x: 0, y: mid - lay.rowH / 2, width: W, height: lay.rowH }),
        el('text', { class: 'lab', x: 0, y: mid + lay.labelY }, r.label),
        segs,
        el('text', { class: 'val', x: total < 0 ? x(0) + 6 : x(total) + 6, y: y + 4 }, text(total)));
    })));
  hover(svg, (k) => describe(rows[Number(k)]));
  return segments.length > 1 ? mountSvg(mount, svg, legendOf(segments)) : mountSvg(mount, svg);
};

export const ideogram = ({
  mount, chromosomes, bands, aria, describe,
}) => {
  const W = widthOf(mount);
  const L = 26;
  const R = 8;
  const T = 4;
  const rowH = W < 560 ? 16 : 20;
  const bh = rowH - 7;
  const max = Math.max(...chromosomes.map((c) => c.length));
  const x = scaleLinear([0, max], [L, W - R]);
  const H = T + chromosomes.length * rowH + 2;
  const flat = [];
  const svg = el('svg', { class: 'ideo', viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    el('defs', {}, el('pattern', {
      id: 'hatch', width: 4, height: 4, patternUnits: 'userSpaceOnUse', patternTransform: 'rotate(45)',
    }, el('rect', { class: 'hatch-bg', width: 4, height: 4 }), el('line', { class: 'hatch-line', x1: 0, y1: 0, x2: 0, y2: 4 }))),
    chromosomes.map((c, i) => el('g', { 'data-chrom': c.chrom, transform: `translate(0 ${T + i * rowH})` },
      el('text', { class: 'lab', x: 0, y: bh / 2 + 4 }, c.chrom),
      bands.filter((b) => b.chrom === c.chrom).map((b) => {
        const k = flat.push(b) - 1;
        const inset = b.stain === 'acen' ? 2 : b.stain === 'stalk' ? 3.5 : 0;
        return el('rect', {
          class: `band ${b.stain}`,
          'data-band': `${b.chrom}${b.name}`,
          'data-k': k,
          'data-rank': k,
          x: x(b.start),
          y: inset,
          width: Math.max(0.6, x(b.end) - x(b.start)),
          height: bh - inset * 2,
          style: `--i:${i}`,
        });
      }),
      el('rect', { class: 'body', x: L, y: 0, width: x(c.length) - L, height: bh, rx: 3 }))));
  hover(svg, (k) => describe(flat[Number(k)]));
  return mountSvg(mount, svg);
};

export const rampLegend = (lo, hi, format, color = 'var(--ink)', steps = 5) => html('ul', { class: 'legend ramp' },
  Array.from({ length: steps }, (_, i) => {
    const v = lo + ((hi - lo) * i) / (steps - 1);
    return html('li', {}, html('i', { style: `background:${color};opacity:${heat(v, lo, hi)}` }), format(v));
  }));

export const strips = ({
  mount, rows, lo, hi, aria, describe, legend,
}) => {
  const W = widthOf(mount);
  const L = 26;
  const R = 50;
  const rowH = W < 560 ? 11 : 14;
  const bh = rowH - 3;
  const most = Math.max(...rows.map((r) => r.bins.length));
  const cw = (W - L - R) / most;
  const H = rows.length * rowH + 2;
  const flat = [];
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    rows.map((r, i) => el('g', { 'data-chrom': r.key, transform: `translate(0 ${i * rowH})` },
      el('text', { class: 'lab', x: 0, y: bh / 2 + 4 }, r.key),
      r.bins.map((v, j) => {
        if (v === null) return null;
        const k = flat.push({ row: r, index: j, value: v }) - 1;
        return el('rect', {
          class: 'cell', 'data-bin': `${r.key}:${j}`, 'data-k': k, 'data-rank': k, x: L + j * cw, y: 0, width: cw + 0.4, height: bh, style: `fill:var(--ink);fill-opacity:${heat(v, lo, hi)};--i:${i}`,
        });
      }),
      el('text', { class: 'val', x: L + r.bins.length * cw + 6, y: bh / 2 + 4 }, r.right))));
  hover(svg, (k) => describe(flat[Number(k)]));
  return mountSvg(mount, svg, legend);
};

export const matrix = ({
  mount, rows, cols, values, aria, describe, format, center = 1,
}) => {
  const W = widthOf(mount);
  const L = 24;
  const T = 22;
  const side = Math.min(64, (Math.min(W, 420) - L) / cols.length);
  const flatValues = values.flat();
  const below = Math.max(...flatValues.map((v) => center - v), 1e-9);
  const above = Math.max(...flatValues.map((v) => v - center), 1e-9);
  const H = T + rows.length * side + 2;
  const svg = el('svg', { viewBox: `0 0 ${L + cols.length * side} ${H}`, role: 'img', 'aria-label': aria, style: `max-width:${L + cols.length * side}px` },
    cols.map((c, j) => el('text', { class: 'tick', x: L + j * side + side / 2, y: 14, 'text-anchor': 'middle' }, c)),
    rows.map((r, i) => el('g', {},
      el('text', { class: 'lab', x: 0, y: T + i * side + side / 2 + 4 }, r),
      cols.map((c, j) => {
        const v = values[i][j];
        const color = v < center ? 'var(--s2)' : 'var(--s1)';
        const shade = v < center ? heat(center - v, 0, below) : heat(v - center, 0, above);
        return el('g', { 'data-cell': `${r}${c}`, 'data-k': i * cols.length + j, 'data-rank': i * cols.length + j },
          el('rect', {
            class: 'cell', x: L + j * side, y: T + i * side, width: side - 2, height: side - 2, style: `fill:${color};fill-opacity:${shade};--i:${i * 4 + j}`,
          }),
          el('text', { class: 'cell-lab', x: L + j * side + (side - 2) / 2, y: T + i * side + (side - 2) / 2 + 4, 'text-anchor': 'middle' }, format(v)));
      }))));
  hover(svg, (k) => describe(Math.floor(Number(k) / cols.length), Number(k) % cols.length));
  return mountSvg(mount, svg);
};

export const stackedArea = ({
  mount, years, series, aria, format = exact, max = null,
}) => {
  const W = widthOf(mount);
  const L = 46;
  const R = 12;
  const T = 12;
  const B = 24;
  const H = 280;
  const rows = years.map((year, i) => Object.fromEntries([['year', year], ...series.map((s) => [s.key, s.values[i]])]));
  const layers = stack(rows, series.map((s) => s.key));
  const totals = rows.map((r) => series.reduce((n, s) => n + r[s.key], 0));
  const ticks = niceTicks(max ?? Math.max(...totals), 4);
  const x = scaleLinear([years[0], years.at(-1)], [L, W - R]);
  const y = scaleLinear([0, ticks.at(-1)], [H - B, T]);
  const area = (points) => {
    const upper = points.map(([, y1], i) => `${x(years[i]).toFixed(1)},${y(y1).toFixed(1)}`);
    const lower = points.map(([y0], i) => `${x(years[i]).toFixed(1)},${y(y0).toFixed(1)}`).reverse();
    return `M${upper.join('L')}L${lower.join('L')}Z`;
  };
  const step = labelStep(years.length, W - L, 40);
  const cross = el('line', { class: 'cross', y1: T, y2: H - B, x1: 0, x2: 0, visibility: 'hidden' });
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: t === 0 ? 'base' : 'grid', x1: L, x2: W - R, y1: y(t), y2: y(t) }),
      el('text', { class: 'tick', x: L - 8, y: y(t) + 4, 'text-anchor': 'end' }, compact(t)))),
    layers.map((layer, i) => el('path', {
      class: 'area', 'data-key': layer.key, d: area(layer.points), style: `fill:${series[i].color};--i:${i * 8}`,
    })),
    years.map((yr) => (yr % step === 0 ? el('text', { class: 'tick', x: x(yr), y: H - 6, 'text-anchor': 'middle' }, yr) : null)),
    cross,
    el('rect', { class: 'hit', 'data-k': 0, x: L, y: T, width: W - L - R, height: H - T - B }));
  hover(svg, (k, event) => {
    const box = svg.getBoundingClientRect();
    const px = ((event.clientX ?? box.left + box.width / 2) - box.left) * (W / box.width);
    const i = event.index ?? Math.max(0, Math.min(years.length - 1, Math.round(((px - L) / (W - L - R)) * (years.length - 1))));
    cross.setAttribute('x1', x(years[i]));
    cross.setAttribute('x2', x(years[i]));
    cross.setAttribute('visibility', 'visible');
    return [
      html('b', {}, `${years[i]}: ${format(totals[i])}`),
      ...[...series].reverse().map((s) => html('span', {}, html('i', { style: `background:${s.color}` }), `${s.label} ${format(rows[i][s.key])}`)),
    ];
  }, { steps: years.length });
  return mountSvg(mount, svg, legendOf(series));
};

export const waffle = ({ mount, parts, aria, format = exact }) => {
  const total = parts.reduce((n, p) => n + p.value, 0);
  const cells = waffleCells(parts.map((p, i) => ({ key: i, value: p.value })));
  const side = Math.min(widthOf(mount), 300) / 10;
  const gap = 2;
  const svg = el('svg', { viewBox: `0 0 ${side * 10} ${side * 10}`, role: 'img', 'aria-label': aria, style: `max-width:${side * 10}px` },
    cells.map((k, n) => el('rect', {
      class: 'cell',
      x: (n % 10) * side,
      y: Math.floor(n / 10) * side,
      width: side - gap,
      height: side - gap,
      'data-k': k,
      style: `fill:${parts[k].color};--i:${n}`,
    })));
  hover(svg, (k) => {
    const p = parts[Number(k)];
    return tipLines(`${((p.value / total) * 100).toFixed(1)}%`, p.label, format(p.value));
  });
  mount.replaceChildren(svg, legendOf(parts.map((p) => ({ label: `${p.label} ${((p.value / total) * 100).toFixed(1)}%`, color: p.color }))));
  return svg;
};

export const treemap = ({ mount, items, height, aria, describe }) => {
  const W = widthOf(mount);
  const H = height ?? Math.round(Math.min(460, Math.max(260, W * 0.5)));
  const rects = squarify(items, W, H);
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    rects.map((r, i) => el('g', { 'data-key': r.key, 'data-k': items.findIndex((it) => it.key === r.key) },
      el('rect', { class: 'cell', x: r.x, y: r.y, width: Math.max(0, r.w - 1), height: Math.max(0, r.h - 1), style: `fill:${r.color ?? 'var(--data-2)'};--i:${i * 3}` }),
      r.w > 46 && r.h > 22 ? el('text', { class: 'tile-lab', x: r.x + 7, y: r.y + 17 }, r.label) : null,
      r.w > 74 && r.h > 40 ? el('text', { class: 'tile-sub', x: r.x + 7, y: r.y + 33 }, r.sub) : null)));
  hover(svg, (k) => describe(items[Number(k)]));
  return mountSvg(mount, svg);
};

export const lines = ({
  mount, series, aria, describe, annotations = [], scale = 'log', yFormat = compact, height = 300, labelW = 92, logStep = 3,
}) => {
  const W = widthOf(mount);
  const L = 46;
  const R = labelW;
  const T = 16 + annotations.length * 18;
  const B = 24;
  const H = height;
  const xs = [...new Set(series.flatMap((s) => s.rows.map((r) => r.x)))].sort((a, b) => a - b);
  const log = scale === 'log';
  const all = series.flatMap((s) => s.rows.map((r) => r.y));
  const ys = all.filter((v) => v > 0);
  const neg = log ? 0 : Math.min(0, ...all);
  const up = niceTicks(Math.max(...ys), 4);
  const unit = up[1] ?? 1;
  const down = neg < 0 ? Array.from({ length: Math.ceil(-neg / unit) }, (_, i) => -(i + 1) * unit).reverse() : [];
  const lo = log ? 10 ** Math.floor(Math.log10(Math.min(...ys))) : (down[0] ?? 0);
  const hi = log ? 10 ** Math.ceil(Math.log10(Math.max(...ys))) : up.at(-1);
  const x = scaleLinear([xs[0], xs.at(-1)], [L, W - R]);
  const y = log ? scaleLog([lo, hi], [H - B, T]) : scaleLinear([lo, hi], [H - B, T]);
  const ticks = log ? logTicks(lo, hi).filter((t) => Math.round(Math.log10(t)) % logStep === 0) : [...down, ...up];
  const dx = xs.length > 1 ? Math.min(...xs.slice(1).map((v, k) => v - xs[k])) : 1;
  const step = labelStep(Math.round((xs.at(-1) - xs[0]) / dx) + 1, W - L - R, 40);
  const labelY = spread(series.map((s) => y(s.rows.at(-1).y)), 14, T, H - B);
  const cross = el('line', { class: 'cross', y1: T, y2: H - B, x1: 0, x2: 0, visibility: 'hidden' });
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: t === 0 && lo < 0 ? 'base' : 'grid', x1: L, x2: W - R, y1: y(t), y2: y(t) }),
      el('text', { class: 'tick', x: L - 8, y: y(t) + 4, 'text-anchor': 'end' }, yFormat(t)))),
    el('line', { class: 'base', x1: L, x2: W - R, y1: H - B, y2: H - B }),
    xs.map((v) => (v % (step * dx) === 0 ? el('text', { class: 'tick', x: x(v), y: H - 6, 'text-anchor': 'middle' }, v) : null)),
    series.map((s, n) => {
      const last = s.rows.at(-1);
      const ly = labelY[n];
      const moved = Math.abs(ly - y(last.y)) > 2;
      return el('g', {},
        el('path', { class: `line${s.alt ? ' alt' : ''}${s.cls ? ` ${s.cls}` : ''}`, d: `M${s.rows.filter((r) => !log || r.y > 0).map((r) => `${x(r.x).toFixed(1)},${y(r.y).toFixed(1)}`).join('L')}`, style: s.color ? `stroke:${s.color}` : null }),
        moved ? el('path', { class: 'leader', d: `M${x(last.x) + 2},${y(last.y)}L${x(last.x) + 7},${ly}` }) : null,
        el('text', { class: 'lab', x: x(last.x) + (moved ? 10 : 8), y: ly + 4 }, s.label));
    }),
    annotations.map((a, n) => {
      const s = series[a.series ?? 0];
      const r = s.rows.find((p) => p.x === a.x);
      if (!r) return null;
      const ty = 14 + n * 18;
      return el('g', {},
        el('path', { class: 'note-line', d: `M${x(r.x)},${ty + 4}V${y(r.y) - 6}` }),
        el('text', { class: 'note', y: ty + 4, ...noteAt(x(r.x), a.text, W) }, a.text));
    }),
    cross,
    el('rect', { class: 'hit', 'data-k': 0, x: L, y: T, width: W - L - R, height: H - T - B }));
  hover(svg, (k, event) => {
    const box = svg.getBoundingClientRect();
    const px = ((event.clientX ?? box.left + box.width / 2) - box.left) * (W / box.width);
    const i = event.index ?? Math.max(0, Math.min(xs.length - 1, Math.round(((px - L) / (W - L - R)) * (xs.length - 1))));
    cross.setAttribute('x1', x(xs[i]));
    cross.setAttribute('x2', x(xs[i]));
    cross.setAttribute('visibility', 'visible');
    return describe(xs[i]);
  }, { steps: xs.length });
  return mountSvg(mount, svg);
};

export const dumbbell = ({
  mount, rows, years, describe, aria, format = compact, tickFormat = compact,
}) => {
  const W = widthOf(mount);
  const R = 60;
  const T = 22;
  const lay = rowLayout(rows, Math.min(170, Math.max(96, W * 0.3)), 24);
  const { L, rowH } = lay;
  const H = T + rows.length * rowH + 6;
  const values = rows.flatMap((r) => [r.a, r.b]).filter((v) => v > 0);
  const lo = 10 ** Math.floor(Math.log10(Math.min(...values)));
  const hi = 10 ** Math.ceil(Math.log10(Math.max(...values)));
  const x = scaleLog([lo, hi], [L, W - R]);
  const ticks = logTicks(lo, hi);
  const svg = el('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': aria },
    ticks.map((t) => el('g', {},
      el('line', { class: 'grid', x1: x(t), x2: x(t), y1: T - 6, y2: H }),
      el('text', { class: 'tick', x: x(t), y: 10, 'text-anchor': 'middle' }, tickFormat(t)))),
    el('g', { class: lay.cls }, rows.map((r, i) => {
      const mid = T + i * rowH + rowH / 2;
      const y = mid + lay.markY;
      return el('g', { 'data-row': r.key, 'data-k': i },
        el('rect', { class: 'row-hit', x: 0, y: mid - rowH / 2, width: W, height: rowH }),
        el('text', { class: 'lab', x: 0, y: mid + lay.labelY }, r.label),
        el('line', { class: 'stem', x1: x(r.a), x2: x(r.b), y1: y, y2: y, style: `--i:${i}` }),
        el('circle', { class: 'dot then', cx: x(r.a), cy: y, r: 4, style: `--i:${i}` }),
        el('circle', { class: 'dot hi', cx: x(r.b), cy: y, r: 4.5, style: `--i:${i}` }),
        el('text', { class: 'val', x: Math.max(x(r.a), x(r.b)) + 9, y: y + 4 }, format(r.b)));
    })));
  hover(svg, (k) => describe(rows[Number(k)]));
  const legend = html('ul', { class: 'legend' },
    html('li', {}, html('i', { class: 'dot-key then' }), String(years[0])),
    html('li', {}, html('i', { class: 'dot-key now' }), String(years[1])));
  mountSvg(mount, svg, legend);
  return svg;
};

export const observeFigures = () => {
  const figs = [...document.querySelectorAll('.fig')];
  if (reduced() || !('IntersectionObserver' in window)) {
    figs.forEach((f) => f.classList.add('in'));
    return;
  }
  const io = new IntersectionObserver((entries) => {
    for (const entry of entries) {
      if (entry.isIntersecting) {
        entry.target.classList.add('in');
        io.unobserve(entry.target);
      }
    }
  }, { threshold: 0.15 });
  figs.forEach((f) => io.observe(f));
};

export const countUp = (node, value, format = exact) => {
  if (reduced() || !('IntersectionObserver' in window)) {
    node.textContent = format(value);
    return;
  }
  node.textContent = format(0);
  const io = new IntersectionObserver((entries) => {
    if (!entries.some((e) => e.isIntersecting)) return;
    io.disconnect();
    const t0 = performance.now();
    const tick = (now) => {
      const t = Math.min(1, (now - t0) / 1500);
      node.textContent = format(Math.round(value * easeOut(t)));
      if (t < 1) requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  });
  io.observe(node);
};

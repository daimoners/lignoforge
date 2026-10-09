import { useEffect, useRef, useState } from "react";
import { PALETTE } from "../lib/format";

function useWidth<T extends HTMLElement>() {
  const ref = useRef<T>(null);
  const [w, setW] = useState(640);
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver(([e]) => setW(Math.max(280, e.contentRect.width)));
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, []);
  return [ref, w] as const;
}

export function niceTicks(min: number, max: number, n = 6): number[] {
  if (!isFinite(min) || !isFinite(max)) return [0, 1];
  if (min === max) { const d = Math.abs(min) || 1; min -= d / 2; max += d / 2; }
  const raw = (max - min) / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) ?? 10 * mag;
  const out: number[] = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) out.push(Number(v.toPrecision(12)));
  return out;
}
const tick = (v: number) => (Math.abs(v) >= 1e4 || (v !== 0 && Math.abs(v) < 1e-2) ? v.toExponential(1) : String(Number(v.toPrecision(4))));

export interface Series { name: string; x: number[]; y: number[]; color?: string }

export function LineChart({ series, xLabel, yLabel, height = 300, meanLine, ariaLabel }: {
  series: Series[]; xLabel: string; yLabel: string; height?: number; meanLine?: boolean; ariaLabel?: string;
}) {
  const [ref, width] = useWidth<HTMLDivElement>();
  const [hover, setHover] = useState<number | null>(null);
  const m = { l: 62, r: 16, t: 12, b: 42 };
  const iw = width - m.l - m.r, ih = height - m.t - m.b;
  const all = series.flatMap(s => s.y), xs = series.flatMap(s => s.x);
  if (!all.length) return <div ref={ref} className="muted">No data</div>;
  let ymin = Math.min(...all), ymax = Math.max(...all);
  const pad = (ymax - ymin) * 0.06 || Math.abs(ymax) * 0.05 || 1;
  ymin -= pad; ymax += pad;
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const sx = (v: number) => m.l + (xmax === xmin ? iw / 2 : ((v - xmin) / (xmax - xmin)) * iw);
  const sy = (v: number) => m.t + ih - ((v - ymin) / (ymax - ymin)) * ih;
  const yt = niceTicks(ymin, ymax, 5), xt = niceTicks(xmin, xmax, Math.max(3, Math.floor(iw / 90)));
  const stride = Math.max(1, Math.ceil(Math.max(...series.map(s => s.x.length)) / 1500));
  const ref0 = series[0];
  const hi = hover !== null ? Math.min(hover, ref0.x.length - 1) : null;

  return (
    <div ref={ref} style={{ position: "relative" }}>
      <svg className="chart" width={width} height={height} role="img" aria-label={ariaLabel ?? `${yLabel} versus ${xLabel}`}
        onMouseLeave={() => setHover(null)}
        onMouseMove={e => {
          const r = e.currentTarget.getBoundingClientRect();
          const x = ((e.clientX - r.left - m.l) / iw) * (xmax - xmin) + xmin;
          let best = 0, bd = Infinity;
          ref0.x.forEach((v, i) => { const d = Math.abs(v - x); if (d < bd) { bd = d; best = i; } });
          setHover(best);
        }}>
        {yt.map(v => <g key={`y${v}`}><line className="grid-line" x1={m.l} x2={width - m.r} y1={sy(v)} y2={sy(v)} />
          <text x={m.l - 8} y={sy(v) + 4} textAnchor="end">{tick(v)}</text></g>)}
        {xt.map(v => <g key={`x${v}`}><line className="axis" x1={sx(v)} x2={sx(v)} y1={m.t + ih} y2={m.t + ih + 4} />
          <text x={sx(v)} y={m.t + ih + 17} textAnchor="middle">{tick(v)}</text></g>)}
        <line className="axis" x1={m.l} x2={width - m.r} y1={m.t + ih} y2={m.t + ih} />
        <line className="axis" x1={m.l} x2={m.l} y1={m.t} y2={m.t + ih} />
        <text x={m.l + iw / 2} y={height - 6} textAnchor="middle">{xLabel}</text>
        <text transform={`translate(14 ${m.t + ih / 2}) rotate(-90)`} textAnchor="middle">{yLabel}</text>
        {series.map((s, k) => {
          const color = s.color ?? PALETTE[k % PALETTE.length];
          const pts = s.x.filter((_, i) => i % stride === 0).map((x, i) => `${sx(x).toFixed(1)},${sy(s.y[i * stride]).toFixed(1)}`).join(" ");
          const mean = s.y.reduce((a, b) => a + b, 0) / s.y.length;
          return <g key={s.name}>
            <polyline points={pts} fill="none" stroke={color} strokeWidth="1.6" strokeLinejoin="round" />
            {meanLine && <line x1={m.l} x2={width - m.r} y1={sy(mean)} y2={sy(mean)} stroke={color} strokeDasharray="5 4" opacity=".6" />}
          </g>;
        })}
        {hi !== null && <line x1={sx(ref0.x[hi])} x2={sx(ref0.x[hi])} y1={m.t} y2={m.t + ih} stroke="var(--faint)" strokeDasharray="3 3" />}
      </svg>
      {hi !== null && <div className="chart-tip" style={{ left: Math.min(sx(ref0.x[hi]) + 12, width - 190), top: m.t + 4 }}>
        <div className="faint">{xLabel}: {tick(ref0.x[hi])}</div>
        {series.map((s, k) => <div key={s.name}><span className="swatch" style={{ background: s.color ?? PALETTE[k % PALETTE.length] }} />
          {series.length > 1 ? `${s.name}: ` : ""}{tick(s.y[Math.min(hi, s.y.length - 1)])}</div>)}
      </div>}
    </div>
  );
}

/** Square matrix heat map with a single-hue sequential scale (0 … max). */
export function Heatmap({ matrix, labels, labelColors, max = 1, valueLabel = "value" }: {
  matrix: number[][]; labels: string[]; labelColors?: string[]; max?: number; valueLabel?: string;
}) {
  const [ref, width] = useWidth<HTMLDivElement>();
  const [tip, setTip] = useState<{ i: number; j: number } | null>(null);
  const n = matrix.length, m = 34;
  const cell = Math.max(14, Math.min(46, (Math.min(width, 620) - m) / n));
  const size = m + cell * n;
  const color = (v: number) => {
    const t = Math.max(0, Math.min(1, v / max));
    return `color-mix(in srgb, var(--primary) ${Math.round(t * 100)}%, var(--surface-2))`;
  };
  return (
    <div ref={ref} style={{ position: "relative" }}>
      <svg className="chart" width={size} height={size} role="img" aria-label="Inter-monomer contact map"
        onMouseLeave={() => setTip(null)}>
        {labels.map((l, i) => <g key={l}>
          <text x={m + cell * i + cell / 2} y={m - 8} textAnchor="middle" style={{ fill: labelColors?.[i] }}>{l}</text>
          <text x={m - 8} y={m + cell * i + cell / 2 + 4} textAnchor="end" style={{ fill: labelColors?.[i] }}>{l}</text></g>)}
        {matrix.map((row, i) => row.map((v, j) =>
          <rect key={`${i}-${j}`} x={m + cell * j} y={m + cell * i} width={cell - 1} height={cell - 1} rx="2"
            fill={i === j ? "var(--surface-2)" : color(v)} onMouseEnter={() => setTip({ i, j })} />))}
      </svg>
      {tip && <div className="chart-tip" style={{ left: m + cell * tip.j + cell + 6, top: m + cell * tip.i }}>
        monomers {labels[tip.i]}–{labels[tip.j]}: {matrix[tip.i][tip.j].toFixed(3)} {valueLabel}</div>}
      <div className="row small muted" style={{ marginTop: 8 }}>
        0 <div style={{ width: 120, height: 8, borderRadius: 3, background: "linear-gradient(90deg, var(--surface-2), var(--primary))", border: "1px solid var(--border)" }} /> {max}
      </div>
    </div>
  );
}

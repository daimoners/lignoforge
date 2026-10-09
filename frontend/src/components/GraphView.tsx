import { useMemo } from "react";
import type { ChainGraph } from "../api/types";
import { MONOMER_COLOR, linkLabel } from "../lib/format";

/** Node-link drawing of the monomer graph: longest path on a line, branches hang off it. */
export function GraphView({ graph }: { graph: ChainGraph }) {
  const layout = useMemo(() => {
    const adj = new Map<number, number[]>();
    graph.nodes.forEach(n => adj.set(n.id, []));
    graph.edges.forEach(e => { adj.get(e.source)!.push(e.target); adj.get(e.target)!.push(e.source); });
    const bfs = (s: number) => {
      const d = new Map([[s, 0]]), par = new Map<number, number>(); const q = [s];
      while (q.length) { const u = q.shift()!; for (const v of adj.get(u)!) if (!d.has(v)) { d.set(v, d.get(u)! + 1); par.set(v, u); q.push(v); } }
      return { d, par };
    };
    const a = [...bfs(graph.nodes[0].id).d.entries()].sort((x, y) => y[1] - x[1])[0][0];
    const { d, par } = bfs(a);
    const b = [...d.entries()].sort((x, y) => y[1] - x[1])[0][0];
    const path: number[] = [b];
    while (path[path.length - 1] !== a) path.push(par.get(path[path.length - 1])!);
    path.reverse();
    const pos = new Map<number, { x: number; y: number }>();
    path.forEach((id, i) => pos.set(id, { x: 50 + i * 86, y: 80 }));
    const placed = new Set(path), used = new Map<number, number>();
    const rest = graph.nodes.map(n => n.id).filter(id => !placed.has(id));
    let guard = 0;
    while (rest.length && guard++ < 1000) {
      const id = rest.shift()!;
      const p = adj.get(id)!.find(n => pos.has(n));
      if (p === undefined) { rest.push(id); continue; }
      const k = (used.get(p) ?? 0); used.set(p, k + 1);
      const pp = pos.get(p)!;
      pos.set(id, { x: pp.x + (k % 2 ? 36 : -36) * (1 + Math.floor(k / 2) * 0.4), y: pp.y + (k % 2 ? 62 : -62) });
    }
    const w = Math.max(...[...pos.values()].map(p => p.x)) + 50;
    const h = Math.max(...[...pos.values()].map(p => p.y)) + 50;
    return { pos, w, h: Math.max(h, 150) };
  }, [graph]);

  return (
    <div style={{ overflow: "auto" }}>
      <svg width={layout.w} height={layout.h} role="img" aria-label="Monomer connectivity graph" className="chart">
        {graph.edges.map((e, i) => {
          const a = layout.pos.get(e.source)!, b = layout.pos.get(e.target)!;
          const cc = !e.linkage.includes("O");
          return <g key={i}>
            <line x1={a.x} y1={a.y} x2={b.x} y2={b.y} stroke="var(--sage)" strokeWidth="2" strokeDasharray={cc ? "5 4" : undefined} />
            <text x={(a.x + b.x) / 2} y={(a.y + b.y) / 2 - 7} textAnchor="middle" style={{ fill: "var(--text)", fontSize: 10.5 }}>{linkLabel(e.linkage)}</text>
          </g>;
        })}
        {graph.nodes.map(n => {
          const p = layout.pos.get(n.id)!;
          return <g key={n.id}><circle cx={p.x} cy={p.y} r="17" fill={MONOMER_COLOR[n.type]} />
            <text x={p.x} y={p.y + 4} textAnchor="middle" style={{ fill: "#fff", fontSize: 12, fontWeight: 600 }}>{n.type}</text>
            <text x={p.x} y={p.y + 31} textAnchor="middle" style={{ fontSize: 10 }}>{n.id}</text></g>;
        })}
      </svg>
    </div>
  );
}

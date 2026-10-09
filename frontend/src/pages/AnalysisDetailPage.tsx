import { Link, useNavigate, useParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import type { AnalysisRecord } from "../api/types";
import { Alert, Card, ConfirmButton, KV, Loading, PageHead, StatusBadge } from "../components/ui";
import { Heatmap, LineChart } from "../components/charts";
import { Topbar } from "../components/Layout";
import { download, fmt } from "../lib/format";
import { useAnalysis, useCatalog, useProjectName, useRun } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";
import { MONOMER_COLOR } from "../lib/format";
import { useGraph } from "../lib/hooks";

export function AnalysisDetailPage() {
  const p = useProjectName();
  const aid = useParams().aid!;
  const base = `/p/${encodeURIComponent(p)}`;
  const a = useAnalysis(p, aid);
  const cat = useCatalog();
  const qc = useQueryClient(), nav = useNavigate(), toast = useToast();
  const del = useMutation({ mutationFn: () => api.deleteAnalysis(p, aid),
    onSuccess: () => { qc.invalidateQueries({ queryKey: ["analyses", p] }); nav(`${base}/analysis`); }, onError: e => toast(errMsg(e), "error") });
  const rec = a.data;
  if (a.isLoading || !rec) return <Loading />;
  const label = cat.data?.observables.find(o => o.id === rec.observable)?.label ?? rec.observable;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/analysis`, label: "Analysis" }, label]} />
    <div className="page">
      <PageHead title={<span className="row" style={{ gap: 12 }}>{label} <StatusBadge status={rec.status} /></span>}
        subtitle={<span>Run <Link to={`${base}/runs/${rec.run_id}`}>{rec.run_id}</Link> · <span className="mono">{rec.id}</span></span>}
        actions={<>
          {rec.result && <button className="btn" onClick={() => download(`${rec.id}_${rec.observable}.json`, JSON.stringify(rec, null, 2), "application/json")}>Export JSON</button>}
          {rec.result && Object.keys(rec.result.series).length > 0 && <button className="btn" onClick={() => download(`${rec.id}_${rec.observable}.csv`, toCsv(rec), "text/csv")}>Export CSV</button>}
          <ConfirmButton label="Delete" onConfirm={() => del.mutate()} /></>} />
      {rec.status === "running" && <Loading text="Computing…" />}
      {rec.status === "failed" && <Alert kind="error" title="The analysis failed.">{rec.error}</Alert>}
      {rec.result && <Results p={p} rec={rec} />}
    </div>
  </>;
}

function toCsv(rec: AnalysisRecord): string {
  const r = rec.result!, keys = Object.keys(r.series);
  return [["time_ps", ...keys].join(","), ...r.time_ps.map((t, i) => [t, ...keys.map(k => r.series[k][i])].join(","))].join("\n");
}

function Results({ p, rec }: { p: string; rec: AnalysisRecord }) {
  const r = rec.result!;
  const keys = Object.keys(r.series);
  const x = r.time_ps;
  return <div className="stack">
    {keys.length > 0 && <div className="grid" style={{ gridTemplateColumns: keys.length > 1 ? "repeat(auto-fit, minmax(420px, 1fr))" : "1fr" }}>
      {keys.map(k => <Card key={k} title={k}>
        <LineChart series={[{ name: k, x, y: r.series[k] }]} xLabel="Time (ps)" yLabel={k} meanLine /></Card>)}</div>}
    {keys.length > 0 && <Card title="Summary statistics" pad={false}><table>
      <thead><tr><th>Quantity</th><th className="num">Mean</th><th className="num">Std. dev.</th><th className="num">Min</th><th className="num">Max</th><th className="num">Mean, 2nd half</th></tr></thead>
      <tbody>{keys.map(k => { const s = r.summary[k]; return <tr key={k}><td>{k}</td>
        <td className="num">{fmt(s.mean)}</td><td className="num">{fmt(s.std)}</td><td className="num">{fmt(s.min)}</td><td className="num">{fmt(s.max)}</td><td className="num">{fmt(s.mean_second_half)}</td></tr>; })}</tbody></table></Card>}
    {rec.observable === "rdf" && <Card title="Radial distribution function">
      <LineChart series={[{ name: "g(r)", x: r.extra.r as number[], y: r.extra.g as number[] }]} xLabel="r (Å)" yLabel="g(r)" />
      <p className="muted small">{String(r.extra.group_a)} ({String(r.extra.n_a)} atoms) – {String(r.extra.group_b)} ({String(r.extra.n_b)} atoms)</p></Card>}
    {rec.observable === "contacts" && <ContactMap p={p} rec={rec} />}
    {rec.observable === "end_to_end" && <Alert kind="info">Measured between the centres of mass of monomers {(r.extra.monomers as number[]).join(" and ")} (the two most distant monomers of the chain).</Alert>}
    {rec.observable === "density" && <Alert kind="info">Density of the whole simulation box. For a solvated system this is the solution density, not the density of the lignin alone.</Alert>}
    {Array.isArray(r.extra.missing) && r.extra.missing.length > 0 && <Alert kind="info">Not present in this run: {(r.extra.missing as string[]).join(", ")}.</Alert>}
    <Card title="Parameters"><KV items={[...Object.entries(rec.params).map(([k, v]) => [k, String(v)] as [string, string]),
      ["frames in trajectory", String(r.extra.n_frames_total ?? "–")]]} /></Card>
  </div>;
}

function ContactMap({ p, rec }: { p: string; rec: AnalysisRecord }) {
  const run = useRun(p, rec.run_id);
  const chainId = run.data?.layout[Number(rec.params.molecule ?? 0) % (run.data?.layout.length || 1)]?.chain_id ?? run.data?.chain_ids[0];
  const graph = useGraph(p, chainId ?? "", !!chainId);
  const ex = rec.result!.extra;
  const colors = graph.data ? Object.fromEntries(graph.data.nodes.map(n => [String(n.id), MONOMER_COLOR[n.type]])) : {};
  const labels = ex.labels as string[];
  return <Card title="Inter-monomer contact map">
    <Heatmap matrix={ex.matrix as number[][]} labels={labels} labelColors={labels.map(l => colors[l])} valueLabel="of frames in contact" />
    <p className="muted small">Fraction of {String(ex.n_frames)} frames in which the minimum heavy-atom distance between two monomers is below {String(ex.cutoff)} Å.
      Labels are coloured by monomer type (H red, G green, S blue).</p>
  </Card>;
}

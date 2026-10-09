import { useMemo, useState } from "react";
import { Link, useNavigate, useSearchParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import { Alert, Card, Field, Loading, PageHead } from "../components/ui";
import { Topbar } from "../components/Layout";
import { useCatalog, useProjectName, useRuns } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

export function NewAnalysisPage() {
  const p = useProjectName();
  const base = `/p/${encodeURIComponent(p)}`;
  const [sp] = useSearchParams();
  const cat = useCatalog(), runs = useRuns(p);
  const nav = useNavigate(), qc = useQueryClient(), toast = useToast();
  const [runId, setRunId] = useState(sp.get("run") ?? "");
  const [obs, setObs] = useState("rg");
  const [vals, setVals] = useState<Record<string, string>>({});
  const [common, setCommon] = useState<Record<string, string>>({});

  const eligible = useMemo(() => (runs.data ?? []).filter(r => (r.status === "done" || r.status === "running") && r.kind !== "em"), [runs.data]);
  const info = cat.data?.observables.find(o => o.id === obs);
  const run = eligible.find(r => r.id === runId);

  const start = useMutation({
    mutationFn: () => {
      const params: Record<string, unknown> = {};
      const put = (k: string, v: string, type: string) => { if (v !== "") params[k] = type === "string" ? v : Number(v); };
      Object.entries(info?.options ?? {}).forEach(([k, o]) => put(k, vals[k] ?? "", o.type));
      Object.entries(cat.data?.common_options ?? {}).forEach(([k, o]) => put(k, common[k] ?? "", o.type));
      if (obs === "energy" && !params.stage && run?.kind === "md") params.stage = "prod";
      return api.startAnalysis(p, { run_id: runId, observable: obs, params });
    },
    onSuccess: a => { qc.invalidateQueries({ queryKey: ["analyses", p] }); nav(`${base}/analysis/${a.id}`); },
    onError: e => toast(errMsg(e), "error"),
  });

  if (cat.isLoading || runs.isLoading || !cat.data) return <Loading />;
  const field = (k: string, o: { type: string; default: number | string | null; unit?: string }, store: Record<string, string>, set: (v: Record<string, string>) => void) =>
    <Field key={k} label={`${k.replace(/_/g, " ")}${o.unit ? ` (${o.unit})` : ""}`}>
      <input type={o.type === "string" ? "text" : "number"} value={store[k] ?? ""} step="any"
        placeholder={o.default === null ? "auto" : String(o.default)} onChange={e => set({ ...store, [k]: e.target.value })} /></Field>;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/analysis`, label: "Analysis" }, "New"]} />
    <div className="page" style={{ maxWidth: 820 }}>
      <PageHead title="New analysis" subtitle="Computed on the saved trajectory of a simulation. Heavy-atom and bead definitions are described in the documentation." />
      {!eligible.length ? <Alert kind="warn">No finished simulations with a trajectory. <Link to={`${base}/runs/new`}>Start one.</Link></Alert> :
      <div className="stack">
        <Card title="1 · Simulation"><Field label="Run"><select value={runId} onChange={e => setRunId(e.target.value)} aria-label="Run">
          <option value="">Select a run…</option>
          {eligible.map(r => <option key={r.id} value={r.id}>{r.name} · {r.id} · {r.kind} · {r.status}</option>)}</select></Field>
          {run?.status === "running" && <div style={{ marginTop: 10 }}><Alert kind="info">This run is still in progress; the analysis uses the frames written so far.</Alert></div>}</Card>
        <Card title="2 · Observable">
          <div className="stack"><Field label="Observable"><select value={obs} onChange={e => { setObs(e.target.value); setVals({}); }} aria-label="Observable">
            {cat.data.observables.map(o => <option key={o.id} value={o.id}>{o.label}</option>)}</select></Field>
            {info && Object.keys(info.options).length > 0 && <div className="form-grid">{Object.entries(info.options).map(([k, o]) => field(k, o, vals, setVals))}</div>}
            {obs === "rdf" && <span className="hint muted small">Groups: <code>solute</code>, <code>solute_heavy</code>, <code>solute_O</code>, <code>water_O</code> or any MDAnalysis selection. Water-based RDFs need a solvated run.</span>}
          </div></Card>
        <Card title="3 · Frames and molecule (optional)"><div className="form-grid">
          {Object.entries(cat.data.common_options).filter(([k]) => obs !== "energy" || k === "stage").map(([k, o]) => field(k, o, common, setCommon))}</div>
          <p className="hint muted small" style={{ marginBottom: 0 }}>Leave empty for defaults: the production stage, the first molecule, all frames. “Molecule” selects among the chains of a multi-chain coarse-grained system (0, 1, 2…).</p></Card>
        <div className="row"><button className="btn primary" disabled={!runId || start.isPending} onClick={() => start.mutate()}>{start.isPending ? "Starting…" : "Run analysis"}</button></div>
      </div>}
    </div>
  </>;
}

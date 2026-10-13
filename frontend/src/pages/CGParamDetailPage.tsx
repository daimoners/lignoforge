import { useEffect, useRef } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import { Alert, Card, ConfirmButton, KV, Loading, PageHead, ProgressBar, StatusBadge } from "../components/ui";
import { LineChart } from "../components/charts";
import { Topbar } from "../components/Layout";
import { fmt, fmtDate, fmtDuration, fmtEta } from "../lib/format";
import { useCgparam, useCgparamLog, useCgparamResults, useProjectName } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

const STAGE_LABEL: Record<string, string> = {
  build: "Build chains", reference: "Atomistic reference", map: "Map & bonded fit",
  iterate: "CG optimisation", validate: "Validation (NPT)", export: "Export",
};
const FILES = ["cg_params.json", "cg_bulk.top", "cg_bulk.gro", "cg_npt.mdp", "report.md"];

export function CGParamDetailPage() {
  const p = useProjectName();
  const jid = useParams().jid!;
  const base = `/p/${encodeURIComponent(p)}`;
  const job = useCgparam(p, jid);
  const running = job.data?.status === "running";
  const res = useCgparamResults(p, jid, running);
  const log = useCgparamLog(p, jid, running);
  const qc = useQueryClient(), nav = useNavigate(), toast = useToast();
  const logRef = useRef<HTMLDivElement>(null);
  useEffect(() => { const el = logRef.current; if (el) el.scrollTop = el.scrollHeight; }, [log.data]);
  const cancel = useMutation({ mutationFn: () => api.cancelCgparam(p, jid),
    onSuccess: () => qc.invalidateQueries({ queryKey: ["cgparam", p, jid] }), onError: e => toast(errMsg(e), "error") });
  const del = useMutation({ mutationFn: () => api.deleteCgparam(p, jid),
    onSuccess: () => { qc.invalidateQueries({ queryKey: ["cgparams", p] }); nav(`${base}/cgparam`); }, onError: e => toast(errMsg(e), "error") });

  const j = job.data;
  if (job.isLoading || !j) return <Loading />;
  const h = res.data?.history ?? [];
  const v = res.data?.validation, ref = res.data?.reference, params = res.data?.params;
  const it = h.map(x => x.iteration);

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/cgparam`, label: "CG parametrisation" }, j.name]} />
    <div className="page">
      <PageHead title={<span className="row" style={{ gap: 12 }}>{j.name} <StatusBadge status={j.status} /></span>}
        subtitle={<span className="mono">{j.id} · started {fmtDate(j.started)} · {fmtDuration(j.started, j.finished)}</span>}
        actions={<>
          {running && <button className="btn" disabled={cancel.isPending} onClick={() => cancel.mutate()}>Cancel</button>}
          {!running && <ConfirmButton label="Delete" onConfirm={() => del.mutate()} />}</>} />
      {running && <div style={{ marginBottom: 16 }}><Card><div className="stack" style={{ gap: 8 }}>
        <div className="row"><b>{j.progress ? `${STAGE_LABEL[j.progress.stage] ?? j.progress.stage}` : "Starting…"}</b><span className="spacer" />
          <span className="muted small">{j.progress?.message}</span></div>
        <ProgressBar fraction={j.live?.current?.fraction ?? j.progress?.fraction} />
        <div className="muted small">
          {j.live?.current && <>Current run <span className="mono">{j.live.current.run}</span>: step {j.live.current.step}/{j.live.current.nsteps}
            {j.live.current.eta_s != null && <> · {fmtEta(j.live.current.eta_s)} left</>}
            {j.live.current.ns_per_day != null && <> · {j.live.current.ns_per_day.toFixed(0)} ns/day</>}<br /></>}
          {j.live?.detail}{j.live?.eta_s != null && <> · <b>whole job: about {fmtEta(j.live.eta_s)} left</b></>}</div></div></Card></div>}
      {j.status === "failed" && <div style={{ marginBottom: 16 }}><Alert kind="error" title="The job failed.">
        {Object.entries(j.stages).filter(([, s]) => s.status === "failed").map(([k, s]) => `${STAGE_LABEL[k]}: ${s.error}`).join(" ") || `Exit code ${j.exit_code}.`} The same command resumes where it stopped.</Alert></div>}
      {j.status === "cancelled" && <div style={{ marginBottom: 16 }}><Alert kind="warn">Cancelled. Finished stages are kept.</Alert></div>}

      <div className="grid cols-2" style={{ marginBottom: 16 }}>
        <Card title="Stages"><KV items={Object.keys(STAGE_LABEL).map(k =>
          [STAGE_LABEL[k], <StatusBadge key={k} status={j.stages[k]?.status ?? "pending"} />] as [string, React.ReactNode])} /></Card>
        <Card title="Setup"><KV items={[
          ["Scheme", j.config.scheme], ["Temperature", `${j.config.temperature} K`],
          ["Ensemble", `${j.config.ensemble.n_distinct} sequences × ${j.config.ensemble.copies} copies, ${j.config.ensemble.spec.n_monomers} monomers`],
          ["Atomistic", `${j.config.atomistic.compress_ps}+${j.config.atomistic.equil_ps}+${j.config.atomistic.prod_ps} ps`],
          ["CG", `${j.iterations}/${j.config.cg.iterations} iterations × ${j.config.cg.iter_ps} ps, cut-off ${j.config.cg.cutoff} nm, ${j.config.cg.mixing}`],
          ...(ref ? [["Reference density", `${fmt(ref.density_g_cm3, 3)} g/cm³`], ["Reference pressure", `${fmt(ref.pressure_bar, 1)} bar`]] as [string, string][] : [])]} /></Card>
      </div>

      {h.length > 0 && <div className="grid cols-2" style={{ marginBottom: 16 }}>
        <Card title="Convergence of g(r)"><LineChart xLabel="iteration" yLabel="RMSE of g(r)" height={240}
          series={[{ name: "overall", x: it, y: h.map(x => x.rdf_rmse.overall) }]} /></Card>
        <Card title="Pressure"><LineChart xLabel="iteration" yLabel="pressure (bar)" height={240}
          series={[{ name: "CG", x: it, y: h.map(x => x.pressure_bar ?? NaN) },
                   { name: "reference", x: it, y: h.map(x => x.pressure_ref_bar), color: "#888" }]} /></Card></div>}

      {v && <div style={{ marginBottom: 16 }}><Card title="Validation (constant pressure)">
        <KV items={[
          ["Density CG / reference", `${fmt(v.density_cg, 3)} / ${fmt(v.density_reference, 3)} g/cm³` + (v.checks.density_rel_err != null ? ` (${(v.checks.density_rel_err * 100).toFixed(1)}%)` : "")],
          ["g(r) RMSE", `${fmt(v.checks.rdf_rmse, 3)} (tolerance ${v.tolerances.rdf_rmse})`],
          ...Object.entries(v.rg).map(([n, r]) => [`Rg ${n}`, `${fmt(r.cg, 3)} vs ${fmt(r.reference, 3)} nm (${(r.rel_err * 100).toFixed(1)}%)`] as [string, string])]} />
        <div className="row" style={{ marginTop: 12 }}>{Object.entries(v.checks.pass).map(([k, ok]) =>
          <span key={k} className={`badge ${ok ? "ok" : "failed"}`}><span className="dot" />{k}: {ok ? "pass" : "fail"}</span>)}</div></Card></div>}

      {params && <div style={{ marginBottom: 16 }}><Card title="Non-bonded parameters (12-6 LJ)" pad={false}>
        <div className="table-wrap" style={{ border: 0 }}><table><thead><tr><th>Bead-type pair</th><th className="num">σ (nm)</th><th className="num">ε (kJ/mol)</th></tr></thead><tbody>
          {Object.entries(params.nonbonded).map(([k, q]) => <tr key={k}><td className="mono">{k}</td><td className="num">{q.sigma.toFixed(4)}</td><td className="num">{q.epsilon.toFixed(3)}</td></tr>)}</tbody></table></div>
        <div className="card-body"><p className="muted small">{params.provenance}</p>
          <div className="row">{FILES.map(f => <a key={f} className="btn small" href={api.cgparamFileUrl(p, jid, f)} download={f}>{f}</a>)}</div></div></Card></div>}

      <Card title="Log" actions={<span className="muted small">{running ? "live · last 400 lines" : "last 400 lines"}</span>}>
        <div className="log" ref={logRef} aria-label="Job log">{log.data || "(no output yet)"}</div></Card>
    </div>
  </>;
}

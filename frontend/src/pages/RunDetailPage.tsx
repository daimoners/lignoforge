import { useEffect, useRef } from "react";
import { Link, useNavigate, useParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import { Alert, Card, ConfirmButton, KV, Loading, PageHead, ProgressBar, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtBytes, fmtDate, fmtDuration } from "../lib/format";
import { useProjectName, useRun, useRunFiles, useRunLog } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

export function RunDetailPage() {
  const p = useProjectName();
  const rid = useParams().rid!;
  const base = `/p/${encodeURIComponent(p)}`;
  const run = useRun(p, rid);
  const running = run.data?.status === "running";
  const log = useRunLog(p, rid, running);
  const files = useRunFiles(p, rid, running);
  const qc = useQueryClient(), nav = useNavigate(), toast = useToast();
  const logRef = useRef<HTMLDivElement>(null);
  useEffect(() => { const el = logRef.current; if (el) el.scrollTop = el.scrollHeight; }, [log.data]);

  const cancel = useMutation({ mutationFn: () => api.cancelRun(p, rid),
    onSuccess: () => qc.invalidateQueries({ queryKey: ["run", p, rid] }), onError: e => toast(errMsg(e), "error") });
  const del = useMutation({ mutationFn: () => api.deleteRun(p, rid),
    onSuccess: () => { qc.invalidateQueries({ queryKey: ["runs", p] }); nav(`${base}/runs`); }, onError: e => toast(errMsg(e), "error") });

  const r = run.data;
  if (run.isLoading || !r) return <Loading />;
  const done = r.status === "done";
  const pr = r.progress;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/runs`, label: "Simulations" }, r.name]} />
    <div className="page">
      <PageHead title={<span className="row" style={{ gap: 12 }}>{r.name} <StatusBadge status={r.status} /></span>}
        subtitle={<span className="mono">{r.id} · started {fmtDate(r.started)} · {fmtDuration(r.started, r.finished)}</span>}
        actions={<>
          {done && r.kind !== "em" && <Link className="btn primary" to={`${base}/analysis/new?run=${r.id}`}>Analyse</Link>}
          {running && <button className="btn" disabled={cancel.isPending} onClick={() => cancel.mutate()}>Cancel run</button>}
          {!running && <ConfirmButton label="Delete" onConfirm={() => del.mutate()} />}
        </>} />
      {running && <div style={{ marginBottom: 16 }}><Card>
        <div className="stack" style={{ gap: 8 }}>
          <div className="row"><b>{pr?.stage ? `Stage: ${pr.stage}` : "Starting…"}</b><span className="spacer" />
            <span className="muted small">{pr?.time_ps != null ? `${pr.time_ps.toFixed(1)} ps` : ""}{pr?.fraction != null ? ` · ${(pr.fraction * 100).toFixed(0)}% of stage` : ""}</span></div>
          <ProgressBar fraction={pr?.fraction} /></div></Card></div>}
      {r.status === "failed" && <div style={{ marginBottom: 16 }}><Alert kind="error" title="The run failed.">
        {r.error ?? `GROMACS exited with code ${r.exit_code}.`} See the end of the log below.</Alert></div>}
      {r.status === "cancelled" && <div style={{ marginBottom: 16 }}><Alert kind="warn">This run was cancelled.</Alert></div>}
      <div className="grid cols-2" style={{ marginBottom: 16 }}>
        <Card title="Setup"><KV items={[
          ["Type", r.kind === "md" ? "Atomistic MD (OPLS-AA)" : r.kind === "em" ? "Energy minimisation" : "Coarse-grained"],
          ["Chains", r.layout.map(l => `${l.chain_id}${l.copies > 1 ? ` ×${l.copies}` : ""} (${l.n_atoms} ${r.kind === "cg" ? "beads" : "atoms"})`).join(", ")],
          ["Threads", r.threads || "all"],
          ...Object.entries(r.params).map(([k, v]) => [k.replace(/_/g, " "), String(v)] as [string, string])]} /></Card>
        <Card title="Files" pad={false}>
          {!files.data?.length ? <div className="empty small">No files yet.</div> : <div className="table-wrap" style={{ border: 0, maxHeight: 280 }}><table><tbody>
            {files.data.map(f => <tr key={f.name}><td className="mono">{f.name}</td><td className="num muted">{fmtBytes(f.size)}</td>
              <td><a href={api.runFileUrl(p, rid, f.name)} download={f.name}>Download</a></td></tr>)}</tbody></table></div>}
        </Card>
      </div>
      <Card title="Log" actions={<span className="muted small">{running ? "live · last 400 lines" : "last 400 lines"}</span>}>
        <div className="log" ref={logRef} aria-label="Simulation log">{log.data || "(no output yet)"}</div></Card>
    </div>
  </>;
}

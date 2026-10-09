import { useMemo, useState } from "react";
import { Link, Route, Routes, useNavigate, useParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import type { ChainRecord, Topology } from "../api/types";
import { Alert, Card, ConfirmButton, Empty, Loading, PageHead, Spinner, Stat, StatusBadge, Tabs } from "../components/ui";
import { DataTable, type Col } from "../components/DataTable";
import { GraphView } from "../components/GraphView";
import { Topbar } from "../components/Layout";
import { Viewer3D } from "../components/Viewer3D";
import { MONOMER_COLOR, download, fmt, fmtSigned, linkLabel, RES_TO_MONOMER } from "../lib/format";
import { useChain, useChargeReport, useGraph, useProjectName, useTopology } from "../lib/hooks";
import { OPLS_INFO } from "../lib/opls";
import { errMsg, useToast } from "../lib/toast";

export function ChainDetailPage() {
  const p = useProjectName();
  const cid = useParams().cid!;
  const base = `/p/${encodeURIComponent(p)}`;
  const chain = useChain(p, cid);
  const nav = useNavigate();
  const qc = useQueryClient();
  const toast = useToast();
  const del = useMutation({
    mutationFn: () => api.deleteChain(p, cid),
    onSuccess: () => { qc.invalidateQueries({ queryKey: ["chains", p] }); nav(`${base}/chains`); },
    onError: e => toast(errMsg(e), "error"),
  });
  const c = chain.data;
  if (chain.isLoading || !c) return <Loading />;
  const t = `${base}/chains/${cid}`;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/chains`, label: "Chains" }, c.name]} />
    <div className="page">
      <PageHead title={<span className="row" style={{ gap: 12 }}>{c.name} <StatusBadge status={c.status} /></span>}
        subtitle={<span className="mono">{c.id}</span>}
        actions={<>
          {c.status === "ready" && <Link className="btn primary" to={`${base}/runs/new?chain=${c.id}`}>Run simulation</Link>}
          <ConfirmButton label="Delete" onConfirm={() => del.mutate()} />
        </>} />
      {c.status === "building" && <Card><div className="empty"><Spinner /> <span style={{ marginLeft: 8 }}>Building the chain and relaxing its geometry… this usually takes seconds to a minute.</span></div></Card>}
      {c.status === "failed" && <Alert kind="error" title="Build failed.">{c.error}</Alert>}
      {c.status === "ready" && <>
        {c.warnings.map(w => <div key={w} style={{ marginBottom: 12 }}><Alert kind="warn">{w}</Alert></div>)}
        <Tabs items={[
          { to: t, label: "Structure", end: true }, { to: `${t}/forcefield`, label: "Force field" },
          { to: `${t}/charges`, label: "Charges" }, { to: `${t}/linkages`, label: "Linkages" }, { to: `${t}/files`, label: "Files" }]} />
        <Routes>
          <Route index element={<StructureTab p={p} c={c} />} />
          <Route path="forcefield" element={<ForceFieldTab p={p} cid={cid} />} />
          <Route path="charges" element={<ChargesTab p={p} cid={cid} />} />
          <Route path="linkages" element={<LinkagesTab p={p} cid={cid} />} />
          <Route path="files" element={<FilesTab p={p} c={c} />} />
        </Routes></>}
    </div>
  </>;
}

function WithTopology({ p, cid, children }: { p: string; cid: string; children: (t: Topology) => React.ReactNode }) {
  const t = useTopology(p, cid);
  if (t.isLoading) return <Loading text="Assigning OPLS-AA types…" />;
  if (t.error) return <Alert kind="error" title="Could not build the topology.">{errMsg(t.error)}</Alert>;
  return <>{children(t.data!)}</>;
}

function StructureTab({ p, c }: { p: string; c: ChainRecord }) {
  return <WithTopology p={p} cid={c.id}>{topo => <div className="split" style={{ gridTemplateColumns: "1fr 280px" }}>
    <Viewer3D topo={topo} />
    <div className="stack">
      <div className="grid" style={{ gridTemplateColumns: "1fr 1fr" }}>
        <Stat label="Monomers" value={c.stats.monomer_count ?? "–"} /><Stat label="Atoms" value={topo.n_atoms} />
        <Stat label="MW (g/mol)" value={fmt(c.stats.MW, 0)} /><Stat label="Net charge" value={fmtSigned(topo.net_charge, 3)} /></div>
      <Card title="Specification">
        <dl className="kv small">
          <dt>Seed</dt><dd>{c.resolved.seed}</dd>
          <dt>Branching</dt><dd>{c.resolved.branching || "linear"}</dd>
          <dt>Max degree</dt><dd>{c.stats.max_degree}</dd></dl>
        <div className="small muted" style={{ marginTop: 8, wordBreak: "break-all" }}>SMILES<br /><code>{c.stats.smiles}</code></div>
        <button className="btn small" style={{ marginTop: 8 }} onClick={() => navigator.clipboard?.writeText(c.stats.smiles ?? "")}>Copy SMILES</button>
      </Card></div></div>}</WithTopology>;
}

function ForceFieldTab({ p, cid }: { p: string; cid: string }) {
  return <WithTopology p={p} cid={cid}>{topo => <ForceFieldView topo={topo} />}</WithTopology>;
}

function ForceFieldView({ topo }: { topo: Topology }) {
  const [q, setQ] = useState("");
  const [type, setType] = useState("");
  const [res, setRes] = useState("");
  const byType = useMemoTypes(topo);
  const rows = topo.atoms.filter(a =>
    (!type || a.opls === type) && (!res || String(a.resid) === res) &&
    (!q || a.name.toLowerCase().includes(q.toLowerCase())));
  const cols: Col<Topology["atoms"][number]>[] = [
    { key: "i", header: "#", num: true, render: a => a.index, sort: a => a.index },
    { key: "name", header: "Atom", render: a => <b>{a.name}</b>, sort: a => a.name },
    { key: "res", header: "Residue", render: a => <span><span className="swatch" style={{ background: MONOMER_COLOR[RES_TO_MONOMER[a.resname] ?? "G"] }} />{a.resid} {a.resname}</span>, sort: a => a.resid },
    { key: "el", header: "El.", render: a => a.element, sort: a => a.element },
    { key: "type", header: "OPLS type", render: a => <code>{a.opls}</code>, sort: a => a.opls },
    { key: "q", header: "q (e)", num: true, render: a => fmtSigned(a.charge), sort: a => a.charge },
    { key: "raw", header: "raw q (e)", num: true, render: a => fmtSigned(a.raw_charge), sort: a => a.raw_charge },
    { key: "d", header: "Δq", num: true, render: a => { const d = a.charge - a.raw_charge; return <span className={Math.abs(d) > 0.01 ? "" : "muted"}>{fmtSigned(d, 4)}</span>; }, sort: a => a.charge - a.raw_charge },
  ];
  return <div className="stack">
    <div className="grid" style={{ gridTemplateColumns: "repeat(6, minmax(0, 1fr))" }}>
      <Stat label="Atoms" value={topo.n_atoms} /><Stat label="Bonds" value={topo.n_bonds} /><Stat label="Angles" value={topo.n_angles} />
      <Stat label="Dihedrals" value={topo.n_dihedrals} /><Stat label="1-4 pairs" value={topo.n_pairs} /><Stat label="Impropers" value={topo.n_impropers} /></div>
    <Card title="OPLS-AA types in this chain" pad={false}>
      <table><thead><tr><th>Type</th><th>Description</th><th className="num">Atoms</th><th className="num">Mean q (e)</th><th className="num">q range (e)</th></tr></thead>
        <tbody>{byType.map(t => <tr key={t.type} className="click" onClick={() => setType(type === t.type ? "" : t.type)} style={type === t.type ? { background: "var(--primary-soft)" } : undefined}>
          <td><code>{t.type}</code></td><td style={{ whiteSpace: "normal" }}>{OPLS_INFO[t.type] ?? "–"}</td><td className="num">{t.n}</td>
          <td className="num">{fmtSigned(t.mean, 3)}</td><td className="num muted">{fmtSigned(t.min, 3)} … {fmtSigned(t.max, 3)}</td></tr>)}</tbody></table>
    </Card>
    <Card title="Atoms" actions={<div className="row">
      <input type="text" placeholder="Filter by atom name" value={q} onChange={e => setQ(e.target.value)} style={{ width: 170 }} aria-label="Filter atoms" />
      <select value={res} onChange={e => setRes(e.target.value)} style={{ width: 130 }} aria-label="Residue"><option value="">All residues</option>
        {[...new Set(topo.atoms.map(a => a.resid))].map(r => <option key={r} value={r}>Residue {r}</option>)}</select>
      {type && <button className="btn small" onClick={() => setType("")}>Type: {type} ✕</button>}</div>} pad={false}>
      <DataTable rows={rows} cols={cols} rowKey={a => a.index} />
    </Card>
    <Alert kind="info">Bonded parameters (bond, angle, dihedral, improper constants) are resolved by GROMACS from <code>oplsaa.ff</code> when the topology is processed; a missing parameter stops the run with an explicit error rather than being guessed.</Alert>
  </div>;
}

function useMemoTypes(topo: Topology) {
  return useMemo(() => {
    const m = new Map<string, number[]>();
    topo.atoms.forEach(a => m.set(a.opls, [...(m.get(a.opls) ?? []), a.charge]));
    return [...m.entries()].sort().map(([type, qs]) => ({
      type, n: qs.length, mean: qs.reduce((a, b) => a + b, 0) / qs.length, min: Math.min(...qs), max: Math.max(...qs) }));
  }, [topo]);
}

function ChargesTab({ p, cid }: { p: string; cid: string }) {
  const r = useChargeReport(p, cid);
  if (r.isLoading) return <Loading />;
  if (r.error || !r.data) return <Alert kind="error">{errMsg(r.error)}</Alert>;
  const d = r.data;
  const maxMean = Math.max(...d.by_type.map(t => Math.abs(t.mean_delta)), 1e-9);
  const top = [...d.atoms].sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta)).slice(0, 10);
  return <div className="stack">
    <Alert kind="info" title="How charges are set.">Atoms receive OPLS-AA fragment charges. A residue that ends up with a small net charge is made exactly neutral by shifting <b>every atom of that residue by the same amount</b> (residual ÷ number of atoms). This report shows what was changed. Total chain charge: <b>{fmtSigned(d.net_charge, 4)} e</b>.</Alert>
    <div className="grid cols-2">
      <Card title="By residue" pad={false}>
        <table><thead><tr><th>Residue</th><th>Linkages</th><th className="num">Raw net (e)</th><th className="num">Max |Δq| (e)</th></tr></thead>
          <tbody>{d.residues.map(r => <tr key={r.residue}><td>{r.residue} <span className="muted">{r.resname}</span></td>
            <td className="muted">{r.linkages.map(linkLabel).join(", ") || "–"}</td>
            <td className="num">{fmtSigned(r.raw_net_charge)}</td><td className="num">{r.max_atom_shift.toFixed(4)}</td></tr>)}</tbody></table>
      </Card>
      <Card title="By OPLS type (mean shift)" pad={false}>
        <table><thead><tr><th>Type</th><th className="num">Atoms</th><th>Mean Δq (e)</th></tr></thead>
          <tbody>{d.by_type.map(t => <tr key={t.opls}><td><code>{t.opls}</code></td><td className="num">{t.n_atoms}</td>
            <td><div className="row" style={{ gap: 8, flexWrap: "nowrap" }}><div className="bar-bg" style={{ width: 90 }}><div className="bar" style={{ width: `${(Math.abs(t.mean_delta) / maxMean) * 100}%`, background: t.mean_delta < 0 ? "var(--info)" : "var(--brown)" }} /></div>
              <span className="mono">{fmtSigned(t.mean_delta, 5)}</span></div></td></tr>)}</tbody></table>
      </Card>
    </div>
    <Card title="Ten largest single-atom changes" pad={false} actions={<button className="btn small" onClick={() => download(`${cid}_charge_report.txt`, d.text)}>Download full report</button>}>
      <table><thead><tr><th>Atom</th><th>Residue</th><th>Type</th><th className="num">Raw (e)</th><th className="num">Final (e)</th><th className="num">Δ (e)</th></tr></thead>
        <tbody>{top.map(a => <tr key={a.atom_index}><td><b>{a.atom}</b> <span className="faint">#{a.atom_index}</span></td><td>{a.residue} {a.resname}</td>
          <td><code>{a.opls}</code></td><td className="num">{fmtSigned(a.raw_charge)}</td><td className="num">{fmtSigned(a.charge)}</td><td className="num">{fmtSigned(a.delta, 5)}</td></tr>)}</tbody></table>
    </Card>
  </div>;
}

function LinkagesTab({ p, cid }: { p: string; cid: string }) {
  const g = useGraph(p, cid);
  const t = useTopology(p, cid);
  if (g.isLoading || t.isLoading) return <Loading />;
  if (!g.data || !t.data) return <Alert kind="error">{errMsg(g.error ?? t.error)}</Alert>;
  const counts = new Map<string, number>();
  g.data.edges.forEach(e => counts.set(e.linkage, (counts.get(e.linkage) ?? 0) + 1));
  return <div className="stack">
    <Card title="Monomer connectivity"><GraphView graph={g.data} />
      <div className="row small muted" style={{ marginTop: 8 }}>
        {(["H", "G", "S"] as const).map(m => <span key={m}><span className="swatch" style={{ background: MONOMER_COLOR[m] }} />{m}</span>)}
        <span>· solid line: C–O–C linkage · dashed: C–C linkage</span></div></Card>
    <div className="grid cols-2">
      <Card title="Linkage counts" pad={false}><table><tbody>{[...counts.entries()].map(([k, v]) =>
        <tr key={k}><td>{linkLabel(k)}</td><td className="num">{v}</td></tr>)}</tbody></table></Card>
      <Card title="Inter-monomer bonds" pad={false}>
        <div className="table-wrap" style={{ border: 0, maxHeight: 280 }}><table>
          <thead><tr><th>Linkage</th><th>Atom 1</th><th>Atom 2</th></tr></thead>
          <tbody>{t.data.linkages.map((l, i) => <tr key={i}><td>{linkLabel(l.linkage)}</td>
            <td>res {l.residue_1} · {l.atom_1}</td><td>res {l.residue_2} · {l.atom_2}</td></tr>)}</tbody></table></div></Card>
    </div>
  </div>;
}

function FilesTab({ p, c }: { p: string; c: ChainRecord }) {
  const toast = useToast();
  const [files, setFiles] = useState<Record<string, string> | null>(null);
  const gen = useMutation({ mutationFn: () => api.writeTopology(p, c.id), onSuccess: setFiles, onError: e => toast(errMsg(e), "error") });
  const names = files ? Object.values(files).map(f => f.split("/").pop()!) : [];
  return <div className="stack" style={{ maxWidth: 720 }}>
    <Card title="Structure"><a className="btn" href={api.structureUrl(p, c.id)} download={`${c.id}.pdb`}>Download PDB</a></Card>
    <Card title="GROMACS input files" actions={<button className="btn primary" onClick={() => gen.mutate()} disabled={gen.isPending}>{gen.isPending ? "Generating…" : files ? "Regenerate" : "Generate files"}</button>}>
      {!files ? <Empty title="Not generated yet">Writes the OPLS-AA topology (<code>.top</code>), coordinates (<code>.gro</code>), a minimisation <code>.mdp</code> and the charge-renormalisation report.</Empty>
        : <table><tbody>{names.map(n => <tr key={n}><td className="mono">{n}</td><td><a href={api.topologyFileUrl(p, c.id, n)} download={n}>Download</a></td></tr>)}</tbody></table>}
    </Card>
  </div>;
}

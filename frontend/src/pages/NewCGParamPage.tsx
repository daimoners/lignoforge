import { useEffect, useState } from "react";
import { Link, useNavigate } from "react-router-dom";
import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import { Alert, Card, Field, Loading, NumberInput, PageHead } from "../components/ui";
import { Topbar } from "../components/Layout";
import { useProjectName, useSystem } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

type Cfg = Record<string, any>;
const clone = (o: Cfg): Cfg => JSON.parse(JSON.stringify(o));

export function NewCGParamPage() {
  const p = useProjectName();
  const base = `/p/${encodeURIComponent(p)}`;
  const nav = useNavigate(), qc = useQueryClient(), toast = useToast();
  const sys = useSystem();
  const defaults = useQuery({ queryKey: ["cgparam-defaults"], queryFn: api.cgparamDefaults, staleTime: Infinity });
  const [cfg, setCfg] = useState<Cfg | null>(null);
  const [name, setName] = useState("");
  const [threads, setThreads] = useState<number | "">(0);
  const [json, setJson] = useState("");
  const [jsonErr, setJsonErr] = useState<string | null>(null);
  useEffect(() => { if (defaults.data && !cfg) setCfg(clone(defaults.data.config)); }, [defaults.data, cfg]);
  useEffect(() => { if (cfg) setJson(JSON.stringify(cfg, null, 2)); }, [cfg]);

  const start = useMutation({
    mutationFn: () => api.startCgparam(p, { config: cfg!, name: name.trim() || undefined, threads: Number(threads) || 0 }),
    onSuccess: j => { qc.invalidateQueries({ queryKey: ["cgparams", p] }); nav(`${base}/cgparam/${j.id}`); },
    onError: e => toast(errMsg(e), "error"),
  });
  if (!defaults.data || !cfg) return <Loading />;

  const set = (path: string[], v: unknown) => { const c = clone(cfg); let o = c; path.slice(0, -1).forEach(k => { o = o[k]; }); o[path[path.length - 1]] = v; setCfg(c); };
  const num = (label: string, path: string[], hint?: string, step: number | "any" = "any") => {
    let v: any = cfg; path.forEach(k => { v = v?.[k]; });
    return <Field label={label} hint={hint}><NumberInput value={v ?? ""} step={step} min={0}
      onChange={x => set(path, x === "" ? 0 : x)} /></Field>;
  };
  const disabled = !sys.data?.features.simulations;
  const ens = cfg.ensemble, mono = ens.spec.monomers;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/cgparam`, label: "CG parametrisation" }, "New"]} />
    <div className="page" style={{ maxWidth: 900 }}>
      <PageHead title="New CG parametrisation" subtitle="Runs locally with GROMACS and continues if you close the browser. For production runs use the command line on the cluster (lignoforge-cg-param)." />
      {disabled && <div style={{ marginBottom: 16 }}><Alert kind="error" title="GROMACS not found.">Install GROMACS to run simulations. <Link to="/system">Details</Link></Alert></div>}
      <div className="stack">
        <Card title="Presets" actions={<span className="muted small">Quick test takes a few minutes and is only meant to check the tools.</span>}>
          <div className="row">
            <button className="btn" onClick={() => setCfg(clone(defaults.data!.quick))}>Quick test</button>
            <button className="btn" onClick={() => setCfg(clone(defaults.data!.config))}>Production defaults</button></div></Card>
        <Card title="1 · Model">
          <div className="grid cols-2">
            <Field label="Name (optional)"><input type="text" value={name} placeholder="e.g. kraft-softwood-500K" onChange={e => setName(e.target.value)} /></Field>
            <Field label="Mapping scheme" hint={defaults.data.schemes[cfg.scheme]}>
              <select value={cfg.scheme} onChange={e => set(["scheme"], e.target.value)}>
                {Object.keys(defaults.data.schemes).map(s => <option key={s} value={s}>{s}</option>)}</select></Field>
            {num("Temperature (K)", ["temperature"], "Reference and CG temperature; above Tg for a melt")}
            {num("Random seed", ["seed"], undefined, 1)}
          </div></Card>
        <Card title="2 · Reference ensemble" actions={<span className="muted small">{ens.n_distinct * ens.copies} chains in the box</span>}>
          <div className="grid cols-2">
            {num("Monomers per chain", ["ensemble", "spec", "n_monomers"], undefined, 1)}
            {num("Distinct chain sequences", ["ensemble", "n_distinct"], "Sequences differ by seed", 1)}
            {num("Copies of each", ["ensemble", "copies"], undefined, 1)}
            {num("Chain seed", ["ensemble", "spec", "seed"], undefined, 1)}
            {(["H", "G", "S"] as const).map(m => <Field key={m} label={`${m} fraction`}><NumberInput value={mono[m] ?? 0} step={0.05} min={0} max={1}
              onChange={x => set(["ensemble", "spec", "monomers", m], x === "" ? 0 : x)} /></Field>)}
          </div>
          <p className="muted small" style={{ marginTop: 8 }}>Linkage mix and chemical modification can be edited in the advanced JSON below.</p></Card>
        <Card title="3 · Atomistic reference run (ps)">
          <div className="grid cols-2">
            {num("Start density (g/cm³)", ["atomistic", "start_density"], "Packing density before compression")}
            {num("Compression (ps)", ["atomistic", "compress_ps"])}
            {num("Compression pressure (bar)", ["atomistic", "compress_pressure"])}
            {num("Equilibration (ps)", ["atomistic", "equil_ps"])}
            {num("Production (ps)", ["atomistic", "prod_ps"], "Sampled for the CG targets")}
            {num("Output interval (ps)", ["atomistic", "out_ps"])}</div></Card>
        <Card title="4 · CG optimisation">
          <div className="grid cols-2">
            {num("Cut-off (nm)", ["cg", "cutoff"], "Box must exceed 2.2 × cut-off")}
            {num("Iterations", ["cg", "iterations"], undefined, 1)}
            {num("CG time per iteration (ps)", ["cg", "iter_ps"])}
            {num("Validation run (ps)", ["cg", "validate_ps"], "Constant-pressure run to check density")}
            <Field label="Non-bonded mixing" hint="free: one ε,σ per bead-type pair. lb: per-type values with Lorentz-Berthelot rules (fewer parameters, better transferable).">
              <select value={cfg.cg.mixing} onChange={e => set(["cg", "mixing"], e.target.value)}><option value="free">free</option><option value="lb">lb</option></select></Field>
            <Field label="Threads (0 = all)"><NumberInput value={threads} onChange={setThreads} min={0} step={1} /></Field>
          </div></Card>
        <Card title="Advanced configuration (JSON)" actions={<span className="muted small">Applies when you leave the field</span>}>
          <textarea className="mono" rows={14} style={{ width: "100%" }} value={json} aria-label="Configuration JSON"
            onChange={e => { setJson(e.target.value); setJsonErr(null); }}
            onBlur={() => { try { setCfg(JSON.parse(json)); setJsonErr(null); } catch (e) { setJsonErr(String(e)); } }} />
          {jsonErr && <Alert kind="error">{jsonErr}</Alert>}</Card>
        <div className="row"><button className="btn primary" disabled={disabled || start.isPending || !!jsonErr} onClick={() => start.mutate()}>Start parametrisation</button>
          <Link className="btn" to={`${base}/cgparam`}>Cancel</Link></div>
      </div></div>
  </>;
}

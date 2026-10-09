import { useMemo, useState } from "react";
import { Link, useNavigate, useSearchParams } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import type { RunKind } from "../api/types";
import { Alert, Card, Field, Loading, NumberInput, PageHead, Segmented } from "../components/ui";
import { Topbar } from "../components/Layout";
import { useCgSets, useChains, useProjectName, useSystem } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

type Num = number | "";
export function NewRunPage() {
  const p = useProjectName();
  const base = `/p/${encodeURIComponent(p)}`;
  const [sp] = useSearchParams();
  const nav = useNavigate(), qc = useQueryClient(), toast = useToast();
  const chains = useChains(p), cgSets = useCgSets(p), sys = useSystem();
  const [kind, setKind] = useState<RunKind>("md");
  const [selected, setSelected] = useState<string[]>(sp.get("chain") ? [sp.get("chain")!] : []);
  const [name, setName] = useState("");
  const [threads, setThreads] = useState<Num>(0);
  const [seed, setSeed] = useState<Num>("");
  // atomistic
  const [temperature, setTemperature] = useState<Num>(300);
  const [solvent, setSolvent] = useState("tip3p");
  const [nvt, setNvt] = useState<Num>(100), [npt, setNpt] = useState<Num>(200), [prod, setProd] = useState<Num>(10);
  const [traj, setTraj] = useState<Num>(10), [box, setBox] = useState<Num>(1.2);
  // coarse-grained
  const [copies, setCopies] = useState<Num>(1), [density, setDensity] = useState<Num>(0.3), [cgNs, setCgNs] = useState<Num>(100);
  const [cgSet, setCgSet] = useState("");

  const ready = useMemo(() => (chains.data ?? []).filter(c => c.status === "ready"), [chains.data]);
  const disabled = !sys.data?.features.simulations;
  const params = useMemo(() => {
    const out: Record<string, unknown> = { threads: Number(threads) || 0 };
    if (seed !== "") out.seed = seed;
    if (kind === "md") Object.assign(out, { temperature, solvent, nvt_ps: nvt, prod_ns: prod, traj_ps: traj, box_distance_nm: box, ...(solvent !== "none" ? { npt_ps: npt } : {}) });
    if (kind === "cg") Object.assign(out, { copies, density, temperature, run_ns: cgNs, ...(cgSet ? { cg_parameters: cgSet } : {}) });
    return out;
  }, [kind, threads, seed, temperature, solvent, nvt, npt, prod, traj, box, copies, density, cgNs, cgSet]);

  const start = useMutation({
    mutationFn: () => api.startRun(p, { chain_ids: selected, kind, params, name: name.trim() || undefined }),
    onSuccess: r => { qc.invalidateQueries({ queryKey: ["runs", p] }); nav(`${base}/runs/${r.id}`); },
    onError: e => toast(errMsg(e), "error"),
  });
  const valid = selected.length >= 1 && (kind === "cg" || selected.length === 1) &&
    Object.entries(params).every(([k, v]) => k === "seed" || k === "cg_parameters" || k === "solvent" || (v !== "" && Number(v) >= 0));

  if (chains.isLoading) return <Loading />;
  const num = (label: string, v: Num, set: (x: Num) => void, hint?: string, min = 0, step: number | "any" = "any") =>
    <Field label={label} hint={hint}><NumberInput value={v} onChange={set} min={min} step={step} /></Field>;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/runs`, label: "Simulations" }, "New"]} />
    <div className="page" style={{ maxWidth: 900 }}>
      <PageHead title="New simulation" subtitle="Runs locally on this computer with GROMACS. You can close the browser; the run continues and is picked up again here." />
      {disabled && <div style={{ marginBottom: 16 }}><Alert kind="error" title="GROMACS not found.">Install GROMACS and restart LignoForge to run simulations. <Link to="/system">Details</Link></Alert></div>}
      <div className="stack">
        <Card title="1 · Type">
          <div className="stack">
            <Segmented<RunKind> label="Simulation type" value={kind} onChange={k => { setKind(k); if (k !== "cg" && selected.length > 1) setSelected(selected.slice(0, 1)); }} options={[
              { value: "md", label: "Atomistic MD" }, { value: "em", label: "Minimisation only" }, { value: "cg", label: "Coarse-grained" }]} />
            <p className="muted small" style={{ margin: 0 }}>{kind === "md"
              ? "OPLS-AA all-atom workflow: box, solvation, minimisation, NVT, NPT and production."
              : kind === "em" ? "Energy minimisation of one chain in a box. Quick check that the topology is sound."
              : "One bead per monomer, stochastic dynamics (implicit solvent). Can combine several chains and copies."}</p>
          </div>
        </Card>
        <Card title={kind === "cg" ? "2 · Chains (one or more)" : "2 · Chain"}>
          {!ready.length ? <Alert kind="warn">No finished chains. <Link to={`${base}/chains/new`}>Build one first.</Link></Alert>
            : <div className="stack" style={{ gap: 8 }}>{ready.map(c => <label key={c.id} className="check">
              <input type={kind === "cg" ? "checkbox" : "radio"} name="chain" checked={selected.includes(c.id)}
                onChange={e => setSelected(kind === "cg" ? (e.target.checked ? [...selected, c.id] : selected.filter(x => x !== c.id)) : [c.id])} />
              <b>{c.name}</b> <span className="muted">{c.stats.monomer_count} units · {c.stats.n_atoms} atoms</span></label>)}</div>}
        </Card>
        <Card title="3 · Parameters">
          <div className="form-grid">
            <Field label="Name (optional)"><input type="text" value={name} onChange={e => setName(e.target.value)} placeholder={`${kind}-${selected[0] ?? ""}`} /></Field>
            {kind !== "em" && num("Temperature (K)", temperature, setTemperature)}
            {kind === "md" && <>
              <Field label="Solvent"><select value={solvent} onChange={e => setSolvent(e.target.value)}>
                <option value="tip3p">Water (TIP3P)</option><option value="spce">Water (SPC/E)</option><option value="none">Vacuum (periodic box)</option></select></Field>
              {num("NVT equilibration (ps)", nvt, setNvt)}
              {solvent !== "none" && num("NPT equilibration (ps)", npt, setNpt)}
              {num("Production (ns)", prod, setProd)}
              {num("Frame interval (ps)", traj, setTraj, "How often coordinates are saved")}
              {num("Box padding (nm)", box, setBox, "Distance from the chain to the box edge")}</>}
            {kind === "cg" && <>
              {num("Copies of each chain", copies, setCopies, undefined, 1, 1)}
              {num("Target density (g/cm³)", density, setDensity, "Sets the box size")}
              {num("Run length (ns)", cgNs, setCgNs)}
              <Field label="CG parameters"><select value={cgSet} onChange={e => setCgSet(e.target.value)}>
                <option value="">Provisional defaults</option>
                {(cgSets.data ?? []).map(s => <option key={s.id} value={s.id}>{s.name} ({s.id})</option>)}</select></Field></>}
            <Field label="CPU threads" hint="0 = use all cores"><NumberInput value={threads} onChange={setThreads} min={0} step={1} /></Field>
            {kind !== "em" && <Field label="Random seed" hint="Leave empty to draw one; it is stored with the run"><NumberInput value={seed} onChange={setSeed} min={0} step={1} /></Field>}
          </div>
          {kind === "cg" && !cgSet && <div style={{ marginTop: 14 }}><Alert kind="warn" title="Provisional parameters.">The default coarse-grained force constants and bead sizes are placeholders. Fit them from atomistic simulations before drawing quantitative conclusions.</Alert></div>}
        </Card>
        <div className="row"><button className="btn primary" disabled={!valid || disabled || start.isPending} onClick={() => start.mutate()}>{start.isPending ? "Preparing…" : "Start simulation"}</button>
          <Link className="btn ghost" to={`${base}/runs`}>Cancel</Link></div>
      </div>
    </div>
  </>;
}

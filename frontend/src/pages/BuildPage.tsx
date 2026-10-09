import { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { keepPreviousData, useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import type { ChainSpec, Monomer, ResolvedSpec } from "../api/types";
import { Alert, Card, Field, Loading, NumberInput, PageHead } from "../components/ui";
import { Topbar } from "../components/Layout";
import { MONOMER_COLOR, MONOMER_NAME, linkLabel } from "../lib/format";
import { useCatalog, useProjectName } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

const BIOMASS = ["softwood", "hardwood", "agricultural_residue", "mixed"];
const PROCESS = ["kraft", "soda", "organosolv", "sulfite", "des", "steam_explosion"];
const DEFAULT_LINKAGES: Record<string, number> = {
  "beta-O-4": 55, "alpha-O-4": 5, "4-O-5": 5, "5-5": 10, "beta-5": 10, "beta-beta": 5,
};

function useDebounced<T>(v: T, ms = 350): T {
  const [d, setD] = useState(v);
  useEffect(() => { const t = setTimeout(() => setD(v), ms); return () => clearTimeout(t); }, [v, ms]);
  return d;
}

function Bars({ items, colors }: { items: [string, number][]; colors?: Record<string, string> }) {
  return <div className="stack" style={{ gap: 8 }}>{items.map(([k, v]) =>
    <div key={k} style={{ display: "grid", gridTemplateColumns: "92px 1fr 52px", gap: 10, alignItems: "center" }}>
      <span>{colors && <span className="swatch" style={{ background: colors[k] }} />}{k}</span>
      <div className="bar-bg"><div className="bar" style={{ width: `${v * 100}%`, background: colors?.[k] }} /></div>
      <span className="num" style={{ textAlign: "right" }}>{(v * 100).toFixed(1)}%</span></div>)}</div>;
}

export function BuildPage() {
  const p = useProjectName();
  const nav = useNavigate();
  const qc = useQueryClient();
  const toast = useToast();
  const cat = useCatalog();
  const base = `/p/${encodeURIComponent(p)}`;

  const [name, setName] = useState("");
  const [n, setN] = useState<number | "">(10);
  const [seed, setSeed] = useState<number | "">(42);
  const [branching, setBranching] = useState(0);
  const [optimize, setOptimize] = useState(true);
  const [mono, setMono] = useState<Record<Monomer, number>>({ H: 10, G: 50, S: 40 });
  const [links, setLinks] = useState<Record<string, number>>(DEFAULT_LINKAGES);
  const [biomass, setBiomass] = useState("hardwood");
  const [process, setProcess] = useState("kraft");
  const [jsonText, setJsonText] = useState("");
  const [showJson, setShowJson] = useState(false);

  const spec: ChainSpec = useMemo(() => ({
    n_monomers: n === "" ? undefined : n, seed: seed === "" ? 42 : seed, branching, optimize,
    monomers: mono, linkages: Object.fromEntries(Object.entries(links).filter(([, v]) => v > 0)),
  }), [n, seed, branching, optimize, mono, links]);
  const dspec = useDebounced(spec);
  const sumM = Object.values(mono).reduce((a, b) => a + b, 0);
  const sumL = Object.values(links).reduce((a, b) => a + b, 0);
  const valid = sumM > 0 && sumL > 0 && n !== "" && n >= 1 && n <= 200;

  const preview = useQuery({
    queryKey: ["resolve", dspec], enabled: sumM > 0 && sumL > 0, placeholderData: keepPreviousData,
    queryFn: () => api.resolveSpec(dspec),
  });

  const estimate = useMutation({
    mutationFn: async (): Promise<ResolvedSpec> => {
      let exp: Record<string, unknown>;
      if (showJson && jsonText.trim()) {
        try { exp = JSON.parse(jsonText); } catch { throw new Error("The JSON input is not valid JSON."); }
      } else exp = { material_origin: { biomass_type: biomass }, extraction_process: { process_type: process } };
      return api.resolveSpec({ experimental_input: exp, seed: seed === "" ? 42 : seed });
    },
    onSuccess: r => {
      setMono({ H: round(r.monomers.H * 100), G: round(r.monomers.G * 100), S: round(r.monomers.S * 100) });
      setLinks(Object.fromEntries(Object.entries(r.linkages).filter(([k]) => cat.data?.linkages.includes(k)).map(([k, v]) => [k, round(v * 100)])));
      toast("Composition and linkages set from literature priors. Adjust them as needed.");
    },
    onError: e => toast(errMsg(e), "error"),
  });

  const build = useMutation({
    mutationFn: () => api.buildChain(p, { ...spec, name: name.trim() || undefined }),
    onSuccess: c => { qc.invalidateQueries({ queryKey: ["chains", p] }); nav(`${base}/chains/${c.id}`); },
    onError: e => toast(errMsg(e), "error"),
  });

  if (cat.isLoading || !cat.data) return <Loading />;
  const r = preview.data;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, { to: `${base}/chains`, label: "Chains" }, "Build"]} />
    <div className="page">
      <PageHead title="Build a chain" subtitle="Define size, monomer composition and inter-unit linkages. The same specification and seed always give the same chain." />
      <div className="split" style={{ gridTemplateColumns: "minmax(340px, 1fr) minmax(300px, 420px)" }}>
        <div className="stack">
          <Card title="Estimate from experimental information">
            <div className="stack">
              <p className="muted small" style={{ margin: 0 }}>Fill the controls below from the literature priors for a biomass and extraction process. Nothing is built until you press “Build chain”.</p>
              {!showJson ? <div className="form-grid">
                <Field label="Biomass"><select value={biomass} onChange={e => setBiomass(e.target.value)}>{BIOMASS.map(b => <option key={b} value={b}>{b.replace(/_/g, " ")}</option>)}</select></Field>
                <Field label="Extraction process"><select value={process} onChange={e => setProcess(e.target.value)}>{PROCESS.map(b => <option key={b} value={b}>{b.replace(/_/g, " ")}</option>)}</select></Field>
              </div> : <Field label="LignoForge input JSON" hint="Same format as the command-line input file (NMR, GPC, elemental data, …).">
                <textarea value={jsonText} onChange={e => setJsonText(e.target.value)} placeholder='{"material_origin": {...}, "extraction_process": {...}}' spellCheck={false} /></Field>}
              <div className="row">
                <button className="btn" disabled={estimate.isPending} onClick={() => estimate.mutate()}>{estimate.isPending ? "Estimating…" : "Estimate composition"}</button>
                <button className="btn ghost small" onClick={() => setShowJson(s => !s)}>{showJson ? "Use simple form" : "Paste input JSON"}</button>
              </div>
            </div>
          </Card>

          <Card title="Size and options">
            <div className="form-grid">
              <Field label="Name (optional)"><input type="text" value={name} placeholder="e.g. HW-kraft-10mer" onChange={e => setName(e.target.value)} /></Field>
              <Field label="Number of monomers" hint="1–200"><NumberInput value={n} onChange={setN} min={1} max={200} /></Field>
              <Field label="Random seed" hint="Same seed → same chain"><NumberInput value={seed} onChange={setSeed} /></Field>
              <Field label={`Branching propensity: ${branching.toFixed(2)}`} hint="0 = strictly linear chain">
                <input type="range" min={0} max={1} step={0.05} value={branching} onChange={e => setBranching(Number(e.target.value))} aria-label="Branching propensity" /></Field>
            </div>
            <label className="check" style={{ marginTop: 12 }}><input type="checkbox" checked={optimize} onChange={e => setOptimize(e.target.checked)} />
              Relax the 3-D geometry (MMFF force field) <span className="muted">— recommended</span></label>
          </Card>

          <Card title="Monomer composition" actions={<span className="muted small">relative weights</span>}>
            <div className="stack">{cat.data.monomers.map(m => <div key={m} style={{ display: "grid", gridTemplateColumns: "150px 1fr 70px", gap: 12, alignItems: "center" }}>
              <span><span className="swatch" style={{ background: MONOMER_COLOR[m] }} />{m} · {MONOMER_NAME[m]}</span>
              <input type="range" min={0} max={100} step={1} value={mono[m]} aria-label={`${m} weight`} onChange={e => setMono({ ...mono, [m]: Number(e.target.value) })} />
              <NumberInput value={mono[m]} min={0} max={100} ariaLabel={`${m} weight value`} onChange={v => setMono({ ...mono, [m]: v === "" ? 0 : v })} /></div>)}</div>
          </Card>

          <Card title="Linkages" actions={<span className="muted small">relative weights</span>}>
            <div className="stack">{cat.data.linkages.map(l => <div key={l} style={{ display: "grid", gridTemplateColumns: "150px 1fr 70px", gap: 12, alignItems: "center" }}>
              <span>{linkLabel(l)}</span>
              <input type="range" min={0} max={100} step={1} value={links[l] ?? 0} aria-label={`${l} weight`} onChange={e => setLinks({ ...links, [l]: Number(e.target.value) })} />
              <NumberInput value={links[l] ?? 0} min={0} max={100} ariaLabel={`${l} weight value`} onChange={v => setLinks({ ...links, [l]: v === "" ? 0 : v })} /></div>)}
              <span className="hint muted small">β-1 is not supported by the structure builder and is therefore not offered.</span></div>
          </Card>
        </div>

        <div className="stack" style={{ position: "sticky", top: 70 }}>
          <Card title="Resulting specification">
            {!r ? <span className="muted">Set non-zero weights to see the specification.</span> : <div className="stack">
              <Bars items={cat.data.monomers.map(m => [m, r.monomers[m]] as [string, number])} colors={MONOMER_COLOR} />
              <hr style={{ border: 0, borderTop: "1px solid var(--border)", width: "100%" }} />
              <Bars items={Object.entries(r.linkages).filter(([, v]) => v > 0).map(([k, v]) => [linkLabel(k), v] as [string, number])} />
              {r.warnings.map(w => <Alert key={w} kind="warn">{w}</Alert>)}
              <span className="small muted">{r.n_monomers} monomers · seed {r.seed} · {r.branching > 0 ? `branching ${r.branching}` : "linear"}</span>
            </div>}
          </Card>
          <button className="btn primary" style={{ justifyContent: "center", padding: "10px 14px" }} disabled={!valid || build.isPending} onClick={() => build.mutate()}>
            {build.isPending ? "Starting…" : "Build chain"}</button>
          {!valid && <Alert kind="warn">Enter 1–200 monomers and at least one non-zero composition and linkage weight.</Alert>}
        </div>
      </div>
    </div>
  </>;
}
const round = (x: number) => Math.round(x * 10) / 10;

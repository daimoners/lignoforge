import { useEffect, useMemo, useRef, useState } from "react";
import * as $3Dmol from "3dmol";
import type { TopoAtom, Topology } from "../api/types";
import { MONOMER_COLOR, MONOMER_NAME, PALETTE, RES_TO_MONOMER, download, fmtSigned } from "../lib/format";
import { useTheme } from "../lib/theme";
import { Segmented } from "./ui";

export type ColorMode = "monomer" | "element" | "opls" | "charge" | "residue";
type Repr = "ballstick" | "stick" | "spacefill" | "cg";

const ELEMENT_COLOR: Record<string, string> = { C: "#8a8f86", H: "#d7d7d0", O: "#d9534f", N: "#4f7bd9" };
const MASS: Record<string, number> = { C: 12.011, H: 1.008, O: 15.999, N: 14.007 };

/** Diverging blue–white–red for partial charges, saturating at ±0.7 e. */
export function chargeColor(q: number): string {
  const t = Math.max(-1, Math.min(1, q / 0.7));
  const mix = (a: number, b: number, k: number) => Math.round(a + (b - a) * k);
  const white = [247, 247, 242], red = [200, 60, 50], blue = [50, 100, 190];
  const tgt = t >= 0 ? red : blue, k = Math.abs(t);
  return `rgb(${mix(white[0], tgt[0], k)},${mix(white[1], tgt[1], k)},${mix(white[2], tgt[2], k)})`;
}

function toPdb(topo: Topology): string {
  const L: string[] = [];
  for (const a of topo.atoms) {
    const nm = a.name.length >= 4 ? a.name.slice(0, 4) : ` ${a.name}`.padEnd(4);
    L.push(`HETATM${String(a.index).padStart(5)} ${nm} ${a.resname.padStart(3)} A${String(a.resid).padStart(4)}    ` +
      `${a.x.toFixed(3).padStart(8)}${a.y.toFixed(3).padStart(8)}${a.z.toFixed(3).padStart(8)}  1.00  0.00          ${a.element.padStart(2)}`);
  }
  for (const [i, j] of topo.bonds) L.push(`CONECT${String(i).padStart(5)}${String(j).padStart(5)}`);
  L.push("END");
  return L.join("\n");
}

export function Viewer3D({ topo, height = 540 }: { topo: Topology; height?: number }) {
  const host = useRef<HTMLDivElement>(null);
  const viewer = useRef<$3Dmol.GLViewer | null>(null);
  const first = useRef(true);
  const marker = useRef<$3Dmol.GLShape | null>(null);
  const { theme } = useTheme();
  const [mode, setMode] = useState<ColorMode>("monomer");
  const [repr, setRepr] = useState<Repr>("ballstick");
  const [showH, setShowH] = useState(true);
  const [labels, setLabels] = useState(false);
  const [spin, setSpin] = useState(false);
  const [picked, setPicked] = useState<TopoAtom | null>(null);

  const oplsTypes = useMemo(() => [...new Set(topo.atoms.map(a => a.opls))].sort(), [topo]);
  const colorOf = useMemo(() => {
    const opls = new Map(oplsTypes.map((t, i) => [t, PALETTE[i % PALETTE.length]]));
    return (a: TopoAtom): string => {
      switch (mode) {
        case "element": return ELEMENT_COLOR[a.element] ?? "#999";
        case "monomer": return MONOMER_COLOR[RES_TO_MONOMER[a.resname] ?? "G"];
        case "opls": return opls.get(a.opls) ?? "#999";
        case "charge": return chargeColor(a.charge);
        case "residue": return PALETTE[(a.resid - 1) % PALETTE.length];
      }
    };
  }, [mode, oplsTypes]);

  useEffect(() => {
    const v = $3Dmol.createViewer(host.current!, { backgroundColor: "white" });
    viewer.current = v;
    first.current = true;
    return () => { v.clear(); if (host.current) host.current.innerHTML = ""; viewer.current = null; };
  }, []);

  useEffect(() => {
    const v = viewer.current;
    if (!v) return;
    const view = first.current ? null : v.getView();
    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
    marker.current = null;
    v.setBackgroundColor(theme === "dark" ? "#1b1e15" : "#ffffff", 1);

    if (repr === "cg") {
      const groups = new Map<number, TopoAtom[]>();
      topo.atoms.forEach(a => groups.set(a.resid, [...(groups.get(a.resid) ?? []), a]));
      const com = new Map<number, { x: number; y: number; z: number; res: string }>();
      groups.forEach((atoms, r) => {
        let m = 0, x = 0, y = 0, z = 0;
        atoms.forEach(a => { const w = MASS[a.element] ?? 12; m += w; x += a.x * w; y += a.y * w; z += a.z * w; });
        com.set(r, { x: x / m, y: y / m, z: z / m, res: atoms[0].resname });
      });
      const seen = new Set<string>();
      topo.linkages.forEach(l => {
        const key = [l.residue_1, l.residue_2].sort().join("-");
        if (seen.has(key)) return;
        seen.add(key);
        v.addCylinder({ start: com.get(l.residue_1)!, end: com.get(l.residue_2)!, radius: 0.28, color: "#8b9083", fromCap: 1, toCap: 1 });
      });
      com.forEach((c, r) => {
        v.addSphere({ center: c, radius: 2.3, color: MONOMER_COLOR[RES_TO_MONOMER[c.res] ?? "G"] });
        if (labels) v.addLabel(String(r), { position: c, fontSize: 13, fontColor: "white", backgroundOpacity: 0, showBackground: false });
      });
    } else {
      v.addModel(toPdb(topo), "pdb");
      const cf = (at: { serial: number }) => colorOf(topo.atoms[at.serial - 1]);
      const style: Record<string, unknown> =
        repr === "stick" ? { stick: { radius: 0.15, colorfunc: cf } }
        : repr === "spacefill" ? { sphere: { scale: 0.95, colorfunc: cf } }
        : { stick: { radius: 0.12, colorfunc: cf }, sphere: { scale: 0.26, colorfunc: cf } };
      v.setStyle({}, style);
      if (!showH) v.setStyle({ elem: "H" }, {});
      if (labels) v.addPropertyLabels("atom", { elem: "H", invert: true } as never,
        { fontSize: 10, fontColor: theme === "dark" ? "#e6e8df" : "#22251d", showBackground: false, alignment: "center" } as never);
      v.setClickable({}, true, (at: { serial: number }) => setPicked(topo.atoms[at.serial - 1]));
    }
    if (view) v.setView(view); else { v.zoomTo(); first.current = false; }
    v.spin(spin);
    v.render();
  }, [topo, repr, showH, labels, colorOf, theme]);   // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => { viewer.current?.spin(spin); }, [spin]);

  useEffect(() => {
    const v = viewer.current;
    if (!v) return;
    if (marker.current) v.removeShape(marker.current);
    marker.current = picked && repr !== "cg"
      ? v.addSphere({ center: picked, radius: 0.62, color: "#f2c230", alpha: 0.55 }) : null;
    v.render();
  }, [picked, repr, topo]);

  const legend = useMemo(() => {
    if (mode === "monomer") return (["H", "G", "S"] as const).filter(m => topo.atoms.some(a => RES_TO_MONOMER[a.resname] === m))
      .map(m => [`${m} · ${MONOMER_NAME[m]}`, MONOMER_COLOR[m]] as const);
    if (mode === "element") return Object.entries(ELEMENT_COLOR).filter(([e]) => topo.atoms.some(a => a.element === e)).map(([e, c]) => [e, c] as const);
    if (mode === "opls") return oplsTypes.slice(0, 16).map((t, i) => [t, PALETTE[i % PALETTE.length]] as const);
    return [];
  }, [mode, topo, oplsTypes]);

  return (
    <div className="stack">
      <div className="row">
        <Segmented<ColorMode> label="Colour by" value={mode} onChange={setMode} options={[
          { value: "monomer", label: "Monomer" }, { value: "element", label: "Element" },
          { value: "opls", label: "OPLS type" }, { value: "charge", label: "Charge" }, { value: "residue", label: "Residue" }]} />
        <Segmented<Repr> label="Representation" value={repr} onChange={setRepr} options={[
          { value: "ballstick", label: "Ball & stick" }, { value: "stick", label: "Sticks" },
          { value: "spacefill", label: "Space-filling" }, { value: "cg", label: "CG beads" }]} />
        <label className="check"><input type="checkbox" checked={showH} disabled={repr === "cg"} onChange={e => setShowH(e.target.checked)} /> Hydrogens</label>
        <label className="check"><input type="checkbox" checked={labels} onChange={e => setLabels(e.target.checked)} /> Labels</label>
        <label className="check"><input type="checkbox" checked={spin} onChange={e => setSpin(e.target.checked)} /> Rotate</label>
        <span className="spacer" />
        <button className="btn small" onClick={() => { viewer.current?.zoomTo(); viewer.current?.render(); }}>Reset view</button>
        <button className="btn small" onClick={() => {
          const u = viewer.current?.pngURI(); if (u) { const a = document.createElement("a"); a.href = u; a.download = "structure.png"; a.click(); }
        }}>Save image</button>
      </div>
      <div className="viewer" style={{ height }} data-testid="viewer3d">
        <div ref={host} style={{ position: "absolute", inset: 0 }} />
        {(legend.length > 0 || mode === "charge") && repr !== "cg" && <div className="overlay legend">
          {mode === "charge" ? <>
            <div style={{ height: 8, borderRadius: 3, background: `linear-gradient(90deg, ${chargeColor(-0.7)}, ${chargeColor(0)}, ${chargeColor(0.7)})`, border: "1px solid var(--border)" }} />
            <div className="row small muted" style={{ justifyContent: "space-between", flexWrap: "nowrap" }}><span>−0.7 e</span><span>0</span><span>+0.7 e</span></div></>
            : legend.map(([l, c]) => <div key={l}><span className="swatch" style={{ background: c }} />{l}</div>)}
        </div>}
        {repr === "cg" && <div className="overlay legend">{(["H", "G", "S"] as const).map(m =>
          <div key={m}><span className="swatch" style={{ background: MONOMER_COLOR[m] }} />{m} · {MONOMER_NAME[m]}</div>)}</div>}
        {picked && repr !== "cg" && <div className="overlay info">
          <div className="row"><b>{picked.name}</b><span className="muted">#{picked.index}</span><span className="spacer" />
            <button className="btn small ghost" aria-label="Close" onClick={() => setPicked(null)}>✕</button></div>
          <div className="small">Residue {picked.resid} ({picked.resname}) · {picked.element}</div>
          <div className="small">OPLS: <code>{picked.opls}</code></div>
          <div className="small">q = {fmtSigned(picked.charge)} e <span className="muted">(raw {fmtSigned(picked.raw_charge)})</span></div>
        </div>}
      </div>
      <button className="btn small" style={{ alignSelf: "flex-start" }}
        onClick={() => download("structure.pdb", toPdb(topo), "chemical/x-pdb")}>Download PDB</button>
    </div>
  );
}

import { Alert, Card, Loading, ModelBadge, PageHead } from "../components/ui";
import { Topbar } from "../components/Layout";
import { useSystem } from "../lib/hooks";

const FEATURE_LABEL: Record<string, string> = {
  build_structures: "Build 3-D structures", topology: "OPLS-AA topologies",
  simulations: "Run simulations (GROMACS)", trajectory_analysis: "Trajectory analysis (MDAnalysis)",
  energy_analysis: "Energy analysis",
};

export function SystemPage() {
  const { data: s, isLoading } = useSystem();
  return <>
    <Topbar crumbs={["System"]} />
    <div className="page">
      <PageHead title="System" subtitle="Installed components and the validation status of each scientific model." />
      {isLoading || !s ? <Loading /> : <div className="stack">
        {s.hints.map(h => <Alert key={h} kind="warn">{h}</Alert>)}
        <div className="grid cols-2">
          <Card title="Environment">
            <table><tbody>
              <tr><td>LignoForge</td><td>{s.lignoforge}</td></tr>
              <tr><td>Python</td><td>{s.python}</td></tr>
              <tr><td>GROMACS</td><td>{s.gromacs.version ?? <span className="badge failed">not found</span>} <span className="faint mono">{s.gromacs.path}</span></td></tr>
              <tr><td>RDKit</td><td>{s.rdkit ?? <span className="badge failed">not installed</span>}</td></tr>
              <tr><td>MDAnalysis</td><td>{s.mdanalysis ?? <span className="badge failed">not installed</span>}</td></tr>
              <tr><td>Workspace</td><td className="mono">{s.workspace}</td></tr>
            </tbody></table>
          </Card>
          <Card title="Available features">
            <table><tbody>{Object.entries(s.features).map(([k, v]) =>
              <tr key={k}><td>{FEATURE_LABEL[k] ?? k}</td><td><span className={`badge ${v ? "ok" : "failed"}`}>{v ? "available" : "unavailable"}</span></td></tr>)}</tbody></table>
          </Card>
        </div>
        <Card title="Model validation status" pad={false}>
          <table><thead><tr><th>Component</th><th>Status</th><th>Details</th></tr></thead>
            <tbody>{Object.entries(s.models).map(([k, v]) =>
              <tr key={k}><td>{k.replace(/_/g, " ")}</td><td><ModelBadge status={v.status} /></td>
                <td style={{ whiteSpace: "normal" }}>{v.detail}</td></tr>)}</tbody></table>
        </Card>
        <Alert kind="info">Results that depend on a <b>provisional</b> component should be treated as exploratory
          until that component has been validated against reference data for your system.</Alert>
      </div>}
    </div>
  </>;
}

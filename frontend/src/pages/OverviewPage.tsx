import { Link } from "react-router-dom";
import { Alert, Card, Empty, Loading, PageHead, Stat, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtDate } from "../lib/format";
import { useChains, useProject, useRuns, useSystem } from "../lib/hooks";
import { useProjectName } from "../lib/hooks";

export function OverviewPage() {
  const p = useProjectName();
  const info = useProject(p), chains = useChains(p), runs = useRuns(p), sys = useSystem();
  const base = `/p/${encodeURIComponent(p)}`;
  if (info.isLoading) return <Loading />;
  const gmxMissing = sys.data && !sys.data.features.simulations;

  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, p]} />
    <div className="page">
      <PageHead title={p} subtitle={info.data && `Created ${fmtDate(info.data.created)} · ${info.data.path}`}
        actions={<Link className="btn primary" to={`${base}/chains/new`}>Build a chain</Link>} />
      {gmxMissing && <div style={{ marginBottom: 16 }}><Alert kind="warn" title="GROMACS not found.">
        You can build and inspect structures and topologies, but simulations are disabled. <Link to="/system">How to fix</Link></Alert></div>}
      <div className="grid cols-4" style={{ marginBottom: 20 }}>
        <Stat label="Chains" value={info.data?.n_chains ?? 0} />
        <Stat label="Simulations" value={info.data?.n_runs ?? 0} />
        <Stat label="Analyses" value={info.data?.n_analyses ?? 0} />
        <Stat label="CG parameter sets" value={info.data?.n_cg ?? 0} />
      </div>
      <div className="grid cols-2">
        <Card title="Recent chains" actions={<Link to={`${base}/chains`}>All</Link>} pad={false}>
          {!chains.data?.length ? <Empty title="No chains yet"><Link to={`${base}/chains/new`}>Build your first chain</Link></Empty>
            : <table><tbody>{[...chains.data].reverse().slice(0, 6).map(c => <tr key={c.id}>
              <td><Link to={`${base}/chains/${c.id}`}>{c.name}</Link></td>
              <td className="muted">{c.stats.monomer_count ?? "–"} units</td><td><StatusBadge status={c.status} /></td></tr>)}</tbody></table>}
        </Card>
        <Card title="Recent simulations" actions={<Link to={`${base}/runs`}>All</Link>} pad={false}>
          {!runs.data?.length ? <Empty title="No simulations yet">Open a chain and choose “Run simulation”.</Empty>
            : <table><tbody>{[...runs.data].reverse().slice(0, 6).map(r => <tr key={r.id}>
              <td><Link to={`${base}/runs/${r.id}`}>{r.name}</Link></td>
              <td className="muted">{r.kind}</td><td><StatusBadge status={r.status} /></td></tr>)}</tbody></table>}
        </Card>
      </div>
    </div>
  </>;
}

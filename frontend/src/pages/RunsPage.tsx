import { Link } from "react-router-dom";
import { Card, Empty, Loading, PageHead, ProgressBar, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtDate, fmtDuration } from "../lib/format";
import { useProjectName, useRuns, useSystem } from "../lib/hooks";

export function RunsPage() {
  const p = useProjectName();
  const { data, isLoading } = useRuns(p);
  const sys = useSystem();
  const base = `/p/${encodeURIComponent(p)}`;
  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, "Simulations"]} />
    <div className="page">
      <PageHead title="Simulations" subtitle="Local GROMACS runs: energy minimisation, full atomistic workflows and coarse-grained systems."
        actions={<Link className={`btn primary ${sys.data && !sys.data.features.simulations ? "disabled" : ""}`} to={`${base}/runs/new`}>New simulation</Link>} />
      {isLoading ? <Loading /> : !data?.length
        ? <Card><Empty title="No simulations yet">Open a chain and choose “Run simulation”, or start one here.</Empty></Card>
        : <div className="table-wrap"><table>
          <thead><tr><th>Name</th><th>Type</th><th>Chains</th><th>Status</th><th style={{ width: 160 }}>Progress</th><th>Started</th><th>Duration</th></tr></thead>
          <tbody>{[...data].reverse().map(r => <tr key={r.id}>
            <td><Link to={`${base}/runs/${r.id}`}><b>{r.name}</b></Link> <span className="faint mono">{r.id}</span></td>
            <td>{r.kind === "md" ? "Atomistic MD" : r.kind === "em" ? "Minimisation" : "Coarse-grained"}</td>
            <td className="muted">{r.chain_ids.join(", ")}</td><td><StatusBadge status={r.status} /></td>
            <td>{r.status === "running" ? <ProgressBar fraction={r.progress?.fraction} /> : ""}</td>
            <td className="muted">{fmtDate(r.started)}</td><td className="muted">{fmtDuration(r.started, r.finished)}</td></tr>)}</tbody></table></div>}
    </div>
  </>;
}

import { Link } from "react-router-dom";
import { Card, Empty, Loading, PageHead, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtDate } from "../lib/format";
import { useAnalyses, useCatalog, useProjectName } from "../lib/hooks";

export function AnalysisListPage() {
  const p = useProjectName();
  const { data, isLoading } = useAnalyses(p);
  const cat = useCatalog();
  const label = (id: string) => cat.data?.observables.find(o => o.id === id)?.label ?? id;
  const base = `/p/${encodeURIComponent(p)}`;
  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, "Analysis"]} />
    <div className="page">
      <PageHead title="Analysis" subtitle="Properties and structural observables computed from finished simulations."
        actions={<Link className="btn primary" to={`${base}/analysis/new`}>New analysis</Link>} />
      {isLoading ? <Loading /> : !data?.length
        ? <Card><Empty title="No analyses yet">Finish a simulation, then compute radius of gyration, end-to-end distance, RDFs, contacts and more.</Empty></Card>
        : <div className="table-wrap"><table>
          <thead><tr><th>Observable</th><th>Run</th><th>Status</th><th>Created</th></tr></thead>
          <tbody>{[...data].reverse().map(a => <tr key={a.id}>
            <td><Link to={`${base}/analysis/${a.id}`}><b>{label(a.observable)}</b></Link> <span className="faint mono">{a.id}</span></td>
            <td><Link to={`${base}/runs/${a.run_id}`}>{a.run_id}</Link></td><td><StatusBadge status={a.status} /></td>
            <td className="muted">{fmtDate(a.created)}</td></tr>)}</tbody></table></div>}
    </div>
  </>;
}

import { Link } from "react-router-dom";
import { Card, Empty, Loading, PageHead, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmt, fmtDate, linkLabel } from "../lib/format";
import { useChains, useProjectName } from "../lib/hooks";

export function ChainsPage() {
  const p = useProjectName();
  const { data, isLoading } = useChains(p);
  const base = `/p/${encodeURIComponent(p)}`;
  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, "Chains"]} />
    <div className="page">
      <PageHead title="Chains" subtitle="Lignin chains generated for this project."
        actions={<Link className="btn primary" to={`${base}/chains/new`}>Build a chain</Link>} />
      {isLoading ? <Loading /> : !data?.length
        ? <Card><Empty title="No chains yet">Specify composition and linkages, or estimate them from experimental data.</Empty></Card>
        : <div className="table-wrap"><table>
          <thead><tr><th>Name</th><th>Status</th><th className="num">Monomers</th><th className="num">Atoms</th><th className="num">MW (g/mol)</th>
            <th>Linkages</th><th>Created</th></tr></thead>
          <tbody>{[...data].reverse().map(c => <tr key={c.id}>
            <td><Link to={`${base}/chains/${c.id}`}><b>{c.name}</b></Link> <span className="faint mono">{c.id}</span></td>
            <td><StatusBadge status={c.status} /></td>
            <td className="num">{c.stats.monomer_count ?? "–"}</td><td className="num">{c.stats.n_atoms ?? "–"}</td>
            <td className="num">{fmt(c.stats.MW, 0)}</td>
            <td className="muted">{Object.entries(c.resolved.linkages).filter(([, v]) => v > 0.005).map(([k, v]) => `${linkLabel(k)} ${(v * 100).toFixed(0)}%`).join(" · ")}</td>
            <td className="muted">{fmtDate(c.created)}</td></tr>)}</tbody></table></div>}
    </div>
  </>;
}

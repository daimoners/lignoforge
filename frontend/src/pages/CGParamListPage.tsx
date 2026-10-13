import { Link } from "react-router-dom";
import { Card, Empty, Loading, PageHead, ProgressBar, StatusBadge } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtDate, fmtDuration } from "../lib/format";
import { useCgparams, useProjectName, useSystem } from "../lib/hooks";

export function CGParamListPage() {
  const p = useProjectName();
  const { data, isLoading } = useCgparams(p);
  const sys = useSystem();
  const base = `/p/${encodeURIComponent(p)}`;
  return <>
    <Topbar crumbs={[{ to: "/", label: "Projects" }, { to: base, label: p }, "CG parametrisation"]} />
    <div className="page">
      <PageHead title="CG parametrisation"
        subtitle="Derive a coarse-grained model of bulk lignin: an atomistic OPLS-AA reference is simulated, then bonded terms (Boltzmann inversion) and non-bonded terms (relative entropy) are fitted by iterated CG runs."
        actions={<Link className={`btn primary ${sys.data && !sys.data.features.simulations ? "disabled" : ""}`} to={`${base}/cgparam/new`}>New parametrisation</Link>} />
      {isLoading ? <Loading /> : !data?.length
        ? <Card><Empty title="No parametrisation jobs yet">Start one to derive CG parameters for a lignin bulk ensemble.</Empty></Card>
        : <div className="table-wrap"><table>
          <thead><tr><th>Name</th><th>Scheme</th><th>T (K)</th><th>Status</th><th>Stage</th><th style={{ width: 160 }}>Progress</th><th>Iterations</th><th>Started</th><th>Duration</th></tr></thead>
          <tbody>{[...data].reverse().map(j => <tr key={j.id}>
            <td><Link to={`${base}/cgparam/${j.id}`}><b>{j.name}</b></Link> <span className="faint mono">{j.id}</span></td>
            <td>{j.config.scheme}</td><td className="num">{j.config.temperature}</td>
            <td><StatusBadge status={j.status} /></td><td className="muted">{j.progress?.stage ?? ""}</td>
            <td>{j.status === "running" ? <ProgressBar fraction={j.progress?.fraction} /> : ""}</td>
            <td className="num">{j.iterations}/{j.config.cg.iterations}</td>
            <td className="muted">{fmtDate(j.started)}</td><td className="muted">{fmtDuration(j.started, j.finished)}</td></tr>)}</tbody></table></div>}
    </div>
  </>;
}

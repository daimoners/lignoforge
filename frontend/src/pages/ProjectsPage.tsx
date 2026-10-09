import { useState } from "react";
import { Link, useNavigate } from "react-router-dom";
import { useMutation, useQueryClient } from "@tanstack/react-query";
import { api } from "../api/client";
import { Alert, Card, Empty, Loading, PageHead } from "../components/ui";
import { Topbar } from "../components/Layout";
import { fmtDate } from "../lib/format";
import { useProjects, useSystem } from "../lib/hooks";
import { errMsg, useToast } from "../lib/toast";

export function ProjectsPage() {
  const projects = useProjects();
  const sys = useSystem();
  const qc = useQueryClient();
  const nav = useNavigate();
  const toast = useToast();
  const [name, setName] = useState("");
  const create = useMutation({
    mutationFn: () => api.createProject(name.trim()),
    onSuccess: p => { qc.invalidateQueries({ queryKey: ["projects"] }); nav(`/p/${encodeURIComponent(p.name)}`); },
    onError: e => toast(errMsg(e), "error"),
  });
  const valid = /^[A-Za-z0-9][A-Za-z0-9_.\- ]{0,63}$/.test(name.trim());

  return <>
    <Topbar crumbs={["Projects"]} />
    <div className="page">
      <PageHead title="Projects" subtitle="A project is a folder holding chains, simulations and analyses. Everything in it is reproducible and can be copied or archived as is." />
      {sys.data && sys.data.hints.length > 0 && <div style={{ marginBottom: 16 }}>
        <Alert kind="warn" title="Setup incomplete.">{sys.data.hints.join(" · ")} <Link to="/system">Details</Link></Alert></div>}
      <div className="split">
        <Card title="New project">
          <form className="stack" onSubmit={e => { e.preventDefault(); if (valid) create.mutate(); }}>
            <div className="field"><label htmlFor="pname">Name</label>
              <input id="pname" type="text" value={name} placeholder="e.g. softwood kraft study" onChange={e => setName(e.target.value)} />
              <span className="hint">Letters, digits, spaces, “_”, “.” and “-”. Stored under {sys.data?.workspace ?? "your workspace"}.</span></div>
            <button className="btn primary" disabled={!valid || create.isPending}>Create project</button>
          </form>
        </Card>
        <div>
          {projects.isLoading ? <Loading /> : !projects.data?.length
            ? <Card><Empty title="No projects yet">Create your first project to start building lignin chains.</Empty></Card>
            : <div className="grid cols-2">{projects.data.map(p =>
                <Link key={p.name} to={`/p/${encodeURIComponent(p.name)}`} style={{ color: "inherit", textDecoration: "none" }}>
                  <div className="card card-body stack" style={{ gap: 6 }}>
                    <h2>{p.name}</h2>
                    <div className="muted small">Created {fmtDate(p.created)}</div>
                    <div className="row small"><span>{p.n_chains} chains</span><span className="faint">·</span>
                      <span>{p.n_runs} simulations</span><span className="faint">·</span><span>{p.n_analyses} analyses</span></div>
                  </div></Link>)}</div>}
        </div>
      </div>
    </div>
  </>;
}

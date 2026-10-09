import { Link, NavLink, Outlet, useMatch } from "react-router-dom";
import { useSystem } from "../lib/hooks";
import { useTheme } from "../lib/theme";
import { IconAnalysis, IconChains, IconMoon, IconOverview, IconProjects, IconRuns, IconSun, IconSystem } from "./icons";

const cls = ({ isActive }: { isActive: boolean }) => (isActive ? "active" : "");

export function Layout() {
  const m = useMatch("/p/:project/*");
  const project = m?.params.project;
  const sys = useSystem();
  const { theme, toggle } = useTheme();
  const base = project ? `/p/${encodeURIComponent(project)}` : "";
  const problems = sys.data ? sys.data.hints.length : 0;

  return (
    <div className="app">
      <aside className="sidebar">
        <Link to="/" className="logo-wrap" aria-label="LignoForge home"><img src="/logo.png" alt="LignoForge — digital wood crafting" /></Link>
        <nav className="nav" aria-label="Main">
          <NavLink to="/" end className={cls}><IconProjects /> Projects</NavLink>
          <NavLink to="/system" className={cls}><IconSystem /> System</NavLink>
          {project && <>
            <div className="nav-title">{decodeURIComponent(project)}</div>
            <NavLink to={base} end className={cls}><IconOverview /> Overview</NavLink>
            <NavLink to={`${base}/chains`} className={cls}><IconChains /> Chains</NavLink>
            <NavLink to={`${base}/runs`} className={cls}><IconRuns /> Simulations</NavLink>
            <NavLink to={`${base}/analysis`} className={cls}><IconAnalysis /> Analysis</NavLink>
          </>}
        </nav>
        <div className="sidebar-foot">
          {sys.data && <Link to="/system" className="row" style={{ color: "inherit", gap: 8 }}>
            <span className={`badge ${problems ? "warn" : "ok"}`}><span className="dot" />{problems ? `${problems} setup issue${problems > 1 ? "s" : ""}` : "Environment OK"}</span></Link>}
          <button className="btn small" onClick={toggle} aria-label="Toggle colour theme">
            {theme === "dark" ? <IconSun /> : <IconMoon />}<span>{theme === "dark" ? "Light mode" : "Dark mode"}</span></button>
          {sys.data && <span className="faint">v{sys.data.lignoforge}</span>}
        </div>
      </aside>
      <div className="main"><Outlet /></div>
    </div>
  );
}

export function Topbar({ crumbs }: { crumbs: (string | { to: string; label: string })[] }) {
  return (
    <header className="topbar">
      <div className="crumbs">{crumbs.map((c, i) => <span key={i} className="row" style={{ gap: 8 }}>
        {i > 0 && <span className="faint">/</span>}
        {typeof c === "string" ? (i === crumbs.length - 1 ? <b>{c}</b> : c) : <Link to={c.to}>{c.label}</Link>}</span>)}</div>
    </header>
  );
}

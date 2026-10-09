import { Route, Routes } from "react-router-dom";
import { Layout } from "./components/Layout";
import { ProjectsPage } from "./pages/ProjectsPage";
import { SystemPage } from "./pages/SystemPage";
import { OverviewPage } from "./pages/OverviewPage";
import { ChainsPage } from "./pages/ChainsPage";
import { BuildPage } from "./pages/BuildPage";
import { ChainDetailPage } from "./pages/ChainDetailPage";
import { RunsPage } from "./pages/RunsPage";
import { NewRunPage } from "./pages/NewRunPage";
import { RunDetailPage } from "./pages/RunDetailPage";
import { AnalysisListPage } from "./pages/AnalysisListPage";
import { NewAnalysisPage } from "./pages/NewAnalysisPage";
import { AnalysisDetailPage } from "./pages/AnalysisDetailPage";
import { Empty } from "./components/ui";

export default function App() {
  return (
    <Routes>
      <Route element={<Layout />}>
        <Route index element={<ProjectsPage />} />
        <Route path="system" element={<SystemPage />} />
        <Route path="p/:project">
          <Route index element={<OverviewPage />} />
          <Route path="chains" element={<ChainsPage />} />
          <Route path="chains/new" element={<BuildPage />} />
          <Route path="chains/:cid/*" element={<ChainDetailPage />} />
          <Route path="runs" element={<RunsPage />} />
          <Route path="runs/new" element={<NewRunPage />} />
          <Route path="runs/:rid" element={<RunDetailPage />} />
          <Route path="analysis" element={<AnalysisListPage />} />
          <Route path="analysis/new" element={<NewAnalysisPage />} />
          <Route path="analysis/:aid" element={<AnalysisDetailPage />} />
        </Route>
        <Route path="*" element={<div className="page"><Empty title="Page not found" /></div>} />
      </Route>
    </Routes>
  );
}

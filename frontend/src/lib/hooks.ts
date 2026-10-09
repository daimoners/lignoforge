import { useQuery } from "@tanstack/react-query";
import { api } from "../api/client";

const live = (s?: string) => s === "building" || s === "running";

export const useSystem = () => useQuery({ queryKey: ["system"], queryFn: api.system, staleTime: 30_000 });
export const useCatalog = () => useQuery({ queryKey: ["catalog"], queryFn: api.catalog, staleTime: Infinity });
export const useProjects = () => useQuery({ queryKey: ["projects"], queryFn: api.projects });
export const useProject = (p: string) => useQuery({ queryKey: ["project", p], queryFn: () => api.project(p) });

export const useChains = (p: string) => useQuery({
  queryKey: ["chains", p], queryFn: () => api.chains(p),
  refetchInterval: q => (q.state.data?.some(c => live(c.status)) ? 1500 : false),
});
export const useChain = (p: string, c: string) => useQuery({
  queryKey: ["chain", p, c], queryFn: () => api.chain(p, c),
  refetchInterval: q => (live(q.state.data?.status) ? 1200 : false),
});
export const useTopology = (p: string, c: string, enabled = true) => useQuery({
  queryKey: ["topology", p, c], queryFn: () => api.topology(p, c), enabled, staleTime: Infinity,
});
export const useGraph = (p: string, c: string, enabled = true) => useQuery({
  queryKey: ["graph", p, c], queryFn: () => api.graph(p, c), enabled, staleTime: Infinity,
});
export const useChargeReport = (p: string, c: string, enabled = true) => useQuery({
  queryKey: ["charges", p, c], queryFn: () => api.chargeReport(p, c), enabled, staleTime: Infinity,
});

export const useRuns = (p: string) => useQuery({
  queryKey: ["runs", p], queryFn: () => api.runs(p),
  refetchInterval: q => (q.state.data?.some(r => live(r.status)) ? 2000 : false),
});
export const useRun = (p: string, r: string) => useQuery({
  queryKey: ["run", p, r], queryFn: () => api.run(p, r),
  refetchInterval: q => (live(q.state.data?.status) ? 1500 : false),
});
export const useRunLog = (p: string, r: string, running: boolean) => useQuery({
  queryKey: ["runlog", p, r], queryFn: () => api.runLog(p, r, 400),
  refetchInterval: running ? 2000 : false,
});
export const useRunFiles = (p: string, r: string, running: boolean) => useQuery({
  queryKey: ["runfiles", p, r], queryFn: () => api.runFiles(p, r),
  refetchInterval: running ? 4000 : false,
});
export const useAnalyses = (p: string) => useQuery({
  queryKey: ["analyses", p], queryFn: () => api.analyses(p),
  refetchInterval: q => (q.state.data?.some(a => a.status === "running") ? 1500 : false),
});
export const useAnalysis = (p: string, a: string) => useQuery({
  queryKey: ["analysis", p, a], queryFn: () => api.analysis(p, a),
  refetchInterval: q => (q.state.data?.status === "running" ? 1200 : false),
});
export const useCgSets = (p: string) => useQuery({ queryKey: ["cg", p], queryFn: () => api.cgSets(p) });

import { useParams } from "react-router-dom";
/** The current project name from the URL (only valid inside /p/:project/*). */
export const useProjectName = () => useParams().project!;

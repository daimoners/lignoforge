import type {
  AnalysisRecord, Catalog, ChainGraph, ChainRecord, ChainSpec, ChargeReport, CGSet,
  ProjectInfo, ResolvedSpec, RunFile, RunRecord, SystemInfo, Topology,
} from "./types";

export class ApiError extends Error {
  constructor(public status: number, message: string) { super(message); }
}

const BASE = "/api/v1";

async function req<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(BASE + path, {
    ...init,
    headers: init?.body ? { "Content-Type": "application/json", ...init?.headers } : init?.headers,
  });
  if (!res.ok) {
    let detail = res.statusText;
    try {
      const j = await res.json();
      detail = typeof j.detail === "string" ? j.detail
        : Array.isArray(j.detail) ? j.detail.map((d: { loc?: string[]; msg: string }) =>
            `${(d.loc ?? []).slice(1).join(".")}: ${d.msg}`).join("; ") : detail;
    } catch { /* non-JSON error body */ }
    throw new ApiError(res.status, detail);
  }
  if (res.status === 204) return undefined as T;
  const ct = res.headers.get("content-type") ?? "";
  return (ct.includes("json") ? res.json() : res.text()) as Promise<T>;
}
const post = <T,>(p: string, body?: unknown) =>
  req<T>(p, { method: "POST", body: body === undefined ? undefined : JSON.stringify(body) });
const del = (p: string) => req<void>(p, { method: "DELETE" });
const enc = encodeURIComponent;
const P = (p: string) => `/projects/${enc(p)}`;

export const api = {
  system: () => req<SystemInfo>("/system"),
  catalog: () => req<Catalog>("/catalog"),
  resolveSpec: (spec: ChainSpec) => post<ResolvedSpec>("/resolve-spec", spec),

  projects: () => req<ProjectInfo[]>("/projects"),
  project: (p: string) => req<ProjectInfo>(P(p)),
  createProject: (name: string) => post<ProjectInfo>("/projects", { name }),

  chains: (p: string) => req<ChainRecord[]>(`${P(p)}/chains`),
  chain: (p: string, c: string) => req<ChainRecord>(`${P(p)}/chains/${c}`),
  buildChain: (p: string, spec: ChainSpec) => post<ChainRecord>(`${P(p)}/chains`, spec),
  deleteChain: (p: string, c: string) => del(`${P(p)}/chains/${c}`),
  structure: (p: string, c: string) => req<string>(`${P(p)}/chains/${c}/structure`),
  graph: (p: string, c: string) => req<ChainGraph>(`${P(p)}/chains/${c}/graph`),
  topology: (p: string, c: string) => req<Topology>(`${P(p)}/chains/${c}/topology`),
  chargeReport: (p: string, c: string) => req<ChargeReport>(`${P(p)}/chains/${c}/charge-report`),
  writeTopology: (p: string, c: string) =>
    post<Record<string, string>>(`${P(p)}/chains/${c}/topology/files`),
  structureUrl: (p: string, c: string) => `${BASE}${P(p)}/chains/${c}/structure`,
  topologyFileUrl: (p: string, c: string, name: string) =>
    `${BASE}${P(p)}/chains/${c}/topology/files/${enc(name)}`,

  runs: (p: string) => req<RunRecord[]>(`${P(p)}/runs`),
  run: (p: string, r: string) => req<RunRecord>(`${P(p)}/runs/${r}`),
  startRun: (p: string, body: { chain_ids: string[]; kind: string; params: Record<string, unknown>; name?: string }) =>
    post<RunRecord>(`${P(p)}/runs`, body),
  runLog: (p: string, r: string, tail = 300) => req<string>(`${P(p)}/runs/${r}/log?tail=${tail}`),
  runFiles: (p: string, r: string) => req<RunFile[]>(`${P(p)}/runs/${r}/files`),
  runFileUrl: (p: string, r: string, name: string) => `${BASE}${P(p)}/runs/${r}/files/${enc(name)}`,
  cancelRun: (p: string, r: string) => post<RunRecord>(`${P(p)}/runs/${r}/cancel`),
  deleteRun: (p: string, r: string) => del(`${P(p)}/runs/${r}`),

  analyses: (p: string) => req<AnalysisRecord[]>(`${P(p)}/analyses`),
  analysis: (p: string, a: string) => req<AnalysisRecord>(`${P(p)}/analyses/${a}`),
  startAnalysis: (p: string, body: { run_id: string; observable: string; params: Record<string, unknown> }) =>
    post<AnalysisRecord>(`${P(p)}/analyses`, body),
  deleteAnalysis: (p: string, a: string) => del(`${P(p)}/analyses/${a}`),

  cgSets: (p: string) => req<CGSet[]>(`${P(p)}/cg`),
};

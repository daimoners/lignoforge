// Hand-maintained mirror of the REST payloads (see lignoforge/web/app.py).

export type Monomer = "H" | "G" | "S";
export type RunKind = "em" | "md" | "cg";
export type Status = "building" | "ready" | "running" | "done" | "failed" | "cancelled";

export interface ModelStatus { status: string; detail: string }
export interface SystemInfo {
  lignoforge: string; python: string; platform: string; workspace: string;
  gromacs: { path: string | null; version: string | null };
  rdkit: string | null; mdanalysis: string | null; networkx: string | null;
  features: Record<string, boolean>;
  hints: string[];
  models: Record<string, ModelStatus>;
}
export interface ObservableInfo {
  id: string; label: string; needs: string;
  options: Record<string, { type: string; default: number | string | null; unit?: string }>;
}
export interface Catalog {
  monomers: Monomer[]; linkages: string[]; observables: ObservableInfo[];
  common_options: Record<string, { type: string; default: number | string | null }>;
  run_kinds: RunKind[]; input_schema: Record<string, unknown>;
}

export interface ProjectInfo {
  name: string; created: string; path: string; lignoforge_version: string;
  n_chains: number; n_runs: number; n_analyses: number; n_cg: number;
}

export interface ChainSpec {
  name?: string; n_monomers?: number;
  monomers?: Partial<Record<Monomer, number>>;
  linkages?: Record<string, number>;
  branching?: number; seed?: number; optimize?: boolean; max_iter?: number;
  experimental_input?: Record<string, unknown>;
}
export interface ResolvedSpec {
  n_monomers: number; monomers: Record<Monomer, number>;
  linkages: Record<string, number>; branching: number; seed: number; warnings: string[];
}
export interface ChainRecord {
  id: string; name: string; status: Status; created: string; error: string | null;
  warnings: string[]; spec: ChainSpec; resolved: ResolvedSpec;
  stats: {
    monomer_count?: number; MW?: number; n_atoms?: number; smiles?: string;
    max_degree?: number; [k: string]: unknown;
  };
}
export interface ChainGraph {
  nodes: { id: number; type: Monomer }[];
  edges: { source: number; target: number; linkage: string }[];
}

export interface TopoAtom {
  index: number; name: string; element: string; resid: number; resname: string;
  opls: string; charge: number; raw_charge: number; x: number; y: number; z: number;
}
export interface Residue {
  residue: number; resname: string; linkages: string[]; raw_net_charge: number; max_atom_shift: number;
}
export interface Topology {
  name: string; net_charge: number; n_atoms: number; n_bonds: number; n_angles: number;
  n_dihedrals: number; n_pairs: number; n_impropers: number;
  atoms: TopoAtom[]; bonds: [number, number][];
  linkages: { linkage: string; residue_1: number; atom_1: string; residue_2: number; atom_2: string }[];
  residues: Residue[];
}
export interface ChargeReport {
  text: string; net_charge: number; residues: Residue[];
  by_type: { opls: string; n_atoms: number; mean_delta: number; min_delta: number; max_delta: number }[];
  atoms: { residue: number; resname: string; atom_index: number; atom: string; opls: string;
           raw_charge: number; charge: number; delta: number }[];
}

export interface RunRecord {
  id: string; name: string; kind: RunKind; chain_ids: string[]; status: Status;
  params: Record<string, unknown>; threads: number; created: string; started: string;
  finished?: string; exit_code: number | null; error: string | null;
  layout: { chain_id: string; n_atoms: number; copies: number }[];
  progress: { stage: string; step: number | null; time_ps: number | null; fraction: number | null } | null;
}
export interface RunFile { name: string; size: number }

export interface SeriesResult {
  observable: string; time_ps: number[]; series: Record<string, number[]>;
  summary: Record<string, { mean: number; std: number; min: number; max: number; mean_second_half: number }>;
  extra: Record<string, unknown>;
}
export interface AnalysisRecord {
  id: string; run_id: string; observable: string; params: Record<string, unknown>;
  status: Status; error: string | null; created: string; result: SeriesResult | null;
}

export interface CGSet { id: string; name: string; source: string; created: string; provenance: string }

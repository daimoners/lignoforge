import type { Monomer } from "../api/types";

export const MONOMER_NAME: Record<Monomer, string> = {
  H: "p-Hydroxyphenyl", G: "Guaiacyl", S: "Syringyl",
};
export const RES_TO_MONOMER: Record<string, Monomer> = { HPU: "H", GYU: "G", SYU: "S" };
export const MONOMER_COLOR: Record<Monomer, string> = {
  H: "#d4574a", G: "#4aa564", S: "#3b82c4",
};
export const LINKAGE_LABEL: Record<string, string> = {
  "beta-O-4": "β-O-4", "alpha-O-4": "α-O-4", "4-O-5": "4-O-5", "5-5": "5-5",
  "beta-5": "β-5", "beta-beta": "β-β", "beta-1": "β-1",
};
export const linkLabel = (l: string) => LINKAGE_LABEL[l] ?? l;

export function fmt(x: number | null | undefined, digits = 3): string {
  if (x === null || x === undefined || Number.isNaN(x)) return "–";
  if (x !== 0 && (Math.abs(x) >= 1e5 || Math.abs(x) < 1e-3)) return x.toExponential(2);
  return x.toFixed(digits);
}
export function fmtSigned(x: number, digits = 4) { return (x >= 0 ? "+" : "") + x.toFixed(digits); }
export function fmtDate(s?: string | null): string {
  if (!s) return "–";
  const d = new Date(s);
  return d.toLocaleString(undefined, { dateStyle: "medium", timeStyle: "short" });
}
export function fmtBytes(n: number): string {
  if (n < 1024) return `${n} B`;
  if (n < 1048576) return `${(n / 1024).toFixed(1)} kB`;
  return `${(n / 1048576).toFixed(1)} MB`;
}
export function fmtDuration(startIso?: string, endIso?: string): string {
  if (!startIso) return "–";
  const s = (new Date(endIso ?? Date.now()).getTime() - new Date(startIso).getTime()) / 1000;
  if (s < 60) return `${Math.max(0, Math.round(s))} s`;
  if (s < 3600) return `${Math.floor(s / 60)} min ${Math.round(s % 60)} s`;
  return `${Math.floor(s / 3600)} h ${Math.floor((s % 3600) / 60)} min`;
}
export function download(filename: string, text: string, mime = "text/plain") {
  const url = URL.createObjectURL(new Blob([text], { type: mime }));
  const a = document.createElement("a");
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}
export const nice = (x: number) => Number(x.toPrecision(12));

/** Categorical palette (colour-blind-friendly, muted). */
export const PALETTE = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948",
  "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#86bcb6", "#d37295", "#8cd17d", "#a0cbe8",
  "#ffbe7d", "#b6992d"];

import { useMemo, useState, type ReactNode } from "react";

export interface Col<T> {
  key: string; header: string; num?: boolean;
  render: (row: T) => ReactNode;
  sort?: (row: T) => number | string;
}

export function DataTable<T>({ rows, cols, rowKey, onRowClick, maxHeight = 520 }: {
  rows: T[]; cols: Col<T>[]; rowKey: (r: T) => string | number;
  onRowClick?: (r: T) => void; maxHeight?: number;
}) {
  const [sortKey, setSortKey] = useState<string | null>(null);
  const [dir, setDir] = useState<1 | -1>(1);
  const sorted = useMemo(() => {
    const c = cols.find(x => x.key === sortKey);
    if (!c?.sort) return rows;
    const f = c.sort;
    return [...rows].sort((a, b) => (f(a) < f(b) ? -dir : f(a) > f(b) ? dir : 0));
  }, [rows, cols, sortKey, dir]);
  return (
    <div className="table-wrap" style={{ maxHeight }}>
      <table>
        <thead><tr>{cols.map(c => <th key={c.key} className={`${c.num ? "num" : ""} ${c.sort ? "sortable" : ""}`}
          aria-sort={sortKey === c.key ? (dir === 1 ? "ascending" : "descending") : undefined}
          onClick={() => { if (!c.sort) return; if (sortKey === c.key) setDir(d => (d === 1 ? -1 : 1)); else { setSortKey(c.key); setDir(1); } }}>
          {c.header}{sortKey === c.key ? (dir === 1 ? " ▲" : " ▼") : ""}</th>)}</tr></thead>
        <tbody>{sorted.map(r => <tr key={rowKey(r)} className={onRowClick ? "click" : undefined} onClick={() => onRowClick?.(r)}>
          {cols.map(c => <td key={c.key} className={c.num ? "num" : ""}>{c.render(r)}</td>)}</tr>)}</tbody>
      </table>
    </div>
  );
}

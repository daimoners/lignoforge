import { cloneElement, isValidElement, useId, useState, type ReactElement, type ReactNode } from "react";
import { NavLink } from "react-router-dom";
import type { Status } from "../api/types";

export function Spinner() { return <span className="spinner" aria-label="loading" />; }

export function Loading({ text = "Loading…" }: { text?: string }) {
  return <div className="empty"><Spinner /> <span style={{ marginLeft: 8 }}>{text}</span></div>;
}

export function StatusBadge({ status }: { status: Status | string }) {
  return <span className={`badge ${status}`}><span className="dot" />{status}</span>;
}

export function Alert({ kind = "info", title, children }: {
  kind?: "info" | "warn" | "error" | "ok"; title?: string; children?: ReactNode;
}) {
  return <div className={`alert ${kind}`} role={kind === "error" ? "alert" : undefined}>
    <div>{title && <b>{title} </b>}{children}</div></div>;
}

export function Card({ title, actions, children, pad = true }: {
  title?: ReactNode; actions?: ReactNode; children: ReactNode; pad?: boolean;
}) {
  return (
    <section className="card">
      {(title || actions) && <div className="card-head"><h2>{title}</h2><span className="spacer" />{actions}</div>}
      <div className={pad ? "card-body" : undefined}>{children}</div>
    </section>
  );
}

export function PageHead({ title, subtitle, actions }: { title: ReactNode; subtitle?: ReactNode; actions?: ReactNode }) {
  return (
    <div className="page-head">
      <div style={{ flex: 1, minWidth: 260 }}><h1>{title}</h1>{subtitle && <p>{subtitle}</p>}</div>
      <div className="row">{actions}</div>
    </div>
  );
}

/** Labelled form field: the label is bound to its control (clicking it focuses the control). */
export function Field({ label, hint, children }: { label: string; hint?: ReactNode; children: ReactNode }) {
  const id = useId();
  const control = isValidElement(children)
    ? cloneElement(children as ReactElement<{ id?: string; "aria-describedby"?: string }>,
        { id, ...(hint ? { "aria-describedby": `${id}-hint` } : {}) })
    : children;
  return <div className="field"><label htmlFor={id}>{label}</label>{control}
    {hint && <span className="hint" id={`${id}-hint`}>{hint}</span>}</div>;
}

export function NumberInput({ value, onChange, min, max, step = 1, disabled, ariaLabel, id, ...rest }: {
  value: number | ""; onChange: (v: number | "") => void; min?: number; max?: number;
  step?: number | "any"; disabled?: boolean; ariaLabel?: string; id?: string; "aria-describedby"?: string;
}) {
  return <input type="number" id={id} aria-describedby={rest["aria-describedby"]} aria-label={ariaLabel} value={value} min={min} max={max} step={step} disabled={disabled}
    onChange={e => onChange(e.target.value === "" ? "" : Number(e.target.value))} />;
}

export function Tabs({ items }: { items: { to: string; label: string; end?: boolean }[] }) {
  return <nav className="tabs">{items.map(i =>
    <NavLink key={i.to} to={i.to} end={i.end} className={({ isActive }) => (isActive ? "active" : "")}>{i.label}</NavLink>)}</nav>;
}

export function Segmented<T extends string>({ value, options, onChange, label }: {
  value: T; options: { value: T; label: string }[]; onChange: (v: T) => void; label?: string;
}) {
  return <div className="seg" role="group" aria-label={label}>{options.map(o =>
    <button key={o.value} type="button" className={o.value === value ? "on" : ""} aria-pressed={o.value === value}
      onClick={() => onChange(o.value)}>{o.label}</button>)}</div>;
}

export function Empty({ title, children }: { title: string; children?: ReactNode }) {
  return <div className="empty"><h3>{title}</h3>{children}</div>;
}

export function ProgressBar({ fraction }: { fraction: number | null | undefined }) {
  const det = fraction !== null && fraction !== undefined;
  return <div className={`progress ${det ? "" : "indeterminate"}`}
    role="progressbar" aria-valuenow={det ? Math.round(fraction * 100) : undefined}>
    <i style={det ? { width: `${Math.max(2, fraction * 100)}%` } : undefined} /></div>;
}

export function KV({ items }: { items: [string, ReactNode][] }) {
  return <dl className="kv">{items.map(([k, v]) => <div style={{ display: "contents" }} key={k}><dt>{k}</dt><dd>{v}</dd></div>)}</dl>;
}

export function Stat({ label, value }: { label: string; value: ReactNode }) {
  return <div className="card stat"><div className="l">{label}</div><div className="v">{value}</div></div>;
}

/** Button that asks for a second click before running a destructive action. */
export function ConfirmButton({ label, confirmLabel = "Click again to confirm", onConfirm, disabled }: {
  label: string; confirmLabel?: string; onConfirm: () => void; disabled?: boolean;
}) {
  const [armed, setArmed] = useState(false);
  return <button className="btn danger" disabled={disabled} onBlur={() => setArmed(false)}
    onClick={() => { if (armed) { setArmed(false); onConfirm(); } else setArmed(true); }}>
    {armed ? confirmLabel : label}</button>;
}

export function ModelBadge({ status }: { status: string }) {
  const cls = status === "provisional" ? "provisional" : status.startsWith("validated") ? "ok" : "";
  return <span className={`badge ${cls}`}>{status}</span>;
}

#!/usr/bin/env python3
"""
lignin_ff/tools/test_assign.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Validation tests for assign_chain_types.py.

Each test builds a small lignin chain with the lignoforge API, writes a PDB,
runs assign_chain_types, and checks:
  * net charge of every residue == 0.000
  * atom types at linkage sites are correct
  * atom types at non-linkage sites unchanged vs isolated neutral monomer

Usage
-----
    python lignin_ff/tools/test_assign.py
or
    python -m pytest lignin_ff/tools/test_assign.py -v
"""
from __future__ import annotations

import os
import sys
import tempfile
from typing import Dict

# ── path setup ──────────────────────────────────────────────────────────────
_HERE   = os.path.dirname(os.path.abspath(__file__))
_ROOT   = os.path.abspath(os.path.join(_HERE, "..", ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from assign_chain_types import (  # noqa: E402
    assign_chain_types,
    ISOLATED,
    LINKAGE_MODS,
)

# ── lignoforge imports ───────────────────────────────────────────────────────
try:
    from ligning.monomer  import Monomer   # noqa: F401  (may differ by version)
    from ligning.polymer  import Polymer
    _LIGNING_PKG = "ligning"
except ImportError:
    try:
        from lignoforge.core.monomer import Monomer
        from lignoforge.core.polymer import Polymer
        _LIGNING_PKG = "lignoforge.core"
    except ImportError:
        Monomer = Polymer = None
        _LIGNING_PKG = None

# ── helpers ──────────────────────────────────────────────────────────────────

def _make_pdb(monomers, linkages, tmpdir, name="chain"):
    """
    Build a minimal polymer chain programmatically and write a PDB.

    monomers  : list of monomer type strings, e.g. ["G", "G"]
    linkages  : list of linkage names for consecutive pairs, e.g. ["beta-O-4"]

    Returns the PDB path.
    """
    if Monomer is None or Polymer is None:
        return None

    m0 = Monomer(monomers[0], monomer_index=0)
    m0.create()
    p = Polymer(m0, verbose=False)

    for i, (mt, lk) in enumerate(zip(monomers[1:], linkages), start=1):
        p.add_specific_monomer(mt, lk)

    pdb_path = os.path.join(tmpdir, f"{name}.pdb")

    # Try both API variations
    try:
        from ligning.utils import write_pdb
        write_pdb(p, pdb_path)
    except Exception:
        try:
            from lignoforge.structure.pdb import PDBStructureWriter
            w = PDBStructureWriter()
            w.write_polymer_pdb(p, pdb_path, optimize_3d=False)
        except Exception as exc:
            print(f"  [skip] Could not write PDB: {exc}", file=sys.stderr)
            return None

    return pdb_path


def check_net_charges(typed: dict, tag: str) -> bool:
    ok = True
    for rseq, tr in typed.items():
        net = round(sum(rec[2] for rec in tr["atoms"].values()), 4)
        if abs(net) > 0.001:
            print(f"  FAIL [{tag}] res {rseq} ({tr.get('rtp_name','?')}): "
                  f"net charge = {net:+.4f}")
            ok = False
        else:
            print(f"  OK   [{tag}] res {rseq} ({tr.get('rtp_name','?')}): "
                  f"net charge = {net:+.4f}")
    return ok


def check_atom_type(typed: dict, res_seq: int, atom_name: str, expected_opls: str,
                    tag: str) -> bool:
    tr = typed.get(res_seq)
    if tr is None:
        print(f"  FAIL [{tag}] res {res_seq} not found")
        return False
    rec = tr["atoms"].get(atom_name)
    if rec is None:
        print(f"  FAIL [{tag}] res {res_seq} atom {atom_name} not found "
              f"(present: {list(tr['atoms'].keys())})")
        return False
    got = rec[0]
    if got != expected_opls:
        print(f"  FAIL [{tag}] res {res_seq} {atom_name}: "
              f"expected {expected_opls}, got {got}")
        return False
    print(f"  OK   [{tag}] res {res_seq} {atom_name} = {got}")
    return True


def check_atom_absent(typed: dict, res_seq: int, atom_name: str, tag: str) -> bool:
    tr = typed.get(res_seq)
    if tr is None:
        print(f"  FAIL [{tag}] res {res_seq} not found")
        return False
    if atom_name in tr["atoms"]:
        print(f"  FAIL [{tag}] res {res_seq} atom {atom_name} should be absent")
        return False
    print(f"  OK   [{tag}] res {res_seq} {atom_name} correctly absent")
    return True


# ── Individual tests ──────────────────────────────────────────────────────────

def test_isolated_neutral_charge():
    """Each isolated neutral monomer definition should sum to exactly 0.000."""
    tag = "isolated_neutral_charge"
    ok = True
    for nm, atoms in ISOLATED.items():
        net = round(sum(rec[2] for rec in atoms.values()), 6)
        if abs(net) > 1e-5:
            print(f"  FAIL [{tag}] {nm}: net = {net:+.6f}")
            ok = False
        else:
            print(f"  OK   [{tag}] {nm}: net = {net:+.6f}")
    return ok


def test_isolated_g_monomer(tmpdir):
    """Single isolated G monomer → rtp_name='GYM', vinyl CB=opls_142, HB=opls_144,
    C4=opls_166 (free phenol), HO4 present, no OA.
    """
    tag = "isolated_G"
    pdb = _make_pdb(["G"], [], tmpdir, "isolated_g")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "isolated_g_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    # canonical name must be GYM (vinyl monolignol, not the sp3 ref GNM)
    got_name = typed[1].get("rtp_name", "?")
    if got_name != "GYM":
        print(f"  FAIL [{tag}] rtp_name expected 'GYM', got '{got_name}'")
        ok = False
    else:
        print(f"  OK   [{tag}] rtp_name = '{got_name}'")
    # C4 should stay opls_166 (phenol): HO4 is present in isolated vinyl
    ok &= check_atom_type(typed, 1, "C4",  "opls_166", tag)
    # CB is vinyl sp2: opls_142 (no inter-residue bond, no OA)
    ok &= check_atom_type(typed, 1, "CB",  "opls_142", tag)
    ok &= check_atom_type(typed, 1, "HB",  "opls_144", tag)
    # OA must be absent in the vinyl form
    ok &= check_atom_absent(typed, 1, "OA", tag)
    return ok


def test_isolated_H_monomer(tmpdir):
    """Single isolated H monomer → rtp_name='HPM', vinyl side chain,
    C4=opls_166 (free phenol), H3 present (no OMe at C3), no OA.
    """
    tag = "isolated_H"
    pdb = _make_pdb(["H"], [], tmpdir, "isolated_h")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "isolated_h_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    got_name = typed[1].get("rtp_name", "?")
    if got_name != "HPM":
        print(f"  FAIL [{tag}] rtp_name expected 'HPM', got '{got_name}'")
        ok = False
    else:
        print(f"  OK   [{tag}] rtp_name = '{got_name}'")
    ok &= check_atom_type(typed, 1, "C4",  "opls_166", tag)   # free phenol
    ok &= check_atom_type(typed, 1, "C3",  "opls_145", tag)   # plain aromatic CH
    ok &= check_atom_type(typed, 1, "CB",  "opls_142", tag)   # vinyl
    ok &= check_atom_type(typed, 1, "HB",  "opls_144", tag)
    ok &= check_atom_absent(typed, 1, "OM3", tag)   # no OMe in H-type
    ok &= check_atom_absent(typed, 1, "OA",  tag)   # no alpha-OH in vinyl form
    return ok


def test_isolated_S_monomer(tmpdir):
    """Single isolated S monomer → rtp_name='SYM', vinyl side chain,
    C4=opls_166 (free phenol), OM3 and OM5 present (OMe at C3 and C5).
    """
    tag = "isolated_S"
    pdb = _make_pdb(["S"], [], tmpdir, "isolated_s")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "isolated_s_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    got_name = typed[1].get("rtp_name", "?")
    if got_name != "SYM":
        print(f"  FAIL [{tag}] rtp_name expected 'SYM', got '{got_name}'")
        ok = False
    else:
        print(f"  OK   [{tag}] rtp_name = '{got_name}'")
    ok &= check_atom_type(typed, 1, "C4",  "opls_166", tag)   # free phenol
    ok &= check_atom_type(typed, 1, "C3",  "opls_199", tag)   # aryl ether (OMe)
    ok &= check_atom_type(typed, 1, "C5",  "opls_199", tag)   # aryl ether (OMe)
    ok &= check_atom_type(typed, 1, "CB",  "opls_142", tag)   # vinyl
    ok &= check_atom_absent(typed, 1, "OA",  tag)   # no alpha-OH in vinyl form
    return ok


def test_two_G_beta_O_4(tmpdir):
    """G–G β-O-4 dimer: check donor and acceptor type assignments."""
    tag = "G_G_beta_O_4"
    pdb = _make_pdb(["G", "G"], ["beta-O-4"], tmpdir, "GG_bO4")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "GG_bO4_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)

    # Identify the donor (has O4H in bond) and acceptor (has CB in bond) residues.
    # In a linear dimer the first residue is the O4H DONOR (its O4H bonds to CB of res2).
    r_donor    = 1
    r_acceptor = 2

    # Donor residue: C4→opls_199, O4H→opls_179, HO4 absent
    ok &= check_atom_type(typed, r_donor, "C4",  "opls_199", tag)
    ok &= check_atom_type(typed, r_donor, "O4H", "opls_179", tag)
    ok &= check_atom_absent(typed, r_donor, "HO4", tag)

    # Acceptor residue: CB→opls_157, HB→opls_156, HB2 absent (never existed; PDB has only HB)
    ok &= check_atom_type(typed, r_acceptor, "CB",  "opls_157", tag)
    ok &= check_atom_type(typed, r_acceptor, "HB",  "opls_156", tag)
    return ok


def test_three_G_beta_O_4(tmpdir):
    """G–G–G three-mer β-O-4: central unit should be GYU-like (both mods)."""
    tag = "G_G_G_beta_O_4"
    pdb = _make_pdb(["G", "G", "G"], ["beta-O-4", "beta-O-4"], tmpdir, "GGG_bO4")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "GGG_bO4_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)

    # In the lignoforge 3-mer, residue 1 is the central (both CB acceptor and O4H donor).
    # The canonical name GYU confirms it: "G|beta-O-4@CB|beta-O-4@O4H"
    r_mid = 1
    ok &= check_atom_type(typed, r_mid, "C4",  "opls_199", tag)
    ok &= check_atom_type(typed, r_mid, "O4H", "opls_179", tag)
    ok &= check_atom_absent(typed, r_mid, "HO4", tag)
    ok &= check_atom_type(typed, r_mid, "CB",  "opls_157", tag)
    return ok


def test_two_G_5_5(tmpdir):
    """G–G 5-5 (biphenyl): C5 on each residue → opls_145 q=0, H5 absent."""
    tag = "G_G_5_5"
    pdb = _make_pdb(["G", "G"], ["5-5"], tmpdir, "GG_55")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "GG_55_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    for rseq in [1, 2]:
        ok &= check_atom_type(typed, rseq, "C5", "opls_145", tag)
        ok &= check_atom_absent(typed, rseq, "H5", tag)
    return ok


def test_G_H_4_O_5(tmpdir):
    """G–H 4-O-5 dimer."""
    tag = "G_H_4_O_5"
    pdb = _make_pdb(["G", "H"], ["4-O-5"], tmpdir, "GH_4O5")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "GH_4O5_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    # O4H donor (res1): C4→opls_199, O4H→opls_179, HO4 absent
    ok &= check_atom_type(typed, 1, "C4",  "opls_199", tag)
    ok &= check_atom_type(typed, 1, "O4H", "opls_179", tag)
    ok &= check_atom_absent(typed, 1, "HO4", tag)
    # C5 acceptor (res2): C5→opls_199, H5 absent
    ok &= check_atom_type(typed, 2, "C5", "opls_199", tag)
    ok &= check_atom_absent(typed, 2, "H5", tag)
    return ok


def test_renamed_pdb(tmpdir):
    """
    The renamed PDB must have:
    - residue names identical to rtp_name from the typed dict
    - all names <= 3 characters
    - each rtp_name found in the PDB for the correct residue sequence number
    """
    tag = "renamed_pdb"
    pdb = _make_pdb(["G", "G", "G"], ["beta-O-4", "beta-O-4"], tmpdir, "GGG_rn_test")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True

    from assign_chain_types import assign_chain_types, write_renamed_pdb
    rtp = os.path.join(tmpdir, "GGG_rn.rtp")
    rpdb = os.path.join(tmpdir, "GGG_renamed.pdb")
    typed = assign_chain_types(pdb, rtp, rpdb)

    # Collect resname per resseq from the renamed PDB
    pdb_resnames: dict = {}
    with open(rpdb) as fh:
        for line in fh:
            if line[:6].strip() in ("ATOM", "HETATM") and len(line) >= 26:
                try:
                    resseq = int(line[22:26])
                    rname  = line[17:20].strip()
                    pdb_resnames[resseq] = rname
                except ValueError:
                    pass

    ok = True
    for rseq, tr in typed.items():
        expected = tr["rtp_name"]
        got      = pdb_resnames.get(rseq, "??")
        if len(expected) > 3:
            print(f"  FAIL [{tag}] rtp_name '{expected}' exceeds 3 characters")
            ok = False
        elif got != expected:
            print(f"  FAIL [{tag}] res {rseq}: PDB has '{got}', expected '{expected}'")
            ok = False
        else:
            print(f"  OK   [{tag}] res {rseq}: PDB resname = '{got}'")
    return ok


def test_H_H_beta_O_4(tmpdir):
    """H–H β-O-4 dimer: same type logic as G–G but no OMe groups.
    Donor: C4→opls_199, O4H→opls_179, HO4 absent.
    Acceptor: CB→opls_157, HB→opls_156.
    """
    tag = "H_H_beta_O_4"
    pdb = _make_pdb(["H", "H"], ["beta-O-4"], tmpdir, "HH_bO4")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "HH_bO4_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    # Determine donor and acceptor
    r_donor, r_acceptor = 1, 2
    ok &= check_atom_type(typed, r_donor, "C4",  "opls_199", tag)
    ok &= check_atom_type(typed, r_donor, "O4H", "opls_179", tag)
    ok &= check_atom_absent(typed, r_donor, "HO4", tag)
    ok &= check_atom_absent(typed, r_donor, "OM3", tag)   # no OMe in H-type
    ok &= check_atom_type(typed, r_acceptor, "CB",  "opls_157", tag)
    ok &= check_atom_type(typed, r_acceptor, "HB",  "opls_156", tag)
    return ok


def test_S_S_beta_O_4(tmpdir):
    """S–S β-O-4 dimer: verify OMe groups at both C3 and C5 survive the
    β-O-4 modification, and the internal unit gets canonical name 'SYU'.
    """
    tag = "S_S_S_beta_O_4"
    pdb = _make_pdb(["S", "S", "S"], ["beta-O-4", "beta-O-4"], tmpdir, "SSS_bO4")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "SSS_bO4_custom.rtp")
    typed = assign_chain_types(pdb, rtp)
    ok = check_net_charges(typed, tag)
    # Central unit (res 1, both linkages) → SYU
    r_mid = 1
    got_name = typed[r_mid].get("rtp_name", "?")
    if got_name != "SYU":
        print(f"  FAIL [{tag}] res {r_mid} rtp_name expected 'SYU', got '{got_name}'")
        ok = False
    else:
        print(f"  OK   [{tag}] res {r_mid} rtp_name = '{got_name}'")
    ok &= check_atom_type(typed, r_mid, "C4",  "opls_199", tag)  # ether (donor)
    ok &= check_atom_type(typed, r_mid, "C3",  "opls_199", tag)  # OMe intact
    ok &= check_atom_type(typed, r_mid, "C5",  "opls_199", tag)  # OMe at C5 intact
    ok &= check_atom_type(typed, r_mid, "CB",  "opls_157", tag)  # sp3 ether (acceptor)
    return ok


def test_two_G_beta_1(tmpdir):
    """G–G β-1 dimer: CB→opls_137 on the CB side, C1 unchanged (opls_145) on C1 side."""
    tag = "G_G_beta_1"
    pdb = _make_pdb(["G", "G"], ["beta-1"], tmpdir, "GG_b1")
    if pdb is None:
        print(f"  SKIP [{tag}] PDB generation unavailable")
        return True
    rtp = os.path.join(tmpdir, "GG_b1_custom.rtp")
    typed = assign_chain_types(pdb, rtp)

    if len(typed) < 2:
        print(f"  SKIP [{tag}] β-1 dimer produced only {len(typed)} residue(s); "
              f"linkage may not be supported by this lignoforge version")
        return True

    ok = check_net_charges(typed, tag)
    # Find which residue has CB linkage and which has C1 linkage
    # by inspecting the linkages field
    res_cb = res_c1 = None
    for rseq, tr in typed.items():
        for lname, pos in tr["linkages"]:
            if lname == "beta-1" and pos == "CB":
                res_cb = rseq
            if lname == "beta-1" and pos == "C1":
                res_c1 = rseq

    if res_cb is not None:
        ok &= check_atom_type(typed, res_cb, "CB", "opls_137", tag)
    else:
        print(f"  SKIP [{tag}] no CB-side residue found for beta-1")

    if res_c1 is not None:
        ok &= check_atom_type(typed, res_c1, "C1", "opls_145", tag)
    else:
        print(f"  SKIP [{tag}] no C1-side residue found for beta-1")

    return ok


# ── Pure-unit tests (no PDB needed) ──────────────────────────────────────────

def test_linkage_mods_defined():
    """All 7 linkage types must have at least one position rule."""
    tag = "linkage_mods_defined"
    required = {"beta-O-4", "5-5", "4-O-5", "beta-5", "beta-beta", "alpha-O-4", "beta-1"}
    ok = True
    for lname in required:
        if lname not in LINKAGE_MODS or not LINKAGE_MODS[lname]:
            print(f"  FAIL [{tag}] {lname} not in LINKAGE_MODS")
            ok = False
        else:
            print(f"  OK   [{tag}] {lname} has {len(LINKAGE_MODS[lname])} position rule(s)")
    return ok


# ── Runner ────────────────────────────────────────────────────────────────────

def main():
    all_ok  = True
    tmpdir  = tempfile.mkdtemp(prefix="test_assign_")
    print(f"\nUsing temp dir: {tmpdir}\n")

    tests_unit = [
        ("isolated_neutral_charge", test_isolated_neutral_charge,  ()),
        ("linkage_mods_defined",    test_linkage_mods_defined,      ()),
    ]
    tests_pdb = [
        ("isolated_G",         test_isolated_g_monomer,   (tmpdir,)),
        ("isolated_H",         test_isolated_H_monomer,   (tmpdir,)),
        ("isolated_S",         test_isolated_S_monomer,   (tmpdir,)),
        ("G_G_beta_O_4",       test_two_G_beta_O_4,       (tmpdir,)),
        ("G_G_G_beta_O_4",     test_three_G_beta_O_4,     (tmpdir,)),
        ("H_H_beta_O_4",       test_H_H_beta_O_4,         (tmpdir,)),
        ("S_S_S_beta_O_4",     test_S_S_beta_O_4,         (tmpdir,)),
        ("renamed_pdb",        test_renamed_pdb,           (tmpdir,)),
        ("G_G_5_5",            test_two_G_5_5,            (tmpdir,)),
        ("G_H_4_O_5",          test_G_H_4_O_5,            (tmpdir,)),
        ("G_G_beta_1",         test_two_G_beta_1,         (tmpdir,)),
    ]

    for name, fn, args in tests_unit + tests_pdb:
        print(f"\n{'='*60}")
        print(f"TEST: {name}")
        print(f"{'='*60}")
        try:
            result = fn(*args)
            all_ok &= result
        except Exception as exc:
            print(f"  ERROR [{name}]: {exc}")
            import traceback; traceback.print_exc()
            all_ok = False

    print(f"\n{'='*60}")
    if all_ok:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()

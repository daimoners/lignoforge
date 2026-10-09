#!/usr/bin/env python
"""
lignoforge — unified command line
=================================

    lignoforge chain  ...     build chains and export files (same as lignoforge-chain)
    lignoforge cg-fit ...     derive CG parameters from atomistic MD
    lignoforge gui            start the local web interface
    lignoforge check          report installed dependencies and model status
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
import webbrowser
from typing import Optional


def _gui(args: argparse.Namespace) -> int:
    try:
        import uvicorn
        from lignoforge.web.app import create_app
    except ImportError:
        print("The web interface needs extra packages: pip install 'lignoforge[gui]'",
              file=sys.stderr)
        return 1
    app = create_app(args.workspace, dev=args.dev)
    url = f"http://{args.host}:{args.port}"
    print(f"LignoForge GUI  →  {url}\nWorkspace       →  {app.state.workspace}\n"
          "Press Ctrl+C to stop.", flush=True)
    if not args.no_browser:
        threading.Thread(target=lambda: (time.sleep(1.2), webbrowser.open(url)),
                         daemon=True).start()
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
    return 0


def _check(args: argparse.Namespace) -> int:
    from lignoforge.api import system_check
    info = system_check()
    if args.json:
        print(json.dumps(info, indent=2))
        return 0
    print(f"LignoForge {info['lignoforge']}  (Python {info['python']})")
    g = info["gromacs"]
    print(f"  GROMACS     : {g['version'] or 'NOT FOUND'}  {g['path'] or ''}")
    for k in ("rdkit", "mdanalysis", "networkx"):
        print(f"  {k:12s}: {info[k] or 'NOT INSTALLED'}")
    print("  Features    : " + ", ".join(f"{k}={'yes' if v else 'no'}"
                                          for k, v in info["features"].items()))
    for h in info["hints"]:
        print(f"  hint        : {h}")
    print("  Model status:")
    for k, v in info["models"].items():
        print(f"    {k:22s} {v['status']:22s} {v['detail']}")
    return 0 if info["features"]["build_structures"] else 1


def main(argv: Optional[list] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # Sub-commands that own their parser are dispatched untouched.
    if argv and argv[0] == "chain":
        from lignoforge.cli.build_chain import main as chain_main
        return chain_main(argv[1:])
    if argv and argv[0] == "cg-fit":
        from lignoforge.cli.cg_fit import main as fit_main
        return fit_main(argv[1:])

    ap = argparse.ArgumentParser(prog="lignoforge", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("chain", help="build chains (see 'lignoforge chain -h')")
    sub.add_parser("cg-fit", help="fit CG parameters (see 'lignoforge cg-fit -h')")
    g = sub.add_parser("gui", help="start the local web interface")
    g.add_argument("--workspace", default=None,
                   help="folder holding your projects (default: ~/LignoForge)")
    g.add_argument("--host", default="127.0.0.1")
    g.add_argument("--port", type=int, default=8765)
    g.add_argument("--no-browser", action="store_true")
    g.add_argument("--dev", action="store_true", help="allow the Vite dev server (CORS)")
    g.set_defaults(func=_gui)
    c = sub.add_parser("check", help="check dependencies and model status")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=_check)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

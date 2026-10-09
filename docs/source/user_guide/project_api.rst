.. _project_api:

=========================================
Project API and Web Backend
=========================================

Everything LignoForge does is available through one UI-independent class,
:class:`lignoforge.api.Project`.  The command line, notebooks and the web
interface all use it, so they behave identically.

.. code-block:: python

   from lignoforge.api import Project

   p = Project.create("my_lignin", "softwood study")
   chain = p.build_chain({"n_monomers": 10, "monomers": {"G": 0.9, "H": 0.1},
                          "linkages": {"beta-O-4": 5, "5-5": 1, "beta-5": 1},
                          "seed": 1})                     # reproducible
   p.topology(chain["id"])               # OPLS types, charges, raw charges
   p.charge_report(chain["id"])          # which charges were renormalised

   run = p.start_run(chain["id"], "md", {"solvent": "tip3p", "prod_ns": 10})
   p.run(run["id"])                      # status, progress; survives restarts
   p.analyze(run["id"], "rg")            # also end_to_end, rmsd, density,
                                         # rdf, contacts, energy
   p.fit_cg_parameters([run["id"]])      # Boltzmann inversion → CG parameter set

A project is a plain folder (``project.json``, ``chains/``, ``runs/``,
``analyses/``, ``cg/``).  Every record stores its specification, the resolved
parameters, the seed and the LignoForge version, so any result can be
reproduced and the folder can be copied or archived as it is.

Chain specification
-------------------

``n_monomers``, ``monomers`` (H/G/S fractions), ``linkages`` (a partial dict
is a full specification: unlisted linkages are 0), ``branching``, ``seed``,
``optimize``, and optionally ``experimental_input`` (the LignoForge input JSON)
to estimate any composition not given explicitly.

Runs
----

Runs execute locally as detached processes; ``run.json`` plus an exit-code file
define the state (``running`` / ``done`` / ``failed`` / ``cancelled``), so
closing the application does not lose them.  Kinds: ``em``, ``md`` (full
atomistic workflow) and ``cg`` (coarse-grained system of one or more chains).
The velocity seed is explicit and recorded in the run parameters.

Web interface
-------------

.. code-block:: bash

   pip install 'lignoforge[gui,analysis]'
   lignoforge check            # dependencies and model-validation status
   lignoforge gui              # http://127.0.0.1:8765, projects in ~/LignoForge

The interface covers the whole workflow: **build** chains (explicit composition
or estimated from biomass/process priors or an input JSON, with a live preview
of the resolved specification), **inspect** them in 3-D (colour by monomer,
element, OPLS type, partial charge or residue; atomistic or coarse-grained
beads), browse the **force field** (types, charges, bonded-term counts), read
the **charge-renormalisation report**, see the **linkage graph**, **run**
minimisation / atomistic / coarse-grained simulations locally with live
progress and log, and **analyse** the trajectories (Rg, end-to-end, RMSD,
density, RDF, contact map, energies) with CSV/JSON export.  The environment
page shows which components are missing and the validation status of every
model; light and dark themes are available.

The server listens on ``127.0.0.1`` only and has no authentication; it is
designed for a single local user.  The REST API is documented interactively at
``/api/v1/docs`` (OpenAPI).

Developing the front-end
------------------------

The interface is a React + TypeScript application in ``frontend/``; the built
bundle is written to ``lignoforge/web/static`` and shipped inside the Python
package, so end users do not need Node.js.  To rebuild it::

   tools/build_frontend.sh                      # type-check + bundle
   cd frontend && npm run dev                   # hot reload (needs `lignoforge gui --dev`)
   cd frontend && npx playwright install chromium && npm run test:e2e

The browser tests start the real backend in a temporary workspace and run a
complete session: project, chain build, 3-D view, force field, charges, a real
GROMACS run and its analyses.

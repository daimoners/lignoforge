.. LignoForge documentation root

======================================================
LignoForge — Lignin Structural Modelling Framework
======================================================

.. rubric:: Version |release|

LignoForge is a **top-down, data-driven Python framework** for the automated
generation of atomistic lignin polymer structural models from experimental and
industrial characterisation data.

Given a compact JSON input that encodes NMR-derived monomer fractions, linkage
distributions, and molecular-weight statistics, LignoForge estimates
literature-calibrated structural priors, runs a hierarchical Metropolis Monte
Carlo simulation to grow a population of polymer chains, and exports the results
in multiple machine-readable formats — JSON atomistic / coarse-grained
topologies, PDB structures, SMILES strings, and interactive HTML visualisations.

----

.. grid:: 2

   .. grid-item-card:: 🚀 Getting Started
      :link: getting_started/index
      :link-type: doc

      Install LignoForge and run your first pipeline in under five minutes.

   .. grid-item-card:: 🔬 Concepts
      :link: concepts/index
      :link-type: doc

      Understand the graph model, the MCMC algorithm, the prior
      estimation strategy, and the OPLS-AA parametrisation of lignin
      residues.

.. grid:: 2

   .. grid-item-card:: 📖 User Guide
      :link: user_guide/index
      :link-type: doc

      Task-oriented how-to guides: input format, output files, manual
      chain construction, custom simulations, and GROMACS force-field
      setup.

   .. grid-item-card:: 🗂 API Reference
      :link: api/index
      :link-type: doc

      Complete docstring reference for every public class, function,
      and constant.

----

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   getting_started/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Concepts

   concepts/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: API Reference

   api/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Project

   project/changelog
   project/about

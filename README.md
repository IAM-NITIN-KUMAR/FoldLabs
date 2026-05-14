# FoldBench AF3 Protein Lab

A static web workbench for exploring protein sequences, planning mutations, preparing AlphaFold 3-style input JSON, and inspecting returned structure files.

## Run



## What it does

- Accepts any amino-acid sequence and cleans unsupported symbols.
- Shows sequence length, approximate mass, hydrophobicity, pI, disorder risk, and residue-class map.
- Lets you queue point mutations and preview likely stability/confidence impact.
- Includes familiar experiment presets: single-chain fold, alanine scan, ligand screen, protein complex, stability rescue, and disorder check.
- Tracks basic protocol conditions such as ligand/cofactor, chain copies, pH, and temperature.
- Generates an AlphaFold 3-style JSON package for a backend or local runner.
- Imports `.pdb`, `.cif`, or `.mmcif` results and renders CA/P traces in the interactive viewer.

## AlphaFold 3 note

The app does not submit jobs directly to the official AlphaFold Server. As of May 2026, the official server is browser-first for non-commercial use, while local AlphaFold 3 execution needs a backend with model assets, databases, and appropriate terms/credentials. FoldBench is built as the front-end lab and handoff layer: export JSON, run AF3 elsewhere, then import the result structure here.

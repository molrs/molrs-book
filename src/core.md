# Core

The core functionality of a cheminformatics toolkit is (in my biased opinion) reading and writing molecules from string representations. Once we can accurately and reliably represent molecules, substructures, and transformations as strings, the rest of cheminformatics can be built on top.

To that end, molrs-core is a very simple toolkit focused on three representations:
- Molecule <-> SMILES
- Substructures <- SMARTS
- Transformations <- SMIRKS

# Parsing SMILES
There have been many SMILES parsers written over the years. The approaches range from c++ black magic ([RDKit](https://github.com/rdkit/rdkit)) to abstract syntax trees ([purr](https://crates.io/crates/purr)).

The approach taken by molrs-core is hopefully a middle ground between lower-level, unreadable magic and higher-level approaches with more abstraction.

We also try to minimize the number of classes we'd have to write. We could have a `SmilesAtom` class and an `Atom` class, one for holding raw data out of a SMILES string and one for holding all the data related to an atom post-perception. However, we can avoid this with deliberate use of Rust's Option enum.

We will implement SMILES parsing as the FromStr trait for our Molecule class.

## The Basic Components of a SMILES String
The basic components are:

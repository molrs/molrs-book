# Functions
The public API for Molecule is pretty bare. There are the usual graph-like utilities such as:
- getting all bonds to an atom.
- getting the bond between two atoms.
- getting the indicies of an atom's neighbors.

There are also chemistry-related utilities such as:
- getting the kekulized form of the molecule.
- getting the delocalized form of the molecule.
- getting the explicit valence of an atom.
- getting an atom's theoretical maximum allowed valence.

Because graph-like utilities are not chemistry related and are honestly boring (which is good!) I don't have much to say about them. However, I have written some extended documentation on:
- how molrs kekulizes molecules.
- how molrs delocalized molecules.
- how molrs computes an atom's maximum allowed valence.

All three are non-trivial problems and kekulization in particular stands out as an algorithmic challenge.

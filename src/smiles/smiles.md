# Reading SMILES

Reading a SMILES into a Molecule has two phases. (1) The SMILES has to be correctly parsed into a raw set of atoms and bonds. (2) The raw set of atoms and bonds needs to be perceived into a molecule with rings, implicit hydrogens, explicit bond orders, and stereochemistry.

As an example, here is a very simple SMILES `CC` representing ethane. A more explicit representation of ethane would look like `H3C-CH3` or maybe `[CH3]-[CH3]`. However, SMILES is a much more implicit representation where things like hydrogens and bond orders are removed from the representation and left up to the cheminformatics package to handle.

First `CC` needs to be parsed from a string to a set of atoms and bonds:
```
Atoms: [
    Atom {
        element: C,
    },
    Atom {
        element: C,
    },
]

Bonds: [
    Bond {
        i: 0,
        j: 0,
        bond_type: Default,
    }
]
```

However, a molecule is more than a list of elements and their connectivity. Thus there is a rich set of perception functions needed to flesh out a molecule and its implicit hydrogens, bond orders, and stereochemistry.

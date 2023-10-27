# Parsing SMILES
There have been many SMILES parsers written over the years. The approaches range from c++ black magic ([RDKit](https://github.com/rdkit/rdkit)) to abstract syntax trees ([purr](https://crates.io/crates/purr)).

The approach taken by molrs-core is hopefully a middle ground between lower-level, unreadable magic and higher-level approaches with more abstraction.

See the OpenSMILES specification for more details on SMILES. Note that molrs intentionally does not implement some aspects of the OpenSMILES specification.

## SMILES Syntax
In SMILES syntax, there are only really five actions that need to be handled:
- Add a new atom.
    - Common elements can be specified without a bracket.
    - Uncommon elements or atoms with attributes need a bracket.
- Add a new bond.
- Edit the bond order of the next bond.
- Manage the root atoms of branches.
- Add ring closures.

Handling all five simultaneously in a single pass through the SMILES string ends up looking pretty complicated. So, we will break it down as much as we can. This will (1) convince me I wrote this function correctly (I've rewritten this function three times already...) and (2) help newcomers understand how molrs parses SMILES strings.

We will start with a parser for a minimal SMILES syntax, and work our way up. You can find all the code here, [smiles-parsing](https://github.com/molrs/smiles-parsing).

### Single Atom, No Brackets
Let's imagine a universe where every molecule can only contain one atom. Let's also imagine that this universe only contains the common elements allowed in SMILES:
- wildcard (*)
- boron
- carbon
- nitrogen
- oxygen
- phosphorus
- sulfur
- fluorine
- chlorine
- bromine
- iodine

In this universe, SMILES syntax is very simple and its related parser is also very simple. All we need to do is read the smiles and find what element it corresponds to. Here's our v1:
```rust
fn smiles_parser_v1(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];

    for c in smi.chars() {
        if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        };
    }

    Ok(Molecule {
        atoms,
        bonds: vec![],
        rings: None,
    })
}
```

Using this function to parse "C" gives us a Molecule with no bonds and one atom:
```
Atom {
    element: C,
    isotope: None,
    charge: 0,
    delocalized: false,
    n_implicit_hydrogens: None,
    n_radical_electrons: None,
    point_chirality: Undefined,
},
```

It works! But what if we want more than this limited set of atoms?

### Single Atom, With Brackets
For elements beyond the small set listed above, it's necessary to introduce bracket syntax. For a deeper dive on bracket syntax, I recommend reading the OpenSMILES standard.

Let's expand our parser to allow bracket atoms. Here's v2:
```rust
fn smiles_parser_v2(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];

    // *** NEW CODE BEGINS ***
    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    // *** NEW CODE ENDS ***

    for c in smi.chars() {
        // *** NEW CODE BEGINS ***
        if c_is_in_bracket {
            let atom: &mut Atom = atoms.last_mut().unwrap();
            if c == ']' {
                c_is_in_bracket = false;
                atom.element = match element_str.parse() {
                    Ok(element) => element,
                    Err(_) => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | invalid element {element_str}"
                        )));
                    }
                };
                if element_str.chars().next().unwrap().is_lowercase() {
                    atom.delocalized = true;
                };
                element_str = String::new();
            } else if c == 'H' {
                if element_str.is_empty() {
                    element_str.push('H');
                } else {
                    atom.n_implicit_hydrogens = Some(1);
                    atom_attribute = AtomAttribute::NImplicitHydrogens;
                };
            } else if c == '+' {
                atom.charge = 1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '-' {
                atom.charge = -1;
                atom_attribute = AtomAttribute::Charge;
            } else if c == '@' {
                match atom.point_chirality {
                    PointChirality::Undefined => {
                        atom.point_chirality = PointChirality::CounterClockwise;
                    }
                    PointChirality::CounterClockwise => {
                        atom.point_chirality = PointChirality::Clockwise;
                    }
                    _ => {
                        return Err(MoleculeError::SmilesParseError(format!(
                            "{smi} | chirality error"
                        )));
                    }
                }
            } else if c.is_alphabetic() {
                element_str.push(c);
            } else if c.is_numeric() {
                let c_as_digit = c.to_digit(10).unwrap();
                match atom_attribute {
                    AtomAttribute::Isotope => match atom.isotope {
                        None => {
                            atom.isotope = Some(c_as_digit as u16);
                        }
                        Some(isotope) => atom.isotope = Some(isotope * 10 + c_as_digit as u16),
                    },
                    AtomAttribute::Charge => {
                        atom.charge *= c_as_digit as i8;
                    }
                    AtomAttribute::NImplicitHydrogens => {
                        atom.n_implicit_hydrogens = Some(c_as_digit as u8);
                    }
                };
            };
        } else if c == '[' {
            c_is_in_bracket = true;
            atom_attribute = AtomAttribute::Isotope;
            atoms.push(Atom::default());
        // *** NEW CODE ENDS ***
        } else if c == 'B'
            || c == 'C'
            || c == 'N'
            || c == 'O'
            || c == 'P'
            || c == 'S'
            || c == 'F'
            || c == 'I'
            || c == '*'
        {
            atoms.push(Atom {
                element: Element::from_str(std::str::from_utf8(&[c as u8]).unwrap()).unwrap(),
                isotope: None,
                charge: 0,
                delocalized: false,
                n_implicit_hydrogens: None,
                n_radical_electrons: None,
                point_chirality: molrs::atom::PointChirality::Undefined,
            });
        } else if c == 'l' && atoms.last().unwrap().element == Element::C {
            atoms.last_mut().unwrap().element = Element::Cl;
        } else if c == 'r' && atoms.last().unwrap().element == Element::B {
            atoms.last_mut().unwrap().element = Element::Br;
        } else {
            return Err(MoleculeError::SmilesParseError(format!(
                "{smi} | invalid char {c}"
            )));
        }
    }

    Ok(Molecule {
        atoms,
        bonds: vec![],
        rings: None,
    })
}
```

As you can see, our code is getting a little more complicated. There are two main points of complication that we had to handle:
- Handle chars inside a bracket and outside a bracket differently.
- Handle all the different atom attributes specifiable by bracket syntax.
    - Element
    - Isotope
    - Charge
    - Number of hydrogens
    - Chirality

Using this function to parse "[18OH-]" gives a Molecule with no bonds, and one atom:
```
Atom {
    element: O,
    isotope: Some(
        18,
    ),
    charge: -1,
    delocalized: false,
    n_implicit_hydrogens: Some(
        1,
    ),
    n_radical_electrons: None,
    point_chirality: Undefined,
}
```

### Linear Connectivity
It's nice to read single atoms but chemistry is all about bonding. SMILES has three kinds of bonds:
- linear bonds.
- branched bonds.
- ring-closing bonds.

Let's start with the simplest, linear bonds. Here's v3:
```rust

```

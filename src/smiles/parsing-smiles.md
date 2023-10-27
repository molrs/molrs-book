# Parsing SMILES
There have been many SMILES parsers written over the years. The approaches range from c++ black magic ([RDKit](https://github.com/rdkit/rdkit)) to abstract syntax trees ([purr](https://crates.io/crates/purr)).

The approach taken by molrs-core is hopefully a middle ground between lower-level, unreadable magic and higher-level approaches with more abstraction.

We also try to minimize the number of classes we'd have to write. We could have a `SmilesAtom` class and an `Atom` class, one for holding raw data out of a SMILES string and one for holding all the data related to an atom post-perception. However, we can avoid this with deliberate use of Rust's Option enum.

We will implement SMILES parsing as the FromStr trait for our Molecule class.

## SMILES Syntax
In SMILES syntax, there are only really 4 actions that need to be handled:
- Add a new atom.
    - Common elements can be specified without a bracket.
    - Uncommon elements or atoms with attributes need a bracket.
- Add a new bond.
- Edit the bond order of the next bond.
- Add ring closures.

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

### Single Atom, With Brackets
While we can read a limited set of elements, we really want to be able to read all elements. For elements beyond the small set listed above, it's necessary to introduce bracket syntax. For a deeper dive on SMILES syntax, I recommend reading the OpenSMILES standard.

Let's expand our parser to allow bracket atoms. Here's v2:
```rust
fn smiles_parser_v2(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();

    for c in smi.chars() {
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

As you can see, our code is getting a little ... spaghetti-like. There are two main points of complication that we had to handle:
- Handle chars inside a bracket and outside a bracket differently.
- Handle all the different atom attributes specifiable by bracket syntax.

However, it works as expected. Using this function to parse "[18OH-]" gives a Molecule with no bonds, and one atom:
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

Let's start with the simplest, linear bonds.



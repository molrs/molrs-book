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
fn smiles_parser_v3(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    // *** NEW CODE BEGINS ***
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];
    // *** NEW CODE ENDS ***

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();

    // *** NEW CODE BEGINS ***
    for (i, c) in smi.chars().enumerate() {
    // *** NEW CODE ENDS ***
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
            // *** NEW CODE BEGINS ***
            if i > 0 {
                bond.i = i - 1;
                bond.j = i;
                bonds.push(bond);
                bond = Bond::default();
            }
            // *** NEW CODE ENDS ***
        // *** NEW CODE BEGINS ***
        } else if c == '-'
            || c == '/'
            || c == '\\'
            || c == ':'
            || c == '='
            || c == '#'
            || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
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
            // *** NEW CODE BEGINS ***
            if i > 0 {
                bond.i = i - 1;
                bond.j = i;
                bonds.push(bond);
                bond = Bond::default();
            }
            // *** NEW CODE ENDS ***
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
        bonds,
        rings: None,
    })
}
```

Because bonds are only linear, connectivity will always be 0-1, 1-2, 2-3, etc. Thus for linear connectivity, we can just use a counter to keep track of atom indices.

Note that atoms are pushed when the first char of the atom is seen. However, bonds are pushed when the latter atom of the bond is seen. This is an important nuance and cadence of pushing atoms and bonds to their respective vectors is an important choice.

Parsing a linear molecule like "CC=C" gives the following bonds:
```
Bond {
    i: 0,
    j: 1,
    bond_type: Default,
},
Bond {
    i: 2,
    j: 3,
    bond_type: Double,
},
```

Let's now move on to branching.

### Branches

I would say branching was the most difficult aspect of SMILES to get right. Instead of just using a simple counter, you need to keep track of where the root atom (the i-th atom) of the next bond should be. We can do this using a stack!

```rust
fn smiles_parser_v4(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    // *** NEW CODE BEGINS ***
    let mut root_atom: Vec<usize> = vec![];
    // *** NEW CODE ENDS ***

    // *** NEW CODE BEGINS ***
    for c in smi.chars() {
    // *** NEW CODE ENDS ***
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
            // *** NEW CODE BEGINS ***
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
            // *** NEW CODE ENDS ***
        } else if c == '-'
            || c == '/'
            || c == '\\'
            || c == ':'
            || c == '='
            || c == '#'
            || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
        // *** NEW CODE BEGINS ***
        } else if c == '(' {
            root_atom.push(*root_atom.last().unwrap());
        } else if c == ')' {
            root_atom.pop();
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
            // *** NEW CODE BEGINS ***
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
            // *** NEW CODE ENDS ***
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
        bonds,
        rings: None,
    })
}
```

To handle branches, we keep a stack of root atoms. Here is how it works on a simple example, "CC(CC)C". Imagine we are iterating through this string char by char.
```
i = 0, c = 'C':
    stack = []
    no bond pushed
    push 0 to stack
    stack = [0]
i = 1, c = 'C':
    stack = [0]
    consume 0 and push bond 0-1
    push 1 to stack
    stack = [1]
i = 2, c = '(':
    stack = [1]
    no bond pushed
    push 1 to stack
    stack = [1, 1]
i = 3, c = 'C':
    stack = [1, 1]
    consume 1 and push bond 1-2
    push 2 to stack
    stack = [1, 2]
i = 4, c = 'C':
    stack = [1, 2]
    consume 2 and push bond 2-3
    push 3 to stack
    stack = [1, 3]
i = 5, c = ')':
    stack = [1, 3]
    pop 3
    stack = [1]
i = 6, c = 'C':
    stack = [1]
    consume 1 and push bond 1-4
    stack = [4]
```

Here are some example molecules:
smi = "CC(=O)C", simple branching
```
Bond {
    i: 0,
    j: 1,
    bond_type: Default,
},
Bond {
    i: 1,
    j: 2,
    bond_type: Double,
},
Bond {
    i: 1,
    j: 3,
    bond_type: Default,
},
```

smi = "CS(=O)(=O)C", serial branching
```
Bond {
    i: 0,
    j: 1,
    bond_type: Default,
},
Bond {
    i: 1,
    j: 2,
    bond_type: Double,
},
Bond {
    i: 1,
    j: 3,
    bond_type: Double,
},
Bond {
    i: 1,
    j: 4,
    bond_type: Default,
},
```

smi = "CC(C(F)F)C", nested branching
```
Bond {
    i: 0,
    j: 1,
    bond_type: Default,
},
Bond {
    i: 1,
    j: 2,
    bond_type: Default,
},
Bond {
    i: 2,
    j: 3,
    bond_type: Default,
},
Bond {
    i: 2,
    j: 4,
    bond_type: Default,
},
Bond {
    i: 1,
    j: 5,
    bond_type: Default,
},
```

And there you have it, we can handle branched SMILES.

### Ring Closures
Finally we come to ring closures. Interestingly, ring closing is the only part of SMILES that isn't linear. When the first atom is marked as a ring-closing atom, you are making a promise to come back later and close it. Then when the second, matching atom is found, you need to connect them together.

Here is an implementation of ring closures. Note that OpenSMILES allows double-digit ring indexes with the '%' char marking double-digit ring indexes. In this minimal version, we leave it out.

```rust
fn smiles_parser_v5(smi: &str) -> Result<Molecule, MoleculeError> {
    let mut atoms: Vec<Atom> = vec![];
    let mut bond: Bond = Bond::default();
    let mut bonds: Vec<Bond> = vec![];
    // *** NEW CODE BEGINS ***
    let mut ring_closures: HashMap<usize, usize> = HashMap::new();
    // *** NEW CODE ENDS ***

    let mut c_is_in_bracket: bool = false;
    let mut atom_attribute: AtomAttribute = AtomAttribute::Isotope;
    let mut element_str: String = String::new();
    let mut root_atom: Vec<usize> = vec![];

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
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
        } else if c == '-' || c == '/' || c == '\\' || c == ':' || c == '=' || c == '#' || c == '$'
        {
            bond.bond_type = BondType::try_from(c).unwrap();
        // *** NEW CODE BEGINS ***
        } else if c.is_numeric() {
            let c_as_digit: usize = c.to_digit(10).unwrap() as usize;
            if let std::collections::hash_map::Entry::Vacant(e) = ring_closures.entry(c_as_digit) {
                e.insert(atoms.len() - 1);
            } else {
                let ring_closure_bond = Bond {
                    i: *ring_closures.get(&c_as_digit).unwrap(),
                    j: atoms.len() - 1,
                    bond_type: bond.bond_type,
                };
                bond.bond_type = BondType::Default;
                bonds.push(ring_closure_bond);
                ring_closures.remove(&c_as_digit);
            }
        // *** NEW CODE ENDS ***
        } else if c == '(' {
            root_atom.push(*root_atom.last().unwrap());
        } else if c == ')' {
            root_atom.pop();
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
            if !root_atom.is_empty() {
                bond.i = root_atom.pop().unwrap();
                bond.j = atoms.len() - 1;
                bonds.push(bond);
                bond = Bond::default();
            }
            root_atom.push(atoms.len() - 1);
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
        bonds,
        rings: None,
    })
}
```

We just need to keep track of ring closure indexes. For this, we use a HashMap (like a dict in Python) where the keys are ring-closure indexes and the values are atom indices.

For "C1CC1", we get:
```
Bond {
    i: 0,
    j: 1,
    bond_type: Default,
},
Bond {
    i: 1,
    j: 2,
    bond_type: Default,
},
Bond {
    i: 0,
    j: 2,
    bond_type: Default,
},
```

## Conclusion
With that, we have finished our SMILES parser. There are a few things we have chosen not to implement:
- Double digit ring closure indexes, prefixed with "%".
- Point chirality beyond tetrahedral centers.
- Allowing multi-molecule molecules (we don't allow molecules to be separated using ".").

However, look at all we did implement! If you really read all the way through that, I commend you. Hopefully you learned somethng about cheminformatics and feel more confident in using molrs yourself.

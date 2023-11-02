# Writing SMILES
Now that we can accurate parse SMILES and perceive implicit properties of molecules, we need to be able to store them back into SMILES.

Our plan will be as follows:
- Convert atoms to separate single-atom SMILES strings.
- Prepend bond type to single-atom SMILES.
- Append ring closure indexes to single-atom SMILES.
- Append or prepend parentheses to single-atom SMILES.

Consider the following SMILES string, "C1(=CCCC1)C". It encodes 1-methyl cyclopentene. If we try to partition this SMILES into single-atom SMILES, we would probably partition something like this:
- "C1"
- "(=C"
- "C"
- "C"
- "C1)"
- "C

As you see, ring closures and bond types are on the inside and parentheses are on the outside.

Here is the code, we will break it down block by block.
```rust
impl ToString for Molecule {
    fn to_string(&self) -> String {
        let mut atom_strs: Vec<String> = self.atoms.iter().map(|atom| atom.to_string()).collect();
        let mut ring_closure_index = 1;

        for (i, atom) in self.atoms.iter().enumerate() {
            if atom.delocalized
                && atom.n_implicit_hydrogens == Some(1)
                && atom.element == Element::N
            {
                atom_strs[i] = "[nH]".to_owned();
            }
        }

        let mut vec_neighbors: Vec<Vec<usize>> = (0..self.atoms.len())
            .map(|i| {
                let mut neighbors = self.atom_neighbor_indicies(i);
                neighbors.retain(|neighbor| *neighbor < i);
                neighbors.sort();

                neighbors
            })
            .collect();

        for (i, neighbors) in vec_neighbors.iter_mut().enumerate() {
            if neighbors.is_empty() {
                continue;
            }

            let last_neighbor = *neighbors.last().unwrap();
            let bond = self.atoms_bond_between(i, last_neighbor).unwrap();
            if bond.bond_type != BondType::Single
                && bond.bond_type != BondType::Delocalized
                && bond.bond_type != BondType::Default
            {
                atom_strs[i].insert(0, bond.bond_type.into());
            }

            for neighbor in &neighbors[..neighbors.len() - 1] {
                let bond = self.atoms_bond_between(i, *neighbor).unwrap();
                if bond.bond_type != BondType::Single
                    && bond.bond_type != BondType::Delocalized
                    && bond.bond_type != BondType::Default
                {
                    atom_strs[i].push(bond.bond_type.into());
                }

                let ring_closure_string: String = if ring_closure_index > 9 {
                    format!("%{ring_closure_index}")
                } else {
                    format!("{ring_closure_index}")
                };
                atom_strs[i].push_str(&ring_closure_string);
                atom_strs[*neighbor].push_str(&ring_closure_string);
                ring_closure_index += 1;
            }
        }

        for (i, neighbors) in vec_neighbors.iter_mut().enumerate() {
            if neighbors.is_empty() {
                continue;
            }

            let last_neighbor = *neighbors.last().unwrap();
            if last_neighbor == i - 1 {
                continue;
            }

            let mut parenthesis_level = 0;
            let mut cursor = (last_neighbor, atom_strs[last_neighbor].chars().count());
            for (j, atom_str) in atom_strs.iter().enumerate().take(i).skip(last_neighbor + 1) {
                for (k, c) in atom_str.chars().enumerate() {
                    if c == '(' {
                        parenthesis_level += 1;
                    }
                    if parenthesis_level == 0 && c == ')' {
                        cursor = (j, k + 1);
                    }
                    if c == ')' {
                        parenthesis_level -= 1;
                    }
                }
            }
            atom_strs[cursor.0].insert(cursor.1, '(');
            atom_strs[i].insert(0, ')');
        }

        atom_strs.join("")
    }
}
```

## Set-up
```rust
        let mut atom_strs: Vec<String> = self.atoms.iter().map(|atom| atom.to_string()).collect();
        let mut ring_closure_index = 1;

        for (i, atom) in self.atoms.iter().enumerate() {
            if atom.delocalized
                && atom.n_implicit_hydrogens == Some(1)
                && atom.element == Element::N
            {
                atom_strs[i] = "[nH]".to_owned();
            }
        }
```

There are two things to set up. First, we need to convert the atoms to a vector of single-atom SMILES. We also create our ring closure counter.

Then we have our manual edits/special rules. So far, we only have one. Aromatic nitrogens with implicit hydrogens need to have their hydrogens explicitly written.

## Collect Neighbors
```rust
        let mut vec_neighbors: Vec<Vec<usize>> = (0..self.atoms.len())
            .map(|i| {
                let mut neighbors = self.atom_neighbor_indicies(i);
                neighbors.retain(|neighbor| *neighbor < i);
                neighbors.sort();

                neighbors
            })
            .collect();
```

Next we obtain a vector of the neighbors for each atom. But we only keep the neighbors that have an index lower than the atom. Ie. for "CCC", `vec_neighbors` will be:
```
[
    [],
    [0],
    [1],
]
```

For "C1CC1", `vec_neighbors` will be:
```
[
    [],
    [0],
    [0, 1],
]
```

Atoms without ring-closing bonds will have one neighbor with an index less than itself. However, atoms with ring-closing bonds will have multiple neighbors with indices less than itself.

Understanding writing SMILES requires being able to reason about the connectivity of a molecule through this `vec_neighbors`.

## Add Bond Types and Ring Closures
```rust
        for (i, neighbors) in vec_neighbors.iter_mut().enumerate() {
            if neighbors.is_empty() {
                continue;
            }

            let last_neighbor = *neighbors.last().unwrap();
            let bond = self.atoms_bond_between(i, last_neighbor).unwrap();
            if bond.bond_type != BondType::Single
                && bond.bond_type != BondType::Delocalized
                && bond.bond_type != BondType::Default
            {
                atom_strs[i].insert(0, bond.bond_type.into());
            }

            for neighbor in &neighbors[..neighbors.len() - 1] {
                let bond = self.atoms_bond_between(i, *neighbor).unwrap();
                if bond.bond_type != BondType::Single
                    && bond.bond_type != BondType::Delocalized
                    && bond.bond_type != BondType::Default
                {
                    atom_strs[i].push(bond.bond_type.into());
                }

                let ring_closure_string: String = if ring_closure_index > 9 {
                    format!("%{ring_closure_index}")
                } else {
                    format!("{ring_closure_index}")
                };
                atom_strs[i].push_str(&ring_closure_string);
                atom_strs[*neighbor].push_str(&ring_closure_string);
                ring_closure_index += 1;
            }
        }
```

There are two phases to this for loop. First, we find the bond between the atom and the neighbor with the highest index. This is because the neighbor with the highest index is the one that is linearly (no rings) connected to the atom. We can then easily extract the bond type and prepend the `atom_str` with the bond type as a char.

Second, we handle neighbors that aren't connected linearly (with rings). We find the relevant ring-closing bond and append the bond type as a char to our atom (ie. "C1CC=1" encodes cyclopropene). Then we format our ring closure as a string and append it to both atoms in the ring-closing bond.

## Handling Branches
```rust
        for (i, neighbors) in vec_neighbors.iter_mut().enumerate() {
            if neighbors.is_empty() {
                continue;
            }

            let last_neighbor = *neighbors.last().unwrap();
            if last_neighbor == i - 1 {
                continue;
            }

            let mut parenthesis_level = 0;
            let mut cursor = (last_neighbor, atom_strs[last_neighbor].chars().count());
            for (j, atom_str) in atom_strs.iter().enumerate().take(i).skip(last_neighbor + 1) {
                for (k, c) in atom_str.chars().enumerate() {
                    if c == '(' {
                        parenthesis_level += 1;
                    }
                    if parenthesis_level == 0 && c == ')' {
                        cursor = (j, k + 1);
                    }
                    if c == ')' {
                        parenthesis_level -= 1;
                    }
                }
            }
            atom_strs[cursor.0].insert(cursor.1, '(');
            atom_strs[i].insert(0, ')');
        }
```

Here we come to the most complicated part of the algorithm. Now that we've prepended bond types and append ring-closures to our single-atom SMILES, we're ready to add the rest of the connectivity information with parentheses.

If the atom's linear bond is with its direct neighbor, (`if last_neighbor == i - 1`) we don't need to do anything. But if the root atom of the bond is an atom that occurs further away, we need to wrap things in parentheses, ie "CC(=O)C" for acetone. The third carbon need to bond with the second carbon so we wrap the oxygen with its bond type in parentheses.

However, this naive approach gets ridiculous once you get serial branching. For SF6, we would want to write it as "FS(F)(F)(F)(F)F". However, naively wrapping the intermediate atoms in parentheses will yield "FS((((F)F)F)F)F". While this is a valid SMILES ... it's not very pretty.

To make serial branching look prettier, we keep track of a cursor. The cursor is initialized to match the naive approach, where the '(' is appended to the root atom of the bond.

However, if we find a ')' in the SMILES in between the root atom and our atom, we will update the cursor to insert the '(' somewhere more aesthetically pleasing but still correct.

Lastly, we need to keep track of the parenthesis depth to handle nested branching, ie. "CC(CC(F)C)(Cl)C". While trying to wrap the "Cl" in parentheses, we need to make sure we ignore the "(F)".

## Conclusion
Because we were able to abstract away the Atom -> String conversion, writing SMILES ends up being much shorter in lines of code. However, as you can see, accurately encoding the branching and cyclized connectivity of a Molecule into a linear String takes some careful reasoning.

# Kekulization
Kekulization is the process of resolving delocalized bonds into specific single or double bonds. For example we can write benzene as a delocalized, 6-fold symmetric SMILES "c1ccccc1" or as a kekulized, 3-fold symmetric SMILES "C1=CC=CC=C1".

While the delocalized form is more convenient to work with, kekulization is necessary for perceiving implicit hydrogens as well as making sure we have a valid molecule. If you never kekulize your molecules, you may be passing around a molecule as absurd as "c1cc1" and never know it.

To kekulize a molecule, we will first iterate through the molecule's rings to find paths in the molecule that need kekulization. Then if we find a path that has an even number of atoms in it, we will kekulize it into a alternating chain of double and single bonds. If the path has an odd number of atoms in it, we will leave it for later and hope that by kekulizing another path, we will be able to kekulize it.

Because molrs keeps track of all rings, the rings of a molecule are overlapping. This removes the need to reiterate through the rings as there is already redundancy in coverage among the rings.

At the end, if there are delocalized atoms remaining, we will consider kekulization to have failed.

## Set-up
```rust
    let mut mol = self.clone();
    let rings = match &self.rings {
        Some(rings) => rings,
        None => {
            return Err(MoleculeError::MissingRingsError(
                "how did you get here?".to_owned(),
            ));
        }
    };
```

We first create a new mol as clone of self. We also borrow rings from self. It should be impossible to not have rings initialized but in case someone is creating a molecule manually (not from SMILES) we need to catch this error.

## Find and Kekulize Paths
```rust
    for ring in rings.iter().rev() {
        // phase 1
        let mut path_breaks = vec![];
        for (i, index) in ring.iter().enumerate() {
            if !mol.atom_needs_kekulization(*index) {
                path_breaks.push(i);
                mol.atoms[*index].delocalized = false;
                for bond in mol.atom_bonds_mut(*index) {
                    if bond.bond_type == BondType::Delocalized {
                        bond.bond_type = BondType::Single;
                    }
                }
            }
        }

        // phase 2
        let mut paths = vec![];
        if path_breaks.is_empty() {
            if ring.len() % 2 == 0 {
                paths.push(ring.clone());
            }
        } else if path_breaks.len() == 1 {
            let path_break = path_breaks[0];
            let mut path = ring[path_break + 1..].to_owned();
            path.extend(&ring[..path_break]);
            paths.push(path);
        } else {
            for window in path_breaks.windows(2) {
                let path = ring[window[0] + 1..window[1]].to_owned();
                if !path.is_empty() {
                    paths.push(path);
                }
            }
            let mut path = ring[path_breaks.last().unwrap() + 1..].to_owned();
            path.extend(&ring[..*path_breaks.first().unwrap()]);
            if !path.is_empty() {
                paths.push(path);
            }
        }

        // phase 3
        for path in paths {
            if path.len() % 2 == 0 {
                for (i, window) in path.windows(2).enumerate() {
                    if i % 2 == 0 {
                        mol.atoms_bond_between_mut(window[0], window[1])
                            .unwrap()
                            .bond_type = BondType::Double;
                    } else {
                        mol.atoms_bond_between_mut(window[0], window[1])
                            .unwrap()
                            .bond_type = BondType::Single;
                    }
                }
                if path.len() > 2 {
                    if let Some(bond) = mol
                        .atoms_bond_between_mut(*path.first().unwrap(), *path.last().unwrap())
                    {
                        bond.bond_type = BondType::Single;
                    }
                }
                for i in path {
                    mol.atoms[i].delocalized = false;
                }
            }
        }
    }
```

There are three phases to this algorithm:
1. We first traverse the ring to find atoms that don't need to be kekulized. For example, pyrrole "c1[nH]ccc1" has a nitrogen with two bonds and an implicit hydrogen. Because its valence is full, it is unable to participate in a double bond. Thus, we could consider this nitrogen as not needing kekulization. The atoms that don't need kekulization are considered as path breaks that split the ring into paths that do need kekulization.
2. Then, we use the path breaks created in phase 1 to find kekulization paths. This is unfortunately a little complicated because we need to handle things slightly differently if `path_breaks.len() == 1` or if `path_breaks.len() > 1`. The other complication is that rings are cyclic but vectors are linear. So we need to manually connect paths that wrap around the vector representing the ring.
3. Lastly, we iterate over the kekulization paths, resolving delocalized bonds into single and double bonds. The one sticky part here is handling the bond that connects the first and last atom of a ring.

## Check for Completion
```rust
    if mol.atoms.iter().any(|atom| atom.delocalized) {
        return Err(MoleculeError::KekulizationError(format!(
            "{} | could not be kekulized",
            mol.to_string()
        )));
    }
    if mol
        .bonds
        .iter()
        .any(|bond| bond.bond_type == BondType::Delocalized)
    {
        return Err(MoleculeError::KekulizationError(format!(
            "{} | could not be kekulized",
            mol.to_string()
        )));
    }

    Ok(mol)
```

Then we iterate over the atoms and bonds to ensure that there are no delocalized labels remaining. Finally we can return the kekulized mol.

## When Does an Atom Need Kekulization?
```rust
fn atom_needs_kekulization(&self, index: usize) -> bool {
    let atom = self.atoms[index];
    atom.delocalized
        && self.atom_n_double_bonds(index) == 0
        && (self.atom_explicit_valence(index) + atom.n_implicit_hydrogens.unwrap_or(0))
            < self.atom_maximum_allowed_valence(index)
}
```

An atom needs kekulization if it is delocalized AND it has not double bonds AND if the explicit valence + number of implicit hydrogens is less than the maximum allowed valence.

Note that delocalized bonds are underestimated to have an explicit valence of 1.

It's also useful to think about it the other way, an atom doesn't need kekulization if it isn't delocalized OR it has double bonds OR its explicit valence + number of implicit hydrogens is equal to or greater than the maximum allowed valence.

## Conclusion
Kekulization took me the longest time out of all the basic molecule functions to get right. It's complicated and somewhat hard to reason about, especially making sure that you've reached completion. Thankfully, after 3 refactors, I've gotten the code to the point where (hopefully) it's exactly as complicated as it needs to be while being correct.

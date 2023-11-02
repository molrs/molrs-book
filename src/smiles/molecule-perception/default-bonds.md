# Default Bonds
SMILES doesn't explicitly provide bond types for single bonds and aromatic bonds. Thankfully, we can simply resolve default bonds. SMILES atoms do explicitly label if they are aromatic or not so we can use these labels to label a bond as single or aromatic.

Note that molrs doesn't like using the label "aromatic". There examples of 4n + 2 pi systems that aren't aromatic. The word "aromatic" is way too loaded and so I prefer the less loaded term "delocalized".

While we will make a different error of assuming 4n pi systems are delocalized (ie. cyclobutadiene), these kinds of systems do shift double bonds around through valence isomerization. For example, this [cyclobutadiene system](https://doi.org/10.1002/anie.199207381) has a valence isomerization barrier of 5.8 kcal/mol.

Ultimately, we will not be able to accurately label molecules as "aromatic" or "delocalized" without QM methods. So we will have to settle for different degrees of "wrongness".

```rust
fn perceive_default_bonds(&mut self) {
    for bond in self.bonds.iter_mut() {
        if bond.bond_type != BondType::Default {
            continue;
        }
        if self.atoms.get(bond.i).unwrap().delocalized
            && self.atoms.get(bond.j).unwrap().delocalized
        {
            bond.bond_type = BondType::Delocalized;
        } else {
            bond.bond_type = BondType::Single;
        }
    }
}
```

As you see, its programmatically quite simple. If both atoms in a default bond are labeled as delocalized, we will set that bond as delocalized. If not, the bond will be set as single.

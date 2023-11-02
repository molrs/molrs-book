# Delocalization


```rust
pub fn delocalized(&self) -> Molecule {
    let mut mol = self.clone();

    for ring in self.rings.as_ref().unwrap().iter().rev() {
        if ring
            .iter()
            .all(|index| self.atom_needs_delocalization(*index))
        {
            for index in ring {
                mol.atoms[*index].delocalized = true;
            }
            mol.atoms_bond_between_mut(*ring.first().unwrap(), *ring.last().unwrap())
                .unwrap()
                .bond_type = BondType::Delocalized;
            for window in ring.windows(2) {
                mol.atoms_bond_between_mut(window[0], window[1])
                    .unwrap()
                    .bond_type = BondType::Delocalized;
            }
        }
    }

    mol
}
```

# Implicit Hydrogens

Implicit hydrogens are probably the most important aspect of molecule perception. While the function below is not too complicated, there are two functions it heavily depends on that are not trivial:
1. `kekulized` is a whole other rabbit hole with its own page.
2. `atom_maximum_allowed_valence` is also very complicated as it has many special rules.

We need a reference molecule that is kekulized because molecules in their delocalized form use delocalized bonds. However, delocalized bonds have an amgiguous bond order. Thus for accurate (and easy) implicit hydrogen perception, it's highly desirable to start from a kekulized form.


```rust
pub fn perceive_implicit_hydrogens(&mut self) -> Result<(), MoleculeError> {
    let mol = self.kekulized()?;

    let bond_orders: Vec<u8> = (0..mol.atoms.len())
        .map(|index| mol.atom_explicit_valence(index))
        .collect();
    let maximum_allowed_valences: Vec<u8> = (0..mol.atoms.len())
        .map(|index| mol.atom_maximum_allowed_valence(index))
        .collect();

    for (i, atom) in self.atoms.iter_mut().enumerate() {
        let bond_order = bond_orders[i];
        let maximum_allowed_valence = maximum_allowed_valences[i];

        if bond_order > maximum_allowed_valence {
            return Err(MoleculeError::BondOrderError(format!(
                "explicit bond order is higher than maximum allowed valence for atom {i}"
            )));
        }

        let n_implicit_hydrogens = maximum_allowed_valence - bond_order;
        if atom.n_implicit_hydrogens.is_none() {
            atom.n_implicit_hydrogens = Some(n_implicit_hydrogens);
            atom.n_radical_electrons = Some(0);
        } else {
            atom.n_radical_electrons =
                Some(n_implicit_hydrogens - atom.n_implicit_hydrogens.unwrap());
        }
    }

    Ok(())
}
```

## Set-up
```rust
    let mol = self.kekulized()?;

    let bond_orders: Vec<u8> = (0..mol.atoms.len())
        .map(|index| mol.atom_explicit_valence(index))
        .collect();
    let maximum_allowed_valences: Vec<u8> = (0..mol.atoms.len())
        .map(|index| mol.atom_maximum_allowed_valence(index))
        .collect();
```

We first create our kekulized mol. Note that we don't modify this kekulized mol nor do we return it. It exists simply for us to reference it.

We also pre-compute our bond orders and maximum allowed valences. We will make use of these to compute the number of implicit hydrogens.

## Simple Math
```rust
    for (i, atom) in self.atoms.iter_mut().enumerate() {
        if bond_order > maximum_allowed_valence {
            return Err(MoleculeError::BondOrderError(format!(
                "explicit bond order is higher than maximum allowed valence for atom {i}"
            )));
        }

        let n_implicit_hydrogens = maximum_allowed_valence - bond_order;
        if atom.n_implicit_hydrogens.is_none() {
            atom.n_implicit_hydrogens = Some(n_implicit_hydrogens);
            atom.n_radical_electrons = Some(0);
        } else {
            atom.n_radical_electrons =
                Some(n_implicit_hydrogens - atom.n_implicit_hydrogens.unwrap());
        }
    }
```

First we check that the bond order doesn't exceed the maximum allowed valence.

After that we can compute the number of implicit hydrogens as the difference between the explicit valence (bond order) and the implicit valence (maximum allowed valence). This of course makes the assumption that implicit valence = maximum allowed valence. This means we have to write our `atom_maximum_allowed_valence` function really well.

If the atom doesn't have `n_implicit_hydrogens` set already, we can set it. If it does already, that means the number of implicit hydrogens was explicitly encoded in the SMILES. In that case, we handle potential radical electrons as expected number of implicit hydrogens - actual number of implicit hydrogens.

## Conclusion
As you can see, perceiving implicit hydrogens is actually pretty trivial. However, it's only trivial if you already have implemented kekulization and computing of maximum allowed valence.

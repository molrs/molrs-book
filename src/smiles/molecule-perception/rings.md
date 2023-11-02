# Rings
Traversing the molecular graph to find rings is a somewhat expensive process. However, it is convenient for downstream function to have access to the rings in a molecule.

We will find rings by traversing the graph and enumerating all paths in the molecule. Then while enumerating paths, we either:
1. Store paths that loop back on themselves and remove them from further path enumeration.
2. Remove paths that reach dead ends from further path enumeration.

Once we've exhaustively enumerated all paths, we deduplicate the rings (we will traverse each ring forwards and backwards) and save them to `self.rings`.

```rust
    pub fn perceive_rings(&mut self) {
        let mut paths = vec![];
        let mut closed_paths = vec![];
        for neighbor_index in self.atom_neighbor_indicies(0) {
            paths.push(Some(vec![0, neighbor_index]));
        }

        while !paths.is_empty() {
            let mut new_paths = vec![];
            for path in paths.iter_mut() {
                if path.is_none() {
                    continue;
                }
                let last_atom_index = *path.as_ref().unwrap().last().unwrap();
                let second_to_last_atom_index = path.as_ref().unwrap().iter().rev().nth(1).unwrap();
                let mut neighbor_indices = self.atom_neighbor_indicies(last_atom_index);
                neighbor_indices
                    .retain(|neighbor_index| neighbor_index != second_to_last_atom_index);
                if neighbor_indices.is_empty() {
                    *path = None;
                } else {
                    for neighbor_index in &neighbor_indices[1..] {
                        let mut new_path = path.clone().unwrap();
                        new_path.push(*neighbor_index);
                        new_paths.push(new_path);
                    }
                    path.as_mut().unwrap().push(neighbor_indices[0]);
                }
            }
            for new_path in new_paths {
                paths.push(Some(new_path));
            }
            for path_option in paths.iter_mut() {
                if let Some(path) = path_option {
                    if let Some(index_of_duplicate) = get_index_of_duplicate(path) {
                        closed_paths.push(path[(index_of_duplicate)..path.len() - 1].to_owned());
                        *path_option = None;
                    }
                }
            }
            for i in (0..paths.len()).rev() {
                if paths[i].is_none() {
                    paths.remove(i);
                }
            }
        }
        closed_paths = deduplicate_vecs(closed_paths);
        closed_paths.sort_by_key(|closed_loop| -(closed_loop.len() as isize));

        self.rings = Some(closed_paths);
    }
```

## Set-up
```rust
        let mut paths = vec![];
        let mut closed_paths = vec![];
        for neighbor_index in self.atom_neighbor_indicies(0) {
            paths.push(Some(vec![0, neighbor_index]));
        }
```

The vector `paths` will be a `Vec<Option<Vec<usize>>>` and will hold all paths while we enumerate them.

The vector `closed_paths` will store paths that close back on themselves.

We initialize paths by creating our initial paths starting from our first atom.

## Traversal
```rust
        while !paths.is_empty() {
            // phase 1
            let mut new_paths = vec![];
            for path in paths.iter_mut() {
                if path.is_none() {
                    continue;
                }
                let last_atom_index = *path.as_ref().unwrap().last().unwrap();
                let second_to_last_atom_index = path.as_ref().unwrap().iter().rev().nth(1).unwrap();
                let mut neighbor_indices = self.atom_neighbor_indicies(last_atom_index);
                neighbor_indices
                    .retain(|neighbor_index| neighbor_index != second_to_last_atom_index);
                if neighbor_indices.is_empty() {
                    *path = None;
                } else {
                    for neighbor_index in &neighbor_indices[1..] {
                        let mut new_path = path.clone().unwrap();
                        new_path.push(*neighbor_index);
                        new_paths.push(new_path);
                    }
                    path.as_mut().unwrap().push(neighbor_indices[0]);
                }
            }
            for new_path in new_paths {
                paths.push(Some(new_path));
            }
            // phase 2
            for path_option in paths.iter_mut() {
                if let Some(path) = path_option {
                    if let Some(index_of_duplicate) = get_index_of_duplicate(path) {
                        closed_paths.push(path[(index_of_duplicate)..path.len() - 1].to_owned());
                        *path_option = None;
                    }
                }
            }
            // phase 3
            for i in (0..paths.len()).rev() {
                if paths[i].is_none() {
                    paths.remove(i);
                }
            }
        }
```

In the core of the algorithm, there are three phases.
1. Find neighbors of the terminal atom of the path and extend the path/add new paths as needed.
2. If there is a duplicated atom index in a path, extract the closed loop from the path.
3. Remove empty paths.

## Clean-up
```rust
        closed_paths = deduplicate_vecs(closed_paths);
        closed_paths.sort_by_key(|closed_loop| -(closed_loop.len() as isize));

        self.rings = Some(closed_paths);
```

At the end, we use a custom util to deduplicate our closed paths. For example, a ring will be traversed as [0, 1, 2] and [0, 2, 1]. Because these both encode the same ring, we only need to keep one.

We then sort the rings from longest to shortest before storing them in `self.rings` as an Option.

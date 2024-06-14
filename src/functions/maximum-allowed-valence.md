# Maximum Allowed Valence
For most atoms, the maximum allowed valence is
`4 - abs(4 - (n_valence_electrons - formal_charge))`.

However, the expanded octet presents an annoying group of exceptions. Elements
like P, S, and Cl can form compounds like PF6- or ClO4-.

The current implementation doesn't do a good job handling these cases and I fear
I will need to hard-code a lot of chemical logic. So I have decided not to show
code yet. But I hope to improve the handling of maximum allowed valence soon.

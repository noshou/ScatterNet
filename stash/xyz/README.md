# stash/xyz

Pulled out of `atomic_radii/` — element-symbol-only xyz scanning doesn't
belong in a radius-lookup library, and coordinates need real parsing
(float triples, not throwaway matches) plus a coordinate-system conversion
(cartesian to polar/spherical) before this is useful for anything beyond
what it already did (find which elements are mentioned in a file).

Not wired into the dune build — no `dune` file here on purpose, so it
doesn't need to compile until it's picked back up.

- `Lexer_xyz.mll` — the three line-structured rules (`atom_count_line`,
  `skip_line`, `xyz_atom_line`) that used to live in
  `atomic_radii/Lexer.mll`, plus the regexes only they used.
- `Xyz.ml` — `fold`/`unique_elements`, single-pass streaming over those
  lexer rules.

When this comes back: `xyz_atom_line`'s coordinate field (currently
`[^ '\n']*`, matched and thrown away) needs to become a real float-triple
grammar, and whatever consumes it needs the cartesian -> polar conversion
this was stashed to make room for.

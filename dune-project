(lang dune 3.20)

(generate_opam_files true)

(package
 (name scatternet)
 (allow_empty)
 (synopsis "OCaml scattering toolkit")
 (depends
  (ocaml (>= 5.1))
  dune
  gsl
  owl
  pyml
  sqlite3
  ppx_inline_test))

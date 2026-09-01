# patches

## `owl-1.2-exponpow-args.patch`

**Applies to**: `owl` version `1.2` (the version pinned in this project's opam
switch as of this writing). Check the current pinned version with
`opam list owl`.

**What it fixes**: `owl` 1.2 does not compile as released. Its
`src/owl/stats/owl_stats_dist_exponpow.c` calls `std_gaussian_rvs` with an
argument in two places, but that function is declared and defined to take
none (`extern double std_gaussian_rvs();`). Every other call site in `owl`'s
own codebase generates a standard gaussian draw and scales it afterward
(e.g. `mu + sigma * std_gaussian_rvs()` in `owl_stats_ziggurat.c`) — this one
file was never updated to match. The patch applies that same pattern.

This is an upstream defect, unrelated to this machine or this project. It
should be reported (or has already been reported/fixed) at
https://github.com/owlbarn/owl — check there first; if a newer `owl` release
fixes it, this patch (and the pin below) can be dropped entirely.

**How to apply**: this project currently gets the fix via an `opam pin`,
not via anything wired into the dune build. To reproduce that pin from
scratch (e.g. on a new machine, or after `opam` state was lost):

```sh
git clone --depth 1 --branch 1.2 https://github.com/owlbarn/owl.git /tmp/owl-src
cd /tmp/owl-src
patch -p1 < /path/to/this/repo/patches/owl-1.2-exponpow-args.patch
opam pin add owl /tmp/owl-src -y
```

Verify with `opam list owl` — it should show `pinned to file:///tmp/owl-src`
(or wherever you applied it).

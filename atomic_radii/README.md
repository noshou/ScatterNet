# atomic_radii

`atomic_radii.db` holds three tables: empirical **ionic** radii (charged
species only), a bare-**element** fallback radius for every element in the
periodic table, and a cache of which charge states each element actually has
data for.

## `ionic_radii`

```sql
CREATE TABLE ionic_radii (
    ion    TEXT PRIMARY KEY,  -- lowercase element symbol + |charge| + sign, e.g. 'fe3+', 'cl1-'
    radius REAL NOT NULL      -- empirical ionic radius, picometers
)
```

212 rows — every ion with an empirically measured radius. No modeled,
estimated, or interpolated values, and no ion without one: this table
answers exactly one question (what's this ion's radius) and nothing else.

### Negative radii: `h1+` (-38.0 pm), `n5+` (-10.4 pm), `c4+` (-8.0 pm)

Three rows have a negative `radius`. These are genuine, correctly-transcribed
Shannon (1976) values, not data errors — verified against the published
table. Shannon's scale is anchored to a reference O^2- radius (140 pm), not
an absolute physical size, so a cation with little or no electron cloud of
its own (a bare proton, or a small, highly-stripped cation like C4+/N5+) can
land below that baseline arithmetically. The number is real and correctly
sourced; it just isn't a usable *physical excluded-volume radius* on its own.

**Consumer guidance**: any code that turns `radius` into a physical
volume/sphere must treat `radius <= 0.0` the same as "no entry found" and
fall through to `atomic_radii` instead — never feed a negative radius
directly into a volume formula. Do **not** fix this by editing
`ionic_radii.db` (e.g. clamping to zero or dropping the rows) — that would
misrepresent Shannon's actual published data. The guard belongs in the
consumer, since it's really "these two physical concepts (crystallographic
ionic radius vs. excluded volume) don't agree for these three ions," not
"the data is wrong."

### Sources

| Citation | DOI / identifier | Ions |
|---|---|---|
| Shannon, R.D. (1976). Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. *Acta Cryst.* A32, 751–767. | `10.1107/S0567739476001551` | 210 ions |
| Feldmann, C. (1995). Zur kristallchemischen Ähnlichkeit von Aurid- und Halogenid-Ionen. *Z. Anorg. Allg. Chem.* 621, 1907–1912. | `10.1002/zaac.19956211113` | `au1-` |
| Pauling, L. (1960). *The Nature of the Chemical Bond*, 3rd ed. Cornell University Press. | ISBN 0-8014-0333-2 (pre-DOI monograph) | `h1-` |

## `atomic_radii`

```sql
CREATE TABLE atomic_radii (
    element     TEXT PRIMARY KEY,  -- lowercase element symbol, no charge, e.g. 'fe', 'rn'
    radius      REAL NOT NULL,     -- Angstrom
    radius_type TEXT NOT NULL,     -- 'vdw' | 'metallic' | 'covalent'
    source      TEXT NOT NULL      -- full citation for this row
)
```

118 rows — every element, Z=1 through Z=118, with no charge state involved.
This is the fallback tier for whatever `ionic_radii` doesn't cover: neutral
atoms, noble gases, and anything the charge-aware lookup misses.

`radius_type` is not decoration — van der Waals, metallic, and covalent radii
are three different physical quantities (covalent radii in particular run
30–70% smaller than van der Waals radii for the same element). They are kept
in one table only because each element gets exactly one row from exactly one
source; nothing here blends quantities of different types for the same atom.
A consumer that needs a single physically-comparable scale across the whole
table should be aware `radius_type` varies by element and treat rows
differently accordingly, rather than dropping every `radius` into one
formula uniformly.

Coverage is exactly partitioned — one source per row, no overlaps:

| `radius_type` | Elements | Count |
|---|---|---|
| `vdw` | Main group: H–Ra/Fr/Rn (groups 1, 2, 13–18) | 44 |
| `metallic` | Transition metals Sc–Hg (minus Rf–Cn), lanthanides La–Lu, actinides Ac–Am | 51 |
| `covalent` | Actinides Cm–Lr, transactinides/superheavies Rf–Og | 23 |

### Sources

| Citation | DOI / identifier | `radius_type` | Elements |
|---|---|---|---|
| Mantina, M.; Chamberlin, A.C.; Valero, R.; Cramer, C.J.; Truhlar, D.G. (2009). Consistent van der Waals Radii for the Whole Main Group. *J. Phys. Chem. A* 113, 5806–5812. | `10.1021/jp8111556` | `vdw` | 43 main-group elements, H–Ra/Fr (all except Rn) |
| Runeberg, N.; Pyykkö, P. (1998). Relativistic pseudopotential calculations on Xe2, RnXe, and Rn2: the van der Waals properties of radon. *Int. J. Quantum Chem.* 66, 131–140. | (pre-DOI; see journal record) | `vdw` | `rn` only — used in place of Mantina's own Rn value (2.20 Å) by deliberate choice; Runeberg & Pyykkö's 2.24 Å comes from a direct relativistic coupled-cluster calculation on Rn2 rather than an interpolation, and is kept as the primary value here rather than a footnote |
| Slater, J.C. (1964). Atomic Radii in Crystals. *J. Chem. Phys.* 41, 3199–3204. | `10.1063/1.1725697` | `metallic` | 51 elements: transition metals Sc–Hg (minus Rf–Cn, undiscovered in 1964), lanthanides La–Lu (full), actinides Ac–Am |
| Pyykkö, P.; Atsumi, M. (2009). Molecular Single-Bond Covalent Radii for Elements 1–118. *Chem. Eur. J.* 15, 186–197. | `10.1002/chem.200800987` | `covalent` | 23 elements: actinides Cm–Lr (Slater 1964 has no data for these) and transactinides/superheavies Rf–Og |

**Note on Slater (1964)**: it's a general bonding-radius compilation derived
from ~1200 measured bond lengths (covalent, metallic, and ionic crystals and
molecules combined), not a strict 12-coordinate metallic-radius table. It was
chosen over splitting the transition-metal/lanthanide/actinide block across
several narrower sources because it's the one source that covers nearly the
whole block under one consistent method — consistency across elements
mattered more here than using the single best number for each element
individually.

**Note on the Cm–Lr values** (`covalent`, Pyykkö & Atsumi 2009): the original
PDF is paywalled, so these 8 values were cross-checked against two
independent database reproductions of the paper's data
([KnowledgeDoor](https://www.knowledgedoor.com/2/elements_handbook/pyykko_covalent_radius.html)
and the [Crystallography Open Database](https://radii.crystallography.net/cgi-bin/cov_radii_table.pl)),
which agree with each other exactly. Wikipedia's "Covalent radius" table
gives a different value for curium specifically (172 pm vs. 166 pm here) and
omits berkelium and californium entirely — that discrepancy is unresolved
against the primary source and is noted here rather than silently picked
one way. The `vdw` values above (Mantina 2009) also could not be checked
against the original ACS PDF directly (403/paywalled) and were reconstructed
from a secondary transcription (the QCElemental dataset) that matches the
commonly-cited figures for this well-known table.

## `element_charges`

```sql
CREATE TABLE element_charges (
    element TEXT NOT NULL,
    charge  INTEGER NOT NULL,             -- signed, e.g. 3 for 'fe3+', -1 for 'cl1-'
    ion     TEXT NOT NULL REFERENCES ionic_radii(ion),
    PRIMARY KEY (element, charge)
)
```

212 rows, one per `ionic_radii` row — a cache, not new data. It exists so a
consumer can answer "what charge states does this element actually have data
for" or "which charge state is closest to the one I want" with one indexed
query instead of re-parsing every `ionic_radii.ion` string on every lookup:

```sql
-- nearest available charge state for element 'fe', target charge +5
SELECT ion, charge FROM element_charges
WHERE element = 'fe' ORDER BY ABS(charge - 5) LIMIT 1;
```

This is the intended way to implement "closest match, or fall back to ground
state" behavior: try `ionic_radii` for the exact ion first, and if that
misses, use `element_charges` to find the nearest charge state for that
element before falling all the way back to `atomic_radii`'s bare-element
value. If `element_charges` regenerated from scratch, it must be rebuilt from
`ionic_radii` (it's a pure derived cache — the parsing regex is
`^([a-z]{1,2})(\d*)([+-])$`, magnitude defaults to 1 when no digits are
present, e.g. `na+` → charge `+1`).

## Everything else is code, not data

Ion-vocabulary indexing (`VOCAB`) and form-factor computation (`f0`, `f1`,
`f2`) are **not** duplicated into any static table here — both are handled
live via [xraydb](https://github.com/xraypy/XrayDB) (public domain,
maintained by Matt Newville — DOI
[10.5281/zenodo.16114067](https://doi.org/10.5281/zenodo.16114067)) at
runtime, in `Preprocess/train_data_accessors/vocab.py` and
`Scattering/formfact.py` respectively (original Python pipeline,
`noshou/APS360`):

- **`VOCAB`** (`vocab.py`): built at import time from
  `xraydb.get_xraydb().f0_ions()`, a fixed 211-entry table (98 neutral
  elements + 111 ionic forms + 2 internal aliases). Source: Waasmaier, D. &
  Kirfel, A. (1995). New analytical scattering-factor functions for free
  atoms and ions. *Acta Cryst.* A51, 416–431. DOI:
  [10.1107/S0108767394013292](https://doi.org/10.1107/S0108767394013292).
- **`f0(q)`** (Thomson/non-resonant term, energy-independent, varies with
  momentum transfer `q`) and **`f1(E)`/`f2(E)`** (Chantler anomalous
  corrections, energy-dependent, `q`-independent): both queried live via
  `xraydb.f0`/`f1_chantler`/`f2_chantler` for whatever `q`/`energy` a given
  job needs. Source for the Chantler values: Chantler, C.T. (1995).
  *J. Phys. Chem. Ref. Data* 24, 71–643. DOI:
  [10.1063/1.555974](https://doi.org/10.1063/1.555974); Chantler, C.T.
  (2000). *J. Phys. Chem. Ref. Data* 29, 597–1048. DOI:
  [10.1063/1.1321055](https://doi.org/10.1063/1.1321055).

For the OCaml port in this repo (`ScatterNet/`), the same split applies:
`Xraydb_bindings` calls the real `xraydb` Python package live via `pyml`
(see `ocaml_migration.md` §1b in `APS360/`), rather than re-implementing
`f0`/`f1`/`f2` from a static snapshot — nothing from either source is cached
to disk. `atomic_radii.db` holds the genuinely static, parameter-free
quantities — ionic and atomic radii — and nothing else.

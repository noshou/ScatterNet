# atomic_radii

`atomic_radii.sqlite3` holds three tables: empirical **ionic** radii (charged species only), a bare-**element** fallback radius for every element in the periodic table, and a cache of which charge states each element actually has data for.

## `ionic_radii`

```sql
CREATE TABLE ionic_radii (
    ion    TEXT PRIMARY KEY,  -- lowercase element symbol + |charge| + sign, e.g. 'fe3+', 'cl1-'
    radius REAL NOT NULL      -- empirical ionic radius, picometers
)
```

### Negative radii: `h1+` (-38.0 pm), `n5+` (-10.4 pm), `c4+` (-8.0 pm)

Three rows have a negative `radius`. These are genuine Shannon (1976) values,  Shannon's scale is anchored to a reference O^2- radius (140 pm), not an absolute physical size, so a cation with little or no electron cloud of its own (a bare proton, or a small, highly-stripped cation like C4+/N5+) can land below that baseline arithmetically.  Any code that turns `radius` into a physical  volume/sphere must treat `radius <= 0.0` the same as "no entry found" and fall through to `atomic_radii`.

### Sources


| Citation                                                                                                                                                      | DOI / identifier                       | Ions     |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------- | ---------- |
| Shannon, R.D. (1976). Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides.*Acta Cryst.* A32, 751–767. | `10.1107/S0567739476001551`            | 210 ions |
| Feldmann, C. (1995). Zur kristallchemischen Ähnlichkeit von Aurid- und Halogenid-Ionen.*Z. Anorg. Allg. Chem.* 621, 1907–1912.                              | `10.1002/zaac.19956211113`             | `au1-`   |
| Pauling, L. (1960).*The Nature of the Chemical Bond*, 3rd ed. Cornell University Press.                                                                       | ISBN 0-8014-0333-2 (pre-DOI monograph) | `h1-`    |

## `atomic_radii`

```sql
CREATE TABLE atomic_radii (
    element     TEXT PRIMARY KEY,  -- lowercase element symbol, no charge, e.g. 'fe', 'rn'
    radius      REAL NOT NULL,     -- Angstrom
    radius_type TEXT NOT NULL,     -- 'vdw' | 'metallic' | 'covalent'
    source      TEXT NOT NULL      -- full citation for this row
)
```

118 rows: every ground state element, Z=1 through Z=118.

`radius_type`: van der Waals, metallic, and covalent radii  are three different physical quantities (covalent radii in particular run  30–70% smaller than van der Waals radii for the same element). 


| `radius_type` | Elements                                                                      | Count |
| --------------- | ------------------------------------------------------------------------------- | ------- |
| `vdw`         | Main group: H–Ra/Fr/Rn (groups 1, 2, 13–18)                                 | 44    |
| `metallic`    | Transition metals Sc–Hg (minus Rf–Cn), lanthanides La–Lu, actinides Ac–Am | 51    |
| `covalent`    | Actinides Cm–Lr, transactinides/superheavies Rf–Og                          | 23    |

### Sources


| Citation                                                                                                                                                                                                                 | DOI / identifier              | `radius_type` | Elements                                                                                                                      |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Mantina, M.; Chamberlin, A.C.; Valero, R.; Cramer, C.J.; Truhlar, D.G. (2009). Consistent van der Waals Radii for the Whole Main Group.*J. Phys. Chem. A* 113, 5806–5812.                                               | `10.1021/jp8111556`           | `vdw`         | 43 main-group elements, H–Ra/Fr (all except Rn)                                                                              |
| Runeberg, N.; Pyykkö, P. (1998). Relativistic pseudopotential calculations on Xe2, RnXe, and Rn2: the van der Waals properties of radon.*Int. J. Quantum Chem.* 66, 131–140.                                           | (pre-DOI; see journal record) | `vdw`         | Runeberg & Pyykkö's 2.24 Å comes from a direct relativistic coupled-cluster calculation on Rn2 rather than an interpolation |
| Teatum, E.T.; Gschneidner, K.A. Jr.; Waber, J.T. (1968). Compilation of Calculated Data Useful in Predicting Metallurgical Behavior of the Elements in Binary Alloy Systems.*LA-4003*, Los Alamos Scientific Laboratory. | (LASL report; no DOI)         | `metallic`    | 51 elements: transition metals Sc–Hg (minus Rf–Cn, undiscovered in 1968), lanthanides La–Lu, actinides Ac–Am              |
| Pyykkö, P.; Atsumi, M. (2009). Molecular Single-Bond Covalent Radii for Elements 1–118.*Chem. Eur. J.* 15, 186–197.                                                                                                   | `10.1002/chem.200800987`      | `covalent`    | 23 elements: actinides Cm–Lr and transactinides/superheavies Rf–Og                                                          |

## **`element_charges`**

```sql
CREATE TABLE element_charges (
    element TEXT NOT NULL,
    charge  INTEGER NOT NULL,             -- signed, e.g. 3 for 'fe3+', -1 for 'cl1-'
    ion     TEXT NOT NULL REFERENCES ionic_radii(ion),
    PRIMARY KEY (element, charge)
)
```

212 rows, one per `ionic_radii` row. Exists to answer "what charge states does this element actually have data  for" or "which charge state is closest to the one I want" with one query instead of re-parsing every `ionic_radii.ion` string on every lookup:

```sql
-- nearest available charge state for element 'fe', target charge +5
SELECT ion, charge FROM element_charges
WHERE element = 'fe' ORDER BY ABS(charge - 5) LIMIT 1;
```

1. Try `ionic_radii` for the exact ion first
2. If that misses, use `element_charges` to find the nearest charge state for that  element
3. Fall back to `atomic_radii`'s bare-element value.
4. If `element_charges` regenerated from scratch, it must be rebuilt from `ionic_radii`

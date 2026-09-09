"""Extract data from xraydb.sqlite into one compact db.

Source: xraydb 4.5.8's xraydb.sqlite. Its LICENSE places xraydb.sqlite and
data_sources/ in the public domain via CC0 1.0.
    - Waasmaier & Kirfel (1995) Acta Cryst A51, 416   -> f0 Gaussian coefficients
    - Chantler FFAST (NIST, fine grid)                -> f1/f2 anomalous terms
Energy/f1/f2 are stored as little-endian Float64 BLOBs: compact, and read back
with one reinterpret rather than a JSON parse per element.
"""
import json, sqlite3, struct, sys, os

SRC = sys.argv[1]; DST = sys.argv[2]
src = sqlite3.connect(SRC)
if os.path.exists(DST): os.remove(DST)
dst = sqlite3.connect(DST)

dst.executescript("""
CREATE TABLE waasmaier (
    ion     TEXT PRIMARY KEY,
    element TEXT NOT NULL,
    z       INTEGER NOT NULL,
    c       REAL NOT NULL,
    a1 REAL NOT NULL, a2 REAL NOT NULL, a3 REAL NOT NULL, a4 REAL NOT NULL, a5 REAL NOT NULL,
    b1 REAL NOT NULL, b2 REAL NOT NULL, b3 REAL NOT NULL, b4 REAL NOT NULL, b5 REAL NOT NULL
);
CREATE INDEX idx_waasmaier_element ON waasmaier(element);
CREATE TABLE chantler (
    element TEXT PRIMARY KEY,
    z       INTEGER NOT NULL,
    npts    INTEGER NOT NULL,
    emin    REAL NOT NULL,
    emax    REAL NOT NULL,
    energy  BLOB NOT NULL,
    f1      BLOB NOT NULL,
    f2      BLOB NOT NULL
);
CREATE TABLE provenance (key TEXT PRIMARY KEY, value TEXT NOT NULL);
""")

def blob(xs): return struct.pack("<%dd" % len(xs), *xs)

n_w = 0
for z, el, ion, off, scale, expo in src.execute(
        "SELECT atomic_number, element, ion, offset, scale, exponents FROM Waasmaier"):
    a = json.loads(scale); b = json.loads(expo)
    assert len(a) == 5 and len(b) == 5, (ion, len(a), len(b))
    dst.execute("INSERT INTO waasmaier VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (ion.lower(), el.lower(), z, off, *a, *b))
    n_w += 1

n_c = 0; tot = 0
for el, en, f1, f2 in src.execute("SELECT element, energy, f1, f2 FROM Chantler"):
    e = json.loads(en); y1 = json.loads(f1); y2 = json.loads(f2)
    assert len(e) == len(y1) == len(y2), el
    # Cs has duplicate grid energies upstream, which makes the s=0 spline throw.
    # Drop exact duplicates, keeping the first occurrence, so the grid is strictly
    # increasing for every element.
    keep = [0] + [i for i in range(1, len(e)) if e[i] > e[i-1]]
    if len(keep) != len(e):
        print("  deduped %s: %d -> %d points" % (el, len(e), len(keep)))
    e = [e[i] for i in keep]; y1 = [y1[i] for i in keep]; y2 = [y2[i] for i in keep]
    z = src.execute("SELECT atomic_number FROM elements WHERE element=?", (el,)).fetchone()[0]
    dst.execute("INSERT INTO chantler VALUES (?,?,?,?,?,?,?,?)",
                (el.lower(), z, len(e), e[0], e[-1], blob(e), blob(y1), blob(y2)))
    n_c += 1; tot += len(e)

for k, v in [
    ("f0_source",      "Waasmaier & Kirfel (1995) Acta Cryst A51, 416-431; doi:10.1107/S0108767394013292"),
    ("f0_form",        "f0(s) = c + sum_{i=1..5} a_i*exp(-b_i*s^2), s = q/(4*pi) [1/Ang], valid 0 <= s <= 6"),
    ("anomalous_source","Chantler FFAST (NIST), fine grid; J. Phys. Chem. Ref. Data 24 71 (1995), 29 597 (2000)"),
    ("anomalous_form", "f1 stored as f1_FFAST - Z + f_rel(3/5 CL) + f_NT (xraydb convention); f = f0 + f1 + i*f2"),
    ("extracted_from", "xraydb 4.5.8 xraydb.sqlite"),
    ("license",        "CC0 1.0 - xraydb LICENSE dedicates xraydb.sqlite and data_sources/ to the public domain"),
]:
    dst.execute("INSERT INTO provenance VALUES (?,?)", (k, v))

dst.commit(); dst.execute("VACUUM"); dst.close()
print("waasmaier: %d species, chantler: %d elements / %d grid points" % (n_w, n_c, tot))
print("size: %d bytes" % os.path.getsize(DST))

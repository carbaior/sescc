#!/usr/bin/env python3
"""
Comprueba las listas de Knobel (1917) contra el catálogo de Ulugh Beg de
Verbunt & van Gent (2012) y genera dos versiones filtradas del CSV de entrada
de SESCC.

Uso (desde sescc-main/catalogs):
    python3 check_knobel.py
    python3 check_knobel.py --dat ulughbeg.dat --csv ulughbeg_verbunt.csv \
                            --knobel knobel1917_derived_latitudes.csv
"""
import argparse, csv, sys

ap = argparse.ArgumentParser()
ap.add_argument("--dat", default="ulughbeg.dat")
ap.add_argument("--csv", default="ulughbeg_verbunt.csv")
ap.add_argument("--knobel", default="knobel1917_derived_latitudes.csv")
ap.add_argument("--maxerr", type=float, default=60.0, help="filtro de error posicional (arcmin)")
a = ap.parse_args()

# ---------- 1. Leer ulughbeg.dat (formato CDS de Verbunt & van Gent) ----------
dat = {}
for line in open(a.dat, encoding="utf-8", errors="replace"):
    t = line.split()
    if not t or not t[0].isdigit():
        continue
    U = int(t[0])
    flag = "c" if len(t) > 1 and t[1] == "c" else ""
    # El signo de latitud es el primer token 'B' o 'A' (la constelación lleva '=')
    k = next((j for j in range(2, len(t)) if t[j] in ("B", "A")), None)
    if k is None:                      # p. ej. U961, sin coordenadas
        dat[U] = dict(flag=flag, lat=None, hip=None, err=None)
        continue
    G, M = int(t[k - 2]), int(t[k - 1])
    sign = 1 if t[k] == "B" else -1
    hip = int(t[k + 2])
    nums = [x.rstrip("*") for x in t[k + 3:] if x.rstrip("*")]
    err = float(nums[-1])
    dat[U] = dict(flag=flag, lat=(sign, G, M), hip=hip, err=err)
print(f"ulughbeg.dat: {len(dat)} entradas leídas")

# ---------- 2. Leer la lista de Knobel ----------
kn = list(csv.DictReader(open(a.knobel, encoding="utf-8")))
IV = {int(r["U"]) for r in kn if r["table"] == "IV"}
V = {int(r["U"]) for r in kn if r["table"] == "V"}
print(f"Knobel: Tabla IV = {len(IV)}, Tabla V = {len(V)} estrellas")

# ---------- 3. Comprobación A: latitudes de la Tabla IV ----------
print("\n[A] Latitudes de la Tabla IV (Knobel) frente a ulughbeg.dat")
bad = 0
for r in kn:
    if r["table"] != "IV":
        continue
    U = int(r["U"]); s = r["knobel_lat"].rstrip("?")
    ks = -1 if s[0] == "-" else 1
    kG, kM = int(s[1:].split()[0]), int(s[1:].split()[1])
    d = dat.get(U)
    if d is None or d["lat"] is None:
        print(f"  U{U}: no está en ulughbeg.dat o sin coordenadas"); bad += 1; continue
    sg, G, M = d["lat"]
    same = (G, M) == (kG, kM) and (sg == ks or (G == 0 and M == 0))
    if not same:
        bad += 1
        print(f"  U{U:4d} {r['knobel_label']:10s} Knobel {r['knobel_lat']:8s}  VvG {'+' if sg>0 else '-'}{G:02d} {M:02d}"
              + ("   (dudosa según Knobel)" if r["doubtful"] == "1" else ""))
print(f"  -> {68-bad} de 68 coinciden" if bad else "  -> las 68 coinciden")

# ---------- 4. Comprobación B: Keid y Arcturo ----------
print("\n[B] Estrellas dominantes")
for name, hip in (("Keid", 19849), ("Arcturo", 69673)):
    Us = [U for U, d in dat.items() if d["hip"] == hip]
    for U in Us:
        where = "Tabla IV" if U in IV else "Tabla V" if U in V else "no está en Knobel"
        print(f"  {name} (HIP {hip}) = U{U}: {where}")
    if not Us:
        print(f"  {name} (HIP {hip}): no aparece")

# ---------- 5. Comprobación C: las 27 de al-Sufi ----------
C = {U for U, d in dat.items() if d["flag"] == "c"}
Cin = {U for U in C if dat[U]["err"] is not None and dat[U]["err"] < a.maxerr}
print(f"\n[C] Marcadas 'c' (de al-Sufi): {len(C)}; pasan el filtro de {a.maxerr:g}': {len(Cin)}")

# ---------- 6. Cuántas excluidas pasan el filtro ----------
def passing(S):
    return {U for U in S if U in dat and dat[U]["err"] is not None and dat[U]["err"] < a.maxerr}
print(f"\n[D] Pasan el filtro de {a.maxerr:g}':  Tabla IV {len(passing(IV))},  "
      f"Tabla IV+V {len(passing(IV | V))}")

# ---------- 7. Escribir los dos CSV filtrados ----------
lines = open(a.csv, encoding="utf-8").read().splitlines()
hdr = next(l for l in lines if l.startswith("#") and "seq" in l)
iseq = hdr.lstrip("# ").split(";").index("seq")
for tag, S in (("sinIV", IV), ("sinIVyV", IV | V)):
    out = a.csv.replace(".csv", f"_{tag}.csv")
    n_out = 0
    with open(out, "w", encoding="utf-8") as f:
        for l in lines:
            if l.startswith("#") or not l.strip():
                f.write(l + "\n"); continue
            if int(l.split(";")[iseq]) in S:
                n_out += 1; continue
            f.write(l + "\n")
        f.write(f"# Excluidas {n_out} estrellas de Knobel (1917) {'Tabla IV' if tag=='sinIV' else 'Tablas IV y V'}\n")
    print(f"\nEscrito {out}: {n_out} filas excluidas")

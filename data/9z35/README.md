# Omit models: 9z35 conformations A and B

This README records the atoms omitted from the original `9z35` PDB structures. The missing-atom identities, residue numbers, and coordinates below were verified by comparing `ATOM` and `HETATM` records in each original/omit pair. The conformation percentages are retained from the supplied README.

| Conformation | Original PDB | Omit PDB | Missing-atom CSV |
|---|---|---|---|
| A (35%) | `9z35_confA.pdb` | `9z35_confA_omitOEXsidechains7HOH.pdb` | [`9z35_missing_atoms.csv`](9z35_missing_atoms.csv) |
| B (65%) | `9z35_confB.pdb` | `9z35_confB_omitOEYsidechains7HOH.pdb` | [`9z35_confB_missing_atoms.csv`](9z35_confB_missing_atoms.csv) |

## Missing-atom summary

| Category | Conf A (OEX) | Conf B (OEY) |
|---|---:|---:|
| Protein atoms | 13 | 13 |
| OEC cofactor atoms | 10 | 11 |
| Water oxygen atoms | 7 | 7 |
| **Total missing atoms** | **30** | **31** |

The same 13 protein atom identities and six water oxygen identities are removed from both conformations. Conformation B additionally omits **O6** from the OEY cofactor, so it has 30 missing atoms compared with 29 in conformation A.

## Removed protein atoms (both conformations)

| Chain | Residue | Removed atoms | Count |
|---|---|---|---:|
| A | ASP 170 | OD1, OD2 | 2 |
| A | GLU 189 | OE1, OE2 | 2 |
| A | HIS 332 | NE2 | 1 |
| A | GLU 333 | OE1, OE2 | 2 |
| A | ASP 342 | OD1, OD2 | 2 |
| A | ALA 344 | O, OXT | 2 |
| C | GLU 354 | OE1, OE2 | 2 |

**Note:** These atoms include carboxylate oxygen atoms, HIS A 332 NE2, and ALA A 344 O/OXT. The supplied README heading “coordinating carboxylate oxygens” therefore does not describe every omitted protein atom.

## Removed oxygen-evolving complex (OEC)

| Conformation | Cofactor residue | Removed manganese atoms | Removed calcium atom | Removed oxygen atoms | Count |
|---|---|---|---|---|---:|
| A (35%) | OEX A 602 | MN1, MN2, MN3, MN4 | CA1 | O1, O2, O3, O4, O5 | 10 |
| B (65%) | OEY A 601 | MN1, MN2, MN3, MN4 | CA1 | O1, O2, O3, O4, O5, O6 | 11 |

## Removed metal-coordinating waters (both conformations)

| Chain | Water residue | Removed atom |
|---|---|---|
| A | HOH 705 | O |
| A | HOH 712 | O |
| A | HOH 738 | O |
| A | HOH 776 | O |
| A | HOH 814 | O |
| A | HOH 817 | O |
| A | HOH 856 | O |

All six waters listed in the supplied README—including HOH A 705 and HOH A 776, previously shown in parentheses—are confirmed as removed by the PDB comparison.

## MN1 center coordinates

MN1 Cartesian coordinates extracted from the original PDB atom records (Å):

| Conformation | Cofactor | X (Å) | Y (Å) | Z (Å) |
|---|---|---:|---:|---:|
| A | OEX A 602 MN1 | -34.828 | 39.221 | 366.933 |
| B | OEY A 601 MN1 | -34.830 | 39.328 | 366.977 |

**Difference from supplied README:** The original README listed Conf A MN1 as `(-34.825, 39.317, 366.906)` Å; the actual atom record in `9z35_confA.pdb` is `(-34.828, 39.221, 366.933)` Å. The README Conf B MN1 value `(-34.830, 39.328, 366.977)` Å agrees with `9z35_confB.pdb`. This document uses the PDB-derived coordinates without assuming the reason for the Conf A discrepancy.

## Comparison integrity checks

| Check | Conformation A | Conformation B |
|---|---:|---:|
| Original `ATOM` + `HETATM` count | 53,199 | 53,201 |
| Omit-file `ATOM` + `HETATM` count | 53,170 | 53,171 |
| Removed records | 30 | 31 |
| Additional atom records in omit file | 0 | 0 |
| Modified retained atom records | 0 | 0 |
| `ANISOU` records in either file | 0 | 0 |

Atom matching was based on record type, chain, residue number, insertion code, residue name, atom name, and alternate-location label. Every retained `ATOM`/`HETATM` record is identical in its original and omit file, including the coordinates. No `ANISOU` records were present in these PDB pairs. Non-coordinate records (such as `REMARK` and `CONECT`) are not part of the atom-by-atom comparison.

## CSV contents

Each CSV contains one row per missing atom, retaining the original PDB record type and atom serial number, chain, residue name/number, insertion code, atom name, alternate-location label, element, and the original Cartesian coordinates (`x_A`, `y_A`, `z_A`, in Å).


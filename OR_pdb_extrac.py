import os
import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────
BASE_DIR      = Path(__file__).parent          # dataset/
RECEPTORS_DIR = BASE_DIR / 'receptors'
DOCKED_CSV    = BASE_DIR / 'docked_coords.csv'
OUTPUT_CSV    = BASE_DIR / 'pocket_coords.csv'

POCKET_CUTOFF = 6.0   # Angstroms — standard binding pocket definition

# ── Load ligand docked coords ─────────────────────────────────────────
print("Loading docked_coords.csv...")
df_lig = pd.read_csv(DOCKED_CSV)
print(f"  {df_lig.shape[0]} ligand atom rows, "
      f"{df_lig[['Receptor','Ligand']].drop_duplicates().shape[0]} pairs")

# Keep only heavy atoms (no hydrogens) for distance calculation
df_lig_heavy = df_lig[~df_lig['element'].isin(['HD','HS','H'])]

# ── PDB parser ────────────────────────────────────────────────────────
def parse_pdb(pdb_path):
    """
    Parse a receptor PDB file.
    Returns DataFrame with columns:
        res_id, res_name, chain, atom_name, x, y, z, element
    Heavy atoms only (no H).
    """
    records = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM','HETATM')):
                continue
            atom_name = line[12:16].strip()
            res_name  = line[17:20].strip()
            chain     = line[21].strip()
            res_id    = int(line[22:26].strip())
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            # Element: try column 76-78 first, fall back to first char of atom_name
            element = line[76:78].strip() if len(line) > 76 else ''
            if not element:
                element = atom_name[0]
            # Skip hydrogens
            if element in ('H', 'D'):
                continue
            records.append({
                'res_id'   : res_id,
                'res_name' : res_name,
                'chain'    : chain,
                'atom_name': atom_name,
                'x': x, 'y': y, 'z': z,
                'element'  : element
            })
    return pd.DataFrame(records)

# ── Pre-load all receptor PDB files ──────────────────────────────────
print("Pre-loading receptor PDB files...")
receptor_pdbs = {}
for receptor_dir in sorted(RECEPTORS_DIR.iterdir()):
    if not receptor_dir.is_dir():
        continue
    receptor = receptor_dir.name
    pdb_path = receptor_dir / f'{receptor}.pdb'
    if not pdb_path.exists():
        print(f"  ✗ Missing PDB: {pdb_path}")
        continue
    df_rec = parse_pdb(pdb_path)
    if df_rec.empty:
        print(f"  ✗ Empty PDB parse: {receptor}")
        continue
    receptor_pdbs[receptor] = df_rec

print(f"  Loaded {len(receptor_pdbs)} receptor PDB files")

# ── Extract pocket residues per pair ──────────────────────────────────
print("Extracting pocket residues...")
all_pocket_rows = []
pairs = df_lig_heavy[['Receptor','Ligand']].drop_duplicates()
n_pairs = len(pairs)

for idx, (_, pair) in enumerate(pairs.iterrows()):
    receptor = pair['Receptor']
    ligand   = pair['Ligand']

    if idx % 500 == 0:
        print(f"  Processing pair {idx}/{n_pairs}: {receptor}_{ligand}")

    # Get ligand atom coords for this pair
    lig_atoms = df_lig_heavy[
        (df_lig_heavy['Receptor'] == receptor) &
        (df_lig_heavy['Ligand']   == ligand)
    ][['x','y','z']].values

    if len(lig_atoms) == 0:
        continue

    # Get receptor atoms
    if receptor not in receptor_pdbs:
        continue
    df_rec = receptor_pdbs[receptor]
    rec_coords = df_rec[['x','y','z']].values

    # ── Distance filter: find pocket residues ─────────────────────────
    # For each receptor atom, compute min distance to any ligand atom.
    # Vectorized: broadcast [n_rec, 3] vs [n_lig, 3] → [n_rec, n_lig]
    # Then take min over ligand atoms → [n_rec] distances
    diffs     = rec_coords[:, np.newaxis, :] - lig_atoms[np.newaxis, :, :]  # [n_rec, n_lig, 3]
    distances = np.sqrt((diffs ** 2).sum(axis=2))                            # [n_rec, n_lig]
    min_dists = distances.min(axis=1)                                        # [n_rec]

    pocket_mask = min_dists <= POCKET_CUTOFF
    df_pocket   = df_rec[pocket_mask].copy()

    if df_pocket.empty:
        print(f"  ✗ Empty pocket: {receptor}_{ligand}")
        continue

    df_pocket.insert(0, 'Ligand',   ligand)
    df_pocket.insert(0, 'Receptor', receptor)
    all_pocket_rows.append(df_pocket)

# ── Write output ──────────────────────────────────────────────────────
df_pocket_all = pd.concat(all_pocket_rows, ignore_index=True)
df_pocket_all.to_csv(OUTPUT_CSV, index=False)

print(f"\nDone.")
print(f"  Total pocket atom rows : {df_pocket_all.shape[0]}")
print(f"  Pairs covered          : {df_pocket_all[['Receptor','Ligand']].drop_duplicates().shape[0]}")
print(f"  Output written to      : {OUTPUT_CSV}")

# ── Sanity check: pocket size distribution ────────────────────────────
pocket_sizes = (
    df_pocket_all
    .groupby(['Receptor','Ligand'])['res_id']
    .nunique()
    .describe()
)
print(f"\nPocket residue count per pair:")
print(pocket_sizes.round(1))

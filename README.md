#!/usr/bin/env python3
"""
Extract pose-1 atom coordinates + metadata from all docked .pdbqt files
under outputs/<Receptor>_<Ligand>/docking_1/output.pdbqt

Reads receptor-ligand pairs from input.csv (columns: Receptor, Ligand)
Writes one row per atom to docked_coords.csv
"""

import os
import re
import pandas as pd

# ── Config ──────────────────────────────────────────────────────────
INPUT_CSV   = 'input.csv'
OUTPUTS_DIR = 'outputs'
DOCKING_RUN = 'docking_1'
OUT_CSV     = 'docked_coords.csv'


def extract_top_pose(pdbqt_path):
    """Return only the ATOM/HETATM lines belonging to MODEL 1."""
    with open(pdbqt_path) as f:
        lines = f.readlines()

    pose_lines = []
    in_model_1 = False
    affinity = None

    for line in lines:
        if line.startswith('MODEL'):
            parts = line.split()
            model_num = parts[1] if len(parts) > 1 else None
            in_model_1 = (model_num == '1')
            continue
        if line.startswith('ENDMDL'):
            if in_model_1:
                break
            continue
        if in_model_1:
            if line.startswith('REMARK VINA RESULT'):
                m = re.search(r'REMARK VINA RESULT:\s*([-\d.]+)', line)
                if m:
                    affinity = float(m.group(1))
            if line.startswith(('ATOM', 'HETATM')):
                pose_lines.append(line)

    return pose_lines, affinity


def parse_atom_line(line):
    """Parse fixed-width PDBQT ATOM/HETATM line into fields."""
    return {
        'atom_serial' : int(line[6:11]),
        'atom_name'   : line[12:16].strip(),
        'res_name'    : line[17:20].strip(),
        'x'           : float(line[30:38]),
        'y'           : float(line[38:46]),
        'z'           : float(line[46:54]),
        'element'     : line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0],
    }


def find_pdbqt_path(receptor, ligand):
    """Build expected path, handling space-in-ligand-name folder variants."""
    ligand_variants = [
        ligand.strip(),
        ligand.strip().replace(' ', '_'),
    ]
    for lig in ligand_variants:
        folder = f"{receptor}_{lig}"
        path = os.path.join(OUTPUTS_DIR, folder, DOCKING_RUN, 'output.pdbqt')
        if os.path.exists(path):
            return path
    return None


def main():
    df_input = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df_input)} receptor-ligand pairs from {INPUT_CSV}")

    records = []
    missing = []
    failed = []

    for _, row in df_input.iterrows():
        receptor = str(row['Receptor']).strip()
        ligand   = str(row['Ligand']).strip()

        path = find_pdbqt_path(receptor, ligand)
        if path is None:
            missing.append((receptor, ligand))
            continue

        try:
            pose_lines, affinity = extract_top_pose(path)
            if not pose_lines:
                failed.append((receptor, ligand, 'no MODEL 1 found'))
                continue

            for line in pose_lines:
                atom = parse_atom_line(line)
                records.append({
                    'Receptor'   : receptor,
                    'Ligand'     : ligand,
                    'affinity'   : affinity,
                    **atom
                })
        except Exception as e:
            failed.append((receptor, ligand, str(e)))

    df_out = pd.DataFrame(records)
    df_out.to_csv(OUT_CSV, index=False)

    print(f"\n✓ Wrote {len(df_out)} atom rows ({df_out[['Receptor','Ligand']].drop_duplicates().shape[0]} pairs) to {OUT_CSV}")

    if missing:
        print(f"\n✗ Missing files for {len(missing)} pairs:")
        for r, l in missing[:10]:
            print(f"    {r}_{l}")
        if len(missing) > 10:
            print(f"    ... and {len(missing) - 10} more")

    if failed:
        print(f"\n⚠ Parse failures for {len(failed)} pairs:")
        for r, l, err in failed[:10]:
            print(f"    {r}_{l}: {err}")


if __name__ == '__main__':
    main()

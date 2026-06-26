# replace the check_files function in check_coverage.py

def check_files(names, directory, extensions=(".pdb", ".pdbqt")):
    present, missing = [], []
    for name in names:
        subdir = directory / name          # e.g. receptors/OR10A2/
        # look for receptors/OR10A2/OR10A2.pdb  OR  OR10A2.pdbqt
        found = any(
            (subdir / f"{name}{ext}").exists() for ext in extensions
        )
        (present if found else missing).append(name)
    return present, missing


































((test1_env) ) dellwsl@DESKTOP-M2LSKB8:~/test1/Automated-Molecular-Docking/dataset$ cat code.py
import os
import time
import requests
import subprocess
import glob
import shutil
import pandas as pd
import numpy as np
from Bio import Entrez
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.action_chains import ActionChains

# ====================== GLOBAL CONFIGURATION ======================
EXCEL_FILE = "/home/dellwsl/test1/Automated-Molecular-Docking/dataset/input.xlsx"
VINA_PATH = "/usr/bin/vina"
OBABEL_PATH = "/usr/bin/obabel"
CHROME_DRIVER_PATH = "/usr/bin/chromedriver"

EMAIL = "atirathpal@gmail.com"
Entrez.email = EMAIL

# Docking Parameters
DOCKING_SEEDS = [100, 200, 300]
NUM_POSES = 10
EXHAUSTIVENESS = 16
ENERGY_RANGE = 4

# Global Directories
BASE_DIR = os.getcwd()
RECEPTOR_WH = os.path.join(BASE_DIR, "receptors")
LIGAND_WH = os.path.join(BASE_DIR, "ligands")
OUTPUT_ROOT = os.path.join(BASE_DIR, "outputs")

for d in [RECEPTOR_WH, LIGAND_WH, OUTPUT_ROOT]:
    os.makedirs(d, exist_ok=True)
# =================================================================

def assess_model(receptor_name, pdb_file):
    """Swiss-Model Structure Assessment (Step_02 style navigation)."""
    print(f"🧪 Running Web-based Validation for {receptor_name}...")

    receptor_dir = os.path.join(RECEPTOR_WH, receptor_name)
    os.makedirs(receptor_dir, exist_ok=True)

    chrome_options = Options()
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": receptor_dir,
        "download.prompt_for_download": False,
        "safebrowsing.enabled": True
    })

    driver = webdriver.Chrome(service=Service(CHROME_DRIVER_PATH), options=chrome_options)

    try:
        driver.get("https://swissmodel.expasy.org/")
        WebDriverWait(driver, 10).until(lambda d: d.execute_script("return document.readyState") == "complete")

        tools_dropdown = WebDriverWait(driver, 30).until(
            EC.element_to_be_clickable((By.XPATH, "//a[contains(@class, 'dropdown-toggle') and text()='Tools']"))
        )
        ActionChains(driver).move_to_element(tools_dropdown).click().perform()

        try:
            structure_assessment_link = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.XPATH, "//a[@href='/assess']"))
            )
        except Exception:
            structure_assessment_link = WebDriverWait(driver, 10).until(
                EC.element_to_be_clickable((By.XPATH, "//a[contains(text(), 'Structure Assessment')]"))
            )

        structure_assessment_link.click()

        file_input = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.XPATH, "//input[@type='file']"))
        )
        file_input.send_keys(os.path.abspath(pdb_file))

        WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.XPATH, f"//*[contains(text(), '{os.path.basename(pdb_file)}')]"))
        )
        print(f"📤 Uploaded: {os.path.basename(pdb_file)}")

        sequence_input = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.NAME, "sequenceSearch"))
        )
        molprobity_input = WebDriverWait(driver, 30).until(
            EC.presence_of_element_located((By.NAME, "molprobity"))
        )

        if not sequence_input.is_selected():
            driver.execute_script("arguments[0].click();", sequence_input)
        if not molprobity_input.is_selected():
            driver.execute_script("arguments[0].click();", molprobity_input)

        start_button = WebDriverWait(driver, 30).until(
            EC.element_to_be_clickable((By.CSS_SELECTOR, "button.css-110xwq5"))
        )
        driver.execute_script("arguments[0].click();", start_button)
        print("🧪 Validation started...")

        download_btn = WebDriverWait(driver, 120).until(
            EC.element_to_be_clickable((By.XPATH, "//button[@title='Save as PNG']"))
        )
        download_btn.click()
        time.sleep(5)

        validation_dir = os.path.join(receptor_dir, f"{receptor_name}_validation")
        os.makedirs(validation_dir, exist_ok=True)

        png_files = glob.glob(os.path.join(receptor_dir, "*.png"))
        if not png_files:
            print("⚠️ PNG not found after assessment.")
        else:
            latest_png = max(png_files, key=os.path.getctime)
            shutil.move(latest_png, os.path.join(validation_dir, f"{receptor_name}_Ramachandran_plot.png"))

        molprobity_table = WebDriverWait(driver, 60).until(
            EC.presence_of_element_located((By.ID, "molprobityResults"))
        )
        soup = BeautifulSoup(molprobity_table.get_attribute("outerHTML"), "html.parser")

        validation_data = {"Metric": [], "Value": [], "Outliers": []}
        for row in soup.find_all("tr"):
            cells = row.find_all("td")
            if len(cells) >= 2:
                metric = cells[0].text.strip()
                value = cells[1].text.strip()
                outliers = metric if row.get("class") == ["outlierContainer"] else ""
                validation_data["Metric"].append(metric)
                validation_data["Value"].append(value)
                validation_data["Outliers"].append(outliers)

        df_val = pd.DataFrame(validation_data)
        df_val.to_csv(os.path.join(validation_dir, f"{receptor_name}_Analysis.csv"), index=False)
        print(f"✅ Validation Complete for {receptor_name}")
        return True

    except Exception as e:
        print(f"⚠️ Assessment failed for {receptor_name}: {e}")
        return False

    finally:
        driver.quit()

def analyze_interactions(receptor_pdbqt, ligand_output_pdbqt):
    """
    Calculates H-bonds and Hydrophobic contacts based on user-defined theory:
    H-Bond: N/O (Rec) to N/O (Lig) < 3.5 A
    Hydrophobic: C (Rec) to C (Lig) < 4.5 A
    """
    def get_atoms_from_lines(lines):
        atoms = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                # Element is at columns 77-78 in PDBQT
                elem = line[76:78].strip()
                coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                atoms.append({'elem': elem, 'coords': coords})
        return atoms

    # 1. Load Receptor Atoms once
    with open(receptor_pdbqt, 'r') as f:
        rec_atoms = get_atoms_from_lines(f.readlines())

    results = []
    current_pose_lines = []

    # 2. Process Ligand Poses
    with open(ligand_output_pdbqt, 'r') as f:
        for line in f:
            if line.startswith("MODEL"):
                current_pose_lines = []
            elif line.startswith(("ATOM", "HETATM")):
                current_pose_lines.append(line)
            elif line.startswith("ENDMDL"):
                lig_atoms = get_atoms_from_lines(current_pose_lines)
                hb, hp = 0, 0

                # --- HYDROGEN BOND LOGIC ---
                # Rec Protein (N/O) to Ligand (N/O) < 3.5 A
                rec_polars = [a for a in rec_atoms if a['elem'] in ['N', 'O']]
                lig_polars = [a for a in lig_atoms if a['elem'] in ['N', 'O']]
                for r in rec_polars:
                    for l in lig_polars:
                        if np.linalg.norm(r['coords'] - l['coords']) < 3.5:
                            hb += 1

                # --- HYDROPHOBIC BOND LOGIC ---
                # Rec Protein (C) to Ligand (C) < 4.5 A (User Theory)
                rec_carbons = [a for a in rec_atoms if a['elem'] == 'C']
                lig_carbons = [a for a in lig_atoms if a['elem'] == 'C']
                for rc in rec_carbons:
                    for lc in lig_carbons:
                        dist = np.linalg.norm(rc['coords'] - lc['coords'])
                        # Updated to strict "Less than 4.5" per your instruction
                        if dist < 4.5:
                            hp += 1

                results.append({'h_bonds': hb, 'hydrophobic': hp})
    return results

def download_fasta(gene_name):
    try:
        print(f"🧬 Searching NCBI: {gene_name}...")
        term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism] AND srcdb_refseq[PROP]"
        handle = Entrez.esearch(db="protein", term=term, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]: return None
        p_id = record["IdList"][0]
        fetch_handle = Entrez.efetch(db="protein", id=p_id, rettype="fasta", retmode="text")
        fasta_data = fetch_handle.read()
        fetch_handle.close()
        return fasta_data
    except Exception as e:
        print(f"❌ NCBI Error: {e}"); return None

def run_swiss_model(fasta_data, receptor_name):
    receptor_dir = os.path.join(RECEPTOR_WH, receptor_name)
    os.makedirs(receptor_dir, exist_ok=True)
    pdb_path = os.path.join(receptor_dir, f"{receptor_name}.pdb")

    driver = webdriver.Chrome(service=Service(CHROME_DRIVER_PATH))
    try:
        driver.get("https://swissmodel.expasy.org/interactive")
        textarea = WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.ID, "id_target")))
        textarea.send_keys(fasta_data)
        build_btn = WebDriverWait(driver, 20).until(EC.element_to_be_clickable((By.ID, "buildButton")))
        driver.execute_script("arguments[0].scrollIntoView(true);", build_btn)
        build_btn.click()
        print(f"⏳ Modeling {receptor_name}...")
        pdb_link = WebDriverWait(driver, 900).until(
            EC.presence_of_element_located((By.XPATH, "//a[contains(@href, '.pdb') and contains(text(), 'PDB format')]"))
        )
        time.sleep(60)
        pdb_url = pdb_link.get_attribute("href")
        res = requests.get(pdb_url)
        with open(pdb_path, "wb") as f:
            f.write(res.content)
        return pdb_path
    except Exception as e:
        print(f"❌ SwissModel Error: {e}")
        return None
    finally:
        driver.quit()

def prepare_receptor(raw_pdb, receptor_name):
    receptor_dir = os.path.join(RECEPTOR_WH, receptor_name)
    os.makedirs(receptor_dir, exist_ok=True)

    chain_a_pdb = os.path.join(receptor_dir, f"{receptor_name}_chainA.pdb")
    final_pdbqt = os.path.join(receptor_dir, f"{receptor_name}.pdbqt")

    with open(raw_pdb, 'r') as infile, open(chain_a_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM', 'HETATM')) and line[21] == 'A':
                outfile.write(line)
            elif line.startswith('END'):
                outfile.write(line)

    subprocess.run(
        [OBABEL_PATH, chain_a_pdb, "-O", final_pdbqt, "-xr", "-h", "--partialcharge", "gasteiger"],
        check=True,
        capture_output=True
    )

    if os.path.exists(chain_a_pdb):
        os.remove(chain_a_pdb)

    return final_pdbqt

def prepare_ligand(ligand_name):
    final_pdbqt = os.path.join(LIGAND_WH, f"{ligand_name}.pdbqt")
    temp_sdf = os.path.join(LIGAND_WH, f"{ligand_name}.sdf")
    try:
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        cid_res = requests.get(f"{base_url}/compound/name/{ligand_name}/cids/TXT")
        cid = cid_res.text.strip().split()[0]
        sdf_res = requests.get(f"{base_url}/compound/cid/{cid}/record/SDF/?record_type=3d")
        with open(temp_sdf, "w") as f: f.write(sdf_res.text)
        subprocess.run([OBABEL_PATH, temp_sdf, "-O", final_pdbqt, "-h", "--partialcharge", "gasteiger"], check=True, capture_output=True)
        if os.path.exists(temp_sdf): os.remove(temp_sdf)
        return final_pdbqt
    except Exception as e:
        print(f"❌ Ligand Error: {e}"); return None

def generate_config(receptor_pdbqt, ligand_pdbqt, output_dir):
    coords = []
    with open(receptor_pdbqt, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = (coords.max(axis=0) - coords.min(axis=0)) + 10.0
    config_path = os.path.join(output_dir, "config.txt")
    with open(config_path, "w") as f:
        f.write(f"receptor = {os.path.abspath(receptor_pdbqt)}\n")
        f.write(f"ligand = {os.path.abspath(ligand_pdbqt)}\n\n")
        f.write(f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n\n")
        f.write(f"size_x = {size[0]:.3f}\nsize_y = {size[1]:.3f}\nsize_z = {size[2]:.3f}\n\n")
        f.write(f"exhaustiveness = {EXHAUSTIVENESS}\nnum_modes = {NUM_POSES}\nenergy_range = {ENERGY_RANGE}\n")
    return config_path

# ======================== MAIN LOOP ========================
def main():
    df = pd.read_excel(EXCEL_FILE)

    failed_pairs = []
    for index, row in df.iterrows():
        rec_name = str(row['Receptor']).strip()
        lig_name = str(row['Ligand']).strip()

        try:
            print(f"\n▶️ Row {index+1}: {rec_name} + {lig_name}")

            pair_folder = os.path.join(OUTPUT_ROOT, f"{rec_name}_{lig_name}")
            os.makedirs(pair_folder, exist_ok=True)

        # 1. RECEPTOR & ASSESSMENT
            rec_dir = os.path.join(RECEPTOR_WH, rec_name)
            rec_pdb = os.path.join(rec_dir, f"{rec_name}.pdb")
            rec_pdbqt = os.path.join(rec_dir, f"{rec_name}.pdbqt")

            if not os.path.exists(rec_pdbqt):
                if not os.path.exists(rec_pdb):
                    fasta = download_fasta(rec_name)
                    if fasta:
                        rec_pdb = run_swiss_model(fasta, rec_name)

                if rec_pdb and os.path.exists(rec_pdb):
                    assess_model(rec_name, rec_pdb)
                    rec_pdbqt = prepare_receptor(rec_pdb, rec_name)

        # 2. LIGAND PREP
            lig_pdbqt = os.path.join(LIGAND_WH, f"{lig_name}.pdbqt")
            if not os.path.exists(lig_pdbqt):
                lig_pdbqt = prepare_ligand(lig_name)

        # 3. DOCKING
            if os.path.exists(rec_pdbqt) and os.path.exists(lig_pdbqt):
                config_file = generate_config(rec_pdbqt, lig_pdbqt, pair_folder)

                for run_num, seed in enumerate(DOCKING_SEEDS, 1):
                    run_dir = os.path.join(pair_folder, f"docking_{run_num}")
                    os.makedirs(run_dir, exist_ok=True)
                    out_pdbqt = os.path.join(run_dir, "output.pdbqt")
                    log_file = os.path.join(run_dir, "log.txt")
                    analysis_csv = os.path.join(run_dir, "pose_analysis.csv")

                    if os.path.exists(analysis_csv):
                        continue

                    print(f"  🚀 Run {run_num} (Seed: {seed})...")
                    result = subprocess.run(
                        [
                            VINA_PATH,
                            "--config", config_file,
                            "--out", out_pdbqt,
                            "--seed", str(seed)
                        ],
                        capture_output=True,
                        text=True
                    )
                    with open(log_file, "w", encoding="utf-8") as f:
                        f.write(result.stdout)

                    if os.path.exists(out_pdbqt):
                        results = analyze_interactions(rec_pdbqt, out_pdbqt)
                        affinities = []
                        if os.path.exists(log_file):
                            with open(log_file, 'r') as lf:
                                for line in lf:
                                    p = line.split()
                                    if len(p) > 3 and p[0].isdigit():
                                        affinities.append(p[1])

                        with open(analysis_csv, 'w') as f:
                            f.write("Pose,Affinity_kcal,H_Bonds,Hydrophobic\n")
                            for i, res in enumerate(results):
                                aff = affinities[i] if i < len(affinities) else "N/A"
                                f.write(f"{i+1},{aff},{res['h_bonds']},{res['hydrophobic']}\n")

                print(f"✅ Row {index+1} Finished.")
            else:
                print(f"🛑 Error preparing inputs for Row {index+1}.")

                failed_pairs.append({
                    "Row": index + 1,
                    "Receptor": rec_name,
                    "Ligand": lig_name,
                    "Error": "Input preparation failed"
                })

                continue

        except Exception as e:
            print(f"❌ Row {index+1} Failed.")
            print(f"❌ {rec_name} + {lig_name}")
            print(f"❌ Error: {e}")

            failed_pairs.append({
                "Row": index + 1,
                "Receptor": rec_name,
                "Ligand": lig_name,
                "Error": str(e)
            })

            continue

    if failed_pairs:
        failed_df = pd.DataFrame(failed_pairs)

        failed_df.to_excel(
            "Failed_Pairs.xlsx",
            index=False
        )

        print(
            f"⚠️ {len(failed_pairs)} failed pairs saved "
            f"to Failed_Pairs.xlsx"
        )

def collect_validation_data():
    validation_rows = []

    if not os.path.exists(RECEPTOR_WH):
        return validation_rows

    receptors = [d for d in os.listdir(RECEPTOR_WH) if os.path.isdir(os.path.join(RECEPTOR_WH, d))]

    for receptor in receptors:
        val_csv = os.path.join(RECEPTOR_WH, receptor, f"{receptor}_validation", f"{receptor}_Analysis.csv")
        if not os.path.exists(val_csv):
            continue

        try:
            vdf = pd.read_csv(val_csv)
            for _, row in vdf.iterrows():
                validation_rows.append({
                    "Receptor": receptor,
                    "Metric": row.get("Metric", ""),
                    "Value": row.get("Value", ""),
                    "Outliers": row.get("Outliers", "")
                })
        except Exception as e:
            print(f"⚠️ Could not read validation CSV for {receptor}: {e}")

    return validation_rows


def master_organizer():
    print("Starting master organizer...")
    final_excel_name = "Final_Docking_Results.xlsx"

    tab1_data = []  # Executive Summary
    tab2_data = []  # Triplicate Consistency
    tab3_data = []  # Raw Deep Dive

    if not os.path.exists(OUTPUT_ROOT):
        print("❌ Output folder not found!")
        return

    pairs = [d for d in os.listdir(OUTPUT_ROOT) if os.path.isdir(os.path.join(OUTPUT_ROOT, d))]
    print("Pairs found:", len(pairs))
    print(pairs[:5])

    for pair in pairs:
        pair_path = os.path.join(OUTPUT_ROOT, pair)
        parts = pair.split('_')
        rec_name = parts[0]
        lig_name = "_".join(parts[1:])

        docking_runs = [
            d for d in os.listdir(pair_path)
            if d.startswith("docking_") and os.path.isdir(os.path.join(pair_path, d))
        ]
        docking_runs.sort()

        pair_pose1_affinities = {}
        pair_pose1_hbonds = []
        pair_pose1_hydrophobic = []

        for run in docking_runs:
            run_id = run.split('_')[1]
            csv_path = os.path.join(pair_path, run, "pose_analysis.csv")

            if not os.path.exists(csv_path):
                continue

            try:
                run_df = pd.read_csv(csv_path)
            except Exception as e:
                print(f"Bad CSV: {csv_path}")
                print(f"Error: {e}")
                continue

            # TAB 3: all poses
            for _, row in run_df.iterrows():
                tab3_data.append({
                    "Receptor": rec_name,
                    "Ligand": lig_name,
                    "Run_ID": run_id,
                    "Pose_Number": row.get("Pose", ""),
                    "Affinity": row.get("Affinity_kcal", ""),
                    "H_Bonds": row.get("H_Bonds", ""),
                    "Hydrophobic": row.get("Hydrophobic", "")
                })

            # TAB 2 + summary stats from Pose 1
            pose1 = run_df[run_df["Pose"] == 1] if "Pose" in run_df.columns else pd.DataFrame()
            if not pose1.empty:
                aff = pd.to_numeric(pose1.iloc[0].get("Affinity_kcal", np.nan), errors="coerce")
                hb = pd.to_numeric(pose1.iloc[0].get("H_Bonds", np.nan), errors="coerce")
                hp = pd.to_numeric(pose1.iloc[0].get("Hydrophobic", np.nan), errors="coerce")

                tab2_data.append({
                    "Receptor": rec_name,
                    "Ligand": lig_name,
                    "Run_ID": run_id,
                    "Pose1_Affinity": aff,
                    "Pose1_HBonds": hb,
                    "Pose1_Hydrophobic": hp
                })

                if not pd.isna(aff):
                    pair_pose1_affinities[f"Run{run_id}_Affinity"] = aff
                if not pd.isna(hb):
                    pair_pose1_hbonds.append(hb)
                if not pd.isna(hp):
                    pair_pose1_hydrophobic.append(hp)

        # TAB 1 summary
        if pair_pose1_affinities:
            aff_values = list(pair_pose1_affinities.values())
            summary_row = {
                "Receptor": rec_name,
                "Ligand": lig_name,
            }
            summary_row.update(pair_pose1_affinities)
            summary_row["Average_Affinity"] = float(np.mean(aff_values))
            summary_row["Std_Deviation"] = float(np.std(aff_values))
            summary_row["Avg_HBonds_Pose1"] = float(np.mean(pair_pose1_hbonds)) if pair_pose1_hbonds else np.nan
            summary_row["Avg_Hydrophobic_Pose1"] = float(np.mean(pair_pose1_hydrophobic)) if pair_pose1_hydrophobic else np.nan
            summary_row["Best_Overall_Affinity"] = float(np.min(aff_values))
            tab1_data.append(summary_row)

    validation_rows = collect_validation_data()

    print(f"📊 Compiling {len(tab1_data)} receptor-ligand pairs into Excel...")
    with pd.ExcelWriter(final_excel_name, engine="xlsxwriter") as writer:
        df_tab1 = pd.DataFrame(tab1_data)
        df_tab2 = pd.DataFrame(tab2_data)
        df_tab3 = pd.DataFrame(tab3_data)
        df_tab4 = pd.DataFrame(validation_rows)

        df_tab1.to_excel(writer, sheet_name="Executive_Summary", index=False)
        df_tab2.to_excel(writer, sheet_name="Run_Consistency", index=False)
        df_tab3.to_excel(writer, sheet_name="Raw_Deep_Dive", index=False)
        df_tab4.to_excel(writer, sheet_name="Validation_All", index=False)

        if not df_tab1.empty:
            workbook = writer.book
            worksheet = writer.sheets["Executive_Summary"]
            header_format = workbook.add_format({"bold": True, "bg_color": "#D7E4BC", "border": 1})
            for col_num, value in enumerate(df_tab1.columns.values):
                worksheet.write(0, col_num, value, header_format)

    print(f"✅ Master Excel created: {os.path.abspath(final_excel_name)}")


if __name__ == "__main__":
    main()
    master_organizer()

import os
import json
import csv

from snakemake.utils import min_version
min_version("9.3.0")

# Load configuration
configfile: "lib/config.yaml"

# working directory
TARGET_DIR = config["target"]
LOG_DIR = os.path.join(TARGET_DIR, "logs")
TMP_DIR = os.path.join(TARGET_DIR, "temp")


INPUT_DIR = os.path.expanduser("~/Documents/boltz2_effector_data/702_boltz2_structures")
OUTPUT_DIR = os.path.join(config["scratch"], config["workdir"])


onstart:
    # Clean up old log directories using full path from config
    shell("rm -rf {LOG_DIR}/*")

# Find input dirs within one or two dir levels.
def find_input_dirs(dir):
    matches = []
    for e in os.listdir(dir):
        tmp_path = os.path.join(dir, e)
        if not os.path.isdir(tmp_path):
            continue
        elif e == config["input"]:
            matches.append(tmp_path)
        else:
            for sub_e in os.listdir(tmp_path):
                sub_tmp_path = os.path.join(tmp_path, sub_e)
                if os.path.isdir(sub_tmp_path) and sub_e == config["input"]:
                    matches.append(sub_tmp_path)
    return matches
DIRS = find_input_dirs(TARGET_DIR):

# Find all the structures within all the input dirs.
def find_samples(dir_list):
    samples = []
    for dir in dir_list:
        files = os.listdir(dir)
        names = [f.split(".")[0] for f in files if f.endswith(".pdb")]
        tmp = [os.path.join(dir, n) for n in names]
        samples.extend(tmp)
    return samples
SAMPLES = find_samples(DIRS)

rule all:
    input: f"{OUTPUT_DIR}/summary.csv"

JSON_KEYS = ["confidence_score","ptm","iptm","complex_plddt","complex_iplddt","complex_pde","complex_ipde"]
TSV_KEYS = ["ipSAE", "pDockQ", "pDockQ2", "LIS"]

rule merge_data:
    input:
        json = expand(f"{INPUT_DIR}/{{sample}}_confidence.json", sample=SAMPLES),
        tsv = expand(f"{TMP_DIR}/{{sample}}_10_10.txt", sample=SAMPLES)
    output: csv = f"{OUTPUT_DIR}/summary.csv"
    run:
        rows = []
        for sample, jPath, tPath in zip(SAMPLES, input.json, input.tsv):
            # Load in data from Boltz-2 confidence JSON - we'll use it later.
            with open(jPath) as f:
                data = json.load(f)
            # Find the correct row in the ipSAE TXT file and extract the data.
            tData = {}
            with open(tPath) as f:
                # Need to skip the first line of ipSAE output files as it's blank.
                next(f)
                reader = csv.DictReader(f, delimiter=" ", skipinitialspace = True)
                for row in reader:
                    if row["Type"] == "max":
                        for k in TSV_KEYS:
                            tData[k] = row.get(k, None)
                        break
            # Now combining the Boltz-2 confidence data with the ipSAE data.
            row = {"sample": sample}
            for k in JSON_KEYS:
                row[k] = data.get(k, None)
            row.update(tData)
            rows.append(row)
        fieldnames = ["sample"] + JSON_KEYS + TSV_KEYS
        with open(output.csv, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


rule run_ispsae:
    input:
        cif = f"{TMP_DIR}/{{sample}}.cif",
        npz = f"{INPUT_DIR}/{{sample}}_pae.npz"
    output: f"{TMP_DIR}/{{sample}}_10_10.txt"
    params:
        extra_args = ""
    threads: 1
    resources: mem_mb = "4000"
    shell: "python src/ipsae.py {input.npz} {input.cif} 10 10"

rule convert_pdb:
    input: f"{INPUT_DIR}/{{sample}}.pdb"
    output: f"{TMP_DIR}/{{sample}}.cif"
    threads: 1
    resources: mem_mb = 2000
    shell:
        """
        python -c "from Bio.PDB import PDBParser, MMCIFIO; parser = PDBParser(QUIET=True); structure = parser.get_structure('structure', '{input}'); io = MMCIFIO(); io.set_structure(structure); io.save('{output}')"
        """

rule download_ipSAE:
    output: "src/ipsae.py"
    params:
        url = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/refs/heads/main/ipsae.py"
    shell: "curl -L {params.url} -o {output}"

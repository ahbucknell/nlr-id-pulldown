import os

from snakemake.utils import min_version
min_version("9.3.0")

# Load configuration
configfile: "lib/config.yaml"

# working directory
INPUT_DIR = os.path.expanduser("~/testing_data")
OUTPUT_DIR = os.path.join(config["scratch"], config["workdir"])
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")
TMP_DIR = os.path.join(OUTPUT_DIR, "temp")

onstart:
    # Clean up old log directories using full path from config
    shell("rm -rf {LOG_DIR}/*")

def find_samples(input_dir):
    files = os.listdir(input_dir)
    samples = [f.split(".")[0] for f in files if f.endswith(".pdb")]
    return samples
SAMPLES = find_samples(INPUT_DIR)

rule all:
    input:
        expand(f"{TMP_DIR}/{{sample}}_10_10.txt",sample=SAMPLES)

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
    input: pdb = f"{INPUT_DIR}/{{sample}}.pdb"
    output: cif = f"{TMP_DIR}/{{sample}}.cif"
    threads: 1
    resources: mem_mb = 2000
    shell:
        """
        python -c "from Bio.PDB import PDBParser, MMCIFIO; parser = PDBParser(QUIET=True); structure = parser.get_structure('structure', '{input.pdb}'); io = MMCIFIO(); io.set_structure(structure); io.save('{output.cif}')"
        """

rule download_ipSAE:
    output: "src/ipsae.py"
    params:
        url = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/refs/heads/main/ipsae.py"
    shell: "curl -L {params.url} -o {output}"

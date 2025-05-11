import os
import subprocess
import sys, yaml, requests
from Bio import AlignIO
from Bio.PDB import PDBList, PDBParser, PDBIO, Select, Structure, Model, Chain

UNIPROT_BASE_URL="https://rest.uniprot.org/uniprotkb/search"
PDBE_MAPPINGS_BASE_URL="https://www.ebi.ac.uk/pdbe/api/mappings/interpro"
HEADERS={"accept": "application/json"}

class ChainSelect(Select):
    def __init__(self, keep_chain):
        self.keep_chain = keep_chain
    def accept_chain(self, chain):
        return (chain.id == self.keep_chain)

# Function to create a new directory, after checking for it's existence
def create_new_directory(directory):
    if os.path.isdir(directory):
        subprocess.run(f"rm -rf {directory}", shell = True)
        subprocess.run(f"mkdir -p {directory}", shell = True)
    else:
        subprocess.run(f"mkdir -p {directory}", shell = True)

def build_query(domain, limit):
    return {
        "query": f"(xref:interpro-{domain}) AND reviewed:true AND (database:pdb)",
        "fields": [
            "accession",
            "xref_pdb"
        ],
        "sort": "accession desc",
        "size": f"{limit}"
    }

def do_request(domain, limit):
    response = requests.get(
        UNIPROT_BASE_URL,
        headers=HEADERS,
        params=build_query(domain, limit),
    )

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    return response.json()

def get_ids(domain, limit):
    data = do_request(domain, limit)

    # ids = []
    pdb_ids = []
    for uniprot in data["results"]:
        # ids.append(uniprot["primaryAccession"])
        if (len(uniprot["uniProtKBCrossReferences"]) != 0 and
            uniprot["uniProtKBCrossReferences"][0]["database"] == "PDB"):
            pdb_ids.append(uniprot["uniProtKBCrossReferences"][0]["id"])

    return pdb_ids

def get_chains(domain, pdb_id):
    resp = requests.get(f"{PDBE_MAPPINGS_BASE_URL}/{pdb_id.lower()}")
    if not resp.ok:
        return set()
    chains = []
    for found_domain, info in resp.json().get(pdb_id.lower(), {}).get("InterPro", {}).items():
        if found_domain != domain:
            continue
        for mapping in info.get("mappings", []):
            chains.append(mapping["chain_id"])
    return set(chains)

# def download_pdb_files(ids, directory):
#     create_new_directory(directory)
#     for id in ids:
#         subprocess.run(f"wget  https://files.rcsb.org/view/{id}.pdb -O {directory}/{id}.pdb", shell = True)

def download_pdb_files(
        pdbl,
        parser,
        io,
        pdb_ids,
        domain,
        output_dir,
        pdb_filter,
    ):
    create_new_directory(output_dir)
    for pdb_id in sorted(pdb_ids):
        chains = get_chains(domain, pdb_id)
        if len(chains) == 0:
            continue

        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=f"{output_dir}/pdb", file_format="pdb")

        # Skip all NMR structures
        is_nmr = False
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("EXPDTA") and "NMR" in line.upper():
                    is_nmr = True
                    break

        if is_nmr:
            continue

        # Filter the PDB file by resolution
        if filter_by_resolution(pdb_file, pdb_filter["resolution_threshold"]):
            continue

        struct = parser.get_structure(pdb_id, pdb_file)

        # Get the first model only
        model = next(struct.get_models())

        for chain_id in chains:
            if chain_id not in model:
                continue

            out_name = os.path.join(f"{output_dir}/pdb", f"{pdb_id}_{chain_id}.pdb")

            # Filter the chain by its length
            chain = model[chain_id]
            length = sum(1 for _ in chain)
            if length > pdb_filter["max_residues_count"]:
                continue

            new_structure = Structure.Structure(pdb_id)
            new_model = Model.Model(0)
            new_chain = Chain.Chain(chain_id)

            for residue in chain:
                new_chain.add(residue.copy())  # ensure clean copy

            new_model.add(new_chain)
            new_structure.add(new_model)

            io.set_structure(new_structure)
            io.save(out_name)

def filter_by_resolution(pdb_file, resolution_threshold):
    resolution_line = subprocess.run(
        f"grep \"^REMARK   2 RESOLUTION.\" {pdb_file}",
        capture_output=True,
        shell=True,
    )

    line = resolution_line.stdout.strip()
    if line:
        parts = line.split()
        if len(parts) > 3:
            try:
                resolution = float(parts[3])
                if resolution > resolution_threshold:
                    return True
            except ValueError:
                return False
    else:
        return False

def align(output_dir):
    create_new_directory(f"{output_dir}/alignment")

    with open(f"{output_dir}/alignment/list.txt", 'w') as list:
        list.write("> .\n")
        for file in os.listdir(f"{output_dir}/pdb"):
            if file.endswith(".pdb"):
                list.write(f"+{output_dir}/pdb/{file}\n")

    subprocess.run(f"mustang -f {output_dir}/alignment/list.txt -o {output_dir}/alignment/alignment -F fasta",
                   shell = True)

def fasta_to_stockholm(output_dir):
    # Read the FASTA alignment
    alignment = AlignIO.read(f"{output_dir}/alignment/alignment.afasta", "fasta")

    # Save in Stockholm format
    AlignIO.write(alignment, f"{output_dir}/alignment/alignment.sto", "stockholm")

def build_model(output_dir, domain):
    # convert to Stockholm format
    fasta_to_stockholm(output_dir)

    create_new_directory(f"{output_dir}/model")

    # build a model
    subprocess.run(
        f"hmmbuild {output_dir}/model/{domain}.hmm {output_dir}/alignment/alignment.sto",
        shell=True,
    )

def annotate_domain(domain):
    pass

if __name__ == "__main__":
    # Init biopython objects
    pdbl = PDBList()
    parser = PDBParser()
    io = PDBIO()

    # Read config
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Get PDB identifiers from Uniprot
    pdb_ids = get_ids(
        config['domain']['interpro_id'],
        config['max_amount_of_proteins'],
    )

    # Download PDB files for chains that contain the domain
    download_pdb_files(
        pdbl,
        parser,
        io,
        pdb_ids,
        config['domain']['interpro_id'],
        config['output_dir'],
        config['pdb_filter'],
    )

    # MSA
    align(config['output_dir'])

    # Build profile HMM
    build_model(config['output_dir'], config['domain']['interpro_id'])



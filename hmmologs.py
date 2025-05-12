import os
import subprocess
import random
import sys, yaml, requests

from Bio import AlignIO
from Bio.PDB import PDBList, PDBParser, PDBIO, Structure, Model, Chain


UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
PDBE_MAPPINGS_BASE_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/interpro"
EBI_FASTA_BASE_URL = "https://www.ebi.ac.uk/pdbe/entry/pdb"
UNIPROT_HEADERS = {"accept": "application/json"}


def create_new_directory(directory):
    """Create a new directory.

    Before create, it checks if the directory already exists and remove it.
    """
    if os.path.isdir(directory):
        subprocess.run(f"rm -rf {directory}", shell=True)
        subprocess.run(f"mkdir -p {directory}", shell=True)
    else:
        subprocess.run(f"mkdir -p {directory}", shell=True)


def build_query(domain, limit):
    """Build a query for UniProt API.
    """
    return {
        "query": f"(xref:interpro-{domain}) AND reviewed:true AND (database:pdb)",
        "fields": [
            "accession",
            "xref_pdb"
        ],
        "sort": "accession desc",
        "size": f"{limit}"
    }


def do_request(url, headers=None, params=None):
    """
    Do HTTP GET request.
    :param url: URL
    :param headers: Dictionary of HTTP Headers
    :param params: Dictionary of query parameters
    :return: Response object
    """
    response = requests.get(
        url,
        headers=headers,
        params=params,
    )

    if not response.ok:
        response.raise_for_status()
        sys.exit()

    return response


def get_ids(domain, limit):
    """
    Obtain PDB IDs that contain the domain.
    :param domain: the domain for a search query
    :param limit: the limit for a search query
    :return: the list of PDB IDs
    """
    data = do_request(
        UNIPROT_BASE_URL,
        UNIPROT_HEADERS,
        build_query(domain, limit),
    ).json()

    pdb_ids = []
    for uniprot in data["results"]:
        if (len(uniprot["uniProtKBCrossReferences"]) != 0 and
                uniprot["uniProtKBCrossReferences"][0]["database"] == "PDB"):
            pdb_ids.append(uniprot["uniProtKBCrossReferences"][0]["id"])

    return pdb_ids


def get_chains(domain, pdb_id):
    """
    Obtain a set of chains that contain the domain.
    :param domain: the domain for a search query
    :param pdb_id: the PDB ID for a search query
    :return: the set of chains of the PDB file that contain this domain
    """
    try:
       resp = do_request(f"{PDBE_MAPPINGS_BASE_URL}/{pdb_id.lower()}")
       chains = []
       for found_domain, info in resp.json().get(pdb_id.lower(), {}).get("InterPro", {}).items():
           if found_domain != domain:
               continue
           for mapping in info.get("mappings", []):
               chains.append(mapping["chain_id"])
       return set(chains)
    except Exception as e:
        return set()


def filter_by_resolution(pdb_file, resolution_threshold):
    """
    Filter PDB file by its resolution.
    :param pdb_file: the PDB file to filter
    :param resolution_threshold: the threshold for the resolution
    :return: Boolean flag that says if PDB should be discarded
    """
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


def download_pdb_files(
        pdbl,
        parser,
        io,
        pdb_ids,
        domain,
        output_dir,
        pdb_filter,
):
    """
    Download chains in the PDB format and filter them according to
    constraints set in the config.
    :param pdbl:
    :param parser:
    :param io:
    :param pdb_ids: the IDs of the PDB files to download
    :param domain: the domain of the interest
    :param output_dir: the output directory
    :param pdb_filter: the dictionary with constraints
    :return:
    """
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


def align(output_dir):
    """
    Align the structures using MUSTANG tool (Multiple Structural Alignment).
    :param output_dir: the output directory
    :return:
    """
    create_new_directory(f"{output_dir}/alignment")

    with open(f"{output_dir}/alignment/list.txt", 'w') as list:
        list.write("> .\n")
        for file in os.listdir(f"{output_dir}/pdb"):
            if file.endswith(".pdb"):
                list.write(f"+{output_dir}/pdb/{file}\n")

    subprocess.run(f"mustang -f {output_dir}/alignment/list.txt -o {output_dir}/alignment/alignment -F fasta",
                   shell=True)


def fetch_fasta_by_pdb_id(pdb_id):
    """
    Fetch FASTA file by PDB ID.
    :param pdb_id: ID of a PDB file
    :return:
    """
    try:
        return do_request(f"{EBI_FASTA_BASE_URL}/{pdb_id.lower()}/fasta").text
    except Exception as e:
        raise Exception(f"Failed to fetch FASTA for {pdb_id}")


def download_validation_set(output_dir, ids, fraction, domain):
    """
    Download a set of PDB structures that are going to be used for the validation of
    the model. For the validation we use the same structures that we used to extract chains,
    but only the subset of them (set validation_set_fraction in the config).
    :param output_dir: the output directory
    :param ids: ids of PDB files to download
    :param fraction: the percentage of the PDB files to download
    :param domain: the domain of the interest
    :return:
    """
    create_new_directory(f"{output_dir}/validation")
    ids = random.sample(ids, int(len(ids) * fraction))
    with open(f"{output_dir}/validation/{domain}.fa", 'w') as file:
        for id in ids:
            file.write(fetch_fasta_by_pdb_id(id))


def fasta_to_stockholm(output_dir):
    """
    Converts FASTA file to the Stockholm alignment file
    :param output_dir: the output directory
    :return:
    """
    # Read the FASTA alignment
    alignment = AlignIO.read(f"{output_dir}/alignment/alignment.afasta", "fasta")

    # Save in Stockholm format
    AlignIO.write(alignment, f"{output_dir}/alignment/alignment.sto", "stockholm")


def build_model(output_dir, domain):
    """
    Build a profile HMM using hmmbuild routine of the HHMER.
    :param output_dir: the output directory
    :param domain: the domain of the interest
    :return:
    """
    # convert to Stockholm format
    fasta_to_stockholm(output_dir)

    create_new_directory(f"{output_dir}/model")

    # build a model
    subprocess.run(
        f"hmmbuild {output_dir}/model/{domain}.hmm {output_dir}/alignment/alignment.sto",
        shell=True,
    )


def validate(output_dir, domain):
    """
    Validate the model using hmmsearch routine of the HHMER.
    :param output_dir: the output directory
    :param domain: the domain of the interest
    :return:
    """
    subprocess.run(
        f"hmmsearch {output_dir}/model/{domain}.hmm {output_dir}/validation/{domain}.fa > {output_dir}/validation/{domain}.output",
        shell=True,
    )


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

    # Download PDB files for validation
    download_validation_set(
        config['output_dir'],
        pdb_ids,
        config['validation_set_fraction'],
        config['domain']['interpro_id'],
    )

    # Validate
    validate(
        config['output_dir'],
        config['domain']['interpro_id'],
    )

    # TODO: download random set of the proteins
    # TODO: validate the model on the random set
    # TODO: draw a confusion matrix.

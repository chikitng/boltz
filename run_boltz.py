import textwrap
import os
import json
import csv
from rdkit import Chem

def generate_protein_ligand_yaml(sequence: str, smiles: str) -> str:
    """
    Generates a YAML configuration string for a protein-ligand pair,
    mimicking the structure of the provided 'prot.yaml' file.

    The structure includes the protein sequence, the ligand SMILES string,
    and defines an affinity property binder.

    Args:
        sequence: The full amino acid sequence of the protein.
        smiles: The SMILES string representation of the ligand molecule.

    Returns:
        A string containing the formatted YAML content.
    """

    # Use a triple-quoted f-string for easy multi-line YAML generation.
    # Note: We use {sequence} and {smiles} to insert the user-provided data.
    yaml_content = textwrap.dedent(f"""\
        version: 1  # Optional, defaults to 1
        sequences:
          - protein:
              id: A
              sequence: {sequence}
          - ligand:
              id: B
              smiles: {smiles}
        properties:
            - affinity:
                binder: B
        """)

    return yaml_content.strip()

# --- Example Usage ---
if __name__ == "__main__":
    # 1. Define the user inputs (using a placeholder sequence and a new SMILES string)
    sequence = "GPLGSPEFMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"

    input_sdf="ligand.sdf"

    molecules = []
    with Chem.SDMolSupplier(input_sdf) as SDF:
        # for each mol
        for mol in SDF:
            if mol is None:
                continue
            molecules.append(mol)


    for m in molecules:
        charge = m.GetPropsAsDict()['charge']
        exp_affinity = m.GetPropsAsDict()['IC50_NM']
        smiles = Chem.MolToSmiles(m)
        header = m.GetProp("_Name")
    
        # 2. Generate the YAML content
        new_yaml_data = generate_protein_ligand_yaml(
            sequence=sequence,
            smiles=smiles
        )

        os.makedirs("input", exist_ok=True)
        os.makedirs("output", exist_ok=True)
        
        with open("input/prot.yaml", "w") as f:
            f.write(new_yaml_data)

        os.system("boltz predict --use_msa_server --output_format=pdb --use_potentials --diffusion_samples=1 input/")

        json_file = open('boltz_results_input/predictions/prot/affinity_prot.json')
        json_str = json_file.read()
        affinity = json.loads(json_str).get("affinity_pred_value")
        affinity_prob = json.loads(json_str).get("affinity_probability_binary")

        affinity = 10**(affinity) # convert to IC50 in uM
        #gnina_score = 1e6 * 10**(-1*gnina_score) # convert to IC50 in uM

        fields=[[smiles, round(affinity,3), round(affinity_prob,3), round(exp_affinity,3)]]

        with open('output/results.csv', 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(fields)
    

        if(os.path.isdir('boltz_results_input')):
            cmd="mv boltz_results_input output/{0}".format(header)
            os.system(cmd)

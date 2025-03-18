# Quick script to figure out for each leadlike if there are any fragment structures it was inspired by using RDKit.

from rdkit import Chem
from rdkit.Chem.rdMolAlign import CalcRMS
import json


def read_molecules(path_to_sdf):
    molecules = {}
    for molecule in Chem.SDMolSupplier(
        path_to_sdf,
        sanitize=False,  # some of these fragments may be wonky
    ):
        molecule = Chem.RemoveHs(molecule, sanitize=False)
        Chem.RemoveStereochemistry(molecule)
        molecules[molecule.GetProp("_Name").replace("_bound", "")] = molecule
    return molecules


def check_overlapping_fragment_2d(molecule, fragments_dict):
    overlap = [
        frag_name
        for frag_name, frag_mol in fragments_dict.items()
        if molecule.HasSubstructMatch(frag_mol)
    ]
    return overlap


def check_rmsd_with_frag(molecule, fragment):
    # checks how much spacial overlap there is between the fragment and the fragment substructure of
    # the lead molecule.

    # in the lead molecule, remove all non-matching atom indices
    matching_atom_idcs = molecule.GetSubstructMatch(fragment)
    non_matching_atom_idcs = [
        atom.GetIdx()
        for atom in molecule.GetAtoms()
        if not atom.GetIdx() in matching_atom_idcs
    ]
    res = Chem.RWMol(molecule)
    res.BeginBatchEdit()
    for aid in non_matching_atom_idcs:
        res.RemoveAtom(aid)
    res.CommitBatchEdit()

    if not Chem.MolToSmiles(res) == Chem.MolToSmiles(fragment):
        return "NaN"

    return CalcRMS(res, fragment)


fragments = read_molecules("../fragments/structures/fragments.sdf")
leadlikes = read_molecules("../leadlikes/leadlikes.sdf")

# now for each leadlike, check if there is substructural overlap with any fragments.
overlap_dict = {}
for lead_name, lead_mol in leadlikes.items():
    # construct the dict for json porting
    overlap_dict[lead_name] = {
        "SMILES": Chem.MolToSmiles(lead_mol),
        "MATCHED_FRAGMENTS": [],
    }

    matched_fragments = check_overlapping_fragment_2d(lead_mol, fragments)
    for frag in matched_fragments:
        match_rmsd = check_rmsd_with_frag(lead_mol, fragments[frag])
        overlap_dict[lead_name]["MATCHED_FRAGMENTS"].append(
            {
                "FRAGMENT_NAME": frag,
                "FRAGMENT_SMILES": Chem.MolToSmiles(fragments[frag]),
                "FRAGMENT_LEAD_OVERLAP_RMSD": match_rmsd,
            }
        )
with open("matched_fragments_to_leadlikes.json", "w") as fp:
    json.dump(overlap_dict, fp, indent=4)

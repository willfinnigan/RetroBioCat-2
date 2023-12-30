from rdkit.Chem import AllChem, rdmolops


def rdkit_smile(smile, warning=False, inchi=False):
    try:
        mol = AllChem.MolFromSmiles(smile)
        if inchi == True:
            inchi_mol = AllChem.inchi.MolToInchi(mol)
            mol = AllChem.inchi.MolFromInchi(inchi_mol)
        smile = AllChem.MolToSmiles(mol)
    except:
        if warning == True:
            print('Warning couldnt convert smiles to rdkit - ' + str(smile))
        return None
    return smile


def get_inchl(smiles):
    mol = AllChem.MolFromSmiles(smiles)
    inchi = AllChem.inchi.MolToInchi(mol)
    return inchi



def split_smi(smi):
    mol = AllChem.MolFromSmiles(smi)
    splitMols = rdmolops.GetMolFrags(mol, asMols=True)
    split_list = []
    for mol in splitMols:
        p_smile = AllChem.MolToSmiles(mol)
        split_list.append(p_smile)
    return split_list


def disable_rdkit_logging():
    """
    Disables RDKit whiny logging.
    """
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')
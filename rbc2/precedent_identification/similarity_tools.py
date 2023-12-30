from typing import List, Tuple, Dict

import pandas as pd
from rdkit import DataStructs, Chem
from rdkit.DataStructs import ExplicitBitVect

def get_single_fp(smi: str) -> ExplicitBitVect:
    mol = Chem.MolFromSmiles(smi)
    return Chem.RDKFingerprint(mol)

def get_fingerprints(smis_to_convert: List[str]) -> Tuple[List[ExplicitBitVect], List[str]]:
    mols = []
    converted_smis = []
    for smi in smis_to_convert:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
                converted_smis.append(smi)
        except:
            pass

    fps = [Chem.RDKFingerprint(mol) for mol in mols]
    return fps, converted_smis

def bulk_similarity(target_fp: ExplicitBitVect, smis_to_compare: List[str], fp_gen_method: callable) -> dict[str: float]:
    fps, smis = fp_gen_method(smis_to_compare)
    sims = DataStructs.BulkTanimotoSimilarity(target_fp, fps)
    return {smi: sim for smi, sim in zip(smis, sims)}

def make_fp_df(df, smi_col):
    """ Given a dataframe of precedent reactions, create a fingerprint dataframe"""

    unique_smis = df[smi_col].unique()
    fps, smis = get_fingerprints(unique_smis)
    fp_df = pd.DataFrame({'smiles': smis, 'fp': fps})
    fp_df.set_index('smiles', inplace=True)
    return fp_df


def get_fingerprints_from_fpdf(smis_to_convert: List[str], fpdf: pd.DataFrame) -> Tuple[List[ExplicitBitVect], List[str]]:
    slice_of_fpdf = fpdf[fpdf.index.isin(smis_to_convert)]
    converted_smis = list(slice_of_fpdf.index)
    fps = list(slice_of_fpdf['fp'])
    return fps, converted_smis


if __name__ == '__main__':
    target_smi = 'CCCC=O'
    compare_with = ['CCCO', 'CCCCN', 'CCC(C)CC(=O)O']
    target_fp = get_single_fp(target_smi)
    sims = bulk_similarity(target_fp, compare_with, get_fingerprints)


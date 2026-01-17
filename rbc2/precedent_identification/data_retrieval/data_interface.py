from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Tuple
import pandas as pd



def get_rxn_smi_vectorized(df, product_column, substrate_columns):
    """
    Generate reaction SMILES strings using a vectorized approach.
    Optimized for cases where all columns exist but may contain null values.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing substrate and product columns
    product_column : str
        Name of the column containing product SMILES
    substrate_columns : list
        List of column names containing substrate SMILES

    Returns:
    --------
    pandas.Series
        Series containing reaction SMILES strings
    """

    if df.shape[0] == 0:
        return

    # Create DataFrame with just substrate columns (all columns exist)
    substrates_df = df[substrate_columns].copy()

    # Convert all columns to string type and handle nulls consistently
    for col in substrate_columns:
        # Convert to string with empty string for nulls
        substrates_df[col] = substrates_df[col].astype(str)
        # Replace various null representations with empty string
        substrates_df[col] = substrates_df[col].replace(['nan', 'None', 'none'], '')

    # Join non-empty substrates with dots for each row
    def join_non_empty(row):
        valid_substrates = [s for s in row if s != '']
        return '.'.join(valid_substrates) if valid_substrates else ''

    joined_substrates = substrates_df.apply(join_non_empty, axis=1)

    # Convert products to string and handle nulls
    products = df[product_column].fillna('').astype(str)
    products = products.replace(['nan', 'None', 'none'], '')

    # Create the reaction SMILES with proper format - without using .loc or .name
    # Simply concatenate the two Series directly with conditional logic
    rxn_smiles = pd.Series('', index=df.index)  # Start with empty strings

    # For rows with substrates, create the full reaction SMILES
    mask = joined_substrates != ''
    rxn_smiles[mask] = joined_substrates[mask] + '>>' + products[mask]

    # For rows without substrates but with products, create rxn_smi as ">>product"
    mask_no_substrates = (joined_substrates == '') & (products != '')
    rxn_smiles[mask_no_substrates] = '>>' + products[mask_no_substrates]

    return rxn_smiles


class PrecedentDataQuery(ABC):

    @abstractmethod
    def __init__(self):
        self.product_column: str = ''
        self.substrate_columns: List[str] = []
        self.enzyme_column: str = ''
        self.id_column: str = ''

    # def get_rxn_smi(self, row):
    #
    #     substrates = []
    #     for col in self.substrate_columns:
    #         if col not in list(row.index):  # if column not in row (which is a Series), skip
    #             continue
    #
    #         if row[col] is not None and row[col] != '':  # if column is not empty, add to substrates
    #             substrates.append(row[col])
    #     substrates = [s for s in substrates if pd.isna(s)==False]
    #
    #     rxn_smi = f"{'.'.join(substrates)}>>{row[self.product_column]}"  # compose rxn_smi string
    #
    #     return rxn_smi


    @abstractmethod
    def query_data(self, **key_atts) -> pd.DataFrame():
        """ Takes a number of keyword arguments, specific for each data source"""
        pass

    @abstractmethod
    def get_fps(self, smis: List[str]) -> Tuple[List, List[str]]:
        """ Takes a list of smiles, and returns a list of fingerprints and a list of smiles that were converted to fingerprints """
        pass



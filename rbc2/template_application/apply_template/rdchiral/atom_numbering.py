
def check_atom_numbering(smarts_str):
    # rchiral doesn't like it if atoms are numbered in the products, but they dont appear in the reactants.
    # therefore, remove numbering in these cases
    if type(smarts_str) != str:
        return smarts_str
    if ">>" not in smarts_str:
        return smarts_str

    reactants, products = smarts_str.split(">>")

    atom_numbers_products = []
    atom_numbers_reactants = []
    for i in range(len(products)-1):
        if products[i] == ':':
            char = products[i:i+2]
            atom_numbers_products.append(char)
    for i in range(len(reactants)-1):
        if reactants[i] == ':':
            char = reactants[i:i+2]
            atom_numbers_reactants.append(char)

    for char in atom_numbers_products:
        if char not in atom_numbers_reactants:
            products = products.replace(char, '')

    new_smarts = reactants + '>>' + products
    return new_smarts


if __name__ == "__main__":
    smarts = "[#6:1][#6H1:2](=[O:3])>>[#6:1][#6H0:2](=[O:3])[OH:4]"
    new_smarts = check_atom_numbering(smarts)
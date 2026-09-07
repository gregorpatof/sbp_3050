from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
import pickle


fingerprint_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def get_mwt(smiles):
    """ Calcul du poids moléculaire à partir de la représentation SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    mwt = Descriptors.MolWt(mol)
    return mwt


def get_tpsa(smiles):
    """ Calcul de la surface polaire topologique (TPSA) à partir de la représentation SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    tpsa = Descriptors.TPSA(mol)
    return tpsa


def get_logp(smiles):
    """ Calcul du coefficient de partition octanol:eau (LogP) à partir de la représentation SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    logp = Descriptors.MolLogP(mol) # utilise la méthode de Crippen
    return logp


def get_n_hbond_donors(smiles):
    """ Calcul du nombre de donneurs de ponts H à partir de la représentation SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    n_hbd = mol.NumHDonors(mol)
    return n_hbd


def mpo_deux_points(valeur, p1, p2):
    """ Calcul d'un score MPO avec deux points d'inflexion, avec plateau à 1 de -inf à p1."""
    valeur, p1, p2 = float(valeur), float(p1), float(p2) # conversion à des nombres à virgules
    if p2 <= p1:
        raise ValueError("p2 doit être supérieur à p1")
    score_mpo = None # initialisation de la variable
    if valeur <= p1:
        score_mpo = 1.0
    elif valeur >= p2:
        score_mpo = 0.0
    else:
        score_mpo = (p2 - valeur) / (p2 - p1)
    return score_mpo


def mpo_quatre_points(valeur, p1, p2, p3, p4):
    """ Calcul d'un score MPO avec quatre points d'inflexion, avec plateau à 0 de -inf à p1,
        plateur à 1 de p2 à p3, et plateau à 0 de p4 à +inf.
    """
    # conversion à des nombres à virgule
    valeur, p1, p2, p3, p4 = float(valeur), float(p1), float(p2), float(p3), float(p4)
    if not (p1 <= p2 <= p3 <= p4):
        raise ValueError("Les quatre points doivent être tels que p1 <= p2 <= p3 <= p4")
    score_mpo = None
    if valeur <= p1 or valeur >= p4:
        return 0.0
    if p2 <= valeur <= p3:
        return 1.0
    if valeur < p2:
        return 1.0 if p1 == p2 else (valeur - p1) / (p2 - p1)
    return 1.0 if p3 == p4 else (p4 - valeur) / (p4 - p3)


def build_fingerprints(smi_file, output_file):
    """ Calcul des empreintes moléculaires pour un fichier .smi complet, et sauvegarde en format .pickle"""
    fingerprints_list = []
    names_list = []
    smiles_list = []
    with open(smi_file) as f:
        lines = f.readlines()
    for line in lines:
        smiles, name = line.split()
        mol = Chem.MolFromSmiles(smiles)
        names_list.append(name)
        smiles_list.append(smiles)
        fingerprints_list.append(fingerprint_generator.GetFingerprint(mol))
    data_dictionary = {'names': names_list,
                       'fingerprints': fingerprints_list,
                       'smiles': smiles_list}
    with open(output_file, 'wb') as f:
        pickle.dump(data_dictionary, f)


def find_closest(query_smiles, pickled_data_path):
    """ Trouve la molécule la plus proche du query_smiles dans un ensemble de molécules aux empreintes moléculaires
        précalculées, et retour le coefficient de Tanimoto (Tc), le nom et la structure de cette molécule.
    """
    with open(pickled_data_path, 'rb') as f:
        data = pickle.load(f)
    mol = Chem.MolFromSmiles(query_smiles)
    query_fp = fingerprint_generator.GetFingerprint(mol)

    # Calucl du Tc entre la molécule "query" et chaque molécule de la base de données
    tcs_list = DataStructs.BulkTanimotoSimilarity(query_fp, data['fingerprints'])

    max_tc = 0 # initialisation à 0
    closest_index = None

    for i, tc in enumerate(tcs_list): # i prendra la valeur de chaque index, tc prendra chaque valeur de la liste
        if tc > max_tc:
            max_tc = tc
            closest_index = i

    closest_name = data['names'][closest_index]
    closest_smiles = data['smiles'][closest_index]
    return [max_tc, closest_name, closest_smiles]

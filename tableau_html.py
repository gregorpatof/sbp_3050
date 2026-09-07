"""
tableau_html.py — Production de tableaux HTML contenant des structures moléculaires.

Ce fichier ne contient rien que vous ayez à comprendre ou à modifier pour la vignette.
Il sert uniquement à produire un tableau propre, lisible dans un navigateur, facilement
copiable par capture d'écran pour être inséré dans Word.

Utilisation typique :

    from tableau_html import tableau_html

    tableau_html(mes_donnees,
                 colonnes_smiles=('smiles1', 'smiles2'),
                 fichier='tableau_medicaments.html',
                 titre='Médicaments')
"""

import base64

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# Facteur de suréchantillonnage : les structures sont dessinées à une résolution
# supérieure à leur taille d'affichage, ce qui les garde nettes à l'écran et à
# l'impression.
_FACTEUR = 2

_STYLE = """
<style>
  body   { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif;
           margin: 30px; color: #222; }
  h1     { font-size: 20px; font-weight: 600; margin-bottom: 16px; }
  table  { border-collapse: collapse; }
  th     { background: #f0f0f0; border: 1px solid #bbb; padding: 8px 12px;
           text-align: left; font-size: 14px; white-space: nowrap; }
  td     { border: 1px solid #ddd; padding: 8px 12px; text-align: center;
           vertical-align: middle; font-size: 14px; }
  tr:nth-child(even) td { background: #fafafa; }
  .invalide { color: #b00; font-style: italic; }
</style>
"""


def _image_molecule(smiles, largeur, hauteur):
    """ Dessine une molécule et retourne une balise <img> contenant l'image encodée
        directement dans le HTML (le fichier reste donc autonome).
    """
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return '<span class="invalide">SMILES invalide</span>'
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    dessinateur = rdMolDraw2D.MolDraw2DCairo(largeur * _FACTEUR, hauteur * _FACTEUR)
    dessinateur.drawOptions().bondLineWidth = 2 * _FACTEUR
    dessinateur.DrawMolecule(mol)
    dessinateur.FinishDrawing()
    png = base64.b64encode(dessinateur.GetDrawingText()).decode()
    return (f'<img src="data:image/png;base64,{png}" '
            f'width="{largeur}" height="{hauteur}">')


def tableau_html(donnees, colonnes_smiles=('SMILES',), fichier='tableau.html',
                 titre=None, largeur=200, hauteur=100, telecharger=True):
    """ Produit un fichier HTML contenant un tableau où chaque colonne de SMILES est
        remplacée par la structure dessinée, et les autres colonnes affichées telles
        quelles.

        donnees         : liste de dictionnaires (une entrée par composé), ou
                          DataFrame pandas
        colonnes_smiles : nom(s) de la ou des colonnes contenant des SMILES
        fichier         : nom du fichier à produire
        titre           : titre affiché au-dessus du tableau (optionnel)
        largeur, hauteur: taille d'affichage des structures, en pixels
        telecharger     : si True et si on est dans Google Colab, télécharge le
                          fichier automatiquement

        L'ordre des colonnes du tableau suit l'ordre des colonnes de vos données.
    """
    if hasattr(donnees, 'to_dict'):          # accepte aussi un DataFrame pandas
        donnees = donnees.to_dict('records')
    if not donnees:
        raise ValueError("Aucune donnée à afficher.")
    if isinstance(colonnes_smiles, str):     # tolère une seule colonne non listée
        colonnes_smiles = (colonnes_smiles,)

    colonnes = list(donnees[0].keys())
    for colonne in colonnes_smiles:
        if colonne not in colonnes:
            raise ValueError(f"La colonne '{colonne}' est absente des données. "
                             f"Colonnes disponibles : {colonnes}")

    html = ['<!DOCTYPE html><html><head><meta charset="utf-8">', _STYLE, '</head><body>']
    if titre:
        html.append(f'<h1>{titre}</h1>')
    html.append('<table>')
    html.append('<tr>' + ''.join(f'<th>{colonne}</th>' for colonne in colonnes) + '</tr>')

    for ligne in donnees:
        cellules = []
        for colonne in colonnes:
            valeur = ligne[colonne]
            if colonne in colonnes_smiles:
                contenu = _image_molecule(valeur, largeur, hauteur)
            elif isinstance(valeur, float):
                contenu = f'{valeur:.2f}'
            else:
                contenu = str(valeur)
            cellules.append(f'<td>{contenu}</td>')
        html.append('<tr>' + ''.join(cellules) + '</tr>')

    html.append('</table></body></html>')

    with open(fichier, 'w', encoding='utf-8') as f:
        f.write(''.join(html))
    print(f"{len(donnees)} composés écrits dans {fichier}")

    if telecharger:
        try:
            from google.colab import files
            files.download(fichier)
        except ImportError:
            pass                              # hors de Colab : le fichier est simplement sauvegardé
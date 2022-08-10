"""
platform: any
env: any with RDKit
name: MEIF.py
MEIF utils
"""
import os

from configs.config import *
from itertools import product
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from rdkit import Chem
import pandas as pd
PROTEIN_ELEMENTS = ["pC", "pN", "pO", "pS", "pOth"]
LIGAND_ELEMENTS = ["lBr", "lC", "lCl", "lF", "lI", "lN", "lO", "lP", "lS", "lOth"]
LIGAND_DESC = ['MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex',
               'qed', 'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons',
               'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ',
               'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n',
               'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Kappa1',
               'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA14', 'SMR_VSA1', 'SMR_VSA10',
               'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7',
               'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12',
               'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6',
               'SlogP_VSA7', 'SlogP_VSA8', 'TPSA', 'EState_VSA1', 'EState_VSA10',
               'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5',
               'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1',
               'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5',
               'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3',
               'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
               'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
               'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors',
               'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
               'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount', 'MolLogP',
               'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_N',
               'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO',
               'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O',
               'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde',
               'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide',
               'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azo', 'fr_barbitur',
               'fr_benzene', 'fr_bicyclic', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester',
               'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone',
               'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
               'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
               'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitroso', 'fr_oxazole',
               'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond',
               'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_pyridine', 'fr_quatN',
               'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
               'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_urea']


class MEIF:
    _ds_version = 2020
    _desc_calculator = MolecularDescriptorCalculator(LIGAND_DESC)

    def _get_elem_head(self):
        """
        get all possible type combinations of meshed elements in proteins and ligands
        :return: header(labels) of meshed element type
        """
        ppl = [i[0] + "-" + i[1] + "-" + i[2] for i in product(PROTEIN_ELEMENTS, PROTEIN_ELEMENTS, LIGAND_ELEMENTS)]
        pll = [i[0] + "-" + i[1] + "-" + i[2] for i in product(PROTEIN_ELEMENTS, LIGAND_ELEMENTS, LIGAND_ELEMENTS)]
        return ppl + pll

    def _load_ligand(self, f_ligd):
        """
        This function takes an SDF for a ligand as input and returns a pandas DataFrame
        with its atom types labeled according to MEIF
        :param f_ligd ligand file in SDF format
        :return: pandas DataFrame
        """
        m = Chem.MolFromMolFile(f_ligd, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        ligd_atoms = []
        for atom in m.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol != "H":
                if "l" + symbol not in LIGAND_ELEMENTS:
                    symbol = "Oth"
                entry = [atom.GetIdx(), "l" + symbol]
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ligd_atoms.append(entry)
        df = pd.DataFrame(ligd_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
        assert len(set(df["ECIF_ATOM_TYPE"]) - set(LIGAND_ELEMENTS)) == 0
        return df

    # callables

    def testMEIF(self):
        """
        MEIF tester
        :return:
        """
        ligand = os.path.join(tmp_dir, "6GGH_ligand.sdf")
        print(self._load_ligand(ligand))


if __name__ == '__main__':
    helper = MEIF()
    helper.testMEIF()

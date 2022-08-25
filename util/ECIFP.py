"""
platform: any
env: any with RDKit
name: ECIFP.py
ECIFP utils. modified from https://github.com/DIFACQUIM/ECIF/blob/master/ecif.py
"""

import os

from configs.config import *
from itertools import product
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from rdkit import Chem
import pandas as pd
from scipy.spatial.distance import cdist

PROTEIN_ELEMENTS = ["C", "N", "O", "S", "POTH"]
LIGAND_ELEMENTS = ["Br", "C", "Cl", "F", "I", "N", "O", "P", "S", "LOTH"]

PROTEIN_ATOMS = ['C;4;1;3;0;0', 'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1',
                 'C;4;3;0;0;0', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                 'C;5;3;0;0;0', 'C;6;3;0;0;0', 'N;3;1;2;0;0', 'N;3;2;0;1;1',
                 'N;3;2;1;0;0', 'N;3;2;1;1;1', 'N;3;3;0;0;1', 'N;4;1;2;0;0',
                 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                 'S;2;1;1;0;0', 'S;2;2;0;0;0', 'POTH']
LIGAND_ATOMS = ['Br;1;1;0;0;0', 'C;3;3;0;1;1', 'C;4;1;1;0;0', 'C;4;1;2;0;0',
                'C;4;1;3;0;0', 'C;4;2;0;0;0', 'C;4;2;1;0;0', 'C;4;2;1;0;1',
                'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1', 'C;4;3;0;0;0',
                'C;4;3;0;0;1', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                'C;4;4;0;0;0', 'C;4;4;0;0;1', 'C;5;3;0;0;0', 'C;5;3;0;1;1',
                'C;6;3;0;0;0', 'Cl;1;1;0;0;0', 'F;1;1;0;0;0', 'I;1;1;0;0;0',
                'N;3;1;0;0;0', 'N;3;1;1;0;0', 'N;3;1;2;0;0', 'N;3;2;0;0;0',
                'N;3;2;0;0;1', 'N;3;2;0;1;1', 'N;3;2;1;0;0', 'N;3;2;1;0;1',
                'N;3;2;1;1;1', 'N;3;3;0;0;0', 'N;3;3;0;0;1', 'N;3;3;0;1;1',
                'N;4;1;2;0;0', 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'N;4;2;2;0;0',
                'N;4;2;2;0;1', 'N;4;3;0;0;0', 'N;4;3;0;0;1', 'N;4;3;1;0;0',
                'N;4;3;1;0;1', 'N;4;4;0;0;0', 'N;4;4;0;0;1', 'N;5;2;0;0;0',
                'N;5;3;0;0;0', 'N;5;3;0;1;1', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                'O;2;2;0;0;0', 'O;2;2;0;0;1', 'O;2;2;0;1;1', 'P;5;4;0;0;0',
                'P;6;4;0;0;0', 'P;6;4;0;0;1', 'P;7;4;0;0;0', 'S;2;1;0;0;0',
                'S;2;1;1;0;0', 'S;2;2;0;0;0', 'S;2;2;0;0;1', 'S;2;2;0;1;1',
                'S;3;3;0;0;0', 'S;3;3;0;0;1', 'S;4;3;0;0;0', 'S;6;4;0;0;0',
                'S;6;4;0;0;1', 'S;7;4;0;0;0', 'LOTH']

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


class ECIFP:
    _ds_version = 2020
    _desc_calculator = MolecularDescriptorCalculator(LIGAND_DESC)
    _possible_pl = [i[0] + "-" + i[1] for i in product(PROTEIN_ATOMS, LIGAND_ATOMS)]
    _cached_protein = None

    def _get_atom_type(self, atom):
        """

        :param atom: RDKit supported atom type
        :return:
          an ECIF::atom_types
        #     Atom symbol;
        #     Explicit valence;
        #     Attached heavy atoms;
        #     Attached hydrogens;
        #     Aromaticity;
        #     Ring membership
        """
        AtomType = [atom.GetSymbol(),
                    str(atom.GetExplicitValence()),
                    str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() != "H"])),
                    str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "H"])),
                    str(int(atom.GetIsAromatic())),
                    str(int(atom.IsInRing())),
                    ]

        return ";".join(AtomType)

    def _load_ligand(self, f_ligd):
        """
        This function takes an SDF for a ligand as input and returns a pandas DataFrame
        with its atom types labeled according to ECIFP
        :param f_ligd ligand file in SDF format
        :return: pandas DataFrame
        """
        m = Chem.MolFromMolFile(f_ligd, sanitize=False)
        m.UpdatePropertyCache(strict=False)
        ligd_atoms = []
        for atom in m.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol != "H":
                if symbol not in LIGAND_ELEMENTS:
                    entry = [atom.GetIdx(), "LOTH"]
                else:
                    entry = [atom.GetIdx(), self._get_atom_type(atom)]
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ligd_atoms.append(entry)
        df = pd.DataFrame(ligd_atoms)
        df.columns = ["ATOM_INDEX", "ECIFP_ATOM_TYPE", "X", "Y", "Z"]
        return df

    def _load_protein(self, f_prot):
        """
        This function takes a PDB for a protein as input and returns a pandas DataFrame
        with its atom types labeled according to ECIFP
        :param f_prot: protein file in PDB format
        :return: pandas DataFrame
        """
        fp = open(f_prot)
        prot_atoms = []
        keys = pd.read_csv(os.path.join(data_dir, "keys.csv"), sep=",")
        for line in fp:
            if line[:4] == "ATOM":
                candidate_symbol = line[12:16].replace(" ", "")
                if len(candidate_symbol) < 4 and candidate_symbol[0] != "H" or (
                        len(candidate_symbol) == 4 and candidate_symbol[0] != "H" and candidate_symbol[1] != "H"):
                    if candidate_symbol[0] not in PROTEIN_ELEMENTS:
                        candidate_symbol = "POTH"
                    else:
                        candidate_symbol = candidate_symbol[0]
                    prot_atoms.append([int(line[6:11]),
                                       line[17:20] + "-" + candidate_symbol,
                                       float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
        fp.close()
        df = pd.DataFrame(prot_atoms, columns=["ATOM_INDEX", "PDB_ATOM", "X", "Y", "Z"])
        df = df.merge(keys, left_on='PDB_ATOM', right_on='PDB_ATOM')[
            ["ATOM_INDEX", "ECIFP_ATOM_TYPE", "X", "Y", "Z"]].sort_values(by="ATOM_INDEX").reset_index(drop=True)
        df["ECIFP_ATOM_TYPE"].fillna("POTH")
        return df

    def _get_pl_pairs(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        This function returns the protein-ligand atom-type pairs for a given distance cutoff
        :param protein_f: pdb file name with dir
        :param ligand_f:  sdf file name with dir
        :param distance_cutoff:
        :return:
        """

        Target = self._load_protein(protein_f)
        Ligand = self._load_ligand(ligand_f)

        for i in ["X", "Y", "Z"]:
            Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
            Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

        # Get all possible pairs
        Pairs = list(product(Target["ECIFP_ATOM_TYPE"], Ligand["ECIFP_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIFP_PAIR"])
        Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
        Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
        Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

        Pairs = pd.concat([Pairs, Distances], axis=1)
        Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
        return Pairs

    def _get_pl_pairs_cached(self, ligand_f, distance_cutoff=6.0):
        """
        This function returns the protein-ligand atom-type pairs for a given distance cutoff, with cached protein
        :param ligand_f:  sdf file name with dir
        :param distance_cutoff:
        :return:
        """

        Target = self._cached_protein
        Ligand = self._load_ligand(ligand_f)

        for i in ["X", "Y", "Z"]:
            Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
            Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

        # Get all possible pairs
        Pairs = list(product(Target["ECIFP_ATOM_TYPE"], Ligand["ECIFP_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIFP_PAIR"])
        Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
        Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
        Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

        Pairs = pd.concat([Pairs, Distances], axis=1)
        Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
        return Pairs

    # callables
    def cache_protein(self, f_prot):
        """
        In scenario where a single protein is frequently used, cache it first to MEIF helper
        :param f_prot: protein file in PDB format
        :return:
        """
        self._cached_protein = self._load_protein(f_prot)

    def get_ligand_features_by_file(self, ligand_f):
        """
        calculate ligand descriptors using RDKit
        :param ligand_f: sdf file name with dir
        :return:
        """
        ligand = Chem.MolFromMolFile(ligand_f, sanitize=False)
        ligand.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(ligand)
        return self._DescCalc.CalcDescriptors(ligand)

    def get_ecifp(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIFP) for a protein-ligand pair
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs(protein_f, ligand_f, distance_cutoff=distance_cutoff)
        ECIFP = [list(Pairs["ECIFP_PAIR"]).count(x) for x in self._possible_pl]
        return ECIFP

    def get_ecifp_cached(self, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIFP) for a protein-ligand pair, with cached protein
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs_cached(ligand_f, distance_cutoff=distance_cutoff)
        ECIFP = [list(Pairs["ECIFP_PAIR"]).count(x) for x in self._possible_pl]
        return ECIFP

    def testLoader(self):
        """
        load protein/ligand tester
        :return:
        """
        protein = os.path.join(tmp_dir, "4gmy_protein.pdb")
        print(self._load_protein(protein))
        ligand = os.path.join(tmp_dir, "4gmy_ligand.sdf")
        print(self._load_ligand(ligand))

    def testECIFP(self):
        """
        ECIFP tester
        :return:
        """
        protein = os.path.join(tmp_dir, "1akt_protein.pdb")
        ligand = os.path.join(tmp_dir, "1akt_ligand.sdf")
        print(self.get_ecifp(protein, ligand, 6.0))


if __name__ == '__main__':
    helper = ECIFP()
    helper.testECIFP()

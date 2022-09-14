"""
platform: any
env: any with RDKit
name: ECIF.py
ECIF utils. modified from https://github.com/DIFACQUIM/ECIF/blob/master/ecif.py
"""

import os

from configs.config import *
from itertools import product
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from rdkit import Chem
import pandas as pd
from scipy.spatial.distance import cdist
from util.RDKitHelper import Mol2MolSupplier

PROTEIN_ELEMENTS = ["C", "N", "O", "S"]
LIGAND_ELEMENTS = ["Br", "C", "Cl", "F", "I", "N", "O", "P", "S"]

PROTEIN_ATOMS = ['C;4;1;3;0;0', 'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1',
                 'C;4;3;0;0;0', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                 'C;5;3;0;0;0', 'C;6;3;0;0;0', 'N;3;1;2;0;0', 'N;3;2;0;1;1',
                 'N;3;2;1;0;0', 'N;3;2;1;1;1', 'N;3;3;0;0;1', 'N;4;1;2;0;0',
                 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                 'S;2;1;1;0;0', 'S;2;2;0;0;0']
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
                'S;6;4;0;0;1', 'S;7;4;0;0;0']

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


class ECIF:
    _ds_version = 2016
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

    def _load_ligands_mol2(self, f_ligds):
        """
        for CASF docking power analysis
        This function takes a mol2 for a list of decoys as input and returns a list of pandas DataFrame
        :param f_ligds: decoys file in mol2 format
        :return: a list of pandas DataFrame
        """
        ret = []
        ligands = Mol2MolSupplier(f_ligds)
        for m in ligands:
            ligd_atoms = []
            for atom in m.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol != "H":
                    if symbol not in LIGAND_ELEMENTS:
                        continue
                    else:
                        entry = [int(atom.GetIdx()), self._get_atom_type(atom)]
                    pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                    entry.append(float("{0:.4f}".format(pos.x)))
                    entry.append(float("{0:.4f}".format(pos.y)))
                    entry.append(float("{0:.4f}".format(pos.z)))
                    ligd_atoms.append(entry)
            df = pd.DataFrame(ligd_atoms)
            df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
            ret.append(df)
        return ret

    def _load_ligands_sdf(self, f_ligds):
        """
        for CASF docking power analysis
        This function takes a sdf for a list of decoys as input and returns a list of pandas DataFrame
        :param f_ligds: decoys file in sdf format
        :return: a list of pandas DataFrame
        """
        ret = []
        ligands = Chem.SDMolSupplier(f_ligds)
        for m in ligands:
            # print(m.GetProp('_Name'))
            ligd_atoms = []
            for atom in m.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol != "H":
                    if symbol not in LIGAND_ELEMENTS:
                        continue
                    else:
                        entry = [int(atom.GetIdx()), self._get_atom_type(atom)]
                    pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                    entry.append(float("{0:.4f}".format(pos.x)))
                    entry.append(float("{0:.4f}".format(pos.y)))
                    entry.append(float("{0:.4f}".format(pos.z)))
                    ligd_atoms.append(entry)
            df = pd.DataFrame(ligd_atoms)
            df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
            ret.append(df)
        return ret

    def _load_ligand(self, f_ligd):
        """
        This function takes an SDF for a ligand as input and returns a pandas DataFrame
        with its atom types labeled according to ECIF
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
                    continue
                else:
                    entry = [int(atom.GetIdx()), self._get_atom_type(atom)]
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ligd_atoms.append(entry)
        df = pd.DataFrame(ligd_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
        return df

    def _load_ligand_mol2(self, f_ligd):
        """
        This function takes a mol2 for a ligand as input and returns a pandas DataFrame
        with its atom types labeled according to ECIF
        :param f_ligd ligand file in mol2 format
        :return: pandas DataFrame
        """
        m = Chem.MolFromMol2File(f_ligd)
        m.UpdatePropertyCache(strict=False)
        ligd_atoms = []
        for atom in m.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol != "H":
                if symbol not in LIGAND_ELEMENTS:
                    continue
                else:
                    entry = [int(atom.GetIdx()), self._get_atom_type(atom)]
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ligd_atoms.append(entry)
        df = pd.DataFrame(ligd_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
        return df

    def _load_ligand_old(self, sdf):
        """
        This function takes an SDF for a ligand as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF
        :param sdf: ligand file
        :return: a pandas DataFrame for the ligand with ECIF::atom_types
        """
        # This function takes an SDF for a ligand as input and returns it as a pandas DataFrame

        m = Chem.MolFromMolFile(sdf, sanitize=False)
        m.UpdatePropertyCache(strict=False)

        ECIF_atoms = []

        for atom in m.GetAtoms():
            if atom.GetSymbol() != "H":  # Include only non-hydrogen atoms
                entry = [int(atom.GetIdx())]
                entry.append(self._get_atom_type(atom))
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float("{0:.4f}".format(pos.x)))
                entry.append(float("{0:.4f}".format(pos.y)))
                entry.append(float("{0:.4f}".format(pos.z)))
                ECIF_atoms.append(entry)

        df = pd.DataFrame(ECIF_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
        return (df)

    def _load_protein(self, f_prot):
        """
        This function takes a PDB for a protein as input and returns a pandas DataFrame
        with its atom types labeled according to ECIF
        :param f_prot: protein file in PDB format
        :return: pandas DataFrame
        """
        fp = open(f_prot)
        prot_atoms = []
        keys = pd.read_csv(os.path.join(data_dir, "keys_ecif.csv"), sep=",")
        for line in fp:
            if line[:4] == "ATOM":
                candidate_symbol = line[12:16].replace(" ", "")
                if len(candidate_symbol) < 4 and candidate_symbol[0] != "H" or (
                        len(candidate_symbol) == 4 and candidate_symbol[0] != "H" and candidate_symbol[1] != "H"):
                    prot_atoms.append([int(line[6:11]),
                                       line[17:20] + "-" + candidate_symbol,
                                       float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
        fp.close()
        df = pd.DataFrame(prot_atoms, columns=["ATOM_INDEX", "PDB_ATOM", "X", "Y", "Z"])
        df = df.merge(keys, left_on='PDB_ATOM', right_on='PDB_ATOM')[
            ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]].sort_values(by="ATOM_INDEX").reset_index(drop=True)
        return df

    def _load_protein_old(self, pdb):
        """
        This function takes a PDB for a protein as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF
        :param pdb: protein file
        :return: a pandas DataFrame for the protein with ECIF::atom_types
        """
        Atom_Keys = pd.read_csv(os.path.join(data_dir, "keys_ecif.csv"), sep=",")
        ECIF_atoms = []

        f = open(pdb)
        for i in f:
            if i[:4] == "ATOM":
                # Include only non-hydrogen atoms
                if (len(i[12:16].replace(" ", "")) < 4 and i[12:16].replace(" ", "")[0] != "H") or (
                        len(i[12:16].replace(" ", "")) == 4 and i[12:16].replace(" ", "")[1] != "H" and
                        i[12:16].replace(" ", "")[0] != "H"):
                    ECIF_atoms.append([int(i[6:11]),
                                       i[17:20] + "-" + i[12:16].replace(" ", ""),
                                       float(i[30:38]),
                                       float(i[38:46]),
                                       float(i[46:54])
                                       ])

        f.close()

        df = pd.DataFrame(ECIF_atoms, columns=["ATOM_INDEX", "PDB_ATOM", "X", "Y", "Z"])
        df = df.merge(Atom_Keys, left_on='PDB_ATOM', right_on='PDB_ATOM')[
            ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]].sort_values(by="ATOM_INDEX").reset_index(drop=True)
        if list(df["ECIF_ATOM_TYPE"].isna()).count(True) > 0:
            print("WARNING: Protein contains unsupported atom types. Only supported atom-type pairs are counted.")
        return (df)

    def _get_pl_pairs_with_decoys_cached_mol2(self, decoy_f, distance_cutoff=6.0):
        """
        for CASF docking power analysis
        This function returns the protein-decoys(multiple ligands) atom-type pairs for a given distance cutoff,
             with cached protein
        :param decoy_f: decoys file in mol2 format
        :param distance_cutoff:
        :return:
        """
        Target = self._cached_protein
        Ligands = self._load_ligands_mol2(decoy_f)
        ret = []
        for Ligand in Ligands:
            for i in ["X", "Y", "Z"]:
                Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
                Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

            # Get all possible pairs
            Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
            Pairs = [x[0] + "-" + x[1] for x in Pairs]
            Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
            Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
            Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
            Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

            Pairs = pd.concat([Pairs, Distances], axis=1)
            Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
            ret.append(Pairs)
        return ret

    def _get_pl_pairs_with_decoys_cached_sdf(self, decoy_f, distance_cutoff=6.0):
        """
        for CASF docking power analysis
        This function returns the protein-decoys(multiple ligands) atom-type pairs for a given distance cutoff,
             with cached protein
        :param decoy_f: decoys file in sdf format
        :param distance_cutoff:
        :return:
        """
        Target = self._cached_protein
        Ligands = self._load_ligands_sdf(decoy_f)
        ret = []
        for Ligand in Ligands:
            for i in ["X", "Y", "Z"]:
                Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
                Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

            # Get all possible pairs
            Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
            Pairs = [x[0] + "-" + x[1] for x in Pairs]
            Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
            Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
            Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
            Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

            Pairs = pd.concat([Pairs, Distances], axis=1)
            Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
            ret.append(Pairs)
        return ret

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
        Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
        Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
        Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
        Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

        Pairs = pd.concat([Pairs, Distances], axis=1)
        Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
        return Pairs

    def _get_pl_pairs_mol2(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        This function returns the protein-ligand atom-type pairs for a given distance cutoff
        :param protein_f: pdb file name with dir
        :param ligand_f:  sdf file name with dir
        :param distance_cutoff:
        :return:
        """

        Target = self._load_protein(protein_f)
        Ligand = self._load_ligand_mol2(ligand_f)

        for i in ["X", "Y", "Z"]:
            Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
            Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

        # Get all possible pairs
        Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
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
        Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
        Pairs = [x[0] + "-" + x[1] for x in Pairs]
        Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
        Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
        Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
        Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

        Pairs = pd.concat([Pairs, Distances], axis=1)
        Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
        return Pairs

    # callables
    def get_pl_pairs_with_decoys_cached(self, decoy_f, distance_cutoff=6.0):
        """
        for CASF docking power analysis
        This function returns the protein-decoys(multiple ligands) atom-type pairs for a given distance cutoff,
             with cached protein
        :param decoy_f: decoys file in mol2 format
        :param distance_cutoff:
        :return:
        """
        Target = self._cached_protein
        Ligands = self._load_ligands_mol2(decoy_f)
        ret = []
        for Ligand in Ligands:
            for i in ["X", "Y", "Z"]:
                Target = Target[Target[i] < float(Ligand[i].max()) + distance_cutoff]
                Target = Target[Target[i] > float(Ligand[i].min()) - distance_cutoff]

            # Get all possible pairs
            Pairs = list(product(Target["ECIF_ATOM_TYPE"], Ligand["ECIF_ATOM_TYPE"]))
            Pairs = [x[0] + "-" + x[1] for x in Pairs]
            Pairs = pd.DataFrame(Pairs, columns=["ECIF_PAIR"])
            Distances = cdist(Target[["X", "Y", "Z"]], Ligand[["X", "Y", "Z"]], metric="euclidean")
            Distances = Distances.reshape(Distances.shape[0] * Distances.shape[1], 1)
            Distances = pd.DataFrame(Distances, columns=["DISTANCE"])

            Pairs = pd.concat([Pairs, Distances], axis=1)
            Pairs = Pairs[Pairs["DISTANCE"] <= distance_cutoff].reset_index(drop=True)
            ret.append(Pairs)
        return ret

    def load_ligands(self, f_ligds):
        """
        for CASF docking power analysis
        This function takes a mol2 for a list of decoys as input and returns a list of pandas DataFrame
        :param f_ligds: decoys file in mol2 format
        :return: a list of pandas DataFrame
        """
        ret = []
        ligands = Mol2MolSupplier(f_ligds)
        for m in ligands:
            ligd_atoms = []
            for atom in m.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol != "H":
                    if symbol not in LIGAND_ELEMENTS:
                        continue
                    else:
                        entry = [int(atom.GetIdx()), self._get_atom_type(atom)]
                    pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                    entry.append(float("{0:.4f}".format(pos.x)))
                    entry.append(float("{0:.4f}".format(pos.y)))
                    entry.append(float("{0:.4f}".format(pos.z)))
                    ligd_atoms.append(entry)
            df = pd.DataFrame(ligd_atoms)
            df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
            ret.append(df)
        return ret

    def cache_protein(self, f_prot):
        """
        In scenario where a single protein is frequently used, cache it first to ECIF helper
        :param f_prot: protein file in PDB format
        :return:
        """
        self._cached_protein = self._load_protein(f_prot)

    def get_ligand_features_by_decoy_mol2(self, f_decoy):
        """
        In scenario where decoys of a ligand is used, the LD of these decoys are the same,
        so ld should only be calculated for once
        :param f_decoy: decoys file in mol2 format
        :return:
        """
        ligand = Mol2MolSupplier(f_decoy)[0]
        return self._desc_calculator.CalcDescriptors(ligand)

    def get_ligand_features_by_decoy_sdf(self, f_decoy):
        """
        In scenario where decoys of a ligand is used, the LD of these decoys are the same,
        so ld should only be calculated for once
        :param f_decoy: decoys file in sdf format
        :return:
        """
        ligand = Chem.SDMolSupplier(f_decoy)[0]
        return self._desc_calculator.CalcDescriptors(ligand)

    def get_cached_protein(self):
        """
        get the cached protein
        :return:
        """
        return self._cached_protein

    def test_ligand_features_decoys(self, decoy_f):
        """
        this is a test to find out whether decoys' RDKit features remain the same
        tester: @ECIFPDockingPowerTest.testDecoysLD
        result: they can be regarded as the same
        :param decoy_f: decoys file in mol2 format
        :return:
        """
        ligands = Mol2MolSupplier(decoy_f)
        tups = []
        for ligand in ligands:
            tups.append(self._desc_calculator.CalcDescriptors(ligand))
        return tups

    def get_ligand_features_by_file_sdf(self, ligand_f):
        """
        calculate ligand descriptors using RDKit
        :param ligand_f: sdf file name with dir
        :return:
        """
        ligand = Chem.MolFromMolFile(ligand_f, sanitize=False)
        if ligand is None:
            return None
        ligand.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(ligand)
        return self._desc_calculator.CalcDescriptors(ligand)

    def get_ligand_features_by_file_mol2(self, ligand_f):
        """
        calculate ligand descriptors using RDKit
        :param ligand_f: mol2 file name with dir
        :return:
        """
        ligand = Chem.MolFromMol2File(ligand_f, sanitize=False)
        if ligand is None:
            return None
        ligand.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(ligand)
        return self._desc_calculator.CalcDescriptors(ligand)

    def get_ecif(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-ligand pair
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs(protein_f, ligand_f, distance_cutoff=distance_cutoff)
        ecif = [list(Pairs["ECIF_PAIR"]).count(x) for x in self._possible_pl]
        count = 0
        for number in ecif:
            if number != 0:
                count += 1

        return ecif

    def get_ecif_mol2(self, protein_f, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-ligand pair
        :param protein_f: pdb file name with dir
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs_mol2(protein_f, ligand_f, distance_cutoff=distance_cutoff)
        ecif = [list(Pairs["ECIF_PAIR"]).count(x) for x in self._possible_pl]
        count = 0
        for number in ecif:
            if number != 0:
                count += 1

        return ecif

    def get_ecif_cached(self, ligand_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-ligand pair, with cached protein
        :param ligand_f: sdf file name with dir
        :param distance_cutoff:
        :return:
        """
        Pairs = self._get_pl_pairs_cached(ligand_f, distance_cutoff=distance_cutoff)
        ecif = [list(Pairs["ECIF_PAIR"]).count(x) for x in self._possible_pl]
        return ecif

    def get_decoys_ecif_cached_mol2(self, decoy_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-decoy pair, with cached protein
        :param decoy_f: decoys file in mol2 format
        :param distance_cutoff:
        :return:
        """
        ret = []
        pairs_list = self._get_pl_pairs_with_decoys_cached_mol2(decoy_f, distance_cutoff)
        for Pairs in pairs_list:
            ret.append([list(Pairs["ECIF_PAIR"]).count(x) for x in self._possible_pl])
        return ret

    def get_decoys_ecif_cached_sdf(self, decoy_f, distance_cutoff=6.0):
        """
        get the fingerprint-like array (ECIF) for a protein-decoy pair, with cached protein
        :param decoy_f: decoys file in sdf format
        :param distance_cutoff:
        :return:
        """
        ret = []
        pairs_list = self._get_pl_pairs_with_decoys_cached_sdf(decoy_f, distance_cutoff)
        for Pairs in pairs_list:
            ret.append([list(Pairs["ECIF_PAIR"]).count(x) for x in self._possible_pl])
        return ret

    def get_possible_pl(self):
        return self._possible_pl

    def testLoader(self):
        """
        load protein/ligand tester
        :return:
        """
        # protein = os.path.join(tmp_dir, "4gmy_protein.pdb")
        # prot_table = self._load_protein(protein)
        # old_prot_table = self._load_protein_old(protein)
        # print(self._load_protein(protein))
        ligand = os.path.join(tmp_dir, "4gmy_ligand.sdf")
        lig_table = self._load_ligand(ligand)
        old_lig_table = self._load_ligand_old(ligand)
        print(self._load_ligand(ligand))
        # print(self._load_ligand_old(ligand))

    def testECIF(self):
        """
        ECIF tester
        :return:
        """
        protein = os.path.join(tmp_dir, "4gmy_protein.pdb")
        ligand = os.path.join(tmp_dir, "4gmy_ligand.sdf")
        ecif = self.get_ecif(protein, ligand, 6.0)
        print(len(ecif))
        print(ecif.count(0)/len(ecif))


if __name__ == '__main__':
    helper = ECIF()
    helper.testECIF()

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
from scipy.spatial.distance import cdist

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
    _possible_ppl = [i[0] + "-" + i[1] + "-" + i[2]
                     for i in product(PROTEIN_ELEMENTS, PROTEIN_ELEMENTS, LIGAND_ELEMENTS)]
    _possible_llp = [i[0] + "-" + i[1] + "-" + i[2]
                     for i in product(LIGAND_ELEMENTS, LIGAND_ELEMENTS, PROTEIN_ELEMENTS)]
    _cached_protein = None

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
        df.columns = ["ATOM_INDEX", "ATOM_TYPE", "X", "Y", "Z"]
        assert len(set(df["ATOM_TYPE"]) - set(LIGAND_ELEMENTS)) == 0
        return df

    def _load_protein(self, f_prot):
        """
        This function takes a PDB for a protein as input and returns a pandas DataFrame
        with its atom types labeled according to MEIF
        :param f_prot: protein file in PDB format
        :return: pandas DataFrame
        """
        fp = open(f_prot)
        prot_atoms = []
        for line in fp:
            if line[:4] == "ATOM":
                candidate_symbol = line[12:16].replace(" ", "")
                if len(candidate_symbol) < 4 and candidate_symbol[0] != "H" or (
                        len(candidate_symbol) == 4 and candidate_symbol[0] != "H" and candidate_symbol[1] != "H"):
                    prot_atoms.append([int(line[6:11]),
                                       "p" + candidate_symbol,
                                       float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
        fp.close()
        df = pd.DataFrame(prot_atoms, columns=["ATOM_INDEX", "ATOM_TYPE", "X", "Y", "Z"])
        df["ATOM_TYPE"].fillna("pOth")
        return df

    def _get_mesh_tag(self, x, y, z):
        """
        create a unique tag of a mesh which serves as an index
        :param x: protein id x in ppl or ligand id x in llp
        :param y: protein id y in ppl or ligand id y in llp
        :param z: protein id in llp or ligand id in ppl
        :return: mesh tag
        """
        return "[{}]-[{}]:[{}]".format(x, y, z) if x < y else "[{}]-[{}]:[{}]".format(y, x, z)

    def _get_pair_tag(self, x, y):
        """
        create a unique tag of a (pp|ll) pair which serves as an index
        :param x: protein id x in ppl or ligand id x in llp
        :param y: protein id y in ppl or ligand id y in llp
        :return: pair tag
        """
        return "[{}]-[{}]".format(x, y) if x < y else "[{}]-[{}]".format(y, x)

    def _get_mesh_tag_prefix(self, x):
        """
        get the in-protein | in-ligand pair tag from a mesh tag
        :param x:
        :return:
        """
        return x.split(":")[0]

    def _get_ppl(self, f_prot, f_ligd, pp_dc=6.0, pl_dc=6.0, maximum_neighbours=None):
        """
        This function returns ppl meshes for given distance cutoffs
        :param f_prot: protein file in PDB format
        :param f_ligd: ligand file in SDF format
        :param pp_dc: the protein-protein atom distance cutoff
        :param pl_dc: the protein-ligand atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return: pandas DataFrame that stands for ppl meshes
        """
        P = self._load_protein(f_prot)
        L = self._load_ligand(f_ligd)

        # filter out protein atoms that are not in "pocket"
        for i in ["X", "Y", "Z"]:
            P = P[P[i] < float(L[i].max()) + pl_dc]
            P = P[P[i] > float(L[i].min()) - pl_dc]

        # create protein-ligand pair table and filter out pairs that exceeds the PL distance threshold
        PL = list(product(P["ATOM_TYPE"], L["ATOM_TYPE"]))
        PL_INDEX = list(product(P["ATOM_INDEX"], L["ATOM_INDEX"]))
        PL = [x[0] + "-" + x[1] for x in PL]
        PL = pd.DataFrame(PL, columns=["PL_PAIR"])
        PL_INDEX = pd.DataFrame(PL_INDEX, columns=["PID", "LID"])
        Dist_PL = cdist(P[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_PL = Dist_PL.reshape(Dist_PL.shape[0] * Dist_PL.shape[1], 1)
        Dist_PL = pd.DataFrame(Dist_PL, columns=["DISTANCE"])
        PL = pd.concat([PL, Dist_PL, PL_INDEX], axis=1)
        PL = PL[PL["DISTANCE"] <= pl_dc].reset_index(drop=True)

        PFilter = pd.DataFrame(PL)
        PFilter.drop_duplicates(subset=["PID"], keep="first", inplace=True)
        PFilter = PFilter["PID"]
        P = P[P["ATOM_INDEX"].isin(PFilter)]

        # create protein-protein-ligand mesh table
        PPL = PL.merge(PL, left_on="LID", right_on="LID")
        PPL["PL_PAIR_x"] = PPL["PL_PAIR_x"].apply(lambda x: x.split("-")[0])
        PPL["PPL_TYPE"] = (PPL["PL_PAIR_x"].astype(str)) + "-" + (PPL["PL_PAIR_y"].astype(str))
        PPL = PPL.drop(columns=["PL_PAIR_x", "DISTANCE_x", "PL_PAIR_y", "DISTANCE_y"], axis=1)

        # create a unique ppl tag, format: [small protein id]-[big protein id]:[ligand id]
        PPL["PPLTag"] = PPL.apply(lambda x: self._get_mesh_tag(x["PID_x"],
                                                               x["PID_y"],
                                                               x["LID"]), axis=1)
        PPL = PPL.drop(columns=["PID_x", "PID_y", "LID"], axis=1)
        PPL.drop_duplicates(subset="PPLTag", keep="first", inplace=True)

        # for test: apply PP distance threshold
        PP = list(product(P["ATOM_INDEX"], P["ATOM_INDEX"]))
        PP = pd.DataFrame(PP, columns=["PID_x", "PID_y"])
        PP["PPTag"] = PP.apply(lambda x: self._get_pair_tag(x["PID_x"], x["PID_y"]), axis=1)
        PP = PP.drop(columns=["PID_x", "PID_y"], axis=1)
        Dist_PP = cdist(P[["X", "Y", "Z"]], P[["X", "Y", "Z"]], metric="euclidean")
        Dist_PP = Dist_PP.reshape(Dist_PP.shape[0] * Dist_PP.shape[1], 1)
        Dist_PP = pd.DataFrame(Dist_PP, columns=["PP_DISTANCE"])
        PP = pd.concat([PP, Dist_PP], axis=1).dropna()
        PP = PP[PP["PP_DISTANCE"] <= pp_dc].reset_index(drop=True)
        PP = PP[PP["PP_DISTANCE"] > 0.0]
        PPL["PPLTagPrefix"] = PPL.apply(lambda x: self._get_mesh_tag_prefix(x["PPLTag"]), axis=1)
        PPL = PPL[PPL["PPLTagPrefix"].isin(PP["PPTag"])]
        return PPL["PPL_TYPE"]

    def _get_ppl_cached(self, f_ligd, pp_dc=6.0, pl_dc=6.0, maximum_neighbours=None):
        """
        This function returns ppl meshes for given distance cutoffs, with cached protein
        :param f_ligd: ligand file in SDF format
        :param pp_dc: the protein-protein atom distance cutoff
        :param pl_dc: the protein-ligand atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return: pandas DataFrame that stands for ppl meshes
        """

        if self._cached_protein is None:
            print("Error: No cached protein found, please cache it first.")
            print("See MEIF.cache_protein()")
            return

        P = pd.DataFrame(self._cached_protein)
        L = self._load_ligand(f_ligd)

        # filter out protein atoms that are not in "pocket"
        for i in ["X", "Y", "Z"]:
            P = P[P[i] < float(L[i].max()) + pl_dc]
            P = P[P[i] > float(L[i].min()) - pl_dc]

        # create protein-ligand pair table and filter out pairs that exceeds the PL distance threshold
        PL = list(product(P["ATOM_TYPE"], L["ATOM_TYPE"]))
        PL_INDEX = list(product(P["ATOM_INDEX"], L["ATOM_INDEX"]))
        PL = [x[0] + "-" + x[1] for x in PL]
        PL = pd.DataFrame(PL, columns=["PL_PAIR"])
        PL_INDEX = pd.DataFrame(PL_INDEX, columns=["PID", "LID"])
        Dist_PL = cdist(P[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_PL = Dist_PL.reshape(Dist_PL.shape[0] * Dist_PL.shape[1], 1)
        Dist_PL = pd.DataFrame(Dist_PL, columns=["DISTANCE"])
        PL = pd.concat([PL, Dist_PL, PL_INDEX], axis=1)
        PL = PL[PL["DISTANCE"] <= pl_dc].reset_index(drop=True)

        PFilter = pd.DataFrame(PL)
        PFilter.drop_duplicates(subset=["PID"], keep="first", inplace=True)
        PFilter = PFilter["PID"]
        P = P[P["ATOM_INDEX"].isin(PFilter)]

        # create protein-protein-ligand mesh table
        PPL = PL.merge(PL, left_on="LID", right_on="LID")
        PPL["PL_PAIR_x"] = PPL["PL_PAIR_x"].apply(lambda x: x.split("-")[0])
        PPL["PPL_TYPE"] = (PPL["PL_PAIR_x"].astype(str)) + "-" + (PPL["PL_PAIR_y"].astype(str))
        PPL = PPL.drop(columns=["PL_PAIR_x", "DISTANCE_x", "PL_PAIR_y", "DISTANCE_y"], axis=1)

        # create a unique ppl tag, format: [small protein id]-[big protein id]:[ligand id]
        PPL["PPLTag"] = PPL.apply(lambda x: self._get_mesh_tag(x["PID_x"],
                                                               x["PID_y"],
                                                               x["LID"]), axis=1)
        PPL = PPL.drop(columns=["PID_x", "PID_y", "LID"], axis=1)
        PPL.drop_duplicates(subset="PPLTag", keep="first", inplace=True)

        # for test: apply PP distance threshold
        PP = list(product(P["ATOM_INDEX"], P["ATOM_INDEX"]))
        PP = pd.DataFrame(PP, columns=["PID_x", "PID_y"])
        PP["PPTag"] = PP.apply(lambda x: self._get_pair_tag(x["PID_x"], x["PID_y"]), axis=1)
        PP = PP.drop(columns=["PID_x", "PID_y"], axis=1)
        Dist_PP = cdist(P[["X", "Y", "Z"]], P[["X", "Y", "Z"]], metric="euclidean")
        Dist_PP = Dist_PP.reshape(Dist_PP.shape[0] * Dist_PP.shape[1], 1)
        Dist_PP = pd.DataFrame(Dist_PP, columns=["PP_DISTANCE"])
        PP = pd.concat([PP, Dist_PP], axis=1).dropna()
        PP = PP[PP["PP_DISTANCE"] <= pp_dc].reset_index(drop=True)
        PP = PP[PP["PP_DISTANCE"] > 0.0]
        PPL["PPLTagPrefix"] = PPL.apply(lambda x: self._get_mesh_tag_prefix(x["PPLTag"]), axis=1)
        PPL = PPL[PPL["PPLTagPrefix"].isin(PP["PPTag"])]
        return PPL["PPL_TYPE"]

    def _get_llp(self, f_prot, f_ligd, ll_dc=6.0, pl_dc=6.0, maximum_neighbours=None):
        """
        This function returns llp meshes for given distance cutoffs
        :param f_prot: protein file in PDB format
        :param f_ligd: ligand file in SDF format
        :param ll_dc: the ligand-ligand atom distance cutoff
        :param pl_dc: the protein-ligand atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return: pandas DataFrame that stands for ppl meshes
        """
        P = self._load_protein(f_prot)
        L = self._load_ligand(f_ligd)

        # filter out protein atoms that are not in "pocket"
        for i in ["X", "Y", "Z"]:
            P = P[P[i] < float(L[i].max()) + pl_dc]
            P = P[P[i] > float(L[i].min()) - pl_dc]

        # create protein-ligand pair table and filter out pairs that exceeds the PL distance threshold
        PL = list(product(P["ATOM_TYPE"], L["ATOM_TYPE"]))
        PL_INDEX = list(product(P["ATOM_INDEX"], L["ATOM_INDEX"]))
        PL = [x[0] + "-" + x[1] for x in PL]
        PL = pd.DataFrame(PL, columns=["PL_PAIR"])
        PL_INDEX = pd.DataFrame(PL_INDEX, columns=["PID", "LID"])
        Dist_PL = cdist(P[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_PL = Dist_PL.reshape(Dist_PL.shape[0] * Dist_PL.shape[1], 1)
        Dist_PL = pd.DataFrame(Dist_PL, columns=["DISTANCE"])
        PL = pd.concat([PL, Dist_PL, PL_INDEX], axis=1)
        PL = PL[PL["DISTANCE"] <= pl_dc].reset_index(drop=True)

        LFilter = pd.DataFrame(PL)
        LFilter.drop_duplicates(subset=["LID"], keep="first", inplace=True)
        LFilter = LFilter["LID"]
        L = L[L["ATOM_INDEX"].isin(LFilter)]

        # create ligand-ligand-protein mesh table
        LLP = PL.merge(PL, left_on="PID", right_on="PID")
        LLP["PL_PAIR_x"] = LLP["PL_PAIR_x"].apply(lambda x: x.split("-")[0])
        LLP["LLP_TYPE"] = (LLP["PL_PAIR_x"].astype(str)) + "-" + (LLP["PL_PAIR_y"].astype(str))
        LLP = LLP.drop(columns=["PL_PAIR_x", "DISTANCE_x", "PL_PAIR_y", "DISTANCE_y"], axis=1)

        # create a unique llp tag, format: [small ligand id]-[big ligand id]:[protein id]
        LLP["LLPTag"] = LLP.apply(lambda x: self._get_mesh_tag(x["LID_x"],
                                                               x["LID_y"],
                                                               x["PID"]), axis=1)
        LLP = LLP.drop(columns=["LID_x", "LID_y", "PID"], axis=1)
        LLP.drop_duplicates(subset="LLPTag", keep="first", inplace=True)

        # for test: apply LL distance threshold
        LL = list(product(L["ATOM_INDEX"], L["ATOM_INDEX"]))
        LL = pd.DataFrame(LL, columns=["LID_x", "LID_y"])
        LL["LLTag"] = LL.apply(lambda x: self._get_pair_tag(x["LID_x"], x["LID_y"]), axis=1)
        LL = LL.drop(columns=["LID_x", "LID_y"], axis=1)
        Dist_LL = cdist(L[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_LL = Dist_LL.reshape(Dist_LL.shape[0] * Dist_LL.shape[1], 1)
        Dist_LL = pd.DataFrame(Dist_LL, columns=["LL_DISTANCE"])
        LL = pd.concat([LL, Dist_LL], axis=1).dropna()
        LL = LL[LL["LL_DISTANCE"] <= ll_dc].reset_index(drop=True)
        LL = LL[LL["LL_DISTANCE"] > 0.0]
        LLP["LLPTagPrefix"] = LLP.apply(lambda x: self._get_mesh_tag_prefix(x["LLPTag"]), axis=1)
        LLP = LLP[LLP["LLPTagPrefix"].isin(LL["LLTag"])]
        return LLP["LLP_TYPE"]

    def _get_llp_cached(self, f_ligd, ll_dc=6.0, pl_dc=6.0, maximum_neighbours=None):
        """
        This function returns llp meshes for given distance cutoffs, with cached protein
        :param f_ligd: ligand file in SDF format
        :param ll_dc: the ligand-ligand atom distance cutoff
        :param pl_dc: the protein-ligand atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return: pandas DataFrame that stands for ppl meshes
        """

        if self._cached_protein is None:
            print("Error: No cached protein found, please cache it first.")
            print("See MEIF.cache_protein()")
            return

        P = pd.DataFrame(self._cached_protein)
        L = self._load_ligand(f_ligd)

        # filter out protein atoms that are not in "pocket"
        for i in ["X", "Y", "Z"]:
            P = P[P[i] < float(L[i].max()) + pl_dc]
            P = P[P[i] > float(L[i].min()) - pl_dc]

        # create protein-ligand pair table and filter out pairs that exceeds the PL distance threshold
        PL = list(product(P["ATOM_TYPE"], L["ATOM_TYPE"]))
        PL_INDEX = list(product(P["ATOM_INDEX"], L["ATOM_INDEX"]))
        PL = [x[0] + "-" + x[1] for x in PL]
        PL = pd.DataFrame(PL, columns=["PL_PAIR"])
        PL_INDEX = pd.DataFrame(PL_INDEX, columns=["PID", "LID"])
        Dist_PL = cdist(P[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_PL = Dist_PL.reshape(Dist_PL.shape[0] * Dist_PL.shape[1], 1)
        Dist_PL = pd.DataFrame(Dist_PL, columns=["DISTANCE"])
        PL = pd.concat([PL, Dist_PL, PL_INDEX], axis=1)
        PL = PL[PL["DISTANCE"] <= pl_dc].reset_index(drop=True)

        LFilter = pd.DataFrame(PL)
        LFilter.drop_duplicates(subset=["LID"], keep="first", inplace=True)
        LFilter = LFilter["LID"]
        L = L[L["ATOM_INDEX"].isin(LFilter)]

        # create ligand-ligand-protein mesh table
        LLP = PL.merge(PL, left_on="PID", right_on="PID")
        LLP["PL_PAIR_x"] = LLP["PL_PAIR_x"].apply(lambda x: x.split("-")[0])
        LLP["LLP_TYPE"] = (LLP["PL_PAIR_x"].astype(str)) + "-" + (LLP["PL_PAIR_y"].astype(str))
        LLP = LLP.drop(columns=["PL_PAIR_x", "DISTANCE_x", "PL_PAIR_y", "DISTANCE_y"], axis=1)

        # create a unique llp tag, format: [small ligand id]-[big ligand id]:[protein id]
        LLP["LLPTag"] = LLP.apply(lambda x: self._get_mesh_tag(x["LID_x"],
                                                               x["LID_y"],
                                                               x["PID"]), axis=1)
        LLP = LLP.drop(columns=["LID_x", "LID_y", "PID"], axis=1)
        LLP.drop_duplicates(subset="LLPTag", keep="first", inplace=True)

        # for test: apply LL distance threshold
        LL = list(product(L["ATOM_INDEX"], L["ATOM_INDEX"]))
        LL = pd.DataFrame(LL, columns=["LID_x", "LID_y"])
        LL["LLTag"] = LL.apply(lambda x: self._get_pair_tag(x["LID_x"], x["LID_y"]), axis=1)
        LL = LL.drop(columns=["LID_x", "LID_y"], axis=1)
        Dist_LL = cdist(L[["X", "Y", "Z"]], L[["X", "Y", "Z"]], metric="euclidean")
        Dist_LL = Dist_LL.reshape(Dist_LL.shape[0] * Dist_LL.shape[1], 1)
        Dist_LL = pd.DataFrame(Dist_LL, columns=["LL_DISTANCE"])
        LL = pd.concat([LL, Dist_LL], axis=1).dropna()
        LL = LL[LL["LL_DISTANCE"] <= ll_dc].reset_index(drop=True)
        LL = LL[LL["LL_DISTANCE"] > 0.0]
        LLP["LLPTagPrefix"] = LLP.apply(lambda x: self._get_mesh_tag_prefix(x["LLPTag"]), axis=1)
        LLP = LLP[LLP["LLPTagPrefix"].isin(LL["LLTag"])]
        return LLP["LLP_TYPE"]

    # callables

    def cache_protein(self, f_prot):
        """
        In scenario where a single protein is frequently used, cache it first to MEIF helper
        :param f_prot: protein file in PDB format
        :return:
        """
        self._cached_protein = self._load_protein(f_prot)

    def get_meif(self, f_prot, f_ligd, ll_dc, pp_dc, lp_dc, maximum_neighbours=None):
        """
        MEIF fingerprint calculation
        :param f_prot: protein file in PDB format
        :param f_ligd: ligand file in SDF format
        :param ll_dc: the ligand-ligand atom distance cutoff
        :param pp_dc: the protein-protein atom distance cutoff
        :param lp_dc: the ligand-protein atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return:
        """
        ppl = self._get_ppl(f_prot, f_ligd, pp_dc, lp_dc, maximum_neighbours)
        llp = self._get_llp(f_prot, f_ligd, ll_dc, lp_dc, maximum_neighbours)
        MEIF = [list(ppl).count(x) for x in self._possible_ppl] + [list(llp).count(x) for x in self._possible_llp]
        return MEIF

    def get_meif_cached(self, f_ligd, ll_dc, pp_dc, lp_dc, maximum_neighbours=None):
        """
        MEIF fingerprint calculation, with cached protein
        :param f_ligd: ligand file in SDF format
        :param ll_dc: the ligand-ligand atom distance cutoff
        :param pp_dc: the protein-protein atom distance cutoff
        :param lp_dc: the ligand-protein atom distance cutoff
        :param maximum_neighbours: maximum number of protein neighbours, for test
        :return:
        """
        if self._cached_protein is None:
            print("Error: No cached protein found, please cache it first.")
            print("See MEIF.cache_protein()")
            return
        ppl = self._get_ppl_cached(f_ligd, pp_dc, lp_dc, maximum_neighbours)
        llp = self._get_llp_cached(f_ligd, ll_dc, lp_dc, maximum_neighbours)
        MEIF = [list(ppl).count(x) for x in self._possible_ppl] + [list(llp).count(x) for x in self._possible_llp]
        return MEIF

    def get_ligand_features_by_file(self, f_ligd):
        """
        calculate ligand descriptors using RDKit
        :param f_ligd: ligand file in SDF format
        :return:
        """
        ligand = Chem.MolFromMolFile(f_ligd, sanitize=False)
        ligand.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(ligand)
        return self._desc_calculator.CalcDescriptors(ligand)

    def get_elem_head(self):
        """
        get all possible type combinations of meshed elements in proteins and ligands
        :return: header(labels) of meshed element type
        """
        return self._possible_ppl + self._possible_llp

    def testLoader(self):
        """
        load protein/ligand tester
        :return:
        """
        protein = os.path.join(tmp_dir, "6GGH_H_protein.pdb")
        count_poth = 0
        for atom in self._load_protein(protein)["ATOM_TYPE"]:
            if len(atom) > 2:
                count_poth += 1
        print(count_poth)

        ligand = os.path.join(tmp_dir, "6GGH_ligand.sdf")
        count_loth = 0
        for atom in self._load_ligand(ligand)["ATOM_TYPE"]:
            if len(atom) > 2:
                count_loth += 1
        print(count_loth)

    def testMEIF(self):
        """
        MEIF tester
        :return:
        """
        protein = os.path.join(tmp_dir, "4gmy_protein.pdb")
        ligand = os.path.join(tmp_dir, "4gmy_ligand.sdf")
        ppl = self._get_ppl(protein, ligand)
        llp = self._get_llp(protein, ligand)
        MEIF = [list(ppl).count(x) for x in self._possible_ppl] + [list(llp).count(x) for x in self._possible_llp]
        print(MEIF.count(0)/len(MEIF))


if __name__ == '__main__':
    helper = MEIF()
    helper.testMEIF()

"""
platform: win
env: any with sklearn
name: ecif_predict.py
do prediction with ECIF models
"""
import os
from configs.config import model_test_dir, casf_dir, tmp_dir, dataset_2016_dir
from configs.config import ecif_gbt, ecif_catboost, ecif_lightgbmxt
import pandas as pd
from util.ECIF import ECIF, LIGAND_DESC
from util.RDKitHelper import get_decoy_names, get_decoy_names_sdf
import pickle
from tqdm import tqdm
import time


class ECIF_Predictor:
    _model = None
    _ecif_helper = None

    def __init__(self, model):
        """

        :param model:
        """
        self._ecif_helper = ECIF()
        self._model = pickle.load(open(model, 'rb'))

    def _predict(self, f_pro, f_lig):
        """
        predict on protein-ligand pairs not in pdbbind dataset
        :param f_pro: protein file, .pdb
        :param f_lig: ligand file, .sdf
        :return:
        """
        ecif = self._ecif_helper.get_ecif(f_pro, f_lig, float(6.0))
        ld = self._ecif_helper.get_ligand_features_by_file(f_lig)
        data = ecif + list(ld)
        cols = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        data_f = pd.DataFrame([data], columns=cols)
        return self._model.predict(data_f)[0]

    def _predict_on_decoy_by_id_mol2(self, lid):
        """
        for CASF docking power analysis
        :param lid: ligand/protein id
        :return:
        """
        f_pro = os.path.join(casf_dir["core"], lid, "{}_protein.pdb".format(lid))
        f_lig = os.path.join(casf_dir["core"], lid, "{}_ligand.sdf".format(lid))
        f_dec = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))

        self._ecif_helper.cache_protein(f_pro)
        ecifs = self._ecif_helper.get_decoys_ecif_cached_mol2(f_dec, float(6.0))
        ld = self._ecif_helper.get_ligand_features_by_decoy_mol2(f_dec)
        cols = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        ret = []
        dec_names = get_decoy_names(f_dec)

        for ecif in tqdm(ecifs):
            data = ecif + list(ld)
            data_f = pd.DataFrame([data], columns=cols)
            ret.append(self._model.predict(data_f)[0])

        result = pd.DataFrame({"#code": dec_names, "score": ret}, columns=["#code", "score"]) \
            .sort_values(by="score", ascending=False)
        result.loc[-1] = ["{}_ligand".format(lid), self._predict(f_pro, f_lig)]
        result.to_csv(os.path.join(tmp_dir, "ecif_decoy_pred", "{}_score.dat".format(lid)), index=False)
        return

    def _predict_on_decoy_by_id_sdf(self, lid):
        """
        for CASF docking power analysis, decoys in sdf format
        :param lid: ligand/protein id
        :return:
        """
        f_pro = os.path.join(casf_dir["core"], lid, "{}_protein.pdb".format(lid))
        f_lig = os.path.join(casf_dir["core"], lid, "{}_ligand.sdf".format(lid))
        f_dec = os.path.join(casf_dir["decoys_sdf"], "{}_decoys.sdf".format(lid))

        self._ecif_helper.cache_protein(f_pro)
        ecifs = self._ecif_helper.get_decoys_ecif_cached_sdf(f_dec, float(6.0))
        ld = self._ecif_helper.get_ligand_features_by_decoy_sdf(f_dec)
        cols = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        ret = []

        dec_names = get_decoy_names_sdf(f_dec)

        for ecif in tqdm(ecifs):
            data = ecif + list(ld)
            data_f = pd.DataFrame([data], columns=cols)
            ret.append(self._model.predict(data_f)[0])

        result = pd.DataFrame({"#code": dec_names, "score": ret}, columns=["#code", "score"]) \
            .sort_values(by="score", ascending=False)
        result.loc[-1] = ["{}_ligand".format(lid), self._predict(f_pro, f_lig)]
        result.to_csv(os.path.join(tmp_dir, "ecif_decoy_pred", "{}_score.dat".format(lid)), index=False)
        return

    # callables
    def load_model(self, model):
        self._model = pickle.load(open(model, 'rb'))

    def multi_ligd_pred(self, test_name, model):
        """
        for extra tests, do estimations on a list of proteins that are not belongs to PDBBind dataset;
        and multiple choices of ligands are matched to one protein
        Caution: file structure of subtest folder should be:
        ├── {test_name}(folder)
            ├── {protein_id}(folder)
                ├── {protein_id}_protein.pdb
                ├── ligs(folder)
                    ├── {ligand_name}.sdf
                    ├── {ligand_name}.sdf
                    ├── ...
            ├── {protein_id}(folder)
                ├── {protein_id}_protein.pdb
                ├── ligs(folder)
                    ├── {ligand_name}.sdf
                    ├── {ligand_name}.sdf
                    ├── ...
            ├── ...
        :param test_name: subtest folder
        :param model: model
        :return:
        """
        print("\n")
        ids = os.listdir(os.path.join(model_test_dir, test_name))
        start = time.perf_counter()
        if model is not None:
            self.load_model(model)
        for i in ids:
            preds = []
            lig_name = []
            f_prot = os.path.join(model_test_dir, test_name, "{}\\{}_protein.pdb".format(i, i))
            f_ligs = os.path.join(model_test_dir, test_name, "{}\\ligs".format(i))
            for j in tqdm(os.listdir(f_ligs)):
                f_lig = os.path.join(f_ligs, j)
                pk = self._predict(f_prot, f_lig)
                preds.append(pk)
                lig_name.append(j[:-4])
            df_protein_name = pd.DataFrame([i] * len(lig_name), columns=["protein"])
            df_lig_name = pd.DataFrame(lig_name, columns=["ligand"])
            df_prediction = pd.DataFrame(preds, columns=["prediction"])
            result = df_protein_name.join(df_lig_name.join(df_prediction))
            result.to_csv(os.path.join(model_test_dir, test_name, "{}\\{}_result.csv".format(i, i)))

        end = time.perf_counter()
        print('\n')
        print('run time: {} seconds'.format(round(end - start)))

    def predict_on_core(self, model):
        """
        for CASF analysis, make predictions on CASF core set
        :param model: specific model generated by sklearn
        :return:
        """
        print("\n")
        ids = os.listdir(casf_dir["core"])
        start = time.perf_counter()
        if model is not None:
            self.load_model(model)
        preds = []
        lig_names = []
        for i in tqdm(ids):
            f_prot = os.path.join(casf_dir["core"], "{}\\{}_protein.pdb".format(i, i))
            f_ligd = os.path.join(casf_dir["core"], "{}\\{}_ligand.sdf".format(i, i))
            pk = self._predict(f_prot, f_ligd)
            preds.append(pk)
            lig_names.append(i)
        result = pd.DataFrame({"#code": lig_names, "score": preds}, columns=["#code", "score"])
        result.to_csv(os.path.join(tmp_dir, "ecif.dat"), index=False)

        end = time.perf_counter()
        print('\n')
        print('run time: {} seconds'.format(round(end - start)))

    def predict_on_decoy_by_id(self, lid, model):
        """
        for CASF docking power analysis
        :param lid: ligand/protein id
        :param model:
        :return:
        """
        start = time.perf_counter()
        if model is not None:
            self.load_model(model)

        f_pro = os.path.join(casf_dir["core"], lid, "{}_protein.pdb".format(lid))
        f_lig = os.path.join(casf_dir["core"], lid, "{}_ligand.sdf".format(lid))
        f_dec = os.path.join(casf_dir["decoys"], "{}_decoys.mol2".format(lid))

        self._ecif_helper.cache_protein(f_pro)
        ecifs = self._ecif_helper.get_decoys_ecif_cached_mol2(f_dec, float(6.0))
        ld = self._ecif_helper.get_ligand_features_by_decoy_mol2(f_dec)
        cols = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        ret = []
        dec_names = get_decoy_names(f_dec)

        for ecif in tqdm(ecifs):
            data = ecif + list(ld)
            data_f = pd.DataFrame([data], columns=cols)
            ret.append(self._model.predict(data_f)[0])

        result = pd.DataFrame({"#code": dec_names, "score": ret}, columns=["#code", "score"]) \
            .sort_values(by="score", ascending=False)
        result.loc[-1] = ["{}_ligand".format(lid), self._predict(f_pro, f_lig)]
        end = time.perf_counter()
        print('\n')
        print('run time: {} seconds'.format(round(end - start)))
        result.to_csv(os.path.join(tmp_dir, "ecif_decoy_pred", "{}_score.dat".format(lid)), index=False)
        return

    def predict_on_decoy_mol2(self, model):
        """
        for CASF docking power analysis
        :param model:
        :return:
        """
        start = time.perf_counter()
        if model is not None:
            self.load_model(model)
        lids = os.listdir(dataset_2016_dir["core"])
        err_lids = []
        for lid in lids:
            try:
                self._predict_on_decoy_by_id_mol2(lid)
            except Exception:
                err_lids.append(lid)
                continue
        end = time.perf_counter()
        print('\n')
        print('run time: {} seconds'.format(round(end - start)))
        print('bad decoy ids:')
        print(err_lids)

    def predict_on_decoy_sdf(self, model):
        """
        for CASF docking power analysis
        :param model:
        :return:
        """
        start = time.perf_counter()
        if model is not None:
            self.load_model(model)
        lids = os.listdir(dataset_2016_dir["core"])
        err_lids = []
        for lid in lids:
            try:
                self._predict_on_decoy_by_id_sdf(lid)
            except Exception:
                err_lids.append(lid)
                continue
        end = time.perf_counter()
        print('\n')
        print('run time: {} seconds'.format(round(end - start)))
        print('bad decoy ids:')
        print(err_lids)

    def predict_on_single_decoy_file(self, f_dec, model):
        """
        do prediction on single decoy file, for those fail in either predict_on_decoy_mol2() or predict_on_decoy_sdf()
        :param f_dec: decoy file, decoys in sdf/mol2 format
        :param model:
        :return:
        """
        if model is not None:
            self.load_model(model)

        lid = f_dec[-16:-12]
        f_pro = os.path.join(casf_dir["core"], lid, "{}_protein.pdb".format(lid))
        f_lig = os.path.join(casf_dir["core"], lid, "{}_ligand.sdf".format(lid))

        self._ecif_helper.cache_protein(f_pro)
        postfix = f_dec.split(".")[-1]
        if postfix in {"sdf", "mol2"}:
            ecifs = self._ecif_helper.__getattribute__("get_decoys_ecif_cached_" + postfix)(f_dec, float(6.0))
            ld = self._ecif_helper.__getattribute__("get_ligand_features_by_decoy_" + postfix)(f_dec)
        else:
            return
        cols = self._ecif_helper.get_possible_pl() + LIGAND_DESC
        ret = []
        if postfix == "sdf":
            dec_names = get_decoy_names_sdf(f_dec)
        else:
            dec_names = get_decoy_names(f_dec)

        for ecif in tqdm(ecifs):
            data = ecif + list(ld)
            data_f = pd.DataFrame([data], columns=cols)
            ret.append(self._model.predict(data_f)[0])

        result = pd.DataFrame({"#code": dec_names, "score": ret}, columns=["#code", "score"]) \
            .sort_values(by="score", ascending=False)
        result.loc[-1] = ["{}_ligand".format(lid), self._predict(f_pro, f_lig)]
        result.to_csv(os.path.join(tmp_dir, "ecif_decoy_pred", "{}_score.dat".format(lid)), index=False)
        return


if __name__ == '__main__':
    predictor = ECIF_Predictor(ecif_catboost)
    # predictor.predict_on_core(None)
    predictor.predict_on_decoy_sdf(None)
    predictor.predict_on_single_decoy_file(os.path.join(casf_dir["decoys"], "4mme_decoys.mol2"), None)



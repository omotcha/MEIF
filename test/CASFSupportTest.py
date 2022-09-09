"""
platform: any
env: any
name: CASFSupportTest.py
CASF Support Tester
"""
import os

from predict.ecif_predict import ECIF_Predictor
from configs.config import ecif_gbt, ecif_wold
from configs.config import casf_dir


def testTargetDecoyPrediction():
    predictor = ECIF_Predictor(ecif_gbt)
    # predictor.predict_on_target_decoy_sdf(
    #     os.path.join(casf_dir["core"], "3ejr", "3ejr_protein.pdb"),
    #     os.path.join(casf_dir["decoys_screening_sdf"], "3ejr", "3ejr_1a30.sdf")
    # )
    # predictor.predict_sdf(
    #     os.path.join(casf_dir["core"], "3ejr", "3ejr_protein.pdb"),
    #     os.path.join(casf_dir["decoys_screening_sdf"], "3ejr", "3ejr_1a30.sdf")
    # )
    pred_with_ld = predictor.predict_sdf(
        os.path.join(casf_dir["core"], "3ejr", "3ejr_protein.pdb"),
        os.path.join(casf_dir["core"], "1a30", "1a30_ligand.sdf")
    )
    predictor = ECIF_Predictor(ecif_wold)
    pred_without_ld = predictor.predict_wold_sdf(
        os.path.join(casf_dir["core"], "3ejr", "3ejr_protein.pdb"),
        os.path.join(casf_dir["core"], "1a30", "1a30_ligand.sdf")
    )
    print(pred_with_ld, pred_without_ld)


if __name__ == '__main__':
    testTargetDecoyPrediction()

"""
platform: win
env: any with sklearn
name: train_with_sklearn.py
train ECIF using sklearn, for comparison tests
modified from https://github.com/DIFACQUIM/ECIF/blob/master/03_Examples(ModelTraining).ipynb
"""
import os
from configs.config import meif_data_dir, data_dir
import pandas as pd
import sklearn.tree
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import pearsonr
from math import sqrt
import pickle


def train_ecif():
    ecif = pd.read_csv(os.path.join(meif_data_dir, "ECIF_6.0.csv"))
    ld = pd.read_csv(os.path.join(meif_data_dir, "LD.csv"))
    binding_data = pd.read_csv(os.path.join(data_dir, "BindingData.csv"))
    ecif = ecif.merge(ld, left_on="PDB", right_on="PDB")
    ecif = ecif.merge(binding_data, left_on="PDB", right_on="PDB")

    x_train = ecif[ecif["SET"] == "Train"][list(ecif.columns)[1:-2]]
    y_train = ecif[ecif["SET"] == "Train"]["pK"]

    x_test = ecif[ecif["SET"] == "Test"][list(ecif.columns)[1:-2]]
    y_test = ecif[ecif["SET"] == "Test"]["pK"]

    GBT = GradientBoostingRegressor(random_state=1206, n_estimators=20000, max_features="sqrt", max_depth=8,
                                    min_samples_split=3, learning_rate=0.005, loss="ls", subsample=0.7)
    GBT.fit(x_train, y_train)

    y_pred_GBT = GBT.predict(x_test)
    print("Pearson correlation coefficient for GBT: ", pearsonr(y_test, y_pred_GBT)[0])
    print("RMSE for GBT:", sqrt(mean_squared_error(y_test, y_pred_GBT)))

    pickle.dump(GBT, open(os.path.join(data_dir, "ECIF6_LD_GBT.pkl"), 'wb'))


def train_ecifp():
    ecifp_test = pd.read_csv(os.path.join(meif_data_dir, "ECIFP_test_6.0.csv"))
    ecifp_train = pd.read_csv(os.path.join(meif_data_dir, "ECIFP_train_6.0.csv"))
    x_test = ecifp_test[list(ecifp_test.columns)[1:-1]]
    y_test = ecifp_test['pk']
    x_train = ecifp_train[list(ecifp_train.columns)[1:-1]]
    y_train = ecifp_train['pk']

    GBT = GradientBoostingRegressor(random_state=1206, n_estimators=20000, max_features="sqrt", max_depth=8,
                                    min_samples_split=3, learning_rate=0.005, loss="ls", subsample=0.7)
    GBT.fit(x_train, y_train)

    y_pred_GBT = GBT.predict(x_test)
    print("Pearson correlation coefficient for GBT: ", pearsonr(y_test, y_pred_GBT)[0])
    print("RMSE for GBT:", sqrt(mean_squared_error(y_test, y_pred_GBT)))

    pickle.dump(GBT, open(os.path.join(data_dir, "ECIFP_LD_GBT.pkl"), 'wb'))


if __name__ == '__main__':
    train_ecifp()

"""
platform: win
env: any with autogluon
name: train_with_ag.py
do training on MEIF tabular data using autogluon model generator
"""
import os
from autogluon.tabular import TabularDataset, TabularPredictor
from configs.config import meif_data_dir, data_dir
from sklearn.model_selection import train_test_split
from autogluon.core.metrics import make_scorer
import sklearn.metrics
import pandas as pd


class MEIF_Trainer:

    _meif_data = None
    _meif_train_data = None
    _meif_test_data = None
    _split_factor = 0
    _distance_cutoff = 6.0

    def __init__(self, distance_cutoff=6.0, splitted=True, tag=""):
        """
        init
        :param distance_cutoff:
        :param splitted: if data is already splitted
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        """
        if splitted:
            if tag is None or tag == "":
                self._meif_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "MEIF_train_{}.csv".format(distance_cutoff)))
                self._ecifp_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "MEIF_test_{}.csv".format(distance_cutoff)))
            else:
                self._meif_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "MEIF_train_{}_{}.csv".format(distance_cutoff, tag)))
                self._meif_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "MEIF_test_{}_{}.csv".format(distance_cutoff, tag)))
            self._split_factor = 0
        else:
            if tag is None or tag == "":
                self._meif_data = TabularDataset(os.path.join(meif_data_dir, "MEIF_{}.csv".format(distance_cutoff)))
            else:
                self._meif_data = TabularDataset(
                    os.path.join(meif_data_dir, "MEIF_{}_{}.csv".format(distance_cutoff, tag)))
            self._split_factor = 0.01465823

        self._distance_cutoff = distance_cutoff

    # callables

    def set_split_factor(self, split_factor):
        """
        set a new split factor
        :param split_factor: size(test) / size(test+train)
        :return:
        """
        self._split_factor = split_factor

    def train(self):
        """
        single-cutoff training using autogluon
        :return:
        """
        if self._split_factor == 0:
            train_x = self._meif_train_data[list(self._meif_train_data.columns)[1:-1]]
            train_y = self._meif_train_data['pk']
            test_x = self._meif_test_data[list(self._meif_train_data.columns)[1:-1]]
            test_y = self._meif_test_data['pk']
        else:
            X = self._meif_data[list(self._meif_data.columns)[1:-1]]
            Y = self._meif_data['pk']
            train_x, test_x, train_y, test_y = train_test_split(X,
                                                                Y,
                                                                test_size=self._split_factor
                                                                if 0 < self._split_factor < 1 else 0.01465823,
                                                                random_state=0)

        print("training...\n")
        ag_r2_scorer = make_scorer(name='r2',
                                   score_func=sklearn.metrics.r2_score,
                                   optimum=1,
                                   greater_is_better=True)
        predictor = TabularPredictor(label='pk',
                                     eval_metric='root_mean_squared_error').fit(train_data=train_x.join(train_y))
        print("\n - finished -\n")
        leaderboard = predictor.leaderboard(test_x.join(test_y), extra_metrics=[ag_r2_scorer], silent=True)
        print(leaderboard[['model', 'score_val', 'r2']])


class ECIFP_Trainer:

    _ecifp_data = None
    _ecifp_train_data = None
    _ecifp_test_data = None
    _split_factor = 0
    _distance_cutoff = 6.0

    def __init__(self, distance_cutoff=6.0, splitted=True, tag=""):
        """
        init
        :param distance_cutoff:
        :param splitted: if data is already splitted
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        """
        if splitted:
            if tag is None or tag == "":
                self._ecifp_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIFP_train_{}.csv".format(distance_cutoff)))
                self._ecifp_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIFP_test_{}.csv".format(distance_cutoff)))
            else:
                self._ecifp_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIFP_train_{}_{}.csv".format(distance_cutoff, tag)))
                self._ecifp_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIFP_test_{}_{}.csv".format(distance_cutoff, tag)))
            self._split_factor = 0
        else:
            if tag is None or tag == "":
                self._ecifp_data = TabularDataset(os.path.join(meif_data_dir, "ECIFP_{}.csv".format(distance_cutoff)))
            else:
                self._ecifp_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIFP_{}_{}.csv".format(distance_cutoff, tag)))
            self._split_factor = 0.01465823

        self._distance_cutoff = distance_cutoff

    # callables

    def set_split_factor(self, split_factor):
        """
        set a new split factor
        :param split_factor: size(test) / size(test+train)
        :return:
        """
        self._split_factor = split_factor

    def train(self):
        """
        single-cutoff training using autogluon
        :return:
        """
        if self._split_factor == 0:
            train_x = self._ecifp_train_data[list(self._ecifp_train_data.columns)[1:-1]]
            train_y = self._ecifp_train_data['pk']
            test_x = self._ecifp_test_data[list(self._ecifp_train_data.columns)[1:-1]]
            test_y = self._ecifp_test_data['pk']
        else:
            X = self._ecifp_data[list(self._ecifp_data.columns)[1:-1]]
            Y = self._ecifp_data['pk']
            train_x, test_x, train_y, test_y = train_test_split(X,
                                                                Y,
                                                                test_size=self._split_factor
                                                                if 0 < self._split_factor < 1 else 0.01465823,
                                                                random_state=0)

        print("training...\n")
        ag_r2_scorer = make_scorer(name='r2',
                                   score_func=sklearn.metrics.r2_score,
                                   optimum=1,
                                   greater_is_better=True)
        predictor = TabularPredictor(label='pk',
                                     eval_metric='root_mean_squared_error').fit(train_data=train_x.join(train_y))
        print("\n - finished -\n")
        leaderboard = predictor.leaderboard(test_x.join(test_y), extra_metrics=[ag_r2_scorer], silent=True)
        print(leaderboard[['model', 'score_val', 'r2']])


class ECIF_Trainer:

    _ecif_data = None
    _ecif_train_data = None
    _ecif_test_data = None
    _split_factor = 0
    _distance_cutoff = 6.0
    _using_binding_data = False

    def __init__(self, distance_cutoff=6.0, splitted=True, using_binding_data=False, tag=""):
        """
        init
        :param distance_cutoff:
        :param splitted: if data is already splitted
        :param  using_binding_data: use original ecif data generation and set-splitting method.
                In this case, splitted should be False.
                Mention that BindingData.csv should be provided
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        """
        self._using_binding_data = using_binding_data
        if splitted:
            if tag is None or tag == "":
                self._ecif_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIF_train_{}.csv".format(distance_cutoff)))
                self._ecif_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIF_test_{}.csv".format(distance_cutoff)))
            else:
                self._ecif_train_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIF_train_{}_{}.csv".format(distance_cutoff, tag)))
                self._ecif_test_data = TabularDataset(
                    os.path.join(meif_data_dir, "ECIF_test_{}_{}.csv".format(distance_cutoff, tag)))
            self._split_factor = 0
        else:
            if self._using_binding_data:
                self._ecif_data = pd.read_csv(os.path.join(meif_data_dir, "ECIF_6.0.csv"))
                ld = pd.read_csv(os.path.join(meif_data_dir, "LD.csv"))
                binding_data = pd.read_csv(os.path.join(data_dir, "BindingData.csv"))
                self._ecif_data = self._ecif_data.merge(ld, left_on="PDB", right_on="PDB")
                self._ecif_data = self._ecif_data.merge(binding_data, left_on="PDB", right_on="PDB")
                self._split_factor = 0
            else:
                if tag is None or tag == "":
                    self._ecif_data = TabularDataset(os.path.join(meif_data_dir, "ECIF_{}.csv".format(distance_cutoff)))
                else:
                    self._ecif_data = TabularDataset(
                        os.path.join(meif_data_dir, "ECIF_{}_{}.csv".format(distance_cutoff, tag)))
                self._split_factor = 0.01465823

        self._distance_cutoff = distance_cutoff

    # callables

    def set_split_factor(self, split_factor):
        """
        set a new split factor
        :param split_factor: size(test) / size(test+train)
        :return:
        """
        self._split_factor = split_factor

    def train(self):
        """
        single-cutoff training using autogluon
        :return:
        """
        if self._split_factor == 0:
            if self._using_binding_data:
                train_x = self._ecif_data[self._ecif_data["SET"] == "Train"][list(self._ecif_data.columns)[1:-2]]
                train_y = self._ecif_data[self._ecif_data["SET"] == "Train"]["pK"]
                test_x = self._ecif_data[self._ecif_data["SET"] == "Test"][list(self._ecif_data.columns)[1:-2]]
                test_y = self._ecif_data[self._ecif_data["SET"] == "Test"]["pK"]
            else:
                train_x = self._ecif_train_data[list(self._ecif_train_data.columns)[1:-1]]
                train_y = self._ecif_train_data['pk']
                test_x = self._ecif_test_data[list(self._ecif_train_data.columns)[1:-1]]
                test_y = self._ecif_test_data['pk']
        else:
            X = self._ecif_data[list(self._ecif_data.columns)[1:-1]]
            Y = self._ecif_data['pk']
            train_x, test_x, train_y, test_y = train_test_split(X,
                                                                Y,
                                                                test_size=self._split_factor
                                                                if 0 < self._split_factor < 1 else 0.01465823,
                                                                random_state=0)

        print("training...\n")
        ag_r2_scorer = make_scorer(name='r2',
                                   score_func=sklearn.metrics.r2_score,
                                   optimum=1,
                                   greater_is_better=True)
        predictor = TabularPredictor(label='pK',
                                     eval_metric='root_mean_squared_error').fit(train_data=train_x.join(train_y))
        print("\n - finished -\n")
        leaderboard = predictor.leaderboard(test_x.join(test_y), extra_metrics=[ag_r2_scorer], silent=True)
        print(leaderboard[['model', 'score_val', 'r2']])


if __name__ == '__main__':
    trainer = ECIF_Trainer(distance_cutoff=6.0, splitted=False, using_binding_data=True, tag="")
    trainer.train()

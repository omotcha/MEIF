"""
platform: win
env: any with autogluon
name: train_with_ag.py
do training on MEIF tabular data using autogluon model generator
"""
import os
from autogluon.tabular import TabularDataset, TabularPredictor
from configs.config import meif_data_dir
from sklearn.model_selection import train_test_split
from autogluon.core.metrics import make_scorer
import sklearn.metrics


class MEIF_Trainer:

    _meif_data = None
    _split_factor = 0.01465823
    _distance_cutoff = 6.0

    def __init__(self, distance_cutoff=6.0, tag=""):
        """

        :param distance_cutoff:
        :param tag: added a postfix to generated csv file name, remember not to add diffusing symbols like ':'
        """
        if tag is None or tag == "":
            self._meif_data = TabularDataset(os.path.join(meif_data_dir, "MEIF_{}.csv".format(distance_cutoff)))
        else:
            self._meif_data = TabularDataset(os.path.join(meif_data_dir, "MEIF_{}_{}.csv".format(distance_cutoff, tag)))
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
        X = self._meif_data[list(self._meif_data.columns)[1:-1]]
        Y = self._meif_data['pk']
        train_x, test_x, train_y, test_y = train_test_split(X,
                                                            Y,
                                                            test_size=self._split_factor if 0 < self._split_factor < 1
                                                            else 0.01465823, random_state=0)
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


if __name__ == '__main__':
    trainer = MEIF_Trainer(6.0, "Test")
    trainer.train()

"""
platform: any
env: any
name: AGTrainTest.py
Autogluon Trainer Tester
"""

from train.train_with_ag import MEIF_Trainer


def testTrainMEIF():
    trainer = MEIF_Trainer(6.0, "Test")
    trainer.train()


if __name__ == '__main__':
    testTrainMEIF()

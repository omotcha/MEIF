## Experiments
#### 1. init
- Tabular data generation time: about 4h
- Models: ag-20220825_020839
- Leaderboard
```angular2html
                  model  score_val        r2
0   WeightedEnsemble_L2  -1.169632  0.545061
1            LightGBMXT  -1.190447  0.530040
2         LightGBMLarge  -1.184534  0.527286
3              CatBoost  -1.193233  0.513538
4         ExtraTreesMSE  -1.229752  0.500113
5              LightGBM  -1.205486  0.494610
6       RandomForestMSE  -1.252995  0.494095
7               XGBoost  -1.226522  0.492893
8       NeuralNetFastAI  -1.298570  0.446091
9        NeuralNetTorch  -1.336706  0.390651
10       KNeighborsDist  -1.420279  0.287897
11       KNeighborsUnif  -1.464124  0.249150
```
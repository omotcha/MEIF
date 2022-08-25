## Experiments
#### 1. init(MEIF with no optims)
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

#### 2. ECIFP(ECIF+escape labels+pdbbind2020)
- Tabular data generation time: about 35min
- Models: ag-20220825_085255
- Leaderboard
```angular2html
                  model  score_val        r2
0   WeightedEnsemble_L2  -1.148280  0.627385
1              CatBoost  -1.172127  0.614037
2         LightGBMLarge  -1.164864  0.601690
3            LightGBMXT  -1.186200  0.592328
4              LightGBM  -1.188535  0.589309
5       NeuralNetFastAI  -1.313653  0.587988
6         ExtraTreesMSE  -1.216072  0.583365
7               XGBoost  -1.195175  0.582069
8       RandomForestMSE  -1.223338  0.572278
9        KNeighborsDist  -1.457172  0.450282
10       NeuralNetTorch  -1.398399  0.449546
11       KNeighborsUnif  -1.497640  0.392478
```
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
- Models: ag-20220902_011835
- Leaderboard
```angular2html
                  model  score_val        r2
0              LightGBM  -1.153731  0.640571
1         LightGBMLarge  -1.134660  0.628050
2              CatBoost  -1.160307  0.625722
3               XGBoost  -1.169528  0.615295
4            LightGBMXT  -1.153279  0.610849
5   WeightedEnsemble_L2  -1.117438  0.599059
6         ExtraTreesMSE  -1.182929  0.561037
7       RandomForestMSE  -1.193287  0.552165
8        NeuralNetTorch  -1.407867  0.473484
9        KNeighborsDist  -1.371300  0.301969
10       KNeighborsUnif  -1.434577  0.290977
11      NeuralNetFastAI  -1.295202 -0.180418
```

#### 3. ECIF
- Tabular data generation time: about 14min
- Models: ag-20220902_034242
- Leaderboard
```angular2html
                  model  score_val        r2
0            LightGBMXT  -1.186731  0.716315
1              CatBoost  -1.179754  0.700723
2              LightGBM  -1.179519  0.699779
3   WeightedEnsemble_L2  -1.148319  0.694330
4               XGBoost  -1.209543  0.675154
5         LightGBMLarge  -1.166375  0.668868
6         ExtraTreesMSE  -1.221901  0.612400
7       NeuralNetFastAI  -1.273477  0.607187
8       RandomForestMSE  -1.233089  0.601882
9        NeuralNetTorch  -1.433347  0.537436
10       KNeighborsDist  -1.494613  0.433635
11       KNeighborsUnif  -1.562092  0.391596
```

#### 4. ECIF AUG
- GBT
```angular2html
Pearson correlation coefficient for GBT:  0.8379730624687073
RMSE for GBT: 1.2764862385011655
```
- Models: ag-20220916_014959
- Leaderboard
```angular2html
                  model  score_val        r2
0         LightGBMLarge  -0.861528  0.667936
1   WeightedEnsemble_L2  -0.853553  0.667059
2              LightGBM  -0.878795  0.657715
3               XGBoost  -0.885987  0.655976
4            LightGBMXT  -0.902830  0.655548
5              CatBoost  -0.937753  0.651557
6         ExtraTreesMSE  -0.956337  0.554002
7       RandomForestMSE  -0.961026  0.540696
8        KNeighborsDist  -1.130164  0.411456
9        NeuralNetTorch  -1.284208  0.389558
10       KNeighborsUnif  -1.160719  0.354650
11      NeuralNetFastAI  -1.005798 -0.339483

```

- Predict Time Used(predict on core2016 #item: 285)
```angular2html
ECIF::GBT: 38 seconds
ECIF::CatBoost: 46 seconds
ECIF::LightGBMXT: 45 seconds
ECIFP::GBT: 39 seconds
ECIFP::CatBoost: 45 seconds
ECIFP::LightGBM: 46 seconds
ECIF_AUG::GBT: 64 seconds
ECIF_AUG::LightGMBLarge: 47 seconds
```

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
- Models: ag-20220901_020745
- Leaderboard
```angular2html
                  model  score_val        r2
0   WeightedEnsemble_L2  -1.105918  0.666506
1              LightGBM  -1.138972  0.650136
2            LightGBMXT  -1.156025  0.647661
3         LightGBMLarge  -1.125751  0.645971
4              CatBoost  -1.148192  0.637515
5               XGBoost  -1.166321  0.617023
6       RandomForestMSE  -1.195796  0.611476
7         ExtraTreesMSE  -1.174478  0.603487
8       NeuralNetFastAI  -1.246362  0.599360
9        NeuralNetTorch  -1.378830  0.499498
10       KNeighborsDist  -1.405895  0.475636
11       KNeighborsUnif  -1.475484  0.426464
```
- Predict Time Used(predict on core2016 #item: 285)
```angular2html
catboost: 52 seconds
lightGBM: 73 seconds
```

#### 3. ECIFP x CASF2016 Result compared with ECIF x CASF 2016 Result

- scoring power

Pearson correlation coefficient (R)

Standard deviation in fitting (SD)
```angular2html
                        R      SD
ECIF                0.863    1.10
ECIFP::catboost     0.808    1.28
ECIFP::lightGBM     0.806    1.29
```

- ranking power

The Spearman correlation coefficient (SP) 

The Kendall correlation coefficient (tau) 

The Predictive index (PI) 

```angular2html
                        SP      tau       PI
ECIF                 0.758    0.674    0.790
ECIFP::catboost      0.661    0.575    0.685
ECIFP::lightGBM      0.644    0.558    0.669
```
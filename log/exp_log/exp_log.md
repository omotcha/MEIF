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

- Predict Time Used(predict on core2016 #item: 285)
```angular2html
ECIF::GBT: 38 seconds
ECIF::CatBoost: 46 seconds
ECIF::LightGBMXT: 45 seconds
ECIFP::GBT: 39 seconds
ECIFP::CatBoost: 45 seconds
ECIFP::LightGBM: 46 seconds
```

#### 4. ECIFP x CASF2016 Result compared with ECIF x CASF 2016 Result

- scoring power

Pearson correlation coefficient (R)

Standard deviation in fitting (SD)
```angular2html
                        R      SD
ECIF::GBT           0.839    1.18
ECIF::CatBoost      0.823    1.24
ECIF::LightGBMXT    0.840    1.18
ECIFP::GBT          0.828    1.22
ECIFP::CatBoost     0.812    1.27
ECIFP::LightGBM     0.826    1.23
```

- ranking power

The Spearman correlation coefficient (SP) 

The Kendall correlation coefficient (tau) 

The Predictive index (PI) 

```angular2html
                        SP      tau       PI
ECIF::GBT            0.711    0.625    0.732
ECIF::CatBoost       0.674    0.565    0.706
ECIF::LightGBMXT     0.691    0.586    0.728
ECIFP::GBT           0.735    0.660    0.755
ECIFP::CatBoost      0.684    0.596    0.705
ECIFP::LightGBM      0.663    0.575    0.679
```

- docking power

number of Top 1/2/3 correct binding poses or success rate(#/SRN) N:1..3

The Spearman correlation coefficient in rmsd range 0..N(SPN) N:2..10


```angular2html
                       SR1      SR2      SR3    |     SP2      SP3      SP4      SP5      SP6      SP7      SP8      SP9     SP10
ECIF::GBT            0.814    0.853    0.877    |   0.222    0.230    0.236    0.246    0.256    0.265    0.281    0.302    0.318
ECIF::CatBoost       0.747    0.821    0.853    |   0.202    0.192    0.196    0.191    0.196    0.203    0.215    0.239    0.255
ECIF::LightGBMXT     0.740    0.825    0.849    |   0.180    0.196    0.195    0.189    0.200    0.206    0.222    0.242    0.256
ECIFP::GBT           0.807    0.856    0.902    |   0.220    0.228    0.239    0.259    0.283    0.301    0.317    0.339    0.358
ECIFP::CatBoost      0.818    0.870    0.912    |   0.224    0.222    0.233    0.252    0.272    0.292    0.312    0.337    0.353
ECIFP::LightGBM      0.698    0.779    0.818    |   0.194    0.189    0.182    0.195    0.215    0.230    0.245    0.265    0.283
```

time used
```angular2html
                   Time(seconds)
ECIF::GBT                   1895
ECIF::CatBoost              2668
ECIF::LightGBMXT            2383
ECIFP::GBT                  2022
ECIFP::CatBoost             2880
ECIFP::LightGBM             2568
```

- screening power

Average enrichment factor amongst Top 1%/5%/10% (EF1/5/10)

Success rate amongst Top 1%/5%/10% (SR1/5/10)
```angular2html
                       EF1      EF5     EF10    |     SR1      SR5     SR10
ECIF::GBT               -        -       -      |      -        -       -
ECIF::CatBoost          -        -       -      |      -        -       -
ECIF::LightGBMXT      1.36     1.77     1.37    |   0.053    0.193    0.281
ECIFP::GBT              -        -       -      |      -        -       -
ECIFP::CatBoost         -        -       -      |      -        -       -
ECIFP::LightGBM         -        -       -      |      -        -       -
```

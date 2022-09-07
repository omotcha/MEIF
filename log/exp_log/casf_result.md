## CASF Results

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

Success rate amongst Top 1%/5%/10% (reverse screening) (RSR1/5/10)
```angular2html
                       EF1      EF5     EF10    |     SR1      SR5     SR10    |    RSR1     RSR5    RSR10
ECIF::GBT             1.36     1.77     1.37    |   0.053    0.193    0.281    |   0.035    0.091    0.200
ECIF::CatBoost          -        -       -      |      -        -       -      |     -        -        -
ECIF::LightGBMXT        -        -       -      |      -        -       -      |     -        -        -
ECIFP::GBT              -        -       -      |      -        -       -      |     -        -        -
ECIFP::CatBoost       1.46     1.29     1.45    |   0.035    0.140    0.298    |   0.046    0.123    0.189
ECIFP::LightGBM         -        -       -      |      -        -       -      |     -        -        -
```

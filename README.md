# MEIF
Meshed Element-Type Interaction Feature. Something like ECIF.

## References

 - ECIF https://github.com/DIFACQUIM/ECIF

## Structure
```text
├── MEIF
    ├── casf                        //CASF analysis
        ├── casf_support.py         //CASF supporter
    ├── configs
        ├── config.py               //project configuration
    ├── data                        //data processing
        ├── AffinityDataMngr.py     //create affinity data file from dataset
        ├── (affinity_data.csv)     //generated by AffinityDataMngr.py, needed for ECIF data preparation
        ├── keys.csv                //keys for exchanging PDB Atom Type with ECIF(P) Atom Type
    ├── log                         //logs
        ├── exp_log                 //for experiments
        ├── screening               //for CASF screening
    ├── predict                     //predictions
        ├── ecif_predict.py         //ECIF predictions
        ├── ecifp_predict.py        //ECIFP predictions
    ├── preprocess
        ├── ecif_collector.py       //collect ECIF fingerprints of all protein-ligand pairs
        ├── ecifp_collector.py      //collect ECIFP fingerprints of all protein-ligand pairs
        ├── meif_collector.py       //collect MEIF fingerprints of all protein-ligand pairs
    ├── test                        //tests
    ├── tmp
    ├── train
        ├── AutogluonModels         //Autogluon saved models
        ├── train_with_ag.py        //train MEIF using autogluon
        ├── train_with_sklearn.py   //train ECIF using sklearn
    ├── util
        ├── ECIF.py                 //ECIF helper class, only for reference and comparative experiments
        ├── ECIFP.py                //ECIFP helper class
        ├── MEIF.py                 //MEIF helper class
        ├── RDKitHelper.py          //extend RDKit functionality
        ├── CASFScreeningHelper.py  //CASF screening analysis supporter
```

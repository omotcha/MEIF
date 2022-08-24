# MEIF
Meshed Element-Type Interaction Feature. Something like ECIF.

## References

 - ECIF https://github.com/DIFACQUIM/ECIF

## Structure
```text
├── MEIF
    ├── configs
        ├── config.py               //project configuration
    ├── data                        //data processing
        ├── AffinityDataMngr.py     //create affinity data file from dataset
        ├── (affinity_data.csv)     //generated by AffinityDataMngr.py, needed for ECIF data preparation
    ├── models
    ├── preprocess
        ├── meif_collector.py       //collect MEIF fingerprints of all protein-ligand pairs
    ├── test                        //tests
    ├── tmp
    ├── train
    ├── util
        ├── MEIF.py                 //MEIF helper class
        ├── ECIF.py                 //ECIF helper class, only for reference and comparative experiments
```
# 2021-NERTHUS-system-dynamics
The Nonproprietary Educative Reactor model with THorium and Uranium fuel Salts, or NERTHUS, is a model of a thermal spectrum molten salt reactor (MSR) written in Serpent 2 with an accompanying python API designed to make testing various fuel cycles simple. This repository contains the neutronics and fuel cycle scripts for the model. The accompanying system dynamics model for NERTHUS can be found [here](https://github.com/ondrejch/2021-NERTHUS-core). 

## deck.py

The deck.py script contains the model for NERTHUS. It is designed to write and run the Serpent input, and easily obtain the results from each run.

#### Quick Start

An example of running the deck, and printing k<sub>eff</sub> with its error. 
```python
from deck import serpDeck

# run the deck, get the results, print k_eff and its error
nert = serpDeck()
nert.full_build_run()
nert.get_results()
print(f"k_eff = {nert.k[0]} ± {nert.k[1]")
```

Running a depletion calculation starting with a FLiBe fuel salt and Refueling with a NaBe fuel salt.

```python
from deck import serpDeck

nert = serpDeck(fuel_salt='flibe', enr=0.02, refuel_salt='nabe', enr_ref=0.1, refuel=True)
nert.full_build_run()
```

### Object Attributes

| Variable     | Default Value | Place Holder                            |
| ------------ | ------------- | ------------ |
| self.deck_name | 'nerthus' | Serpent input file name |
| self.qsub_name | 'run.sh' | Shell file name to run the deck |
| self.nuc_libs | 'ENDF7' | Nuclear Data Library |


```python
self.fs_lib = '09c'     # Cross section temperature selection for fuel salt
```







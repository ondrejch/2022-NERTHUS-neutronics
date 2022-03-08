# 2021-NERTHUS-system-dynamics
The Nonproprietary Educative Reactor model with THorium and Uranium fuel Salts, or NERTHUS, is a model of a thermal spectrum molten salt reactor (MSR) written in Serpent 2 with an accompanying python API designed to make testing various fuel cycles simple. This repository contains the neutronics and fuel cycle scripts for the model. The accompanying system dynamics model for NERTHUS can be found [here](https://github.com/ondrejch/2021-NERTHUS-core). 

## deck.py

The deck.py script contains the model for NERTHUS. It is designed to write and run the Serpent input, and easily obtain the results from each run.

### Quick Start

```python
from deck import serpDeck

nert = serpDeck()
nert.full_build_run()
nert.get_results()
print(f"k_eff = {nert.k[0]} ± {nert.k[1]")
```








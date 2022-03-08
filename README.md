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

The object attributes for the NERTHUS model

```python
self.deck_name         = 'nerthus'                          # Serpent input file name
self.qsub_name         = 'run.sh'                           # Shell file name which runs SERPENT
self.nuc_libs          = 'ENDF7'                            # Nuclear data library
self.fs_lib            = '09c'                              # Cross section temperature selection for fuel salt
self.gr_lib            = '09c'                              # XS temp. selection for graphite
self.lib               = '09c'                              # XS temp. selection for other materials
self.histories         = 20000                              # Number of histories to run per generation
self.ngen              = 200                                # Number of active generations
self.nskip             = 60                                 # Number of inactive generations
self.queue             = 'fill'                             # NECluster torque queue ('local' to run on your machine)
self.ompcores          = 8                                  # OMP cores used when running SERPENT
self.memory            = 20                                 # Memory in GB requested for node
self.thermal_expansion = True                               # Bool to include thermal expansion; if False, reactor is modeled at 900K
self.refuel            = refuel                             # Bool to run burnup calculation
self.feedback          = False                              # Bool to use materials card or restart file
self.restart_file      = self.deck_name+".wrk"              # Name of restart file
self.feedback_index    = 0                                  # index of burnstep to read material definitions from
self.do_plots          = False                              # Bool to plot core
self.deck_path         = os.getcwd() + f'/{self.deck_name}' # Directory where SERPENT is ran
self.add_to_deck       = ""                                 # Additional Serpent inputs you want to add to the deck
self.burn_steps        = [[2, 0.0208], [1, 0.9584], [1, 2], [1, 4], [22, 7], [44, 30]] # depletion steps  


```







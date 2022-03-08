# 2021-NERTHUS-system-dynamics
The Non-proprietary Educative Reactor model with THorium and Uranium fuel Salts, or NERTHUS, is a model of a thermal spectrum molten salt reactor (MSR) written in Serpent 2 with an accompanying python API designed to make testing various fuel cycles simple. This repository contains the neutronics and fuel cycle scripts for the model. The accompanying system dynamics model for NERTHUS can be found [here](https://github.com/ondrejch/2021-NERTHUS-core). 

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

```python
self.deck_name         = 'nerthus'                          # Serpent input file name
self.qsub_name         = 'run.sh'                           # Shell file name which runs SERPENT
self.nuc_libs          = 'ENDF7'                            # Nuclear data library
self.fs_lib            = '09c'                              # Cross section temperature selection for fuel salt
self.gr_lib            = '09c'                              # Cross section temperature selection for graphite
self.lib               = '09c'                              # Cross section temperature selection for other materials
self.histories         = 20000                              # Number of histories to run per generation
self.ngen              = 200                                # Number of active generations
self.nskip             = 60                                 # Number of inactive generations
self.queue             = 'fill'                             # Torque queue ('local' to run on your machine)
self.ompcores          = 8                                  # OMP cores used when running SERPENT
self.memory            = 20                                 # Memory in GB requested for node
self.thermal_expansion = True                               # Bool to include thermal expansion; if False, reactor is modeled at 900K
self.refuel            = refuel                             # Bool to run burnup calculation
self.feedback          = False                              # Bool to use restart file
self.restart_file      = self.deck_name+".wrk"              # Name of restart file
self.feedback_index    = 0                                  # Index of burnup step to read material definitions from restart file
self.do_plots          = False                              # Bool to plot core
self.deck_path         = os.getcwd() + f'/{self.deck_name}' # Directory where SERPENT is ran
self.add_to_deck       = ""                                 # Additional Serpent inputs you want to add to the deck
self.burn_steps        = [[2, 0.0208], [1, 0.9584], [1, 2], [1, 4], [22, 7], [44, 30]] # depletion steps  
self.control_rods      = {0:0, 1:0, 2:0, 3:0}               # 0:center, 1:top, 2:bottom left, 3:bottom right; 0:fully removed, 1:fully inserted
self.fs_vol            = 13670000 if self.refuel else None  # Fuel salt volume if refueling
self.fs_dens_tempK     = 900.0                              # Fuel salt temp. used for density calc. [K]
self.fs_mat_tempK      = 900.0                              # Fuel salt temp. used for material XS [K]
self.mod_tempK         = 950.0                              # Graphite temp.
self.mod_boron         = 2e-6                               # Boron in graphite (2ppm default)
self.refuel_rate       = 1e-9                               # Refuel rate of the reactor
```

While most variable are self explanatory, some need more detail to be used.
`self.burn_steps` is a list of lists where the first value of the internal lists are the number of times that step will be run, and the second value is the time-step of each run in days.
For example `self.burn_steps = [[10, 1], [5, 30]]` will run 10 one day burn-up steps, and 5 thirty day burn-up steps totaling a 160 day depletion.

`self.control_rods` is a dictionary with each key, value pair corresponding to a control rod.
When the value is 0, the control rod is removed.
When the value is 1, the control rod is inserted.
There is not yet support for partial insertion of control rods in the NERTHUS model.

### Object Methods


| Method             | Purpose                                                                           |
| ------------------ | --------------------------------------------------------------------------------- |
| `save_deck()`      | Writes the Serpent Input to the specified directory                               |
| `save_qsub_file()` | Writes the shell file to run the model to the specified directory                 |
| `run_deck()`       | Runs the model                                                                    |
| `full_build_run()` | Saves the deck and shell file, and runs the deck (a combination of the previous 3)|
| `get_results()`    | stores the results of a run into object attributes                                |
| `cleanup()`        | removes the directory and all files in the directory the model was run            | 

Note: the `get_results()` method retrieves the criticality constant, neutron generation time and the delayed neutron fractions which are stored in `self.k`, `self.ngt`, `self.beta_tot`, and `self.betas`, respectfully. When running the model without depletion, the variables are a list with the first value being the value returned from Serpent, and the second value being the associated error. When run with depletion, each value is a list of lists where each index stores the value and error at each depletion step.

## burn.py
The burn.py script takes the NERTHUS model and uses it to calculate useful data pertaining to the fuel cycle.

### Quick Start

Calculating the critical enrichment, refuel rate and fuel salt feedback coefficients for the FLiBe-Nabe example from above.

```python
from burn import burn

nert = burn('flibe', 'nabe')
nert.get_enrichment()
nert.get_refuel_rate()
nert.get_feedbacks('fs.tot')
```
Note: The refuel rate and feedback calculations can take a while to run and submit lots of jobs, so if you are playing with the model, adjust `self.ngen`, `self.nskip` and `self.histories` to shorten the time it takes to run the model.


### Object Attributes

```python
self.fuel_salt:str = fuel_salt
self.refuel_salt:str = refuel_salt
self.queue:str = 'fill'
self.memory:int = 30
self.ompcores:int = 8
self.histories:int = 20000
self.ngen:int = 200
self.nskip:int = 60

 # Enrichment search varibles
self.enr_path:str = os.getcwd() + '/enr_search'
self.enr_min:float = 0.01
self.enr_max:float = 0.2
self.enr_eps:float = 1e-9
self.rho_tgt:float = 100.0
self.rho_eps:float = 100.0
self.conv_enr:float = None
self.RhoData = namedtuple("rhoData", 'enr rho rho_err')
self.rholist:list = []
self.iter_max:int  = 20

# refuel rate variales
self.refuel_path:str = os.getcwd() + '/refuel'
self.refuel_enr:float = .1
self.refuel_min:float = 1e-10
self.refuel_max:float = 1e-5
self.refuel_eps:float = 1e-9
self.k_diff_tgt:float = 0.003
self.k_diff_eps:float = 0.003
self.refuelData = namedtuple("refuelData", 'rate k k_err')
self.refuel_list:list = []
self.refuel_iter:int = 20

# feedback coefficient variables
self.feedback_path = os.getcwd() + '/feedback'
self.feedback_temps:list = [800.0, 850.0, 900.0, 950.0, 1000.0]
self.base_temp:float = 900.0
self.feedback_runs:dict = {}
self.burnup_steps:int = 72
self.smoothing_window:int = 11
```

### Object Method

















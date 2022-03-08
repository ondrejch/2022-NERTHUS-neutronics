#!/usr/bin/env python2
'''
Module for writing and running NERTHUS reactor using the SERPENT Monte Carlo code.

visit https://thorconpower.com/docs/exec_summary2.pdf to view TMSR-500 reactor concept
NERTHUS is inspired by.
'''

from salts import Salt
from textwrap import dedent
import math
import os
import serpentTools
import shutil



serpentTools.settings.rc['verbosity'] = 'error'

# Graphite Constants
GRAPHITE_CTE:float = 3.5e-6                    # Graphite linear expansion coefficient [m/m per K]
GRAPHITE_RHO:float = 1.80                      # Graphite density at 950 K [g/cm3]

# Dictionary of fuel salts and their composition
SALTS = {
    'thorConSalt'   : '76%NaF + 12%BeF2 + 9.5%ThF4 + 2.5%UF4',        #NaFBeTh12
    'thorCons_ref'  : '76%NaF + 12%BeF2 + 10.2%ThF4 + 1.8%UF4',       #NaFBeTh12
    'flibe'         : '72%LiF + 16%BeF2 + 12%UF4',                    #flibe
    'nabe'          : '76%NaF + 12%BeF2 + 12%UF4'
}


class serpDeck(object):
    '''
    class which writes NERTHUS core for SERPENT
    '''

    def __init__(self, fuel_salt:str='flibe', enr:float=0.17,
                    refuel_salt:str='flibe', enr_ref:float=0.20, refuel:bool=False) -> None:


        # Check is salt is defined
        try:
            self.salt_formula = SALTS[fuel_salt]
            self.salt_formula_r = SALTS[refuel_salt]
        except ValueError:
           ValueError("Salt "+fuel_salt+" is undefined.")

        # Initialize Reactor Parameters

        # Run parameters
        self.deck_name:str              = 'nerthus'             # SERPENT input file name
        self.qsub_name:str              = 'run.sh'              # Shell file name which runs SERPENT
        self.nuc_libs:str               = 'ENDF7'               # Nuclear data library
        self.fs_lib:str                 = '09c'                 # XS temp. selection for fuel salt
        self.gr_lib:str                 = '09c'                 # XS temp. selection for graphite
        self.lib:str                    = '09c'                 # XS temp. selection for other materials
        self.histories:int              = 20000                 # Number of histories to run per generation
        self.ngen:int                   = 200                   # Number of active generations
        self.nskip:int                  = 60                    # Number of inactive generations
        self.queue:str                  = 'fill'                # NECluster torque queue ('local' to run on your machine)
        self.ompcores:int               = 8                     # OMP cores used when running SERPENT
        self.memory:int                 = 20                    # Memory in GB requested for node
        self.thermal_expansion:bool     = True                  # Bool to include thermal expansion; if False, reactor is modeled at 900K
        self.refuel:bool                = refuel                # Bool to run burnup calculation
        self.feedback:bool              = False                 # Bool to use materials card or restart file
        self.restart_file:str           = self.deck_name+".wrk" # Name of restart file
        self.feedback_index:int         = 0                     # index of burnstep to read material definitions from
        self.do_plots:bool              = False                 # Bool to plot core
        self.deck_path:str              = os.getcwd() + f'/{self.deck_name}' # Directory where SERPENT is ran
        self.add_to_deck:str            = ""                    # Additional Serpent inputs you want to add to the deck
        self.burn_steps:list            = [[2, 0.0208], [1, 0.9584], [1, 2], [1, 4], [22, 7], [44, 30]]


        # Control rods: 0=removed, 1=fully inserted
        self.control_rods:dict              = {0:0, 1:0, 2:0, 3:0} # 0:center, 1:top, 2:bottom left, 3:bottom right; 0:fully removed, 1:fully inserted

        self.fs_vol:int = 13670000 if self.refuel else None     # Fuel salt volume if refueling
        self.fs_dens_tempK:float            = 900.0             # Fuel salt temp. used for density calc. [K]
        self.fs_mat_tempK:float             = 900.0             # Fuel salt temp. used for material XS [K]
        self.mod_tempK:float                = 950.0             # Graphite temp.
        self.mod_boron:float                = 2e-6              # boron in graphite (2ppm default)

        # Make fuel salt
        self.e                  = enr                           # Fuel salt enrichment
        self.salt_name          = fuel_salt                     # Fuel salt name
        self.fuel_salt          = Salt(self.salt_formula, self.e) # Fuel salt object (see salts.py)


        self.e_ref              = enr_ref                       # Refuel salt enrichment
        self.salt_name_ref      = refuel_salt                   # Refuel salt name
        self.refuel_salt        = Salt(self.salt_formula_r, self.e_ref) # Refuel salt object

        self.refuel_rate:float  = 1e-9

    def _make_ellipsoid(self, pos:list, axes:list, name:str=None):
        '''creates A B C D E F G H I J values for ellipsoid surface in SERPENT'''
        x, y, z = pos
        a, b, c = axes
        A = (1/(a**2))
        B = (1/(b**2))
        C = (1/(c**2))
        surface = f'surf {name} quadratic {A:.8E} {B:.8E} {C:.8E} 0 0 0 0 0 0 -1'
        translation = f'trans S {name} {x} {y} {z}'
        return surface, translation

    def _rotate(self, point:list, rotation:float):
        '''rotates a 2D point around 0,0'''
        x, y = point[0], point[1]
        x_rot = x * math.cos(math.radians(rotation)) - y * math.sin(math.radians(rotation))
        y_rot = x * math.sin(math.radians(rotation)) + y * math.cos(math.radians(rotation))
        return [x_rot, y_rot]
    
    def _translate(self, point:list, pos:list):
        '''Moves a point to a new location'''
        x, y = point[0], point[1]
        x_tran = pos[0]
        y_tran = pos[1]

        x = x + x_tran
        y = y + y_tran
        return [x,y]

    def _make_plane(self, point1:list, point2:list, name:str) -> str:
        x1, y1 = point1[0], point1[1]
        x2, y2 = point2[0], point2[1]
        plane = f'\nsurf {name} plane {x1:.8f} {y1:.8f} 0.0 {x2:.8f} {y2:.8f} 0.0 {x2:.8f} {y2:.8f} -1.0'
        return plane

    def _GLE(self, point = None) -> list:
        '''Method for calculating change in distance due to thermal expansion
        of graphite. If thermal expansion is excluded, the geometry is modeled at 900k '''
        if self.thermal_expansion == True:
            if type(point) is list:
                x0, y0 = point[0], point[1]
                xf = x0 * (1.0 + GRAPHITE_CTE * (self.mod_tempK - 293.0))
                yf = y0 * (1.0 + GRAPHITE_CTE * (self.mod_tempK - 293.0))
                result = [xf, yf]
            else:
                pf = point * (1.0 + GRAPHITE_CTE * (self.mod_tempK - 293.0))
                result = pf
            return result
        else:
            if type(point) is list:
                x0, y0 = point[0], point[1]
                xf = x0 * (1.0 + GRAPHITE_CTE * (606.0))
                yf = y0 * (1.0 + GRAPHITE_CTE * (607.0))
                result = [xf, yf]
            else:
                pf = point * (1.0 + GRAPHITE_CTE * (607.0))
                result = pf
            return result

    def _GDE(self) -> float:
        '''Return new density based on graphite thermal expansion'''
        if self.thermal_expansion:
            unit_f  = (1.0 + GRAPHITE_CTE * (self.mod_tempK - 950.0))
            rho_f   = GRAPHITE_RHO / unit_f**3
            return rho_f
        else:
            unit_f  = (1.0 + GRAPHITE_CTE * (-950.0))
            rho_f   = GRAPHITE_RHO / unit_f**3
            return rho_f

    def _make_surfs_and_cells(self) -> str:
        '''method for writing the surfaces and cells for the NERTHUS model'''

        surfs_and_cells_cards:str = dedent('''
            % ====================================== %
            % ===== NERTHUS Serpent input file ===== %
            % ====================================== %

            set title "NERTHUS"

            % =================================== %
            % ===== SURFACES AND CELL CARDS ===== %
            % =================================== %
            ''')

        pot_wall = dedent('''
            % --- STAINLESS STEEL WALL OF REACTOR --- %
            surf pot_wall_outer cyl 0.0 0.0 248.05
            surf pot_wall_inner cyl 0.0 0.0 243.05
            surf pot_wall_top pz  234.7
            surf pot_wall_bot pz -234.7

            % --- CELL FOR POT WALL --- %
            cell pot_wall 0 sus316_stainless_steel
            -pot_wall_outer pot_wall_inner pot_wall_bot -pot_wall_top
            ''')

        surfs_and_cells_cards += pot_wall

        pot_top_outer, pot_top_outer_trans = self._make_ellipsoid([0,0,239.7], [250.74, 250.74, 47.9], 'pot_top_outer')
        pot_top_inner, pot_top_inner_trans = self._make_ellipsoid([0,0,239.7], [243.05, 243.05, 42.9], 'pot_top_inner')

        pot_top = dedent(f'''
            % --- TOP WALL OF POT --- %
            % - wall of the top - %
            {pot_top_outer}
            {pot_top_outer_trans}
            {pot_top_inner}
            {pot_top_inner_trans}

            % - brim at the top of the reactor - %
            surf pot_brim_top pz 239.7
            surf pot_brim_outer cyl 0.0 0.0 254

            % - pipes for control rods and fuel salt - %
            surf pot_top_tube_out cyl 0.0 -50.5 20.5
            surf pot_top_tube_in  cyl 0.0 -50.5 15.5
            surf pot_top_ctrl_out cyl 0.0 0.0 24.5
            surf pot_top_ctrl_in  cyl 0.0 0.0 19.5
            surf pot_top_plane    pz 295

            % --- CELL FOR POT TOP WALL --- %
            cell pot_top 0 sus316_stainless_steel
            (-pot_top_outer pot_top_inner pot_brim_top pot_top_tube_out pot_top_ctrl_out):
            (-pot_brim_outer pot_wall_inner pot_wall_top -pot_brim_top):
            (-pot_top_tube_out pot_top_tube_in pot_top_inner -pot_top_plane pot_brim_top):
            (-pot_top_ctrl_out pot_top_ctrl_in pot_top_inner -pot_top_plane pot_brim_top)
            ''')

        surfs_and_cells_cards += pot_top

        pot_bot_outer, pot_bot_outer_trans = self._make_ellipsoid([0,0,-234.7], [248.05, 248.05, 49.4], 'pot_bot_outer')
        pot_bot_inner, pot_bot_inner_trans = self._make_ellipsoid([0,0,-234.7], [243.05, 243.05, 44.4], 'pot_bot_inner')

        pot_bot = dedent(f'''
            % --- BOTTOM WALL OF REACTOR --- %
            % - wall of the bottom - %
            {pot_bot_outer}
            {pot_bot_outer_trans}
            {pot_bot_inner}
            {pot_bot_inner_trans}

            % - tube for fuel salt - %
            surf pot_bot_tube_out cyl 0.0 0.0 24.8
            surf pot_bot_tube_in  cyl 0.0 0.0 19.8
            surf pot_bot_plane    pz -295

            % --- CELL FOR POT BOTTOM WALL --- %
            cell pot_bot 0 sus316_stainless_steel
            (-pot_bot_outer pot_bot_inner -pot_wall_bot pot_bot_tube_out):
            (-pot_bot_tube_out pot_bot_tube_in pot_bot_plane -pot_wall_bot pot_bot_inner)
            ''')

        surfs_and_cells_cards += pot_bot


        void = dedent('''
            % --- VOID CELLS --- %
            cell void_mid 0 outside
            pot_wall_outer -pot_wall_top pot_wall_bot

            cell void_up 0 outside
            pot_top_plane

            cell void_down 0 outside
            -pot_bot_plane

            cell void_brim 0 outside
            pot_brim_outer pot_wall_top -pot_brim_top

            cell void_top 0 outside
            pot_brim_top -pot_top_plane pot_top_tube_out pot_top_ctrl_out pot_top_outer

            cell void_bot 0 outside
            -pot_wall_bot pot_bot_plane pot_bot_tube_out pot_bot_outer
            ''')

        surfs_and_cells_cards += void


        shield_inner = self._GLE(230.0)
        shield_outer = self._GLE(240.0)
        shield_top = self._GLE(189.0)
        shield_bot = self._GLE(-189.0)

        shield = dedent(f'''
            % --- BORON CARBIDE SHIELD --- %
            % - shield surfaces - %
            surf shield_outer cyl 0.0 0.0 {shield_outer:.8f}
            surf shield_inner cyl 0.0 0.0 {shield_inner:.8f}
            surf shield_top   pz {shield_top:.8f}
            surf shield_bot   pz {shield_bot:.8f}

            % --- CELL FOR B4C SHIELD --- %
            cell shield 0 B4C_shield
            -shield_outer shield_inner -shield_top shield_bot
            ''')

        surfs_and_cells_cards += shield


        # Plug values
        plug_out = self._GLE(236.15)
        plug_top = self._GLE(-191.0)
        plug_mid = self._GLE(-217.36)
        plug_bot = self._GLE(-266.0)
        # Pretty rounded coreners
        plug_corner_cyl = self._GLE(226.15)
        plug_corner_plane = self._GLE(-201.0)
        plug_corner_rad = self._GLE(10.0)
        # Graphite dome
        plug_z = self._GLE(-217.36)
        plug_AB = self._GLE(236.15)
        plug_C  = self._GLE(53.3)
        plug_dome, plug_dome_trans = self._make_ellipsoid([0,0,plug_z], [plug_AB, plug_AB, plug_C], 'plug_dome')


        plug = dedent(f'''
            % --- GRAPHITE PLUG REFLECTOR --- %
            % - plug surfaces - %
            surf plug_out cyl 0.0 0.0 {plug_out:.8f}
            surf plug_top   pz  {plug_top:.8f}
            surf plug_bot   pz  {plug_bot:.8f}
            surf plug_mid   pz  {plug_mid:.8f}
            surf plug_corner_plane pz {plug_corner_plane:.8f}
            surf plug_corner_cyl cyl 0.0 0.0 {plug_corner_cyl:.8f}
            surf plug_corner torz 0.0 0.0 {plug_corner_plane:.8f} {plug_corner_cyl:.8f} {plug_corner_rad:.8f} {plug_corner_rad:.8f}
            {plug_dome}
            {plug_dome_trans}


            % --- CELL FOR PLUG GRAPHITE--- %
            cell plug plug graphite
            (-plug_top plug_corner_plane -plug_corner_cyl):
            (-plug_corner_plane plug_mid -plug_out):
            (-plug_dome -plug_mid plug_bot):
            (-plug_corner plug_corner_cyl plug_corner_plane)

            % --- CELL FOR PLUG FUEL SALT --- %
            cell plug_fuel_salt plug fuelsalt
            -plug_bot: plug_top:
            (plug_dome -plug_mid plug_bot):
            (plug_out plug_mid -plug_corner_plane):
            (plug_corner plug_corner_cyl plug_corner_plane)


            % --- CELL TO FILL PLUG --- %
            cell plug_fill 0 fill plug
            (-plug_top -pot_wall_inner pot_wall_bot):
            (-pot_bot_inner -pot_wall_bot)
            ''')

        surfs_and_cells_cards += plug

        hat_out = self._GLE(233.0)
        hat_bot = self._GLE(191.0)
        hat_mid = self._GLE(239.0)
        hat_top = self._GLE(268.2)
        # Pretty rounded corners
        hat_corner_cyl = self._GLE(223.0)
        hat_corner_plane = self._GLE(201.0)
        hat_corner_rad = self._GLE(10.0)
        # Graphite dome
        hat_z = self._GLE(239.0)
        hat_AB = self._GLE(233.0)
        hat_C = self._GLE(35.7)
        hat_dome, hat_dome_trans = self._make_ellipsoid([0.0,0.0,hat_z], [hat_AB, hat_AB, hat_C], 'hat_dome')

        hat = dedent(f'''
            % --- GRAPHITE HAT REFLECTOR --- %
            % - hat surfaces - %
            surf hat_out cyl 0.0 0.0 {hat_out:.8f}
            surf hat_top pz {hat_top:.8f}
            surf hat_bot pz {hat_bot:.8f}
            surf hat_mid pz {hat_mid:.8f}
            surf hat_corner_plane pz {hat_corner_plane:.8f}
            surf hat_corner_cyl cyl 0.0 0.0 {hat_corner_cyl:.8f}
            surf hat_corner torz 0.0 0.0 {hat_corner_plane} {hat_corner_cyl} {hat_corner_rad} {hat_corner_rad}
            {hat_dome}
            {hat_dome_trans}

            % --- CELL FOR HAT GRAPHITE --- %
            cell hat hat graphite
            (hat_bot -hat_corner_plane -hat_corner_cyl rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (hat_corner_plane -hat_mid -hat_out rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (hat_mid -hat_dome -hat_top rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (-hat_corner hat_corner_cyl -hat_corner_plane)

            % --- CELL FOR HAT FUELSALT --- %
            cell hat_fuel_salt hat fuelsalt
            (hat_top rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (-hat_bot rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (hat_dome hat_mid -hat_top):
            (hat_out -hat_mid hat_corner_plane):
            (hat_corner hat_corner_cyl -hat_corner_plane)
            

            % --- CELL TO FILL HAT --- %
            cell hat_fill 0 fill hat
            (hat_bot -pot_wall_inner -pot_brim_top):
            (pot_brim_top -pot_top_inner)
            ''')


        surfs_and_cells_cards += hat

        salt_cell = dedent('''
            % --- FUEL SALT CELL --- %
            cell fuel_salt 0 fuelsalt
            (shield_top shield_inner -hat_bot -pot_wall_inner):             % Above shield
            (shield_outer -shield_top shield_bot -pot_wall_inner):          % Outside shield
            (-shield_bot shield_inner plug_top -pot_wall_inner):            % Below shield
            (pot_bot_inner -pot_wall_bot -pot_bot_tube_in pot_bot_plane):   % In bottom pipe
            (pot_top_inner pot_brim_top -pot_top_tube_in -pot_top_plane):   % In top pipe
            (pot_top_inner pot_brim_top -pot_top_ctrl_in -pot_top_plane rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3) % In top control pipe
            ''')

        surfs_and_cells_cards += salt_cell

        # Logs and Guide rods

        log_top = self._GLE(189)
        log_bot = self._GLE(-189)

        log = dedent(f'''
            % --- CORE DEFINITION --- %
            % - top and bottom of logs
            surf log_top pz {log_top}
            surf log_bot pz {log_bot}
            ''')

        def make_slab(trans:list, name:str):
            '''Creates input for NERTHUS slab'''

            slab = dedent(f'''
                % --- {name.upper()} DEFINITION
                % - SURFACES''')

            slab_points = {
                0   : [ 0.0  ,  0.0      ],
                1   : [ 0.0  , -4.77     ],
                2   : [ 0.38 , -4.77     ],
                3   : [ 0.38 , -5.37     ],
                4   : [ 0.0  , -5.37     ],
                5   : [ 0.0  , -14.76    ],
                6   : [ 0.38 , -14.76    ],
                7   : [ 0.38 , -15.36    ],
                8   : [ 0.0  , -15.36    ],
                9   : [ 0.0  , -19.7     ],
                10  : [-3.384, -17.74624 ],
                11  : [-3.384, -17.836247],
                12  : [-3.984, -17.836247],
                13  : [-3.984,   1.864375],
                14  : [-0.6  ,  -0.09    ],
                15  : [-0.6  ,   0.0     ],
            }
            # Expand points from themal expansion then translate them to new location
            for p in slab_points:
                slab_points[p] = self._translate(self._GLE(slab_points[p]), trans)

            plane_names = []
            # List of points which share a common plane, so they are not made twice
            common_plane = [4, 6, 8]

            count = 1
            for i in slab_points:
                if i in common_plane: # Check to see if plane exists
                    pass
                elif i == 15: # Make the last plane
                    plane_name = name + '_plane' + str(count)
                    plane_names.append(plane_name)
                    slab += self._make_plane(slab_points[15], slab_points[0], plane_name)
                else:
                    # Make plane and keep track of names
                    plane_name = name + '_plane' + str(count)
                    plane_names.append(plane_name)
                    slab += self._make_plane(slab_points[i], slab_points[i+1], plane_name)
                    count += 1

            # Make cells for log, in universe 2
            gr_cell = dedent(f'''\n
                % - CELLS FOR {name.upper()}
                cell {name} 2 graphite
                (-{plane_names[0]} -{plane_names[6]} -{plane_names[9]} -{plane_names[10]} -log_top log_bot):
                ( {plane_names[0]} -{plane_names[1]} -{plane_names[2]} -{plane_names[3]} -log_top log_bot):
                ( {plane_names[0]} -{plane_names[4]} -{plane_names[2]} -{plane_names[5]} -log_top log_bot):
                ( {plane_names[6]} -{plane_names[7]} -{plane_names[8]} -{plane_names[9]} -log_top log_bot):
                ( {plane_names[10]} -{plane_names[11]} -{plane_names[12]} -{plane_names[0]} -log_top log_bot)

                ''')
            slab += gr_cell

            cell = dedent(f'''
                #(-{plane_names[0]} -{plane_names[6]} -{plane_names[9]} -{plane_names[10]} -log_top log_bot:
                 {plane_names[0]} -{plane_names[1]} -{plane_names[2]} -{plane_names[3]} -log_top log_bot:
                 {plane_names[0]} -{plane_names[4]} -{plane_names[2]} -{plane_names[5]} -log_top log_bot:
                 {plane_names[6]} -{plane_names[7]} -{plane_names[8]} -{plane_names[9]} -log_top log_bot:
                 {plane_names[10]} -{plane_names[11]} -{plane_names[12]} -{plane_names[0]} -log_top log_bot)''')

            return slab, cell

        def make_yoke(trans:list, rot:float, name:str):
            '''Writes yoke input for NERTHUS'''

            yoke = dedent(f'''
                % --- {name.upper()} DEFINITION
                % - SURFACES''')
            yoke_points = {
                0 : [ 0.0 ,   0.0],
                1 : [ 1.6 ,  -0.79],
                2 : [ 1.6 , -20.96],
                3 : [ 0.34, -21.59],
                4 : [-1.6 , -20.64],
                5 : [-1.6 ,  -0.79]
            }

            # Expand points from themal expansion then rotate and translate
            for i in yoke_points:
                yoke_points[i] = self._GLE(yoke_points[i])
                if rot != None:
                    yoke_points[i] = self._rotate(yoke_points[i], rot)
                if trans != None:
                    yoke_points[i] = self._translate(yoke_points[i], trans)

            plane_names = []

            for i in yoke_points:
                if i == 5:
                    plane_name = name + '_plane5'
                    plane_names.append(plane_name)
                    yoke += self._make_plane(yoke_points[5], yoke_points[0], plane_name)
                else:
                    plane_name = name + '_plane' + str(i)
                    plane_names.append(plane_name)
                    yoke += self._make_plane(yoke_points[i], yoke_points[i+1], plane_name)

            # Make cell for yoke
            cell = dedent(f'''\n
                % - CELLS FOR {name.upper()}
                cell {name} 2 graphite
                -{plane_names[0]} -{plane_names[1]} -{plane_names[2]}
                -{plane_names[3]} -{plane_names[4]} -{plane_names[5]}
                -log_top log_bot
                ''')
            yoke += cell

            cell = dedent(f'''
                #(-{plane_names[0]} -{plane_names[1]} -{plane_names[2]}
                -{plane_names[3]} -{plane_names[4]} -{plane_names[5]}
                -log_top log_bot)''')

            return yoke, cell

        # Make the log, yokes first
        yoke1_trans = self._GLE([0, -0.24])
        yoke2_trans = self._rotate(yoke1_trans, -120)

        yoke1, yoke_fs1 = make_yoke(yoke1_trans, 0.0, 'yoke1')
        yoke2, yoke_fs2 = make_yoke(yoke2_trans,-120.0, 'yoke2')

        log += yoke1 + yoke2

        # Make slabs
        slab1, slab_fs1 = make_slab(self._GLE([-1.99,-0.7]), 'slab1')
        slab2, slab_fs2 = make_slab(self._GLE([-6.36,1.8185]), 'slab2')
        slab3, slab_fs3 = make_slab(self._GLE([-10.74,4.3428]), 'slab3')
        slab4, slab_fs4 = make_slab(self._GLE([-15.125, 6.87]), 'slab4')


        log += slab1 + slab2 + slab3 + slab4

        # Make fuelsalt cell for log
        log += dedent('''
            % - LOG FUELSALT CELL
            cell log_salt 2 fuelsalt
            #(-guide_cyl log_top: -guide_cyl -log_bot)
            ''')

        log_fuelsalt_cells = yoke_fs1 + yoke_fs2 + slab_fs1 + slab_fs2 + slab_fs3 + slab_fs4

        log += log_fuelsalt_cells

        guide_rod_rad = self._GLE(4.45)

        # Make guide rods for log
        log += dedent(f'''
            % --- GUIDE ROD AT TOP AND BOTTOM OF LOG --- %
            % - SURFACES
            surf guide_cyl cyl 0 0 {guide_rod_rad}

            % - CELLS
            cell guide_rod 2 graphite
            -guide_cyl log_top: -guide_cyl -log_bot

            % - SET SYMMETRY
            set usym 2   3   3  0  0 150 120
            ''')
        surfs_and_cells_cards += log

        # Make the control rod log
        ctrl_half_width = self._GLE(19.055)

        ctrl_log = dedent(f'''
            % --- CONTROL LOG AND ROD DEFINITION
            % - SURFACES
            surf ctrl_hex hexxc 0 0 {ctrl_half_width:.8f}
            ''')

        gr_cell = dedent(f'''
            % - CELLS
            cell ctrl_log 7 graphite
            -ctrl_hex -log_top log_bot
            ''')

        fs_cell = dedent(f'''
            % - FUELSALT CELL
            cell ctrl_fs 7 fuelsalt
            (log_top -ctrl_hex rod_chan_0 rod_chan_1 rod_chan_2 rod_chan_3):
            (-log_bot -ctrl_hex):
            (log_bot -log_top ctrl_hex):
            (log_top ctrl_hex):
            (-log_bot ctrl_hex):
            #(''')

        # Values for control rod
        channel_radius:float = 0.525  # [cm] radius of small fuelsalt channels
        center_rod_radius:float = 5.3 # [cm] radius of center control rod
        outer_rod_radius:float = 6.0  # [cm] radius of outer control rod radius
        outer_distance:float = 13.3   # [cm] distance from center of log to center of outer control rods
        channel_distance:float = 2.5  # [cm] distance from center fuelsalt channel to center of outer fuelsalt channels

        # Make the small channels for cooling the log
        # 3 groups of channels with 7 holes per channel (center and 6 surrounding)
        for group in [1,2,3]:
            theta = group * 120.0 + 30.0
            x1 = self._GLE(outer_distance * math.cos(math.radians(theta)))
            y1 = self._GLE(outer_distance * math.sin(math.radians(theta)))
            name = f'group{group}cyl0'
            gr_cell += f' {name}'
            fs_cell += f' {name}'
            ctrl_log += f'surf {name} cyl {x1:.8f} {y1:.8f} {self._GLE(channel_radius):.8f}\n'
            for hole in [1,2,3,4,5,6]:
                phi = hole * 60.0 - 60.0
                x2 = self._GLE((channel_distance * math.cos(math.radians(phi))) + x1)
                y2 = self._GLE((channel_distance * math.sin(math.radians(phi))) + y1)
                name = f'group{group}cyl{hole}'
                gr_cell += f' {name}\n' if hole == 6 else f' {name}'
                fs_cell += f' {name}\n' if hole == 6 else f' {name}'
                ctrl_log += f'surf {name} cyl {x2:.8f} {y2:.8f} {self._GLE(channel_radius):.8f}\n'


        # Make the large channels where the control rods are inserted
        ctrl_channels = dedent('''
        % --- FILL CONTROL ROD CHANNELS --- %
        ''')
        name = 'rod_chan_0'
        gr_cell += f' {name}'
        ctrl_log += f'surf {name} cyl 0 0 {self._GLE(center_rod_radius):.8f}\n'

        if self.control_rods[0] == 0:
            ctrl_channels += dedent(f'''
                % - control log universe
                cell {name}_ctrl_log 7 fuelsalt
                log_bot -{name}

                % - hat universe
                cell {name}_hat hat fuelsalt
                -{name}

                % - base universe
                cell {name}_base 0 fuelsalt
                -pot_top_plane pot_brim_top pot_top_inner -{name}
                ''')
        else:
            ctrl_channels += dedent(f'''
                % - surfaces for control rod - %
                surf ctrl_rod0 cyl 0.0 0.0 4.7
                surf ctrl_arm0 cyl 0.0 0.0 3.0

                % - control log boron cell - %
                cell ctrl0_log 7 B4C_natural shield_bot -shield_top -ctrl_rod0
                % - control log steel cell - %
                cell ctrl0_log_arm 7 sus316_stainless_steel shield_top -ctrl_arm0
                % - control log salt cell - %
                cell ctrl0_log_salt 7 fuelsalt
                (-{name} ctrl_rod0 -shield_top shield_bot):
                (-shield_bot -{name} log_bot):
                (shield_top ctrl_arm0 -rod_chan_0)

                % - hat steel cell - %
                cell ctrl0_hat_arm hat sus316_stainless_steel -ctrl_arm0
                % - hat salt cell - %
                cell ctrl0_hat_salt hat fuelsalt -{name} ctrl_arm0 shield_bot

                % - base universe steel cell - %
                cell ctrl0_arm 0 sus316_stainless_steel -pot_top_plane -ctrl_arm0 pot_top_inner pot_brim_top
                % - base universe salt cell - %
                cell ctrl0_salt 0 fuelsalt -pot_top_plane ctrl_arm0 -rod_chan_0 pot_top_inner pot_brim_top

                ''')

        for hole in [1,2,3]:
            theta = hole * 120.0 - 30
            x = self._GLE(outer_distance * math.cos(math.radians(theta)))
            y = self._GLE(outer_distance * math.sin(math.radians(theta)))
            name = f'rod_chan_{hole}'
            gr_cell += f' {name}'
            ctrl_log += f'surf {name} cyl {x:.8f} {y:.8f} {self._GLE(outer_rod_radius):.8f}\n'
            # Add surfaces and cells for slat or control rods
            if self.control_rods[hole] == 0:
                ctrl_channels += dedent(f'''
                    % - control log universe
                    cell {name}_ctrl_log 7 fuelsalt
                    log_bot -{name}

                    % - hat universe
                    cell {name}_hat hat fuelsalt
                    -{name}

                    % - base universe
                    cell {name}_base 0 fuelsalt
                    -pot_top_plane pot_brim_top pot_top_inner -{name}
                    ''')
            else:
                ctrl_channels += dedent(f'''
                    % - surfaces for control rod - %
                    surf ctrl_rod{hole} cyl {x:.8f} {y:.8f} 4.7
                    surf ctrl_arm{hole} cyl {x:.8f} {y:.8f} 3.0
    
                    % - control log boron cell - %
                    cell ctrl{hole}_log 7 B4C_natural shield_bot -shield_top -ctrl_rod{hole}
                    % - control log steel cell - %
                    cell ctrl{hole}_log_arm 7 sus316_stainless_steel shield_top -ctrl_arm{hole}
                    % - control log salt cell - %
                    cell ctrl{hole}_log_salt 7 fuelsalt
                    (-{name} ctrl_rod{hole} -shield_top shield_bot):
                    (-shield_bot -{name} log_bot):
                    (shield_top ctrl_arm{hole} -rod_chan_{hole})
    
                    % - hat steel cell - %
                    cell ctrl{hole}_hat_arm hat sus316_stainless_steel -ctrl_arm{hole}
                    % - hat salt cell - %
                    cell ctrl{hole}_hat_salt hat fuelsalt -{name} ctrl_arm{hole} shield_bot
    
                    % - base universe steel cell - %
                    cell ctrl{hole}_arm 0 sus316_stainless_steel -pot_top_plane -ctrl_arm{hole} pot_top_inner pot_brim_top
                    % - base universe salt cell - %
                    cell ctrl{hole}_salt 0 fuelsalt -pot_top_plane ctrl_arm{hole} -rod_chan_{hole} pot_top_inner pot_brim_top

                ''')

        ctrl_log += gr_cell + '\n'
        ctrl_log += fs_cell + ')\n'

        surfs_and_cells_cards += ctrl_log
        surfs_and_cells_cards += ctrl_channels

        # Make lattice for full core
        ref_hex_half_width = self._GLE(20.055)

        lattice = dedent(f'''
            % --- SOLID GRAPHITE HEXAGON FOR OUTER REFLECTOR
            % - SURFACES
            surf ref_hex hexxc 0 0 {ref_hex_half_width:.8f}

            % - CELLS
            cell ref_hex 3 graphite -ref_hex -log_top log_bot

            % - FUELSALT
            cell ref_hex_fs 3 fuelsalt #(-ref_hex -log_top log_bot)
            ''')


        lattice_pitch = self._GLE(38.1104)

        lattice += dedent(f'''
            % --- LATTICE DEFINITION FOR WHOLE CORE
            %lat <uni> <type> 0 0 <nx> <ny> <p>
            lat lattice 2 0 0 17 17 {lattice_pitch:.8f}
            %0 1 2 3 4 5 6 7 c 9 0 1 2 3 4 5 6
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 2 2 2 2 3 3 3 3 %1
             3 3 3 3 3 3 3 2 2 2 2 2 2 2 3 3 3 %2
             3 3 3 3 3 3 2 2 2 2 2 2 2 2 3 3 3 %3
             3 3 3 3 3 2 2 2 2 2 2 2 2 2 3 3 3 %4
             3 3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 %5
             3 3 3 3 2 2 2 2 7 2 2 2 2 3 3 3 3 %c
             3 3 3 2 2 2 2 2 2 2 2 2 2 3 3 3 3 %7
             3 3 3 2 2 2 2 2 2 2 2 2 3 3 3 3 3 %8
             3 3 3 2 2 2 2 2 2 2 2 3 3 3 3 3 3 %9
             3 3 3 2 2 2 2 2 2 2 3 3 3 3 3 3 3 %10
             3 3 3 3 2 2 2 2 3 3 3 3 3 3 3 3 3 %11
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
            ''')

        surfs_and_cells_cards += lattice

        # Fill core with Lattice
        core = dedent(f'''
            % - CORE
            cell core 0 fill lattice -shield_inner -hat_bot plug_top
            ''')

        surfs_and_cells_cards += core

        return surfs_and_cells_cards

        # Wall of reactor
    def _make_mat_cards(self) -> str:
        '''Creates material definitions for SERPENT input'''
        material_cards = dedent('''
        % ===== Material Cards ===== %
        ''')

        # Stainless steel
        sus316_stainless_steel = dedent(f'''
            % sus316 stainless steel for outer walls
            mat sus316_stainless_steel -7.9 rgb 222 58 74
             6000.{self.lib} -0.0008   % Natural Carbon
             7014.{self.lib} -0.001    % Natural Nitrogen
            14000.{self.lib} -0.0075   % Natural Silicon
            15031.{self.lib} -0.00045  % Natural Phosphorous
            16000.{self.lib} -0.0003   % Natural Sulfur
            24000.{self.lib} -0.17     % Natural Chromium
            25055.{self.lib} -0.02     % Natural Manganese
            26000.{self.lib} -0.65495  % Natural Iron
            28000.{self.lib} -0.12     % Natural Nickel
            42000.{self.lib} -0.025    % Natural Molybdenum
            ''')
        material_cards += sus316_stainless_steel

        # Boron Carbide
        B4C_natural = dedent(f'''
            % Boron Carbide with naturally enriched boron
            mat B4C_natural -2.52 rgb 43 77 227
             5010.{self.lib} 0.16 % Boron 10
             5011.{self.lib} 0.64 % Boron 11
             6000.{self.lib} 0.2  % Natural Carbon
            ''')
        material_cards += B4C_natural

        # 10% Boron Carbide 90% graphite for shield
        B4C_shield = dedent(f'''
            % 10% Boron Carbide 90% graphite, used for shield
            mat B4C_shield -{0.9 * self._GDE() + 0.8 * 2.52} rgb 99 73 214
             5010.{self.lib} 0.16   % Boron 10
             5011.{self.lib} 0.64   % Boron 11
             6000.{self.lib} 0.92   % Natural Carbon
            ''')
        material_cards += B4C_shield

        # Graphite
        gr_frac  = 1 - self.mod_boron
        b10_frac = 0.2 * self.mod_boron
        b11_frac = 0.8 * self.mod_boron
        gr_dens = self._GDE()

        graphite = dedent(f'''
            % Graphite moderator
            mat graphite -{gr_dens:.8f} moder graph 6000
            tms {self.mod_tempK} rgb 59 59 59
             6000.{self.gr_lib} {gr_frac}   % Natural Carbon
             5010.{self.gr_lib} {b10_frac}  % Boron 10
             5011.{self.gr_lib} {b11_frac}  % Boron 11
            % Thermal scattering libraries for graphite
            therm graph 0 gre7.04t gre7.08t gre7.12t gre7.16t gre7.18t gre7.22t
            ''')
        material_cards += graphite

        fuelsalt = '\n' + self.fuel_salt.serpent_mat(self.fs_dens_tempK, self.fs_mat_tempK, \
                                            'fuelsalt', self.fs_lib, self.fs_vol, '54 227 167')
        material_cards += fuelsalt

        if self.refuel:
            refuel_salt = '\n' + self.refuel_salt.serpent_mat(self.fs_dens_tempK, self.fs_mat_tempK, \
                                          'refuelsalt', self.fs_lib, self.fs_vol, '54 227 167')
            material_cards += refuel_salt

        if self.feedback:
            material_cards += f"\nset rfr idx {self.feedback_index} {self.restart_file}\n"

        return material_cards

    def _make_data_cards(self) -> str:
        data_cards = dedent(f'''
            % ===== Data Cards ===== %
            set power 557000000 % Watts
            set pop {self.histories} {self.ngen} {self.nskip} % {self.histories} neutrons, {self.ngen} active cycles, {self.nskip} inactive cycles
            ''')

        if self.nuc_libs == 'ENDF7':
            data_cards += dedent('''
                % Data Libraries
                set acelib "sss_endfb7u.sssdir"
                set declib "sss_endfb7.dec"
                set nfylib "sss_endfb7.nfy"
                ''')
        else:
            print('Use ENDF7, or edit the source code to include other libraries')
            exit()

        if self.do_plots:
            data_cards += dedent('''
                % --- PLOTS
                plot 1 5000 5000 0 -300 300 -300 300
                %plot 2 5000 5000 0 -300 300 -300 300
                %plot 3 5000 5000 0 -300 300 -300 300

                ''')

        if self.refuel:
            data_cards += dedent(f'''
                % --- REPROCESSING CARDS

                set rfw 1 % write restart file

                % - OFFGAS TANK
                mat offgas -0.001 burn 1 vol 1e9 tmp {self.fs_mat_tempK}
                 2004.{self.fs_lib} 1

                % - OVERFLOW TANK
                mat overflow -0.001 burn 1 vol 1e9 tmp {self.fs_mat_tempK}
                2004.{self.fs_lib} 1

                % --- MASS FLOW DEFINITIONS

                set pcc 0 % predictor-corrector turned off for depletion

                mflow U_in
                all {self.refuel_rate}

                mflow offgasratecore
                Ne 1e-2
                Ar 1e-2
                He 1e-2
                Kr 1e-2
                Xe 1e-2
                Rn 1e-2

                % Account for increase in volume with refueling
                mflow over
                all {self.refuel_rate}

                rep source_rep
                rc refuelsalt fuelsalt U_in 0
                rc fuelsalt offgas offgasratecore 1
                rc fuelsalt overflow over 1

                % 4 years of burnup
                dep
                pro source_rep
                daystep
                ''')

            step_string = ""
            size = 0

            for (num, val) in self.burn_steps:
                for _ in range(num):
                    size += len(str(val)) + 1
                    if size > 45:
                        step_string += f"{val}\n"
                        size = 0
                    else:
                        step_string += f"{val} "
                    

            data_cards += dedent(step_string)
        return data_cards

    def _get_deck(self) -> str:
        deck  = self._make_surfs_and_cells()
        deck += self._make_mat_cards()
        deck += self._make_data_cards()
        if self.add_to_deck != None:
            deck += dedent(self.add_to_deck)
        return deck

    def save_deck(self) -> None:
        if self.fs_mat_tempK < 600:
            print("Error: Fuel Salt temperature is below 600 (your core is frozen solid)\nPlease increase \"self.fs_mat_tempK\"")
            quit()

        if self.mod_tempK > 3000:
            print("Error: Graphite temperature is above 3000K\n Please lower \"self.mod_tempK\"")
            quit()

        try:
            os.makedirs(self.deck_path, exist_ok=True)
            with open(self.deck_path + '/' + self.deck_name, 'w') as outfile:
                outfile.write(self._get_deck())
        except IOError as e:
            print('[ERROR] Unable to write file: ',
                self.deck_path + '/' + self.deck_name)
            print(e)

    def save_qsub_file(self) -> None:
        'Writes run file for TORQUE.'
        qsub_content = dedent(f'''
            #!/bin/bash
            #PBS -V
            #PBS -N NERTHUS
            #PBS -q {self.queue}
            #PBS -l nodes=1:ppn={self.ompcores}
            #PBS -l mem={self.memory}GB

            hostname
            rm -f done.dat
            cd ${{PBS_O_WORKDIR}}
            module load mpi
            module load serpent
            sss2 -omp {self.ompcores} {self.deck_name} > myout.out
            awk 'BEGIN{{ORS="\\t"}} /ANA_KEFF/ || /CONVERSION/ {{print $7" "$8;}}' {self.deck_name}_res.m > done.out
            rm {self.deck_name}.out
            ''')
        try:                # Write the deck
            f = open(self.deck_path + f'/{self.qsub_name}', 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to qsub file")
            print(e)

    def run_deck(self) -> None:
        'Runs the deck using qsub_path script'
        if self.queue == 'local':    # Run the deck locally
            os.chdir(self.deck_path)
            os.system(f'sss2 -omp {self.ompcores} {self.deck_name} > done.out')
            #os.system(f'rm {self.deck_name}.out')
            os.chdir('/..')
        else:               # Submit the job on the cluster
            os.system('cd ' + self.deck_path + f' && qsub {self.qsub_name}')

    def full_build_run(self) -> None:
        self.save_deck()
        self.save_qsub_file()
        self.run_deck()

    def get_results(self) -> bool:
        if os.path.exists(self.deck_path+'/done.out') and \
            os.path.getsize(self.deck_path+'/done.out') > 30:
            pass
        else:                   # Calculation not done yet
            return False
        # Get results
        if not self.refuel: # No refueling
            results = serpentTools.read(self.deck_path + '/' + self.deck_name + "_res.m")
            # Get k-eff value and error
            k = results.resdata['anaKeff'][0]
            k_err = results.resdata['anaKeff'][1] * k
            self.k  = [k, k_err]
            # Get neutron generation time and error
            ngt = results.resdata["adjNauchiGenTime"][0]
            ngt_err = results.resdata["adjNauchiGenTime"][1] * ngt
            self.ngt   = [ngt, ngt_err]
            # Get betas and error
            betas = results.resdata["adjNauchiBetaEff"]
            self.betas = []
            for i in range(len(betas)//2):
                if i == 0:
                    beta_tot = betas[i]
                    beta_tot_err = betas[i+1] * beta_tot
                    self.beta_tot = [beta_tot, beta_tot_err]
                else:
                    beta = betas[i*2]
                    beta_err = betas[i*2+1] * beta
                    self.betas.append([beta, beta_err])
            return True

        if self.refuel: # refueling
            results = serpentTools.read(self.deck_path + '/' + self.deck_name + "_res.m")
            burn_results = serpentTools.read(self.deck_path + '/' + self.deck_name + "_dep.m")
            # daysteps for burnup calc
            self.days = burn_results.days
            # K values
            self.k = [[k[0], k[1]*k[0]] for k in results.resdata['anaKeff']]
            # Neutron generation times for burnup calc
            self.ngt = [[n[0],n[1]*n[0]] for n in results.resdata['adjNauchiGenTime']]
            # Beta total for burnup calc
            self.beta_tot = [[b[0], b[1]*b[0]] for b in results.resdata['adjNauchiBetaEff']]
            # Betas for burnup calc
            betas = results.resdata['adjNauchiBetaEff']
            self.betas = []
            for i in range(len(self.days)):
                beta = []
                for j in range(len(betas[i])//2):
                    if j == 0:
                        pass
                    else:
                        b = betas[i][j*2]
                        b_err = betas[i][j*2 + 1] * b
                        beta.append([b,b_err])
                self.betas.append(beta)
            return True

    def cleanup(self, purge:bool=True):
        'Delete the run directories'
        if os.path.isdir(self.deck_path):
            if purge:
                shutil.rmtree(self.deck_path)
            else:
                with os.scandir(self.deck_path) as it:
                    for entry in it:
                        if entry.is_file():
                            os.remove(entry)



if __name__ == '__main__':
    test = serpDeck(fuel_salt='flibe', enr=0.2, refuel_salt='flibe', enr_ref=0.1, refuel=False)
    test.ngen = 80
    test.nskip = 20
    test.histories = 100
    test.queue = 'local'
    test.ompcores = 16
    some_string = """
            %Hello, this is a test
            %This is another test
            """
    test.add_to_deck = some_string
    test.save_deck()


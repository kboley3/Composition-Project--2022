
# This file is part of ExoPlex - a self consistent planet builder
# Copyright (C) 2017 - by the ExoPlex team, released under the GNU
# GPL v2 or later.

import os
import sys
import time
import matplotlib.font_manager as font_manager
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import FunctionsbyMass as fbm

start_time = time.time()

# hack to allow scripts to be placed in subdirectories next to exoplex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))


combine_phases = True
grids = True
known=False
Fname= 'CMB_P19000GPa_T20000K_FeMg_1p0_SiMg_1p0_v2'
fullOutput= False

Mass = 1.3 # in Earth masses
#create filename to store values

Output_filename = 'Test'
#Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
CaMg = 0.03
SiMg = 0.87
AlMg = 0.1
FeMg = 0.76

mlayers=1000
clayers=1200

pltname='MagmaPlot'

Planet, num_core_layers, num_mantle_layers, wt_frac_water=fbm.RunExoPlex(Mass,CaMg, SiMg, AlMg,FeMg,mlayers, clayers,output=Output_filename,file=Fname, grids=grids, known=known, fullOutput=fullOutput)
fbm.PlotProfile(Planet,num_core_layers, num_mantle_layers, wt_frac_water,'Save', verbose=True, melt=False)

print("--- %s seconds ---" % (time.time() - start_time))


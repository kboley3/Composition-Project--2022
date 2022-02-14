
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
from pdb import set_trace as b
start_time = time.time()

# hack to allow scripts to be placed in subdirectories next to exoplex:
if not os.path.exists('ExoPlex') and os.path.exists('../ExoPlex'):
    sys.path.insert(1, os.path.abspath('..'))

import ExoPlex as exo

def RunExoPlex( Mass,CaMg, SiMg, AlMg,FeMg,water, wlayers, mlayers, clayers, output='Test',file='CMB_P19000GPa_T20000K_FeMg_1p0_SiMg_1p0_v2', grids=True, known=False, fullOutput=False, Dry=True, magma=False):
    
    Pressure_range_mantle_UM = '1000 1400000'
    Temperature_range_mantle_UM = '1400 3000'
    
    Pressure_range_mantle_LM = '1000000 7500000'
    Temperature_range_mantle_LM = '2000 5000'
    
    core_rad_frac_guess = .3
    water_rad_frac_guess = 0.1
    water_potential_temp = 300
    
    combine_phases = True
    use_grids = grids
    knownfile=known
    Fname= file

    Mass_planet = Mass # in Earth masses
    #create filename to store values

    Output_filename = output
    #Next user must input the ratios by mole (Earth is Ca/Mg = .07, Si.Mg = 0.90, Al/Mg = 0.09, Fe/Mg = 0.9)
    


    #How much water do you want in your planet? By mass fraction.
    wt_frac_water = water

    #Don't forget that if you have water you need to add water layers
    number_h2o_layers = wlayers
    
    
    #Now we can mix various elements into the core or mantle
    wt_frac_Si_core = 0. #by mass <1
    wt_frac_O_core = 0. #by mass
    wt_frac_S_core = 0. #by mass
    mol_frac_Fe_mantle =0#by mole

    #What potential temperature (in K) do you want to start your mantle adiabat?
    Mantle_potential_temp = 2000.

    #Input the resolution of your upper mantle and lower mantle composition, density grids
    #These are input as number of T, P points. 50 50 = 2500 grid points, which takes about
    #5 minutes to calculate. Lower mantle resolution does not need to be higher since it's
    #mostly ppv.
    resolution_UM = '50 50'
    resolution_LM = '20 20'

    #lastly we need to decide how many layers to put in the planet. This is the resolution of
    #the mass-radius sampling.
    num_mantle_layers = mlayers
    num_core_layers = clayers



    number_of_runs = 100
    Output_radii = []
    Output_mass = []


    ######### Initalize and run ExoPlex


    compositional_params = [wt_frac_water,FeMg,SiMg,CaMg,AlMg,mol_frac_Fe_mantle,wt_frac_Si_core, \
                          wt_frac_O_core,wt_frac_S_core,combine_phases,use_grids]

    if use_grids == True and knownfile==False:
        filename = exo.functions.find_filename(compositional_params)
    elif knownfile==True: 
        filename=Fname
    else:
        filename=''

    structure_params =  [Pressure_range_mantle_UM,Temperature_range_mantle_UM,resolution_UM,
                         Pressure_range_mantle_LM, Temperature_range_mantle_LM, resolution_LM,
                         core_rad_frac_guess,Mantle_potential_temp,water_rad_frac_guess,water_potential_temp]

    #print(structure_params)
    layers = [num_mantle_layers,num_core_layers,number_h2o_layers]

    #This is where we actually run the planet. First PerPlex grids of mineralogy, density,
    #Cp and alpha are calculated and stored in the Solutions folder. If the file already exists
    #(in name, not necessarily in composition), then PerPlex is not run again.

    Planet = exo.run_planet_mass(Mass_planet,compositional_params,structure_params,layers,filename, Dry, magma )

    #Planet is a dictionary containing many parameters of interest:
    #Planet.get('radius') = list of the radial points from calculation (m)
    #Planet.get('mass') = list of the cumulative mass at each radius point from calculation (kg)
    #Planet.get('density') = list of densities from calculation (kg/m^3)
    #Planet.get('temperature') = list of temperature points from calculation (K)
    #Planet.get('gravity') = list of gravity points from calculation (SI)
    #Planet.get('pressure') = list of pressure points from calculation (bar)
    #Planet.get('alpha') = list of values of thermal expansivity points from calculation (1/K)
    #Planet.get('cp') = list of values of specific heat points from calculation (SI)
    #Planet.get('phases') = list of phases and their molar fractions
    
    if fullOutput:
        #If you'd like the full output, uncomment out these lines!
        Output_filename = Output_filename + '_Radius_'+ str('%.2f'%(Planet['radius'][-1]/6371e3))
        exo.functions.write(Planet,Output_filename)
        
    return Planet, num_core_layers, num_mantle_layers, wt_frac_water
    

def PlotProfile(Planet,num_core_layers, num_mantle_layers, wt_frac_water, filename, verbose=True, melt=True):
    if verbose:
        
        print()
        print("Mass = ", '%.4f'%(Planet['mass'][-1]/5.97e24), "Earth masses")
        print("Radius = ", '%.3f'%(Planet['radius'][-1]/6371e3), "Earth radii")
        print("Core Mass Fraction = ", '%.2f'%(Planet['mass'][num_core_layers]/Planet['mass'][-1]))
        print("Core Radius Fraction = ", '%.2f'%(Planet['radius'][num_core_layers]/Planet['radius'][-1]))
        print("CMB Pressure = " ,'%.2f' % (Planet['pressure'][num_core_layers]/10000), "GPa")
        print("number of oceans:",'%.2f' % (wt_frac_water*Planet['mass'][-1]/1.4e21))
       
    
    if melt:
        
        Pres=Planet['pressure'][num_core_layers: (num_mantle_layers+num_core_layers)]/10000
        Temp=Planet['temperature'][num_core_layers: (num_mantle_layers+num_core_layers)]
        meltpoint, meltpointP, ind, MT=MeltCurve(Pres, Temp)
    
    
    Teal='#11BAC8'
    dPurp='#7E04C6'
    Green='#BF076C'
    Purp='#0D6592'
    
    
    figure =plt.figure(figsize = (20,15))
    
    ax1 = plt.subplot2grid((4, 3), (0, 0), colspan=3, rowspan=1)
    ax2 = plt.subplot2grid((4, 3), (1, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((4, 3), (2, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((4, 3), (3, 0), colspan=3, rowspan=1)
    
    #Density profile
    ax1.plot(Planet['radius'] / (1.e3*6378.14), Planet['density'] / 1.e3, linewidth=2., color=Green)
    ax1.set_ylim(0., (max(Planet['density']) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Planet['radius']) /(1.e3*6378.14))
    ax1.set_ylabel("Density ( g/cm$^3$)",fontsize=20,fontfamily='serif')
    
    
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    plt.rcParams["legend.markerscale"] = 1
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    
    
    # Make a subplot showing the calculated pressure profile
    ax2.plot(Planet['radius'] / (1.e3*6378.14), Planet['pressure'] / 1.e4, linewidth=2., color=Teal)
    ax2.set_ylim(0., (max(Planet['pressure']) / 1e4) + 10.)
    ax2.set_xlim(0., max(Planet['radius']) / (1.e3*6378.14))
    ax2.set_ylabel("Pressure (GPa)",fontsize=20,fontfamily='serif')
    ax2.axvline(Planet['radius'][num_core_layers]/(1.e3*6378.14) , color = 'k', linestyle = '--')
    if melt: 
        ax2.axvline(Planet['radius'][num_core_layers+ind]/(1.e3*6378.14) , color = Purp, linestyle = '--')
        ax2.text(Planet['radius'][num_core_layers+ind]/(1.e3*6378.14) + 0.01, 0.5*max(Planet['pressure']/1.e4), 'Melt Radius: ' + str(round(Planet['radius'][num_core_layers+ind]/(1.e3*6378.14),2)) + r'$R_\oplus$',fontsize=15,fontfamily='serif', color=Purp)
        
    ax2.text(Planet['radius'][num_core_layers]/(1.e3*6378.14) +0.1, 0.5*max(Planet['pressure']/1.e4), 'CMB Pressure: ' + str(round(Planet['pressure'][num_core_layers]/1.e4,2)) + ' GPa',fontsize=15,fontfamily='serif')
    ax2.text(Planet['radius'][0]/(1.e3*6378.14) + 0.1, 0.5*max(Planet['pressure']/1.e4), 'Core Pressure: ' + str(round(Planet['pressure'][0]/1.e4,2)) + ' GPa',fontsize=15,fontfamily='serif')
    
    
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    plt.rcParams["legend.markerscale"] = 1
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    
    
    # Make a subplot showing the calculated gravity profile
    ax3.plot(Planet['radius'] / (1.e3*6378.14), Planet['gravity'], linewidth=2., color=dPurp)
    ax3.set_ylabel("Gravity (m/s$^2)$",fontsize=20,fontfamily='serif')
    ax3.set_xlim(0., max(Planet['radius']) / (1.e3*6378.14))
    ax3.set_ylim(0., max(Planet['gravity']) + 0.5)
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    plt.rcParams["legend.markerscale"] = 1
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    
    
    # Make a subplot showing the calculated temperature profile
    ax4.plot(Planet['radius'] / (1.e3*6378.14), Planet['temperature'], 'g', linewidth=2.)
    if melt: 
        ax4.plot(Planet['radius'][num_core_layers: (num_mantle_layers+num_core_layers)] / (1.e3*6378.14), MT,'--', color=Purp, linewidth=2.)
        
    ax4.set_ylabel("Temperature ($K$)",fontsize=20,fontfamily='serif')
    ax4.set_xlabel(r"Radius ($R_\oplus$)",fontsize=20,fontfamily='serif')
    ax4.set_xlim(0., max(Planet['radius']) / (1.e3*6378.14))
    ax4.axvline(Planet['radius'][num_core_layers]/(1.e3*6378.14) , color = 'k', linestyle = '--')
   # ax4.set_ylim(0., max(Planet['temperature']) + 100)
    ax4.axvline(Planet['radius'][num_core_layers]/1.e3, color = 'k', linestyle = '--')
    ax4.text(Planet['radius'][num_core_layers]/(1.e3*6378.14) + 0.1, 0.7*max(Planet['temperature']), 'CMB Temperature: ' + str(round(Planet['temperature'][num_core_layers],2)) + ' K',fontsize=15,fontfamily='serif')
    ax4.text(Planet['radius'][0]/(1.e3*6378.14) + 0.1, 0.7*max(Planet['temperature']), 'Core Temperature: ' + str(round(Planet['temperature'][0],2)) + ' K',fontsize=15,fontfamily='serif')
    #ax4.text(Planet['radius'][num_core_layers]/(1.e3*6378.14) + 0.3, 0.2*max(Planet['temperature']), 'Melting Curve', color= Purp, fontsize=15,fontfamily='serif')
    
    start, end = ax4.get_xlim()
    #ax4.xaxis.set_ticks(np.arange(start, end, 0.712123))
    #ax4.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    
    font = font_manager.FontProperties(family='serif',style='normal', size=20)
    plt.legend(prop=font )
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    plt.savefig(filename+'.png', dpi=200)
    plt.show()
    
def MeltCurve(pres, T, Dry=True, plot=True): 
    
    [idx]=np.where(pres<189.75)
    
    
    
    P1=pres[idx[0]: len(pres)]
    P=pres[0:idx[0]]
    
    
    
    if Dry: 
        a=np.array([1831,4.6,0.33]) # Units a1=K, a2=GPa, a3= none (ratio)
        MgMelt=a[0]*(1+P/a[1])**a[2]
        
        b=np.array([5400,140,0.48]) # Units b1=K, b2=GPa, b3= none (ratio)
        MgMelt1=b[0]*(P1/b[1])**b[2]
        
        MgMeltT=np.concatenate((MgMelt,MgMelt1))
        Compare= MgMeltT-T
        i= np.argwhere(np.isclose(MgMeltT,T, atol=10)).reshape(-1)
        
        
        
        
        if plot: 
            plt.plot(pres, T, label='ExoPlex')
            plt.plot(pres, MgMeltT, label='Melting Curve')
            
            plt.xlabel('Pressure (GPa)')
            plt.ylabel('Temperature (K) ')
            plt.legend()
            plt.savefig('plt.png', dpi=200)
        
        if len(i)==0: 
            meltpoint=T[0]
            meltpointP=pres[0]
            i=0
        
        else:
            meltpoint=T[i[len(i)-1]]
            meltpointP=pres[i[len(i)-1]]
            i=i[len(i)-1]
            
        return meltpoint, meltpointP , i, MgMeltT
        
    
    elif Dry=='iron': 

        xFeO=.0818+1e-15
        c=np.array([360,0.0818,102,64.1,3.62]) # Units a1=K, a2=GPa, a3= none (ratio)
        FeMelt=MgMelt+c[0]*(c[1]-xFeO)
        
        FeMelt1=MgMelt1+(c[2]+c[3]*P-c[4]*P**2)*(c[1]-xFeO)
        
        FeMeltT=np.concatenate((FeMelt,FeMelt1))
    
    elif Dry==False:

        xH2O=2
        d=np.array([43,0.75]) # Units a1=K, a2=GPa, a3= none (ratio)
        WetMelt=FeMelt-d[0]*xH2O**d[1]
        
        WetMelt1=FeMelt1-d[0]*xH2O**d[1]
        
        WetMeltT=np.concatenate((WetMelt,WetMelt1))
    P=np.concatenate((P,P1))
    
    return T
    
    
def PlotMR(Mdry,Rdry, M2,R2,M8,R8, M15,R15,M25,R25,M48,R48, filename='MR'):
    
    
    
    Teal='#11BAC8'
    dPurp='#7E04C6'
    Pink='#BF076C'
    Green='#259A00'
    Ora='#E39605'
    BG='#709B9A'
    
    
    
    fig, ax = plt.subplots(figsize = (12,9))
    
    #linestyle=(0, (3, 10, 1, 10))
    #Density profile
    ax.plot(Mdry, Rdry, linewidth=2., color=Pink, label='Anhyrdous')
    ax.plot(M2, R2, linewidth=2., color=Teal, label='Hydrous 2wt%')
    ax.plot(M8, R8, linewidth=2., color=Green, label='Hydrous 8wt%')
    
    ax.plot(M15, R15, linewidth=2., color=dPurp, label='Anhyrdous (Solomatova)')
    ax.plot(M25, R25, linewidth=2.,color=BG, label='Hydrous 2.2wt% (Solomatova)')
    ax.plot(M48, R48, linewidth=2.,color=Ora, label='Carbonated 5.2wt% (Solomatova)')
    ax.set_ylabel(r"Radius (R$_\oplus$)",fontsize=20,fontfamily='serif')
    ax.set_xlabel(r"Mass (M$_\oplus$)",fontsize=20,fontfamily='serif')
    
    yticks = ticker.MaxNLocator(10)

    # Set the yaxis major locator using your ticker object. You can also choose the minor
    # tick positions with set_minor_locator.
    #ax.yaxis.set_major_locator(yticks)
    ax.xaxis.set_major_locator(yticks)
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    plt.rcParams["legend.markerscale"] = 1
    
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    
    plt.minorticks_on()
    
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.001f'))
    
    font = font_manager.FontProperties(family='serif',style='normal', size=20)
    plt.legend(prop=font )
    plt.savefig(filename+'.png', dpi=200)
    
    plt.show()
   
def PlotSing(Mdry,Rdry,x,y,filename='Single'):
    
    
    
    Teal='#11BAC8'
    dPurp='#7E04C6'
    Green='#BF076C'
    
    
    
    fig, ax = plt.subplots(figsize = (12,9))
    
    
    #Density profile
    ax.plot(Mdry, Rdry, linewidth=3., color=Green, label='Anhyrdous')
    ax.set_ylabel(y,fontsize=20,fontfamily='serif')
    ax.set_xlabel(x,fontsize=20,fontfamily='serif')
    
    yticks = ticker.MaxNLocator(10)

    # Set the yaxis major locator using your ticker object. You can also choose the minor
    # tick positions with set_minor_locator.
    #ax.yaxis.set_major_locator(yticks)
    ax.xaxis.set_major_locator(yticks)
    
    plt.tick_params(which='both',direction='in',labelsize=16);
    plt.rcParams["legend.markerscale"] = 1
    
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    
    plt.xticks(fontsize=15,fontfamily='serif')
    plt.yticks(fontsize=15,fontfamily='serif')
    
    plt.minorticks_on()
    #ax.set_ylim(.31, .32)
    
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.001f'))
    
    font = font_manager.FontProperties(family='serif',style='normal', size=20)
    plt.legend(prop=font )
    plt.savefig(filename+'.png', dpi=200)
    
    plt.show()
    
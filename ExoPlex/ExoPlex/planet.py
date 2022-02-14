import sys
import numpy as np
Earth_radius = 6.371e6
Earth_mass = 5.97e24

from ExoPlex import minphys as minphys

def initialize_by_radius(*args):
    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the radius
   input by user.

    Parameters
    ----------
    radius_planet: float
        input radius of planet in Earth radii

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
        keys = 'radius','density','temperature','gravity','pressure', 'alpha','cp','Vphi''Vp','Vs','K'
    """
    radius_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids = compositional_params

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = structural_params[8]
    water_potential_temp = structural_params[9]



    # if there is a water layer, the imput temperature is lowered because that temperature is for the crustal layer
    # also 50 shells are used for the water layer hence the nh20 vaiable

    if wt_frac_water == 0. and number_h2o_layers > 0:
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    # used for compressiofh2on funciton

    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    # 15 mineral phases + 2ice + liquid water #phasechange
    planet_radius_guess = radius_planet*Earth_radius
    water_thickness_guess = water_rad_frac*planet_radius_guess
    core_thickness_guess = core_rad_frac * (planet_radius_guess-water_thickness_guess)
    mantle_thickness_guess = planet_radius_guess - water_thickness_guess - core_thickness_guess

    CMB_T_guess = (4180.*radius_planet-2764.*pow(radius_planet,2.)+1219.*pow(radius_planet,3.)) + (Mantle_potential_temp-1600)*(0.82+pow(radius_planet,1.81))
    #print("CMBT guess",CMB_T_guess)
    CMB_P_guess = 10000*(262.*radius_planet-550.*pow(radius_planet,2.) + 432.*pow(radius_planet,3.))

    dP_dr = (CMB_P_guess-5000)/(num_mantle_layers)
    dT_dr = (CMB_T_guess-Mantle_potential_temp)/(num_mantle_layers)

    for i in range(num_layers):

        if i <= num_core_layers:
            radius_layers[i]=((float(i)/num_core_layers)*core_thickness_guess)

        elif i < (num_core_layers+num_mantle_layers):
            radius_layers[i]=(core_thickness_guess+((float(i-num_core_layers)/num_mantle_layers)*mantle_thickness_guess))
            #density_layers[i]=3100.

        else:

            radius_layers[i]=core_thickness_guess+mantle_thickness_guess+\
                             ((float(i-num_core_layers-num_mantle_layers)/number_h2o_layers)*water_thickness_guess)
            #density_layers[i]=1100.

    for i in range(num_layers):

            if i<number_h2o_layers:
                Pressure_layers[i] = 5000.
                Temperature_layers[i] = 300.

            elif i <= number_h2o_layers+num_mantle_layers-1:
                Pressure_layers[i] = 5000. + dP_dr*(i-number_h2o_layers)
                Temperature_layers[i] = Mantle_potential_temp+dT_dr*(i-number_h2o_layers)
                if Temperature_layers[i]>7000:
                    Temperature_layers[i] = 7000
            else:
                Pressure_layers[i] = Pressure_layers[number_h2o_layers+num_mantle_layers]  + 3.5*(dP_dr*num_mantle_layers/num_core_layers)*(i-number_h2o_layers-num_mantle_layers)
                Temperature_layers[i] = 1900.

    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers= Temperature_layers[::-1]
    #Pressure_layers[-1] = 3000.

    #print(Temperature_layers[num_core_layers+1])
    radius_layers[-1] = planet_radius_guess
    #initial temperature guess of 0.5 K per km
    keys = ['radius','density','temperature','gravity','pressure',\
            'alpha','cp','Vphi''Vp','Vs','K']


    return dict(zip(keys,[radius_layers, density_layers,Temperature_layers,gravity_layers, Pressure_layers,
                          alpha, cp,Vphi,Vp,Vs,K]))

def initialize_by_mass(*args):
    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the mass
   input by user.

    Parameters
    ----------
    mass_planet: float
        input radius of planet in Earth mass

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    densityEarthUnits=0.23896328515
    mass_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]
    core_mass_frac = args[4]

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids = compositional_params

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = structural_params[8]
    water_potential_temp = structural_params[9]

    if wt_frac_water == 0. and number_h2o_layers > 0:
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    radius_layers = np.zeros(num_layers)

    # used for compressiofh2on funciton
    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    # 15 mineral phases + 2ice + liquid water #phasechange
    water_mass = (wt_frac_water*mass_planet)*Earth_mass
    core_mass = (core_mass_frac * (mass_planet*Earth_mass-water_mass))
    mantle_mass = (mass_planet*Earth_mass)-water_mass-core_mass


    Radius_planet_guess = 1.

    mass_layers[0] = 0
    #Update to use BurnMan
    if num_mantle_layers > 0:
        if num_core_layers >0:
            for i in range(num_layers):

                if i <num_core_layers:

                    radius_layers[i] = (float(i)/float(num_layers))*(Radius_planet_guess*Earth_radius)
                    mass_layers[i]  = core_mass/num_core_layers

                    Temperature_layers[i] = 10000.

                elif i < (num_core_layers+num_mantle_layers):

                    radius_layers[i] = (float(i)/float(num_layers))*(Radius_planet_guess*Earth_radius)
                    mass_layers[i] = (mantle_mass/num_mantle_layers)
                    R=((3/(4*np.pi))*((mass_planet/densityEarthUnits)))**(1/3)
                    
                    T_CMB=4180*R-2764*R**2+1219*R**3
    
                    Temperature_layers[i] =Temperature_layers[i-1]-(T_CMB-300)/num_mantle_layers # changed 1600.
                    
                    if i== num_layers-1:
                        print(Temperature_layers[i], R)
                        import matplotlib.pyplot as plt
                        plt.plot(radius_layers,Temperature_layers)
                        plt.savefig("Temp.png")
                   
                else:
                    radius_layers[i] = (float(i)/float(num_layers))*(Radius_planet_guess*Earth_radius)
                    mass_layers[i] = (water_mass/number_h2o_layers)


                    Temperature_layers[i] = 300.
                    

            for i in range(num_layers):
                #in the water
                if i >= num_core_layers+num_mantle_layers:
                    
                    Pressure_layers[i] = 1
                else:
                    #in the core
                    P_CMB=262*R-550*R**2+432*R**3
                    Pressure_layers[i] =(float((5000.-(P_CMB*10000))/float(num_core_layers+num_mantle_layers))*float(i)
                                          + P_CMB*10000)#Pressure_layers[i-1]-(P_CMB-3000)/num_mantle_layers
                    ''' (float((5000.-(300.*10000))/float(num_core_layers+num_mantle_layers))*float(i)
                                          + 300.*10000)
                    Pressure_layers[i-1]+(P_CMB-1)/num_mantle_layers
                    
                   
                    '''
                    '''if i== num_layers-1:
                        print(Pressure_layers[i])
                        plt.figure()
                        plt.plot(radius_layers,Pressure_layers)
                        plt.savefig("Pres2.png")'''
        else:
            for i in range(num_mantle_layers):
                mantle_mass = mass_planet * Earth_mass

                radius_layers[i] = (float(i) / float(num_layers)) * (Radius_planet_guess * Earth_radius)
                mass_layers[i] = (mantle_mass / num_mantle_layers)
                Pressure_layers[i] = (1e8 / 10000)
                
                T_CMB=4180*radius_layers[i]-2764*radius_layers[i]**2+1219*radius_layers[i]**3
    
                Temperature_layers[i] = Temperature_layers[i-1]+(T_CMB-2000)/num_mantle_layers # changed 1600.
                print(Temperature_layers[i])
                


    elif number_h2o_layers > 0:
        water_mass = mass_planet*Earth_mass
        for i in range((number_h2o_layers)):
            radius_layers[i] = (float(i) / float(num_layers)) * (Radius_planet_guess * Earth_radius)
            mass_layers[i] = (water_mass / number_h2o_layers)
            Temperature_layers[i] = 300.
            Pressure_layers[i] = 1e3
    else:
        core_mass = mass_planet*Earth_mass
        for i in range((num_core_layers)):
            radius_layers[i] = (float(i) / float(num_layers)) * (Radius_planet_guess * Earth_radius)
            mass_layers[i] = core_mass / num_core_layers
            Temperature_layers[i] = 0.

    mass_update = np.zeros(num_layers)

    for i in range(len(mass_layers)):
        mass_update[i]=(sum(mass_layers[:i+1]))

    mass_layers= mass_update

    keys = ['mass', 'density', 'temperature', 'gravity', 'pressure', \
            'alpha', 'cp', 'Vphi''Vp', 'Vs', 'K']

    return dict(zip(keys, [mass_layers, density_layers, Temperature_layers, gravity_layers, Pressure_layers,
                           alpha, cp, Vphi, Vp, Vs, K]))



def initialize_by_mass(*args):
    """
   This module creates the dictionary of lists for each planetary parameter (e.g., density) for a planet of the mass
   input by user.

    Parameters
    ----------
    mass_planet: float
        input radius of planet in Earth radii

    structural_params: list
        Structural parameters of the planet; See example for description

    compositional_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    mass_planet = args[0]
    structural_params = args[1]
    compositional_params = args[2]
    num_mantle_layers, num_core_layers, number_h2o_layers = args[3]
    core_mass_frac = args[4]
    
 

    wt_frac_water, FeMg, SiMg, CaMg, AlMg, mol_frac_Fe_mantle, wt_frac_Si_core, \
     wt_frac_O_core, wt_frac_S_core,combine_phases,use_grids = compositional_params

    core_rad_frac = structural_params[6]
    Mantle_potential_temp = structural_params[7]
    water_rad_frac = structural_params[8]
    water_potential_temp = structural_params[9]

    Radius_planet_guess = pow(mass_planet*5.97e24 / 5500. / (4*np.pi/3.),1/3.)/6371e3


    water_thickness_guess = water_rad_frac*Radius_planet_guess
    core_thickness_guess = core_rad_frac * (Radius_planet_guess-water_thickness_guess)
    mantle_thickness_guess = Radius_planet_guess - water_thickness_guess - core_thickness_guess

    if wt_frac_water == 0. and number_h2o_layers > 0:
       print ("You have layers of water but no water!")
       number_h2o_layers = 0

    num_layers = num_core_layers+num_mantle_layers + number_h2o_layers # add 50 shells if there is an h2O layer
    # arrays to be used

    # these grids are initialized in this function and are then a passed around through the other routines and edited each iteration
    radius_layers = np.zeros(num_layers)
    density_layers = np.zeros(num_layers)
    volume_layers = np.zeros(num_layers)
    mass_layers = np.zeros(num_layers)
    cumulative_mass = np.zeros(num_layers)
    Temperature_layers = np.zeros(num_layers)
    radius_layers = np.zeros(num_layers)

    # used for compressiofh2on funciton
    gravity_layers = np.zeros(num_layers)
    Pressure_layers = np.zeros(num_layers)
    alpha = np.zeros(num_layers)
    cp = np.zeros(num_layers)
    Vphi = np.zeros(num_layers )

    Vs = np.zeros(num_layers )
    Vp = np.zeros(num_layers)
    K =  np.zeros(num_layers)

    # 15 mineral phases + 2ice + liquid water #phasechange
    water_mass = (wt_frac_water*mass_planet)*Earth_mass
    core_mass = (core_mass_frac * (mass_planet*Earth_mass-water_mass))
    mantle_mass = (mass_planet*Earth_mass)-water_mass-core_mass


    #update radius guess to scale with mass


    mass_layers[0] = 0
    #Update to use BurnMan
    CMB_T_guess = (4180. * Radius_planet_guess - 2764. * pow(Radius_planet_guess, 2.) + 1219. * pow(Radius_planet_guess,
                                                                                                    3.)) + (
                              Mantle_potential_temp - 1600) * (0.82 + pow(Radius_planet_guess, 1.81))
    # print("CMBT guess",CMB_T_guess)
    CMB_P_guess = 10000 * (
                262. * Radius_planet_guess - 550. * pow(Radius_planet_guess, 2.) + 432. * pow(Radius_planet_guess, 3.))

    dP_dr = (CMB_P_guess - 5000.) / (num_mantle_layers)
    dT_dr = (CMB_T_guess - Mantle_potential_temp) / (num_mantle_layers)

    for i in range(num_layers):

        if i <= num_core_layers:
            radius_layers[i] = ((float(i) / num_core_layers) * core_thickness_guess)
            mass_layers[i] = core_mass / num_core_layers

        elif i < (num_core_layers + num_mantle_layers):
            radius_layers[i] = (core_thickness_guess + ((float(i - num_core_layers) / num_mantle_layers) * mantle_thickness_guess))
            mass_layers[i] = (mantle_mass / num_mantle_layers)

        else:

            radius_layers[i] = core_thickness_guess + mantle_thickness_guess + \
                               ((float(i - num_core_layers - num_mantle_layers) / number_h2o_layers) * water_thickness_guess)
            mass_layers[i] = (mantle_mass / num_mantle_layers)

    for i in range(num_layers):

        if  i < number_h2o_layers:
            Pressure_layers[i] = 5000.
            Temperature_layers[i] = 300.

        elif i <= number_h2o_layers + num_mantle_layers - 1:
            Pressure_layers[i] = 5000. + dP_dr * (i - number_h2o_layers)
            Temperature_layers[i] = Mantle_potential_temp + dT_dr * (i - number_h2o_layers)
            if Temperature_layers[i] > 7000:
                Temperature_layers[i] = 7000
        else:
            Pressure_layers[i] = Pressure_layers[number_h2o_layers + num_mantle_layers] + 3.5 * (
                    dP_dr * num_mantle_layers / num_core_layers) * (i - number_h2o_layers - num_mantle_layers)
            Temperature_layers[i] =  4500+ (float((Temperature_layers[num_core_layers]+1000-(Temperature_layers[num_core_layers])))/float(num_core_layers))\
                                    *float((i-num_layers))
       
        
                                    
    Pressure_layers = Pressure_layers[::-1]
    Temperature_layers = Temperature_layers[::-1]

    mass_update = np.zeros(num_layers)

    for i in range(len(mass_layers)):
        mass_update[i] = (sum(mass_layers[:i + 1]))


    mass_layers = mass_update
    mass_layers[-1] = mass_planet * Earth_mass

    keys = ['mass', 'density', 'temperature', 'gravity', 'pressure', \
            'alpha', 'cp', 'Vphi''Vp', 'Vs', 'K']

    return dict(zip(keys, [mass_layers, density_layers, Temperature_layers, gravity_layers, Pressure_layers,
                           alpha, cp, Vphi, Vp, Vs, K]))

def compress_radius(*args):
    """
   This module iterates the density, mass within a sphere, adiabatic temperature and gravity integrals for a planet of radius R
   until convergence is reached. Convergence is defined as the change from the previous run to the current is the
   difference in the density of all layers is <1e-6.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet

     grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases

    Core_wt_per: float
        Composition of the Core

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """

    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]
    structural_params= args[3]
    layers= args[4]
    Dry=args[5]
    Magma=args[6]
    
    n_iterations = 1
    max_iterations = 100
    
    if Magma: 
        number=10
    else:
        number=1


    old_rho = [10  for i in range(len(Planet['density']))]
    converge = False
    while n_iterations <= max_iterations and converge == False:
        #print "iteration #",n_iterations


        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print ("Density has a nan")
                print (i, Planet['pressure'][i],Planet['temperature'][i])
                print

                sys.exit()

        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers , Magma,Dry)
        Planet['gravity'] = minphys.get_gravity(Planet,layers)
        Planet['pressure'] = minphys.get_pressure(Planet,layers)

        if n_iterations >number:
            Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers, Magma)
            converge, old_rho = minphys.check_convergence(Planet['density'], old_rho)

        n_iterations+=1

    return Planet

def compress_mass(*args):
    """
   This module iterates the density, mass within a sphere, adiabatic temperature and gravity integrals for a planet of Mass M
   until convergence is reached. Convergence is defined as the change from the previous run to the current is the
   difference in the density of all layers is <1e-6.

    Parameters
    ----------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet

     grids: list of lists
        UM and LM grids containing pressure, temperature, density, expansivity, specific heat and phases

    Core_wt_per: float
        Composition of the Core

    structural_params: list
        Structural parameters of the planet; See example for description

    layers: list
        Number of layers for core, mantle and water
    Returns
    -------
    Planet: dictionary
        Dictionary of initial guess of pressure, temperature, expansivity, specific heat and phases for modeled planet
    """
    Planet = args[0]
    grids = args[1]
    Core_wt_per = args[2]
    structural_params= args[3]
    layers= args[4]
    dry=args[5]
    magma=args[6]
    n_iterations = 1
    max_iterations = 20
    
    if magma:
        number=10
    else:
        number=1


    old_r = [10  for i in range(len(Planet['mass']))]
    converge = False

    while n_iterations <= max_iterations and converge == False:
        print ("iteration #",n_iterations)
        if n_iterations>number:
            converge,old_r = minphys.check_convergence(Planet['density'],old_r)

        for i in range(len(Planet['density'])):
            if np.isnan(Planet['density'][i]) == True:
                print ("Density has a nan")
                print (i, Planet['pressure'][i],Planet['temperature'][i])
                print
                sys.exit()


        Planet['density'] = minphys.get_rho(Planet,grids,Core_wt_per,layers,dry,magma )
        Planet['radius'] = minphys.get_radius(Planet, layers)
        Planet['gravity'] = minphys.get_gravity(Planet,layers)
        Planet['temperature'] = minphys.get_temperature(Planet, grids, structural_params, layers,magma)
        Planet['pressure'] = minphys.get_pressure(Planet,layers)
        n_iterations+=1

    return Planet




# coding: utf-8

# # Installing Cantera
# For this notebook you will need [Cantera](http://www.cantera.org/), an open source suite of object-oriented software tools for problems involving chemical kinetics, thermodynamics, and/or transport processes.
# Fortunately there are now Anaconda packages (in the `cantera` channel), so to install you can simply type 
# ```
# conda install -c cantera cantera
# ```
# at your terminal (if you can remember back to when you installed Anaconda!).
# 
# If you are on Windows *and* you get an error like `ImportError: DLL load failed: The specified module could not be found.` then you probably also need to install the Visual C++ Redistributable which you can get [from Microsoft here](https://www.microsoft.com/en-us/download/details.aspx?id=48145).  If you don't have this problem, let me know and I can updated these instructions!
# 
# There are other, more difficult, ways to install it in [the instructions](http://www.cantera.org/docs/sphinx/html/install.html) if you can't get the Anaconda packages to work. It is also already on the COE computer lab 274 Snell (though there you will have to `pip install jupyter` to get this notebook working).

# In[1]:

# First, import cantera, with the nickname `ct` to save us some typing later.
import cantera as ct

# Then the usual suspects:
import numpy as np
get_ipython().magic('matplotlib inline')
from matplotlib import pyplot as plt


# # Heptane combustion

# Download the reduced n-heptane model from LLNL https://combustion.llnl.gov/archived-mechanisms/alkanes/heptane-reduced-mechanism. Save the files alongside this python notebook. These files are in "CHEMKIN" format. First, we have to convert them into a format that is usable by Cantera. 
# This may take a while and issue some warnings, but then end by saying `Validating mechanism...PASSED`:

# In[ ]:

from cantera import ck2cti
ck2cti.main(['--input=heptanesymp159_mec.txt',
             '--thermo=heptanesymp_therm.txt',
             '--permissive',
             '--output=heptanesymp159.cti'])


# Clearly, with 160 species and 1540 reactions, this mechanism is more detailed than any we have considered before!
# Now, let's create a 'Solution' phase in Cantera called `gas` from the Cantera mechanism file we just created.

# In[ ]:

gas = ct.Solution('heptanesymp159.cti')


# Let's examine some of the reactions and species in the mechanism.
# 
# This will return the first 10 reactions:

# In[ ]:

gas.reaction_equations(np.arange(10))


# And this will return a list of the chemical species names, joined by spaces:

# In[ ]:

print(" ".join(gas.species_names))


# Knowing what all those species names mean is a [formidable challenge](http://www.northeastern.edu/comocheng/2014/04/nsf-grant-to-identify-discrepancies/) but we are [making headway](http://www.northeastern.edu/comocheng/2015/05/uscombustionmeeting/) (and more help is welcome).
# 
# For now, lets loop through all the species looking for ones with 7 carbons and 16 hydrogen atoms, which should be all the isomers of heptane.

# In[ ]:

for species in gas.species():
    if species.composition == {'C':7, 'H':16}:
        print(species.name)


# There is only one!
# Based on the name beginning with 'n' let's assume it represents normal-heptane (all 7 carbons in a single chain with no branching), which is the fuel that we want to simulate. Now we need to find the index number for this species.

# In[ ]:

i_heptane = gas.species_names.index('nc7h16')
# alternatively, this shortcut:
i_heptane = gas.species_index('nc7h16')
print("heptane is species index {0}".format(i_heptane))


# To specify the state of a system we must supply two intensive variables (temperature, pressure, density, specific entropy, specific enthalpy, specific volume) and the composition (mass or mole fractions). We will set the temperature, pressure, and mole fractions. In cantera, mole fractions are `X` and mass fractions are `Y`. We can then print some properties of our gas system by typing  `gas()`.

# In[ ]:

gas.TPX = 1000, 10e5, 'nc7h16:1.0'
gas()


# To find equilbrium you must specify which two intensive variables to hold constant. We'll find the equilbrium at constant Temperature and Pressure, then print the properties again.

# In[ ]:

gas.equilibrate('TP')
gas()


# You will recall from Thermodynamics II that a system going to equilibrium at constant T and P should minimize the specific Gibbs free energy of the system. Sure enough, it has gone down (compare the "Gibbs function" in the "1 kg" columns above.  To check that number represents what we expect (this will be returned in Cantera's default SI units, a combination of K, m<sup>3</sup>, Pa, J, kg, kmol; in this case J/kg)

# In[ ]:

print(gas.h - gas.T * gas.s)
print(gas.g)


# Now lets find the equilibrium composition at 1 bar pressure and a range of temperatures between 100 and 2000 K

# In[ ]:

temperatures = np.arange(100,2000,20)
# make a big array to store the results in
equilibrium_mass_fractions = np.zeros((len(temperatures), gas.n_species))
for i, T in enumerate(temperatures):
    gas.TP = T, 1e5
    gas.equilibrate('TP')
    print(T,end=" ")
    equilibrium_mass_fractions[i,:] = gas.Y


# Now plot the equilibrium mass fractions as a function of temperature. With 160 lines, let's forgo the legend and instead label the biggest peaks directly.

# In[ ]:

plt.plot(temperatures,equilibrium_mass_fractions)
plt.xlabel("Temperature (K)")
plt.ylabel("Equilibrium mole fraction")
for i, name in enumerate(gas.species_names):
    Y = equilibrium_mass_fractions[:,i]
    if max(Y)> 0.08:
        peakT = temperatures[Y.argmax()]
        peakY = max(Y)
        plt.text(peakT,peakY, name)
plt.show()


# ## Question (a) 
# What do you notice about the species that peaks at 100K, and the ones that peak at 2000K? Can you explain or justify this?

# To see some of the complexity hidden at low concentrations, let's plot the y axis on a logarithmic scale:

# In[ ]:

plt.semilogy(temperatures,equilibrium_mass_fractions)
plt.ylim(1e-30,1)
plt.xlabel("Temperature (K)")
plt.ylabel("Equilibrium mole fraction")
plt.show()


# If you think about how many reactions are equilibrated, it was remarkably quick!
# 
# Now we'll add some air, which is mostly nitrogen and oxygen. First of all, find the names of anything with just 2 oxygen atoms or just 2 nitrogen atoms.

# In[ ]:

for species in gas.species():
    if species.composition == {'O':2} or species.composition == {'N':2}:
        print(species.name)


# Now look up and store the species indices

# In[ ]:

i_oxygen = gas.species_names.index('o2')
print("oxygen is species index {0}".format(i_oxygen))
i_nitrogen = gas.species_names.index('n2')
print("nitrogen is species index {0}".format(i_nitrogen))


# ## Question (b)
# 
# For a "stoichiometric" mixture of n-heptane and air (enough oxygen to reach complete combustion) how many moles of heptane and how many moles of nitrogen should you have for one mole of oxygen?  Assume air is 80% nitrogen and 20% oxygen.

# In[ ]:

oxygen_mole = 1. # moles oxygen
## ANSWER:


# In[ ]:

## Some checks
assert nitrogen_mole / oxygen_mole == 4, "Assume air is 80% nitrogen and 20% oxygen"
assert oxygen_mole / heptane_mole == 3+1+3//5*3+8-5//3, "C7H16 + ?? O2 => 8 H2O + 7 CO2"


# Now use those to make a string for the '`X`' when we set `gas.TPX`. Although we call it a mole fraction, they don't need to add up to one: Cantera will normalize it, preserving the ratios.  Then print it, use it, and check it.

# In[ ]:

X_string = 'nc7h16:{0},o2:{1},n2:{2}'.format(heptane_mole, oxygen_mole, nitrogen_mole)
print("The 'X' will be set to {0!r}".format(X_string))
gas.TPX = 1000, 10e5, X_string
gas()
assert round(gas.concentrations[i_oxygen] / gas.concentrations[i_heptane], 2) == 11


# ## Question (c)
# We can do an equilibrium analysis like before, but before you do,
# starting with a stoichiometric mixture of fuel and air
# what do you expect the equilibrium composition to mostly consist of?
# (Imagine all reactions are fast with no barriers)

# In[ ]:

temperatures = np.arange(100,2000,20)
# make a big array to store the results in
equilibrium_mass_fractions = np.zeros((len(temperatures), gas.n_species))
for i, T in enumerate(temperatures):
    gas.TP = T, 1e5
    gas.equilibrate('TP')
    print(T, end=" ")
    equilibrium_mass_fractions[i,:] = gas.Y
plt.plot(temperatures,equilibrium_mass_fractions)
for i, name in enumerate(gas.species_names):
    Y = equilibrium_mass_fractions[:,i]
    if max(Y)> 0.08:
        peakT = temperatures[Y.argmax()]
        peakY = max(Y)
        plt.text(peakT,peakY, name)
plt.show()


# ## Kinetics
# Now we are done with equilbria, let's calculate some kinetics!
# 
# Cantera can do complex networks of reactors with valves, flow controllers, etc.
# but we will make a simple "reactor network" with just one constant volume ideal gas batch reactor.

# In[ ]:

gas.TPX = 800, 10e5, X_string

reactor = ct.IdealGasReactor(gas)
reactor_network = ct.ReactorNet([reactor])

start_time = 0.0  #starting time
end_time = 4e-3 # seconds
n_steps = 251
times = np.linspace(start_time, end_time, n_steps)
concentrations = np.zeros((n_steps, gas.n_species))
pressures = np.zeros(n_steps)
temperatures = np.zeros(n_steps)

print_data = True
if print_data:
    #this just gives headings
    print('{0:>10s} {1:>10s} {2:>10s} {3:>14s}'.format(
            't [s]','T [K]','P [Pa]','u [J/kg]')) 

for n, time in enumerate(times):
    if time > 0:
        reactor_network.advance(time)
    temperatures[n] = reactor.T
    pressures[n] = reactor.thermo.P
    concentrations[n,:] = reactor.thermo.concentrations
    if print_data:
        print('{0:10.3e} {1:10.3f} {2:10.3f} {3:14.6e}'.format(
                 reactor_network.time, reactor.T, reactor.thermo.P, reactor.thermo.u))


# Now let's plot some graphs to see how things look

# In[ ]:

plt.plot(times*1e3, concentrations[:,i_heptane])
plt.ylabel("Heptane concentration (kmol/m3)")
plt.xlabel("Time (ms)")
plt.ylim(0,)
plt.show()
plt.plot(times*1e3, pressures/1e5)
plt.xlabel("Time (ms)")
plt.ylabel("Pressure (bar)")
plt.show()
plt.plot(times*1e3, temperatures)
plt.xlabel("Time (ms)")
plt.ylabel("Temperature (K)")
plt.show()


# Although the timescale is milliseconds instead of hours, that looks remarkably like the thermal runaway reaction that caused the T2 laboratory explosion that we studied last lecture. This time, however, it's not just a thermal runaway but a chemical runaway - it's the gradual accumulation of reactive radical species like `OH` that is auto-catalytic.
# 
# Let's look at some of the other species:

# In[ ]:

# skip the zeroth species which is nitrogen
plt.plot(times*1e3, concentrations[:,1:])
plt.ylim(0,)
plt.ylabel("Concentration")
plt.xlabel("Time (ms)")
for i, name in enumerate(gas.species_names):
    if i==0: continue
    concentration = concentrations[:,i]
    peak_concentration = max(concentration)
    if peak_concentration > 0.001:
        peak_time = times[concentration.argmax()]
        plt.text(peak_time*1e3, peak_concentration, name)
plt.show()


# Let's zoom in on the y axis by making it logarithmic:

# In[ ]:

plt.semilogy(times*1e3, concentrations)
plt.ylim(1e-15,1)
plt.ylabel("Concentration")
plt.xlabel("Time (ms)")
plt.show()


# What a mess! Let's zoom in a little and see if we can pick out any significant intermediates

# In[ ]:

plt.semilogy(times*1e3, concentrations)
plt.ylim(1e-4,1)

# Add some labels
for t in [1.5, 3]:
    i = (times*1e3>t).nonzero()[0][0]
    time = times[i]*1e3
    for j, name in enumerate(gas.species_names):
        concentration = concentrations[i,j]
        if concentration > 1e-4:
            plt.text(time, concentration, name)
plt.ylabel("Concentration")
plt.xlabel("Time (ms)")
plt.show()


# Not really!  We would have to do a flux analysis and [reaction path diagram](http://www.cantera.org/docs/sphinx/html/cython/examples/kinetics_reaction_path.html) to see what is going on. 
# 
# ## Defining ignition delay time.
# We want to identify when the ignition occurs, so that we could compare our simulation with an experiment.
# Some experiments measure pressure rise; some monitor the concentration of an intermediate like `OH` via laser absorption; but other studies monitor the luminescence of excited `OH*` decaying to ground state `OH` (which it does by emitting a photon). This process is proportional to the rate of formation (not concentration) of `OH*`, which is predominantly made by reaction of `CH` with `O2`, so it is pretty closely proportional to the product `[CH][O2]`, i.e. "brightest flash of light" is propontional to “peak `OH*` emission” which can be modeled as “peak in the product of `[CH]` and `[O2]`”.  Likewise photoemission from creation of excited `CH*` can be modeled reasonably as the product `[C2H][O]`. When modeling an experiment it's important to know precisely what the experimenter measurend and how they defined their derived parameters. For now we'll look for the peak in `OH*` emission:

# In[ ]:

i_ch = gas.species_index('ch')
i_o2 = gas.species_index('o2')
excited_oh_generation = concentrations[:,i_ch] * concentrations[:,i_o2]
plt.plot(times*1e3, excited_oh_generation)
plt.xlabel("Time (ms)")
plt.ylabel("Excited OH* emission (arbitrary units)")
plt.show()
ignition_time = times[excited_oh_generation.argmax()]
print("Ignition delay time is {0} ms".format(ignition_time * 1e3))


# Now let's put it all together, into a function that takes temperature, pressure, and stoichiometry, and predicts ignition delay time for n-heptane. It's a bit different from before - now we let the ODE solver choose the array of times, which means we don't know how long it will be when we begin, so we have to use lists (which can grow as we add to them) and convert to arrays when we've finished. 

# In[ ]:

def get_ignition_delay(temperature, pressure = 10.,
                       stoichiometry = 1.0, plot = False):
    """
    Get the ignition delay time in miliseconds, at the specified
    temperature (K), pressure (bar), and stoichiometry 
    (stoichiometric = 1.0, fuel-rich > 1.0, oxygen-rich < 1.0).
    Default pressure is 10.0 bar, default stoichoimetry is 1.0.
    If plot=True then it draws a plot (default is False).
    """
    oxygen_mole = 1. 
    nitrogen_mole = 4*oxygen_mole
    heptane_mole = stoichiometry/11
    X_string = 'nc7h16:{0},o2:{1},n2:{2}'.format(heptane_mole, oxygen_mole, nitrogen_mole)

    gas.TPX = temperature, pressure*1e5, X_string
    reactor = ct.IdealGasReactor(gas)
    reactor_network = ct.ReactorNet([reactor])

    time = 0.0
    end_time = 10e-3 
    
    # Use lists instead of arrays, so they can be any length
    times = []
    concentrations = []
    pressures = []
    temperatures = []
    
    print_data = True
    while time < end_time:
        time = reactor_network.time
        times.append(time)
        temperatures.append(reactor.T)
        pressures.append(reactor.thermo.P)
        concentrations.append(reactor.thermo.concentrations)
        # take a timestep towards the end_time.
        # the size of the step will be determined by the ODE solver
        # depending on how quickly things are changing.
        reactor_network.step(end_time)
    
    print("Reached end time {0:.2f} ms in {1} steps".format(times[-1]*1e3, len(times)))
    # convert the lists into arrays
    concentrations = np.array(concentrations)
    times = np.array(times)
    pressures = np.array(pressures)
    temperatures = np.array(temperatures)

    if plot:
        plt.subplot(2,1,1)
        plt.plot(times*1e3, pressures/1e5)
        plt.ylabel("Pressure (bar)", color='b')
        ax2 = plt.gca().twinx()
        ax2.set_ylabel('Temperature (K)', color='r')
        ax2.plot(times*1e3, temperatures, 'r')
    i_ch = gas.species_index('ch')
    i_o2 = gas.species_index('o2')
    excited_oh_generation = concentrations[:,i_o2] * concentrations[:,i_ch]
    if plot:
        plt.subplot(2,1,2)
        plt.plot(times*1e3, excited_oh_generation, 'g')
        plt.ylabel("OH* emission")
        plt.ylim(0,max(1e-8,1.1*max(excited_oh_generation)))
        plt.xlabel("Time (ms)")
        plt.tight_layout()
        plt.show()
    step_with_highest_oh_gen = excited_oh_generation.argmax()
    if step_with_highest_oh_gen > 1 and excited_oh_generation.max()>1e-20:
        ignition_time_ms = 1e3 * times[step_with_highest_oh_gen]
        print("At {0} K {1} bar, ignition delay time is {2} ms".format(temperature, pressure, ignition_time_ms))
        return ignition_time_ms
    else:
        print("At {0} K {1} bar, no ignition detected".format(temperature, pressure))
        return np.infty


# Let's test it at 1000 K, 10 bar.

# In[ ]:

get_ignition_delay(1000, 10, plot=True)


# Now let's repeat it at a range of temperatures and pressures, and plot all the delay times on one graph

# In[ ]:

temperatures = np.linspace(1000,1500.,25)
ignition_delay_times = np.zeros_like(temperatures)
for P in [10,50]:
    for i,T in enumerate(temperatures):
        ignition_delay_times[i] = get_ignition_delay(T, P)
    plt.semilogy(1000./temperatures, ignition_delay_times, 'o-', label='{0} bar'.format(P))

plt.legend(loc='best')
plt.xlabel("1000K / temperature")
plt.ylabel("Ignition delay time (ms)")
plt.ylim(1e-2,)
plt.show()


# ## Question (d)
# Explain why this look as you would expect from Arrhenius behaviour.

# ## Question (e)
# Repeat the analysis but going down to 650K (i.e. cover the range 650-1500K).
# Describe and try to explain what you find.

# In[ ]:

temperatures = np.linspace(650,1500.,25)
ignition_delay_times = np.zeros_like(temperatures)
for P in [10,50]:
    for i,T in enumerate(temperatures):
        ignition_delay_times[i] = get_ignition_delay(T, P)
    plt.semilogy(1000./temperatures, ignition_delay_times, 'o-', label='{0} bar'.format(P))

plt.legend(loc='best')
plt.xlabel("1000K / temperature")
plt.ylabel("Ignition delay time (ms)")
plt.ylim(1e-2,)
plt.show()


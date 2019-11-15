# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import os	# To change directory
#os.chdir('/home/aleksandar/Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

#constants



# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_fwf(file, header=None)
	return df

#
def expectation_E(temp):
    """
    calculating the expectation value for the energy with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
    
    returns
    -------
    the expected energy value
    """
    z = 4*(3 + np.cosh(8/temp))
    return -32*np.sinh(8/temp)/z

def expectation_E_squered(temp):
    """
    calculating the expectation value for the energy squered with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
    
    returns
    -------
    the expected energy^2 value
    """
    z = 4*(3 + np.cosh(8/temp))
    return 4*64*np.cosh(8/temp) / z

def expected_magnetization(temp):
    """
    calculating the expected value for the magnetization with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp

    returns
    -------
    the expected magnetization value
    """
    z = 4*(3 + np.cosh(8/temp))
    return 8*(2 + np.exp(8/temp)) / z

def expected_magnetization_squered(temp):
    """
    calculating the expected value for the magnetization squered with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
    
    returns
    -------
    the expected magnetization^2 value
    """
    z = 4*(3 + np.cosh(8/temp))
    return 32*np.exp(8/temp) / z #weird comment: I actually find 32, oh well, this is per particle at least

def criticalheat(temp):
    """
    calculating the expected value for the critical heat with coupling constant, J = 1, k = 1 and beta = 1/temperature.
    since for our purpose both J, k and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
    
    returns
    -------
    the critical heat value
    """
    return 1 / (temp*temp)*(expectation_E_squered(temp) - expectation_E(temp)**2)

def suceptibility(temp):
    """
    calculating the expected value for the suceptibility with coupling constant, J = 1, k = 1 and beta = 1/temperature.
    since for our purpose both J, k and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
    
    returns
    -------
    the value for suceptibility
    """
    return 1 / temp*(expected_magnetization_squered(temp) - expected_magnetization(temp)**2)


temperature = np.arange(0.5, 4.0, 0.1)
plt.gca().set_title("gjkdg")
plt.subplot(2,2,1)
plt.plot(temperature, expectation_E(temperature))
plt.ylabel(r" Mean Energy, $\langle E \rangle$", fontsize=12)
#plt.title(r"Mean values of E and M, Cv and $\chi$ as functions of T")

plt.subplot(2,2,2)
plt.plot(temperature, expected_magnetization(temperature))
plt.ylabel(r"Mean Magnetization, $\langle M \rangle$", fontsize=12)

plt.subplot(2,2,3)
plt.plot(temperature, criticalheat(temperature))
plt.ylabel(r"critical heat, Cv", fontsize=12)

plt.subplot(2,2,4)
plt.plot(temperature, suceptibility(temperature))
plt.ylabel(r"Suceptibility, $\chi$", fontsize=12)

plt.subplot(224)
plt.xlabel(r"Temperature, $[k_b / J]$", fontsize=12)
plt.subplot(223)
plt.xlabel(r"Temperature, $[k_b / J]$", fontsize=12)
plt.show()




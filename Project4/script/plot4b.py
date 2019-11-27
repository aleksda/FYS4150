# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import os	# To change directory
#os.chdir('../Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")


# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, names = ['MCs','E','E2','M','M2','absM'], 
		engine='python', sep=" ")
	return df


def expectation_E(temp):
    """
    calculating the expectation value for the energy with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the system
    returns
    -------
    the expected energy value
    """
    z = 4 * (3 + np.cosh(8 / temp))
    return -32 * np.sinh(8 / temp) / z

def expectation_E_squered(temp):
    """
    calculating the expectation value for the energy squered with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the system
    returns
    -------
    the expected energy^2 value
    """
    z = 4 * (3 + np.cosh(8 / temp))
    return 4 * 64 * np.cosh(8 / temp) / z

def expected_magnetization(temp):
    """
    calculating the expected value for the magnetization with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the sustem
    returns
    -------
    the expected magnetization value
    """
    z = 4 * (3 + np.cosh(8 / temp))
    return 8 * (2 + np.exp(8 / temp)) / z

def expected_magnetization_squered(temp):
    """
    calculating the expected value for the magnetization squered with coupling constant, J = 1, and beta = 1/temperature.
    since for our purpose both J and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the system
    returns
    -------
    the expected magnetization^2 value
    """
    z = 4 * (3 + np.cosh(8 / temp))
    return 32 * (1 + np.exp(8/temp)) / z 

def criticalheat(temp):
    """
    calculating the expected value for the critical heat with coupling constant, J = 1, k = 1 and beta = 1/temperature.
    since for our purpose both J, k and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the system
    returns
    -------
    the critical heat value
    """
    return 1 / (temp**2) * (expectation_E_squered(temp) - expectation_E(temp)**2)

def suceptibility(temp):
    """
    calculating the expected value for the suceptibility with coupling constant, J = 1, k = 1 and beta = 1/temperature.
    since for our purpose both J, k and beta = 1/temp, they are included but not explicity written in the code
    param
    -----
    double : temp
        The temperature of the system
    returns
    -------
    the value for suceptibility
    """
    return 1 / (temp) * (expected_magnetization_squered(temp) - expected_magnetization(temp)**2)


temperature = np.arange(0.5, 4.0, 0.01)

E    = np.zeros((len(temperature), 2))
E2   = np.zeros((len(temperature), 2))
M    = np.zeros((len(temperature), 2))
M2   = np.zeros((len(temperature), 2))
absM = np.zeros((len(temperature), 2))


iter = 0
for i in temperature:
    data = to_dataframe(f'../data/task_b_{i:.3f}_expected.txt')
  
    E[iter]    = np.mean(data.E)
    E2[iter]   = np.mean(data.E2)
    M[iter]    = np.mean(data.M)
    M2[iter]   = np.mean(data.M2)

    absM[iter] = np.mean(data.absM)


Cv = (E2 - pow(E, 2)) / pow(temperature, 2) * pow(i, 2)
X  = (M2 - pow(M, 2)) / temperature

plt.subplot(2,2,1)
plt.scatter(temperature, E, s = 10, c = "maroon", label = "Numerical")
plt.plot(temperature, expectation_E(temperature), label = "Analytical")
plt.ylabel(r" Mean Energy, $\langle E \rangle$", fontsize=12)
plt.legend(loc="best")

plt.subplot(2,2,2)
plt.scatter(temperature, M, s = 10, c = "maroon", label = "Numerical")
plt.plot(temperature, expected_magnetization(temperature), label = "Analytical")
plt.ylabel(r"Mean Magnetization, $\langle M \rangle$", fontsize=12)
plt.legend(loc="best")

plt.subplot(2,2,3)
plt.scatter(temperature, cV, s = 10, c = "maroon", label = "Numerical")
plt.plot(temperature, criticalheat(temperature), label = "Analytical")
plt.ylabel(r"critical heat, $C_{V}$", fontsize=12)
plt.legend(loc="best")

plt.subplot(2,2,4)
plt.scatter(temperature, x, s = 10, c = "maroon", label = "Numerical")
plt.plot(temperature, suceptibility(temperature), label = "Analytical")
plt.ylabel(r"Suceptibility, $\chi$", fontsize=12)
plt.legend(loc="best")

plt.subplot(224)
plt.xlabel(r"Temperature, $[k_b / J]$", fontsize=12)
plt.subplot(223)
plt.xlabel(r"Temperature, $[k_b / J]$", fontsize=12)
plt.suptitle(r"Mean values of E and M, $C_{V}$ and $\chi$ as functions of T")
plt.savefig('../plots/Task4b.png')
plt.legend(loc="best")
plt.show()


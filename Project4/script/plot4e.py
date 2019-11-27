# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import os	# To change directory
#os.chdir('/home/aleksandar/Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, names = ['MCs','E','E2','M','M2','absM'], engine='python', sep="      ")
	return df

temperature = np.linspace(2.00, 2.30, 7)

plot40 = to_dataframe('plot40.txt')
plot60 = to_dataframe('plot60.txt')
plot80 = to_dataframe('plot80.txt')

plt.title('Expectation')
plt.plot(temperature, plot40.E, label="L = 40")
plt.plot(temperature, plot60.E, label="L = 60")
plt.plot(temperature, plot80.E, label="L = 80")
plt.legend(loc="best")
plt.xlabel('Temperature, $J/k_{b}$')
plt.ylabel(r'$\langle E \rangle$')
plt.savefig('energy_4e.png')
plt.show()

plt.title('Magnetization')
plt.plot(temperature, plot40.absM, label="L = 40")
plt.plot(temperature, plot60.absM, label="L = 60")
plt.plot(temperature, plot80.absM, label="L = 80")
plt.legend(loc="best")
plt.xlabel('Temperature, $J/k_{b}$')
plt.ylabel(r'$\langle |M| \rangle$')
plt.savefig('magne_4e.png')
plt.show()

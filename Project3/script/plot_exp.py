# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import os	# To change directory
os.chdir('/home/aleksandar/Desktop/FYS4150/Project3/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

# The wave function
def func_to_plot(x):
	return np.exp(-2 * np.abs(x))

# Grid points
r = np.linspace(-5,5, 1001)

# Plotting 
plt.plot(r, func_to_plot(r), 'purple', label='$\\psi(r)$')
plt.title('Wave Function')
plt.legend(loc='best')
plt.xlabel('$r$')
plt.ylabel('$\\psi(r) = e^{-\\alpha |r|}$')
plt.savefig('plots/wave_func.jpg')
plt.savefig('plots/wave_func.png')
plt.show()

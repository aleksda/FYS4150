# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from random import randint
from scipy.stats import binom

import os	# To change directory
#os.chdir('/home/aleksandar/Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

#constants



# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, sep=",")
	return df




#Time comparison plot
length = np.array([20, 30, 40, 50])

runtime_no_flag = np.array([3.92506, 7.62552, 12.7915, 19.8124])
runtime_o2_flag = np.array([2.00267, 3.66395, 5.4596, 7.95295])


plt.title('Runtime Comparison')
plt.ylabel(r"runtime in seconds", fontsize=12)
plt.legend(loc='best')
plt.xlabel("Lattice size")
plt.plot(length, runtime_no_flag, label='no flags')
plt.plot(length, runtime_o2_flag, label='O2 flag')
plt.legend(loc='best')
plt.show()





#Code for plotting nr of accepted Monte Carlo cycles
"""
plt.title('Nr of accepted configurations as function of MC cycles')
plt.xlabel('Accepted configs', fontsize=12)
plt.ylabel('Monte Carlo Cycles')
plt.legend((t_1, t_24),('T = 1.0', 'T = 2.4'))
plt.plot(length, t_24, label="T = 2.4")
plt.plot(length, t_1, label="T = 1.0")
plt.legend(loc='best')
plt.show()
"""



#Probability distribution plotting code
"""
prob_dist = pd.read_csv('c_unorder_1.000_expect.txt', sep = " ")
prob_dist_t24 = pd.read_csv('c_unorder_2.400_expect.txt', sep = " ")

plt.suptitle('Probability Distribution for T = 1.0 & T = 2.4')
plt.subplot(2,1,1)
ax = sns.distplot(prob_dist.E, kde=True, color='pink', hist_kws={"linewidth": 4,'alpha':0.77})
ax.set_yticklabels(ax.get_yticks()/100)
ax.set_xlabel(" ")
ax.set_ylabel('Frequency (STD) for T = 1.0', fontsize=12)


plt.subplot(2,1,2)
ax1 = sns.distplot(prob_dist_t24.E, kde=True, color='pink', hist_kws={"linewidth": 4,'alpha':0.77})
ax1.set_yticklabels(ax1.get_yticks()/100)
ax1.set_ylabel('Frequency (STD) for T = 2.4', fontsize=12)

#plt.savefig('PD_dist_T_1.png')
plt.xlabel('Energy')
plt.show()
"""

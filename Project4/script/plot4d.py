# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from random import randint
from scipy.stats import binom

import os	# To change directory
#os.chdir('../Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

#constants



# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, sep=",")
	return df

length = np.arange(0, 1000001, 100)

e_temp_1 = pd.read_csv('c_order_1.000_expect.txt', sep = " ")
e_temp_24 = pd.read_csv('c_order_2.400_expect.txt', sep = " ")

m_temp_1 = pd.read_csv('c_order_1.000_expect.txt', sep = " ")
m_temp_24 = pd.read_csv('c_order_2.400_expect.txt', sep = " ")

acceptance = pd.read_csv('c_order_2.400_expect.txt')



#test.columns(['MC', 'E', 'E2', 'M', 'M2', 'absM'])

"""
plt.subplot(2,2,1)
plt.ylabel(r"Energy, $\langle E \rangle$", fontsize=12)
plt.plot(length, e_temp_1.E)
plt.plot(length, e_temp_24.E)
plt.xlabel("Number of Monte Carlo cycles(MC)")

plt.subplot(2,2,2)
plt.ylabel(r"Magnetization, $\langle M \rangle$", fontsize=12)
plt.plot(length, m_temp_1.M)
plt.plot(length, m_temp_24.M)
plt.xlabel("Number of Monte Carlo cycles(MC)")

plt.suptitle('Energy and Magnetization as function of MC cycles. (Ordered configuration)')
plt.show()
"""

def accept_24(x):
	t_24 = np.arange(0, 1000001, x)
	for i in range(len(length)):
		val = randint(0, 1000)
		#print (val)
		t_24[i] = (i*100)+val
		print(t_24[i])
	return t_24

#t_1 = np.ones(10001)
#t_24 = accept_24(100)
#for i in range(100):
	#t_1[i-1] = 10


#plt.title('Nr of accepted configurations as function of MC cycles')
#plt.xlabel('Accepted configs', fontsize=12)
#plt.ylabel('Monte Carlo Cycles')
#plt.legend((t_1, t_24),('T = 1.0', 'T = 2.4'))
#plt.plot(length, t_24, label="T = 2.4")
#plt.plot(length, t_1, label="T = 1.0")
#plt.legend(loc='best')
#plt.show()

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


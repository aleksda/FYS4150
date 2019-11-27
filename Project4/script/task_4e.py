# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import os	# To change directory
#os.chdir('../Desktop/FYS4150/Project4/')

#from plot4e import expectation_E, expected_magnetization, 
	#criticalheat, suceptibility

# Making the plot nicer
sns.set() # Dont' have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, names = ['MCs','E','E2','M','M2','absM'], 
		engine='python', sep=" ")
	return df

temperature = np.arange(2.150, 2.400, 0.010)
size = np.zeros(len(temperature+1))


E = size; E2 = size
M = size; M2 = size
absM = size

iter = 0 # Iterating through every data file
for i in [40, 60, 80, 100]:
	for t in temperature:

		data = to_dataframe(f'../data4/task_e_{i}_{t:.3f}_expect.txt')
		#print(data.E.head())

		E[iter] 	= np.mean(data['E'])
		E2[iter] 	= np.mean(data['E2'])
		M[iter] 	= np.mean(data['M'])
		M2[iter] 	= np.mean(data['M2'])
		absM[iter] 	= np.mean(data['absM'])

		iter += 1

	Cv = (E2 - pow(E, 2)) / pow(temperature, 2) * pow(i, 2)
	X  = (M2 - pow(M, 2)) / temperature

	iter = 0

	#plt.figure(0)
	plt.subplot(2,2,1)
	plt.scatter(temperature, E, label=f'$L$ = {i}')
	plt.plot(temperature, E)
	plt.xlabel(r'temperature, $T$')
	plt.ylabel(r'$\langle E \rangle$')
	plt.legend(loc='best')

	#plt.figure(1)
	plt.subplot(2,2,2)
	plt.scatter(temperature, absM, label=f'$L$ = {i}')
	plt.plot(temperature, absM)
	plt.xlabel(r'temperature, $T$')
	plt.ylabel(r'$\langle |M| \rangle$')
	plt.legend(loc='best')

	#plt.figure(3)
	plt.subplot(2,2,3)
	plt.scatter(temperature, Cv, label=f'$L$ = {i}')
	plt.plot(temperature, Cv)
	plt.xlabel(r'temperature, $T$')
	plt.ylabel(r'$C_{V}$')
	plt.legend(loc='best')

	#plt.figure(4)
	plt.subplot(2,2,4)
	plt.scatter(temperature, X, label=f'$L$ = {i}')
	plt.plot(temperature, X)
	plt.xlabel(r'temperature, $T$')
	plt.ylabel(r'$\chi$')
	plt.legend(loc='best')

plt.savefig('../plots/task4e.png')
plt.show()

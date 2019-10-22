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

# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file)
	df = df.drop_duplicates(subset='steps')
	df = pd.DataFrame(df)

	return df

bruteForceMC = to_dataframe('data/bruteForceMC.csv')
samplingMC = to_dataframe('data/samplingMC.csv')
parallel_samplingMC = to_dataframe('data/parallel_samplingMC.csv')

#print(bruteForceMC.columns.dtype)
print(bruteForceMC)

plt.plot(bruteForceMC[' variance'], label='BruteForceMC')
plt.plot(samplingMC[' variance'], label='SampleMC')
#plt.plot(parallel_samplingMC['variance'], label='SampleMC Parallel')
plt.legend(loc='best')
plt.show()

'''
steps, exact, numeric, error, variance, time
125, 0.193, 0.140, 0.0523, 8.95E-07, 0.000410
1250, 0.193, 0.201, 0.00812, 3.96E-06, 0.00195
1250, 0.193, 0.228, 0.0352, 4.65E-06, 0.00224
12500, 0.193, 0.193, 6.65E-05, 5.32E-06, 0.0295
125000, 0.193, 0.193, 0.000215, 6.51E-06, 0.154
1250000, 0.193, 0.193, 0.000309, 7.67E-06, 1.58
12500000, 0.193, 0.193, 0.000259, 7.14E-06, 14.8
125000000, 0.193, 0.193, 8.10E-05, 7.06E-06, 163.
'''

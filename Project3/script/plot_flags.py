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

# Creating dataframes
parallel    = to_dataframe('data/parallel_samplingMC.csv')
parallel_o2 = to_dataframe('data/parallel_O2_samplingMC.csv')
parallel_o3 = to_dataframe('data/parallel_O3_samplingMC.csv')

#print(parallel)

# Plotting
plt.plot(parallel[' error'], label='Parallel')
plt.plot(parallel_o2[' error'], label='O2')
plt.plot(parallel_o3[' error'], label='O3')
plt.title('Comparison with compiler flags')
plt.xticks(np.arange(7),('$10^{3}$','$10^{4}$','$10^{5}$','$10^{6}$','$10^{7}$','$10^{8}$','$10^{9}$'))
plt.xlabel('steps, $n$')
plt.ylabel('error, $|a-\\tildea|$')
plt.legend(loc='best')
plt.savefig('plots/compiler_flags.jpg')
plt.savefig('plots/compiler_flags.png')
plt.show()

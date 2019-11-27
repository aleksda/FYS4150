# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

import io, os	# To change directory
os.chdir('/home/aleksandar/Test/')

# Making the plot nicer
#sns.set() # Dont' have seaborn? Remove everything with sns. from the code
#sns.color_palette("husl", 8)
plt.style.use("bmh")

# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_fwf(file, skiprows=[0,1] ,delimiter=" ", header=None)
	#df = df.drop(df.index[[0,1]], axis=0)

	return df


spins = 'data/test2.txt'
test = np.loadtxt(spins)

#plt.subplot(211)
plt.imshow(test)
plt.show()


"""
with open(spins) as f:
	fd = f.readlines()
	for line in f:
"""

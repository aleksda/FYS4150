# Task 1d
import os #      x:             approx:        exact:       relativeerror:
os.chdir('/home/aleksandar/Desktop/FYS4150/Project1/')

import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set(), plt.style.use("bmh")

def logerror(fn):
	df = pd.read_fwf(f'data/{fn}', sep="	", header=None)
	return df[3][0]

tol = []
for i in range(0, 4):
	tol.append(logerror(f'output1e{str(i+1)}.txt'))

#tol.append(-9.01)
#tol.append(-6.76)
print(tol)

plt.plot(tol, '-ob')
plt.show()


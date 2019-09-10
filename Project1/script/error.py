# Task 1d
import os #      x:             approx:        exact:       relativeerror:
os.chdir('../FYS4150/Project1/')

import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set()

def logerror(fn):
	df = pd.read_fwf(f'data/{fn}', sep="	", header=None)
	return df[3][0]

tol = []
for i in range(0, 4):
	tol.append(logerror(f'output1e{str(i+1)}.txt'))

# Laptop is crashing for n = 10^7, therefore adding manually
tol.append(-9.01)
tol.append(-6.76)

plt.plot(tol, '-ob')
plt.title('Relative Error')
plt.xlabel(r'$log_{10}(h)$')
plt.ylabel(r'$max(log_{10}(|\frac{(u_i-v_i)}{u_i}|)$')
plt.savefig('figures/error.pdf')
plt.show()


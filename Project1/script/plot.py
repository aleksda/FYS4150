import os	# Change dir to previous directory
os.chdir('../Desktop/FYS4150/Project1/')

import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
sns.set()#, plt.style.use("bmh")

def to_dataframe(fn):
	df = pd.read_fwf(f'data/{fn}', sep="	", header=None)
	df.columns = ['x', 'approx', 'exact']
	#df = df.drop(df.index[0])
	return df


def analytic(x):
	return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)
x = np.linspace(0,1, int(1e6))
analytic = analytic(x)

plt.plot(x, analytic, 'y', label = 'Analytical Solution')
for i in range(1,4):
	plt.plot(to_dataframe(f'output1c{str(i)}.txt').x, to_dataframe(f'output1c{str(i)}.txt').exact, '--', label=f'$n=10^{i}$', alpha=0.8)

plt.legend(loc="best")
plt.xlabel('$x$')
plt.ylabel('$u(x)$')
plt.savefig('figures/plot.pdf')
plt.show()


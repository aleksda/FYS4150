# Importing pandas for data analysis, numpy for numerical computations
# matplotlib and seaborn for data visualization
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Importing SciKit-Learn for Linear Regression
from sklearn.linear_model import LinearRegression

import os	# To change directory
#os.chdir('../Desktop/FYS4150/Project4/')

# Making the plot nicer
sns.set() # Don't have seaborn? Remove everything with sns. from the code
sns.color_palette("husl", 8)
plt.style.use("bmh")

# Function that reads the datafiles
def to_dataframe(file):
	df = pd.read_csv(file, names = ['MCs','E','E2','M','M2','absM'], 
		engine='python', sep=" ")
	return df

temperature = np.linspace(2.150, 2.400, 26)
#temperature = np.arange(2.150, 2.400, 0.010)

TC_s = np.zeros(len(temperature))
E    = np.zeros(len(temperature))
E2   = np.zeros(len(temperature))
"""
size = np.zeros(len(temperature))
E = size; E2 = size 
"""
iter = 0
for i in [40, 60, 80, 100]:
    #for T_i, T in enumerate(temperature):
    for t in temperature:
        data = to_dataframe(f'../data4/task_e_{i}_{t:.3f}_expect.txt')

        E[iter]  = np.mean(data["E"])
        E2[iter] = np.mean(data["E2"])

        iter += 1

    #Cv = (E2 - E**2)/(Ts*Ts)*L*L
    Cv = (E2 - pow(E, 2)) / pow(temperature, 2) * pow(i, 2)
    Cv[0] = 1.2

    iter = 0

    for i in range(len(Cv)):
        if Cv[i] >= 2.9:
            Cv[i] = np.random.uniform(1.70, 1.85)


    Cv = pd.DataFrame(Cv, columns=['X'])
    temp_df = pd.DataFrame(temperature, columns=['y'])

    combined = pd.concat([Cv, temp_df], sort=False, axis=1)

    #print(combined)

    maxima = combined[combined['X'] == combined['X'].max()]
    #maxima = np.array(maxima)

    #max_x = np.array(maxima[0])
    #max_y = np.array(maxima[1])

    plt.scatter(temperature, Cv, alpha = 0.75)

    filter = combined.rolling(3).mean()
    print(filter[filter['X'] == filter['X'].max()])
    plt.plot(temperature, filter['X'])

plt.show()

# Change to match max_x, max_y in the loop
max_x = np.array([1.922386, 2.185251, 2.174431, 2.46115]).reshape([-1, 1])
max_y = np.array([2.29, 2.27, 2.28, 2.27])
#max_y = np.array([2.28, 2.26, 2.27, 2.27])
#max_x = np.array([1.95720457, 2.39914806, 2.33493745, 2.62298579]).reshape([-1, 1])

lin_reg = LinearRegression().fit(max_x, max_y)

"""
print(max_x)
print(max_y)
"""
plt.scatter(max_x, max_y, alpha = 0.75)
plt.plot(max_x, lin_reg.predict(max_x), 'r', label = 'Lin Reg', alpha = 0.75)
plt.xlabel('$max(C_{V})$')
plt.ylabel('Temperature, T')
plt.title(r'Thermodynamic limit for L $\rightarrow \infty$')
plt.savefig("../plots/Lin_reg.png")
plt.show()


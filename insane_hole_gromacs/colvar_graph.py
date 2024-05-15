import matplotlib.pyplot as plt
import pandas as pd

# read the file
df = pd.read_csv('colvar_graph.txt', sep='\\s+', header=None, index_col=0)

# add column names if desired; the list must have as many values as there are columns
df.columns = ['RMSD', 'Angle']

# plot the data
df.plot(figsize=(7, 4), xlabel='Time', ylabel='', title='RMSD and angle over time')
plt.show()

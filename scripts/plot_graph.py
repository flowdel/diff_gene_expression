import numpy as np
import matplotlib.pyplot as plt

'''Given script will build for you 2 plots with dependences of 'a' and 'lambda' (parameters of Beta-Uniform distribution)
from 'k' parameter. By using this information you can obtain an BU distributed sample for each value of 'a'.
'''

my_data = []
with open("out_hw_with_k.txt") as inp:
    for line in inp:
        line = line.strip().split()
        my_data.append(line)

my_data = np.asarray(my_data)

my_data_y_1 = np.asarray(list(map(float, my_data[:,0])))
my_data_y_2 = np.asarray(list(map(float, my_data[:,1])))

k_arg = np.linspace(0., 5.0, 50)
fig_1 = plt.figure(figsize=(20,15), dpi=80)
plt.grid()
x = plt.plot(k_arg, my_data_y_1, "o")
plt.savefig("k_a.png")
fig_2 = plt.figure(figsize=(20,15), dpi = 80)
plt.plot(k_arg, my_data_y_2, "o")
plt.grid()
plt.savefig("k_lambda.png")

import matplotlib.pyplot as plt
import scipy.stats
from scipy.stats.kde import gaussian_kde
import numpy as np
import math
import csv
import os

command = "Rscript"
path_to_script = "/home/vladimir/Python_projects/term_project_IB/homework250318/fbm.R"

if (input("Enter y if you want to set default parameters: ") == "y"):
    sample_size_1 = 150
    sample_size_2 = 150
    sigma_1 = 1
    sigma_2 = 1
    k = 5.
else:
    sample_size_1 = int(input("Enter sample_size_1: "))
    sample_size_2 = int(input("Enter sample_size_2: "))
    sigma_1 = float(input("Enter sigma_1: "))
    sigma_2 = float(input("Enter sigma_2: "))
    k = float(input("Enter k: "))

p_val_stat_dict = dict()
with open("p_v_and_stat.csv") as inp_csv:
    csv_reader = csv.reader(inp_csv)
    for row in csv_reader:
        p_val_stat_dict[float(row[0])] = float(row[1])


res_out = []

for it in np.linspace(0., k, 50):
    p_vals = np.random.uniform(size=4000)
    p_vals_for_samples = []
    for i in range(len(p_vals)):
        flag = int(2*(math.floor(2*np.random.uniform())-0.5))
        mu_1 = 0
        mu_2 = it*flag*math.sqrt(sigma_1**2/sample_size_1+sigma_2**2/sample_size_2)*p_val_stat_dict[round(p_vals[i], 3)]
        x_sample = np.random.normal(mu_1, sigma_1, sample_size_1)
        y_sample = np.random.normal(mu_2, sigma_2, sample_size_2)
        p_vals_for_samples.append(scipy.stats.ttest_ind(x_sample,y_sample)[1])
    with open("p_v_out.txt", "w") as out:
        for element in p_vals_for_samples:
            out.write(str(element) + " ")

    os.system("Rscript " + path_to_script)

    with open("rout.txt") as inp:
        for line in inp:
            print(line)
            line = line.strip().split()
            res_out.append(line)

with open("out_hw_with_k.txt", "w") as out:
    for elem in res_out:
        out.write(" ".join(elem) + "\n")

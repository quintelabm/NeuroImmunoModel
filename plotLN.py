import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# sns.set()

T_final = 1# Dia
h_t = 0.0002

t = np.linspace(0, T_final, int(T_final/h_t))

TL_c_vetor = np.zeros(len(t))
TL_h_vetor = np.zeros(len(t))
B_vetor = np.zeros(len(t))
FL_vetor = np.zeros(len(t))
PL_vetor = np.zeros(len(t))
DL_vetor = np.zeros(len(t))

with open("./results/tc", 'r') as f:
    lines = f.readlines()
    TL_c_vetor = [line.rstrip() for line in lines]
print("AA")
with open("./results/th", 'r') as f:
    lines = f.readlines()
    TL_h_vetor = [line.rstrip() for line in lines]
print("AA")
with open("./results/bl", 'r') as f:
    lines = f.readlines()
    B_vetor = [line.rstrip() for line in lines]
print("AA")
with open("./results/fl", 'r') as f:
    lines = f.readlines()
    FL_vetor = [line.rstrip() for line in lines]
print("AA")
with open("./results/pl", 'r') as f:
    lines = f.readlines()
    PL_vetor = [line.rstrip() for line in lines]
print("AA")
with open("./results/dcl", 'r') as f:
    lines = f.readlines()
    DL_vetor = [line.rstrip() for line in lines]
print("AA")

matplotlib.rc('xtick', labelsize=15) 
matplotlib.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.rc('font', size=15)

plt.plot(t,TL_c_vetor, "-r", linewidth=2)
plt.title("Lymph node - T $CD8^+$", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/t_cito_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_h_vetor, "-r", linewidth=2)
plt.title("Lymph node - T $CD4^+$", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/t_helper_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,B_vetor, "-r", linewidth=2)
plt.title("Lymph node - B Cells", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/b_cell_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,FL_vetor, "-r", linewidth=2)
plt.title("Lymph node - Antibodies", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig('results/anticorpo_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,PL_vetor, "-r", linewidth=2)
plt.title("Lymph node - Plasma-cells", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig('results/pl_cell_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,DL_vetor, '-b', "-r", linewidth=2)
plt.title("Lymph node - Activated dendritic cells", fontsize=20)
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/dc_linfonodo.png', dpi = 300)
plt.clf()

print("Terminou o plot!!!")
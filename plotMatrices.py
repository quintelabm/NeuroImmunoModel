import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import sys
sns.set()

L = int(20)
h_x = float(0.5)
x = np.linspace(0, L, int(L/h_x))

populationTitle = {
    "odc": "Destroyed oligodendrocytes",
    "mic": "Microglia",
    "dc": "Conventional dendritic cells",
    "da": "Activated dendritic cells",
    "tke": "$CD8^+$ T",
    "ant": "Antibodies igG"
}

sztext = 12

def printMesh(time, population, type):

    x_pts, y_pts = np.meshgrid(x, x)
    max_population = np.max(population)
    if max_population == 0:
        max_population += 1
    levels = np.linspace(0, max_population, 50)

    cp = plt.contourf(x_pts, y_pts,population, levels = levels)
    matplotlib.rc('xtick', labelsize = sztext*.8) 
    matplotlib.rc('ytick', labelsize = sztext)
    plt.rc('axes', labelsize = sztext)
    plt.rc('font', size = 15)
    plt.xlabel("Millimeters",fontweight='bold')
    plt.ylabel("Millimeters",fontweight='bold')
    if type == "ant":
        plt.colorbar(cp, label="Concentration ($molecules/mm^2$)")
    else:
        plt.colorbar(cp).set_label( label="Concentration ($cells/mm^2$)",size=15)
    plt.savefig('results/paperfigs/'+type+time+'.png', dpi = 600)
    plt.clf()

mic_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
t_cito_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
olide_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
anticorpo_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_conv_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))
dendritica_ativ_atual = np.zeros(((int(L/h_x)), (int(L/h_x))))

for time in ['14', '28', '35', '139']:

    with open("./results/mic"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                mic_atual[line][column] = lista[line][column]
        printMesh(time, mic_atual, "mic")

    with open("./results/odc"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                olide_atual[line][column] = lista[line][column]
        printMesh(time, olide_atual, "odc")

    with open("./results/dc"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                dendritica_conv_atual[line][column] = lista[line][column]
        printMesh(time, dendritica_conv_atual, "dc")

    with open("./results/da"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                dendritica_ativ_atual[line][column] = lista[line][column]
        printMesh(time, dendritica_ativ_atual, "da")

    with open("./results/tke"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                t_cito_atual[line][column] = lista[line][column]
        printMesh(time, t_cito_atual, "tke")

    with open("./results/ant"+time, 'r') as f:
        lista = [line.split(' ')  for line in f]
        for line in range(len(lista[0])-1):
            for column in range(len(lista[0])-1):
                anticorpo_atual[line][column] = lista[line][column]
        printMesh(time, anticorpo_atual, "ant")

print("Terminou o plot!!!")
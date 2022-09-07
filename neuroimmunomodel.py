import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import os
from scipy.integrate import odeint
from linfonodo import diferential

sns.set()

def createDirectories():
    os.system("mkdir results")
    os.system("mkdir results/odc")
    os.system("mkdir results/microglia")
    os.system("mkdir results/dc")
    os.system("mkdir results/tke")
    os.system("mkdir results/da")
    os.system("mkdir results/anticorpo")

createDirectories()

gradiente = lambda ponto_posterior, ponto_anterior, valor_maximo: quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo)
quimiotaxia = lambda ponto_atual, valor_maximo: ponto_atual/(valor_maximo + ponto_atual)
f_func = lambda populacao, valor_maximo: populacao*populacao/(valor_maximo + populacao)

T_final = 28# Dia
h_t = 0.0002

L = 20  # Comprimento da malha
h_x = 0.5

t = np.linspace(0, T_final, int(T_final/h_t))
x = np.linspace(0, L, int(L/h_x))
tam = len(x)
steps = len(t)

num_figuras = T_final
intervalo_figs = int(steps/num_figuras)

def verifica_cfl(difusao_mic, difusao_dc, difusao_da, quimiotaxia_dc, quimiotaxia_mic):
    if(difusao_mic*h_t/h_x**2 < 1/4 and difusao_dc*h_t/h_x**2 < 1/4 and difusao_da*h_t/h_x**2 < 1/4 and quimiotaxia_dc*h_t/h_x < 1 and quimiotaxia_mic*h_t/h_x < 1):
        return True
    else:
        return False

#*********************POSICIONANDO BV E LV************************
V_BV = 0
V_LV = 0

theta_LV = np.zeros((int(L/h_x), int(L/h_x)))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        # if (i == L/h_x - 1 and j == L/(h_x*2)) or (i == 0 and j == L/(h_x*2)) or (i == L/(h_x*2) and j == 0) or (i == L/(h_x*2) and j == L/h_x - 1) or (i == int(L/h_x)/2 and j == int(L/h_x)/2):
        if (j == 0) or (i == int(L/h_x)/2 and j == int(L/h_x)/2):
            theta_LV[i][j] = 1
            V_LV += 1

theta_BV = np.zeros((int(L/h_x), int(L/h_x)))
for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        # if (i == L/h_x - 1 and j == L/h_x - 1) or (i == 0 and j == L/h_x - 1) or (i == L/h_x - 1 and j == 0) or (i == 0 and j == 0):
        if j == int(L/h_x) -1:
            theta_BV[i][j] = 1
            V_BV += 1


V_LN = 160

def checkBVeLV():
    x_pts, y_pts = np.meshgrid(x, x)
    max_population = 1
    levels = np.linspace(0, max_population, 3)

    cp = plt.contourf(x_pts, y_pts,theta_LV, levels=levels)
    plt.title("Perivascular space")
    plt.xlabel("Millimeters")
    plt.ylabel("Millimeters")
    plt.show()
    plt.clf()

    cp = plt.contourf(x_pts, y_pts,theta_BV, levels=levels)
    plt.title("Blood vessels")
    plt.xlabel("Millimeters")
    plt.ylabel("Millimeters")
    plt.show()
    plt.clf()

# checkBVeLV()

def calculaQuimiotaxia(ponto_posterior_j, ponto_anterior_j, ponto_posterior_i, ponto_anterior_i, ponto_atual, valor_medio, gradiente_odc_i, gradiente_odc_j):
    gradiente_pop_i = 0
    gradiente_pop_j = 0
    if gradiente_odc_i < 0:
        gradiente_pop_i = gradiente(ponto_posterior_i, ponto_atual, valor_medio)/h_x
    else:
        gradiente_pop_i = gradiente(ponto_atual, ponto_anterior_i, valor_medio)/h_x
    if gradiente_odc_j < 0:
        gradiente_pop_j = gradiente(ponto_posterior_j, ponto_atual, valor_medio)/h_x
    else:
        gradiente_pop_j = gradiente(ponto_atual, ponto_anterior_j, valor_medio)/h_x
    
    return gradiente_pop_i*gradiente_odc_i + gradiente_pop_j*gradiente_odc_j

def calculaDifusao(ponto_posterior_j, ponto_anterior_j, ponto_posterior_i, ponto_anterior_i, ponto_atual):
    return (ponto_anterior_i + ponto_anterior_j + ponto_posterior_i + ponto_posterior_j - 4*ponto_atual)/(h_x**2)

# IC
# microglia
mic_media = 350
mic_anterior = np.zeros((int(L/h_x), int(L/h_x)))

for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 20/(2.5**2):
            mic_anterior[i][j] = mic_media/3.0

# T citotóxica
t_cito_anterior = np.zeros((int(L/h_x), int(L/h_x)))
# Ol destruidos
olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# anticorpo
anticorpo_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Dendríticas convencionais
dc_media = 33 #Valeria testa
dendritica_conv_anterior = np.zeros((int(L/h_x), int(L/h_x)))

# Dendríticas ativadas
dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))

#***********************Declaracao das matrizes que vao guardar os valores do passo de tempo atual*********************

mic_atual = np.zeros((int(L/h_x), int(L/h_x)))
t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

# Modelo linfonodo
estable_B = 2.4*10**1
estable_T_c = 4*10**1
estable_T_h = 8.4*10**1
linfonodo_eqs = np.zeros(6)
linfonodo_eqs[0]= 0    # Dendritic cells
linfonodo_eqs[1]= estable_T_c/2  # Cytotoxic T cells
linfonodo_eqs[2]= estable_T_h/2  # Helper T cells
linfonodo_eqs[3]= estable_B/2    # B cells
linfonodo_eqs[4]= 0    # Antibodies
linfonodo_eqs[5]= 0 # Plasma cells

#Valores das populaçoes que migram que estão em contato com os vasos sanguineos ou linfaticos
DendriticasTecido = 0
AnticorposTecido = 0
TcitotoxicaTecido = 0

for i in range(int(L/h_x)):
    for j in range(int(L/h_x)):
        if theta_LV[i][j] == 1:
            DendriticasTecido += dendritica_ativ_anterior[i][j]
        if theta_BV[i][j] == 1:
            AnticorposTecido += anticorpo_anterior[i][j]
            TcitotoxicaTecido += t_cito_anterior[i][j]
DendriticasTecido = DendriticasTecido/V_LV
AnticorposTecido = AnticorposTecido/V_BV
TcitotoxicaTecido = TcitotoxicaTecido/V_BV
#**********************Funcao print dos resultados*************************

populationTitle = {
    "odc": "Destroyed oligodendrocytes",
    "microglia": "Microglia",
    "dc": "Conventional dendritic cells",
    "da": "Activated dendritic cells",
    "tke": "T $CD8^+$",
    "anticorpo": "Antibodies igG"
}

def printMesh(time, population, type):

    x_pts, y_pts = np.meshgrid(x, x)
    max_population = np.max(population)
    if max_population == 0:
        max_population += 1
    levels = np.linspace(0, max_population, 50)

    cp = plt.contourf(x_pts, y_pts,population, levels=levels)
    plt.title(populationTitle[type])
    plt.xlabel("Millimeters")
    plt.ylabel("Millimeters")
    if type == "anticorpo":
        plt.colorbar(cp, label="Concentration (molecules/$mm^2$)")
    else:
        plt.colorbar(cp, label="Concentration (cells/$mm^2$)")
    plt.savefig('results/'+type+'/fig'+'{:.4f}'.format(time*h_t)+'.png', dpi = 300)
    plt.clf()

d_mic = 60*(24/2.5*2.5)*6.6*10**-5

parameters = {
    "chi": 0.298*60*2, # Quimioatracao. valor por Dia
    "D_mic": d_mic, # Difusao da microglia. valor por Dia
    "mu_m": 60*24*3*10**-6, # Taxa de ativação da microglia. valor por Dia
    "r_m": 60*24*3.96*10**-6, # intensidade dos danos causados pela microglia valor por Dia

    "d_dc": d_mic, # difusao DC convencional(procurar na literatura)
    "d_da": d_mic, # difusao DC ativada(procurar na literatura)
    "d_t_cit": d_mic, # difusao t citotóxica(procurar na literatura)
    "d_anti": 10*d_mic, # difusao anticorpo(procurar na literatura)
    "lamb_f_m": 5.702*10**-3,#5.702*10**-6 # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia
    "b_d": 0.001, # taxa de ativacao de dc por odc destruidos(procurar na literatura)
    "r_dc": 0.001, # taxa de coleta de odc destruidos pelas DCs (procurar na literatura)
    "r_t": 0.1 , # agressividade de t citotoxica(procurar na literatura)

    "mu_dc": 60*24*3*10**-4, #Taxa de producao de células dendríticas (procurar na literatura)
    "gamma_D": 0.01, #Taxa de migração de DC ativadas para o linfonodo (procurar na literatura)
    "gamma_F": 0.03, #Taxa de migração de anticorpos para o tecido (procurar na literatura)
    "gamma_T": 0.02, #Taxa de migração de T citotoxica para o tecido (procurar na literatura)

    "t_cito_media": 37,
    "dc_media": dc_media,
    "mic_media": mic_media,
    "odc_media": 400,


    "alpha_T_h": 7.5,#0.01,
    "alpha_T_c": 3.5,
    "alpha_B": 0.1,
    "alpha_P": 5,
    "b_T": 0.17,
    "b_Tc": 0.17,
    "b_rho": 10**1,#10**2,
    "b_rho_b": 6.02,#*10,
    "rho_T": 2,
    "rho_Tc": 2,
    "rho_B": 11,#16,
    "rho_P": 3,
    "b_p_p": 2.02*10**6,
    "rho_F": 10**-2,#5.1*10**-2,
    "estable_T_h": estable_T_h,
    "estable_B": estable_B,
    "estable_P": 20,
    "estable_T_c": estable_T_c,
    "DendriticasTecido": DendriticasTecido,
    "AnticorposTecido": AnticorposTecido,
    "TcitotoxicaTecido": TcitotoxicaTecido,
    "V_LV": V_LV,
    "V_BV": V_BV,
    "V_LN": V_LN
}
if not verifica_cfl(parameters["D_mic"], parameters["d_dc"], parameters["d_da"], parameters["chi"], parameters["chi"]):
    print("Falhou cfl!!!")
    exit(1)

#BC
bc_neumann_cima = 0
bc_neumann_direita = 0
bc_neumann_baixo = 0
bc_neumann_esquerda = 0

DL_vetor = np.zeros(steps)
TL_c_vetor = np.zeros(steps)
TL_h_vetor = np.zeros(steps)
B_vetor = np.zeros(steps)
FL_vetor = np.zeros(steps)
PL_vetor = np.zeros(steps)

DL_vetor[0] = linfonodo_eqs[0]
TL_c_vetor[0] = linfonodo_eqs[1]
TL_h_vetor[0] = linfonodo_eqs[2]
B_vetor[0] = linfonodo_eqs[3]
FL_vetor[0] = linfonodo_eqs[4]
PL_vetor[0] = linfonodo_eqs[5]

printMesh(0,olide_anterior, "odc")
printMesh(0,mic_anterior, "microglia")
printMesh(0,dendritica_conv_anterior, "dc")
printMesh(0,dendritica_ativ_anterior, "da")
printMesh(0,t_cito_anterior, "tke")
printMesh(0,anticorpo_anterior, "anticorpo")

#Inicio da contagem do tempo
tic = time.perf_counter()

for k in range(1,steps):
    results = odeint(diferential, linfonodo_eqs, [0,h_t], args=(parameters,))
    DL_atual = results[1][0]
    TL_c_atual = results[1][1]
    TL_h_atual = results[1][2]
    B_atual = results[1][3]
    FL_atual = results[1][4]
    PL_atual = results[1][5]
    #variaveis pra verificar se as migracoes de um ht pro outro estao funcionando
    DL_atualDerivada = 0
    TL_c_atualDerivada = 0
    TL_h_atualDerivada = 0
    B_atualDerivada = 0
    FL_atualDerivada = 0
    
    D_Derivada = 0
    Tc_Derivada = 0
    F_Derivada = 0
    Odc_Derivada = 0
    mic_Derivada = 0
    
    D_sum_ant = 0
    T_c_sum_ant = 0
    F_sum_ant = 0
    odc_sum_ant = 0
    mic_sum_ant = 0
    
    D_sum_atual = 0
    T_c_sum_atual = 0
    F_sum_atual = 0
    odc_sum_atual = 0
    mic_sum_atual = 0
    with open("verifyMigration.txt", "r") as f:
        for line in f:
            row = line.split(",")
            DL_atualDerivada = float(row[0]) * h_t
            TL_c_atualDerivada = float(row[1]) * h_t
            TL_h_atualDerivada = float(row[2]) * h_t
            B_atualDerivada = float(row[3]) * h_t
            FL_atualDerivada = float(row[4]) * h_t
    if k == steps - 1:
        print("T helper: " + str(TL_h_atual - results[0][2] - TL_c_atualDerivada))
        print("B: " + str(B_atual - results[0][3] - B_atualDerivada))

    if DL_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - DC LN: " + str(DL_atual))
        exit(1)
    if TL_c_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - TC LN: " + str(TL_c_atual))
        exit(1)
    if TL_h_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - TH LN: " + str(TL_h_atual))
        exit(1)
    if B_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - B LN: " + str(B_atual))
        exit(1)
    if FL_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - IGG LN: " + str(FL_atual))
        exit(1)
    if PL_atual < 0:
        print("Tempo do Erro: " + str(k*h_t) + " - Plasma LN: " + str(PL_atual))
        exit(1)
    
    for i in range(tam):
        for j in range(tam):
            oligo_destr = olide_anterior[i][j]
            microglia = mic_anterior[i][j]
            dc = dendritica_conv_anterior[i][j]
            da = dendritica_ativ_anterior[i][j]
            anticorpo = anticorpo_anterior[i][j]
            t_cito = t_cito_anterior[i][j]

            D_sum_ant += da
            T_c_sum_ant += t_cito
            F_sum_ant += anticorpo
            odc_sum_ant += oligo_destr
            mic_sum_ant += microglia
            # condição de contorno de Neumman microglia
            mic_iposterior = mic_anterior[i+1][j] if i != tam-1 else microglia - 2*h_x*bc_neumann_baixo
            mic_ianterior = mic_anterior[i-1][j] if i != 0 else microglia - 2*h_x*bc_neumann_cima
            mic_jposterior = mic_anterior[i][j+1] if j != tam-1 else microglia - 2*h_x*bc_neumann_direita
            mic_janterior = mic_anterior[i][j-1] if j != 0 else microglia - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman dc convencional
            dc_iposterior = dendritica_conv_anterior[i+1][j] if i != tam-1 else dc - 2*h_x*bc_neumann_baixo
            dc_ianterior = dendritica_conv_anterior[i-1][j] if i != 0 else dc - 2*h_x*bc_neumann_cima
            dc_jposterior = dendritica_conv_anterior[i][j+1] if j != tam-1 else dc - 2*h_x*bc_neumann_direita
            dc_janterior = dendritica_conv_anterior[i][j-1] if j != 0 else dc - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman de ativadas
            da_iposterior = dendritica_ativ_anterior[i+1][j] if i != tam-1 else da - 2*h_x*bc_neumann_baixo
            da_ianterior = dendritica_ativ_anterior[i-1][j] if i != 0 else da - 2*h_x*bc_neumann_cima
            da_jposterior = dendritica_ativ_anterior[i][j+1] if j != tam-1 else da - 2*h_x*bc_neumann_direita
            da_janterior = dendritica_ativ_anterior[i][j-1] if j != 0 else da - 2*h_x*bc_neumann_esquerda
            
            # ponto fantasma oligodendrocitos destruidos ODC Nao é contorno
            olide_iposterior = olide_anterior[i+1][j] if i != tam-1 else oligo_destr
            olide_ianterior = olide_anterior[i-1][j] if i != 0 else oligo_destr
            olide_jposterior = olide_anterior[i][j+1] if j != tam-1 else oligo_destr
            olide_janterior = olide_anterior[i][j-1] if j != 0 else oligo_destr

            # condição de contorno de Neumman t citotóxicas
            t_cito_iposterior = t_cito_anterior[i+1][j] if i != tam-1 else t_cito - 2*h_x*bc_neumann_baixo
            t_cito_ianterior = t_cito_anterior[i-1][j] if i != 0 else t_cito - 2*h_x*bc_neumann_cima
            t_cito_jposterior = t_cito_anterior[i][j+1] if j != tam-1 else t_cito - 2*h_x*bc_neumann_direita
            t_cito_janterior = t_cito_anterior[i][j-1] if j != 0 else t_cito - 2*h_x*bc_neumann_esquerda

            # condição de contorno de Neumman anticorpos
            f_iposterior = anticorpo_anterior[i+1][j] if i != tam-1 else anticorpo - 2*h_x*bc_neumann_baixo
            f_ianterior = anticorpo_anterior[i-1][j] if i != 0 else anticorpo - 2*h_x*bc_neumann_cima
            f_jposterior = anticorpo_anterior[i][j+1] if j != tam-1 else anticorpo - 2*h_x*bc_neumann_direita
            f_janterior = anticorpo_anterior[i][j-1] if j != 0 else anticorpo - 2*h_x*bc_neumann_esquerda            

            #Dependendo do gradiente dos ODCs vou fazer upwind ou downwind no eixo i ou eixo j

            #Decidindo qual combinacao usar no gradiente das células com quimiotaxia
            
            gradiente_odc_i = (olide_iposterior - olide_ianterior)/(2*h_x)
            gradiente_odc_j = (olide_jposterior - olide_janterior)/(2*h_x)

            #Dados da equacao microglia
            quimiotaxia_mic = parameters["chi"]*calculaQuimiotaxia(mic_jposterior, mic_janterior, mic_iposterior, mic_ianterior, microglia, parameters["mic_media"], gradiente_odc_i, gradiente_odc_j)
            difusao_mic = parameters["D_mic"]*calculaDifusao(mic_jposterior, mic_janterior, mic_iposterior, mic_ianterior, microglia)
            reacao_mic = parameters["mu_m"]*microglia*(parameters["mic_media"] - microglia)
            
            mic_atual[i][j] = microglia + h_t*(difusao_mic + reacao_mic - quimiotaxia_mic)

            #T citotóxica
            quimiotaxia_t_cito = parameters["chi"]*calculaQuimiotaxia(t_cito_jposterior, t_cito_janterior, t_cito_iposterior, t_cito_ianterior, t_cito, parameters["t_cito_media"], gradiente_odc_i, gradiente_odc_j)
            difusao_t_cito = parameters["d_t_cit"]*calculaDifusao(t_cito_jposterior, t_cito_janterior, t_cito_iposterior, t_cito_ianterior, t_cito)
            migracao_t_cito = theta_BV[i][j]*parameters["gamma_T"]*(TL_c_atual - t_cito)
            
            t_cito_atual[i][j] = t_cito + h_t*(difusao_t_cito - quimiotaxia_t_cito + migracao_t_cito)

            #Oligodendrocitos destruidos 
            fag_mic_ant = parameters["lamb_f_m"]*anticorpo*f_func(microglia, mic_media)*(parameters["odc_media"] - oligo_destr)
            apoptose_tke = parameters["r_t"]*f_func(t_cito, parameters["t_cito_media"])*(parameters["odc_media"] - oligo_destr)

            olide_atual[i][j] = oligo_destr + h_t*(parameters["r_m"]*f_func(microglia, mic_media)*(parameters["odc_media"] - oligo_destr) + fag_mic_ant + apoptose_tke)

            #Anticorpo
            difusao_anticorpo = parameters["d_anti"]*calculaDifusao(f_jposterior, f_janterior, f_iposterior, f_ianterior, anticorpo)
            reacao_anticorpo = fag_mic_ant #Mesmo termo que soma na equacao das ODCs
            migracao_anticorpo = theta_BV[i][j]*parameters["gamma_F"]*(FL_atual - anticorpo)

            anticorpo_atual[i][j] = anticorpo + h_t*(difusao_anticorpo - reacao_anticorpo + migracao_anticorpo)

            #DC convencional
            quimiotaxia_dc = parameters["chi"]*calculaQuimiotaxia(dc_jposterior, dc_janterior, dc_iposterior, dc_ianterior, dc, parameters["dc_media"], gradiente_odc_i, gradiente_odc_j)
            difusao_dc = parameters["d_dc"]*calculaDifusao(dc_jposterior, dc_janterior, dc_iposterior, dc_ianterior, dc)
            reacao_dc = parameters["mu_dc"]*oligo_destr*(parameters["dc_media"] - dc)
            ativacao_dc_da = parameters["b_d"]*oligo_destr*dc

            dendritica_conv_atual[i][j] = dc + h_t*(reacao_dc + difusao_dc - quimiotaxia_dc - ativacao_dc_da)
            
            #DA ativada
            difusao_da = parameters["d_da"]*calculaDifusao(da_jposterior, da_janterior, da_iposterior, da_ianterior, da)
            migracao_da = theta_LV[i][j]*parameters["gamma_D"]*(DL_atual - da)

            dendritica_ativ_atual[i][j] = da + h_t*(difusao_da + ativacao_dc_da + migracao_da)

            Tc_Derivada += h_t*(difusao_t_cito - quimiotaxia_t_cito + migracao_t_cito)
            D_Derivada += h_t*(difusao_da + ativacao_dc_da + migracao_da)
            F_Derivada += h_t*(difusao_anticorpo - reacao_anticorpo + migracao_anticorpo)
            Odc_Derivada += h_t*(parameters["r_m"]*f_func(microglia, mic_media)*(parameters["odc_media"] - oligo_destr) + fag_mic_ant + apoptose_tke)
            mic_Derivada += h_t*(difusao_mic + reacao_mic - quimiotaxia_mic)

            D_sum_atual += dendritica_ativ_atual[i][j]
            T_c_sum_atual += t_cito_atual[i][j]
            F_sum_atual += anticorpo_atual[i][j]
            odc_sum_atual += olide_atual[i][j]
            mic_sum_atual += mic_atual[i][j]

            if microglia < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel microglia: " + str(microglia))
                exit(1)
            if da < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel DA: " + str(da))
                exit(1)
            if dc < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel dc: " + str(dc))
                exit(1)
            if t_cito < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel t_cito: " + str(t_cito))
                exit(1)
            if anticorpo < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel anticorpo: " + str(anticorpo))
                exit(1)
            if oligo_destr < 0:
                print("Tempo do Erro: " + str(k*h_t) + " - Variavel oligo_destr: " + str(oligo_destr))
                exit(1)
    if k == steps - 1:
        print("Dendriticas: " + str(results[1][0] + D_sum_atual - results[0][0] - D_sum_ant - DL_atualDerivada - D_Derivada))
        print("T citotoxicas: " + str(results[1][1] + T_c_sum_atual - results[0][1] - T_c_sum_ant - TL_c_atualDerivada - Tc_Derivada))
        print("Anticorpos: " + str(results[1][4] + F_sum_atual - results[0][4] - F_sum_ant - FL_atualDerivada - F_Derivada))
        print("Odcs: " + str(odc_sum_atual - odc_sum_ant - Odc_Derivada))
        print("Microglia: " + str(mic_sum_atual - mic_sum_ant - mic_Derivada))
    olide_anterior = np.copy(olide_atual)
    dendritica_conv_anterior = np.copy(dendritica_conv_atual)
    dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
    t_cito_anterior = np.copy(t_cito_atual)
    anticorpo_anterior = np.copy(anticorpo_atual)
    mic_anterior = np.copy(mic_atual)
    
    #Atualização da concentração das populações que migram.
    #Valores das populaçoes que migram que estão em contato com os vasos sanguineos ou linfaticos
    DendriticasTecido = 0
    AnticorposTecido = 0
    TcitotoxicaTecido = 0
    
    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if theta_LV[i][j] == 1:
                DendriticasTecido += dendritica_ativ_atual[i][j]
            if theta_BV[i][j] == 1:
                AnticorposTecido += anticorpo_atual[i][j]
                TcitotoxicaTecido += t_cito_atual[i][j]

    parameters["TcitotoxicaTecido"] = TcitotoxicaTecido/V_BV
    parameters["DendriticasTecido"] = DendriticasTecido/V_LV
    parameters["AnticorposTecido"] = AnticorposTecido/V_BV
    
    linfonodo_eqs = [DL_atual, TL_c_atual, TL_h_atual, B_atual, FL_atual, PL_atual]
    DL_vetor[k] = DL_atual
    TL_c_vetor[k] = TL_c_atual
    TL_h_vetor[k] = TL_h_atual
    B_vetor[k] = B_atual
    FL_vetor[k] = FL_atual
    PL_vetor[k] = PL_atual

    if k%intervalo_figs ==0 or k == steps-1:
        printMesh(k,olide_anterior, "odc")
        printMesh(k,mic_anterior, "microglia")
        printMesh(k,dendritica_conv_anterior, "dc")
        printMesh(k,dendritica_ativ_anterior, "da")
        printMesh(k,t_cito_anterior, "tke")
        printMesh(k,anticorpo_anterior, "anticorpo")
        print("Tempo: "+ str(k*h_t))
        print("IgG-T: " + str(AnticorposTecido))

#Fim da contagem do tempo
toc = time.perf_counter()

final_time = (toc - tic)/60

print("Tempo de execução: " + str(final_time) + " min")

#Transforma de dias para horas no plot
# t = np.multiply(t,24)

plt.plot(t,DL_vetor)
plt.title("Lymph node - Activated dendritic cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/dc_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_c_vetor)
plt.title("Lymph node - T $CD8^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/t_cito_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,TL_h_vetor)
plt.title("Lymph node - T $CD4^+$")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/t_helper_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,B_vetor)
plt.title("Lymph node - B Cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Cells/$mm^2$)")
plt.savefig('results/b_cell_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,FL_vetor)
plt.title("Lymph node - Antibodies")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig('results/anticorpo_linfonodo.png', dpi = 300)
plt.clf()

plt.plot(t,PL_vetor)
plt.title("Lymph node - Plasma-cells")
plt.xlabel("Time (days)")
plt.ylabel("Concentration (Molecules/$mm^2$)")
plt.savefig('results/pl_cell_linfonodo.png', dpi = 300)
plt.clf()
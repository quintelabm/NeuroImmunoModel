import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import saltelli
from SALib.analyze import sobol
from linfonodo import diferential

gradiente = lambda ponto_posterior, ponto_anterior, valor_maximo: quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo)
quimiotaxia = lambda ponto_atual, valor_maximo: ponto_atual/(valor_maximo + ponto_atual)
f_func = lambda populacao, valor_maximo: populacao*populacao/(valor_maximo + populacao)

class model: 

    def __init__(self, ht, hx, Lenght, Tf, params):
        self.ht = ht
        self.hx = hx
        self.Lenght = Lenght
        self.Tf = Tf
        self.parameters = params

        self.tvet = np.linspace(0, Tf, int(Tf/ht))
        self.xvet = np.linspace(0, Lenght, int(Lenght/hx))
        self.tamMesh = len(self.xvet)
        self.steps = len(self.tvet)
        self.num_figuras = Tf
        self.intervalo_figs = int(self.steps/self.num_figuras)
        self.theta_BV = np.zeros((int(Lenght/hx), int(Lenght/hx)))
        self.theta_LV = np.zeros((int(Lenght/hx), int(Lenght/hx)))
        self.V_BV = 0
        self.V_LV = 0
        
        
    def verifica_cfl(self, difusao_mic, difusao_dc, difusao_da, quimiotaxia_dc, quimiotaxia_mic):
        if(difusao_mic*self.ht/self.hx**2 < 1/4 and difusao_dc*self.ht/self.hx**2 < 1/4 and difusao_da*self.ht/self.hx**2 < 1/4 and quimiotaxia_dc*self.ht/self.hx < 1 and quimiotaxia_mic*self.ht/self.hx < 1):
            return True
        else:
            return False

    def define_BV_PVS(self):
        self.theta_LV = np.zeros((int(self.Lenght/self.hx), int(self.Lenght/self.hx)))
        for i in range(int(self.Lenght/self.hx)):
            for j in range(int(self.Lenght/self.hx)):
                if (i == self.Lenght/self.hx - 1 and j == self.Lenght/(self.hx*2)) or (i == 0 and j == self.Lenght/(self.hx*2)) or (i == self.Lenght/(self.hx*2) and j == 0) or (i == self.Lenght/(self.hx*2) and j == self.Lenght/self.hx - 1) or (i == int(self.Lenght/self.hx)/2 and j == int(self.Lenght/self.hx)/2):
                    self.LV[i][j] = 1
                    self.V_LV += 1

        theta_BV = np.zeros((int(self.Lenght/self.hx), int(self.Lenght/self.hx)))
        for i in range(int(self.Lenght/self.hx)):
            for j in range(int(self.Lenght/self.hx)):
                if (i == self.Lenght/self.hx - 1 and j == self.Lenght/self.hx - 1) or (i == 0 and j == self.Lenght/self.hx - 1) or (i == self.Lenght/self.hx - 1 and j == 0) or (i == 0 and j == 0):
                    self.theta_BV[i][j] = 1
                    self.V_BV += 1
        self.parameters["V_LV"] = self.V_LV
        self.parameters["V_BV"] = self.V_BV
    
    def checkBVeLV(self,xvet):
        x_pts, y_pts = np.meshgrid(xvet, xvet)
        max_population = 1
        levels = np.linspace(0, max_population, 3)

        cp = plt.contourf(x_pts, y_pts,self.theta_LV, levels=levels)
        plt.title("Perivascular space")
        plt.xlabel("Millimeters")
        plt.ylabel("Millimeters")
        plt.show()
        plt.clf()

        cp = plt.contourf(x_pts, y_pts,self.theta_BV, levels=levels)
        plt.title("Blood vessels")
        plt.xlabel("Millimeters")
        plt.ylabel("Millimeters")
        plt.show()
        plt.clf()
    
modelo = model(.0002, .5, 20, 28, {"V_LV": 0, "V_BV": 0})

# checkBVeLV()

def calculaQuimiotaxia(ponto_posterior_j, ponto_anterior_j, ponto_posterior_i, ponto_anterior_i, ponto_atual, valor_medio, gradiente_odc_i, gradiente_odc_j, h_x):
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

def calculaDifusao(ponto_posterior_j, ponto_anterior_j, ponto_posterior_i, ponto_anterior_i, ponto_atual, h_x):
    return (ponto_anterior_i + ponto_anterior_j + ponto_posterior_i + ponto_posterior_j - 4*ponto_atual)/(h_x**2)

# IC
# microglia
def init_tissue(L, h_x):
    mic_media = 350
    mic_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if (i-int(L/h_x)/2)**2 + (j-int(L/h_x)/2)**2 < 20/(2.5**2):
                mic_anterior[i][j] = mic_media/3.0

    # T citotóxica
    t_cito_anterior = np.zeros((int(L/h_x), int(L/h_x)))
    # Ol destruidos
    odc_media = 400
    olide_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # anticorpo
    anticorpo_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # Dendríticas convencionais
    dc_media = 33 #Valeria testa
    dendritica_conv_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # Dendríticas ativadas
    dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    mic_atual = np.zeros((int(L/h_x), int(L/h_x)))
    t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
    olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
    anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

# Modelo linfonodo
def init_lymph(linfonodo_eqs, estable_B, estable_T_c, estable_T_h):
    linfonodo_eqs[0]= 0    # Dendritic cells
    linfonodo_eqs[1]= 0.2  # Cytotoxic T cells
    linfonodo_eqs[2]= 0.4  # Helper T cells
    linfonodo_eqs[3]= estable_B    # B cells
    linfonodo_eqs[4]= 0    # Antibodies

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
    plt.savefig('../results/'+type+'/fig'+'{:.4f}'.format(time*h_t)+'.png', dpi = 300)
    plt.clf()

d_mic = (60*24*6.6/(2.5**2))*10**-5

def solve(chi = 0.298*60*2, d_mic = 1520*10**-5, mu_m = 60*24*3*10**-6, r_m = 60*24*3.96*10**-6, d_dc = 1520*10**-5, d_da = 1520*10**-5, d_t_cit = 1520*10**-5, d_anti = 1520*10**-4, lamb_f_m = 60*24*3.96*10**-6, b_d = 0.001, r_dc = 0.001, r_t = 0.1, mu_dc = 60*24*3*10**-4, gamma_D = 0.01, gamma_F = 0.03, gamma_T = 0.2, alpha_T_h = 0.01 , alpha_T_c = 0.5, alpha_B = 1, b_T = 0.017, b_Tc = 0.005, b_rho = 10**5, b_rho_b = 6.02*10**3, rho_T = 2, rho_Tc = 2, rho_B = 16, rho_F = 5.1*10**2, estable_T_h = 8.4*10**-3, estable_B = 8.4*10**-4, estable_T_c = 8.4*10**-3):
    linfonodo_eqs = np.zeros(5)
    
    V_BV = 0
    V_LV = 0

    theta_LV = np.zeros((int(L/h_x), int(L/h_x)))
    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if (i == L/h_x - 1 and j == L/(h_x*2)) or (i == 0 and j == L/(h_x*2)) or (i == L/(h_x*2) and j == 0) or (i == L/(h_x*2) and j == L/h_x - 1) or (i == int(L/h_x)/2 and j == int(L/h_x)/2):
                theta_LV[i][j] = 1
                V_LV += 1

    theta_BV = np.zeros((int(L/h_x), int(L/h_x)))
    for i in range(int(L/h_x)):
        for j in range(int(L/h_x)):
            if (i == L/h_x - 1 and j == L/h_x - 1) or (i == 0 and j == L/h_x - 1) or (i == L/h_x - 1 and j == 0) or (i == 0 and j == 0):
                theta_BV[i][j] = 1
                V_BV += 1

    V_LN = 160

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
    
    dendritica_conv_anterior = np.zeros((int(L/h_x), int(L/h_x)))

    # Dendríticas ativadas
    dendritica_ativ_anterior = np.zeros((int(L/h_x), int(L/h_x)))
    mic_atual = np.zeros((int(L/h_x), int(L/h_x)))
    t_cito_atual = np.zeros((int(L/h_x), int(L/h_x)))
    olide_atual = np.zeros((int(L/h_x), int(L/h_x)))
    anticorpo_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_conv_atual = np.zeros((int(L/h_x), int(L/h_x)))
    dendritica_ativ_atual = np.zeros((int(L/h_x), int(L/h_x)))

    init_lymph(linfonodo_eqs, estable_B, estable_T_c, estable_T_h)

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
    qoi = np.zeros(len(t))
    parameters = {
        "chi": chi, # Quimioatracao. valor por Dia
        "D_mic": d_mic, # Difusao da microglia. valor por Dia
        "mu_m": mu_m, # Taxa de ativação da microglia. valor por Dia
        "r_m": r_m, # intensidade dos danos causados pela microglia valor por Dia

        "d_dc": d_dc, # difusao DC convencional(procurar na literatura)
        "d_da": d_da, # difusao DC ativada(procurar na literatura)
        "d_t_cit": d_t_cit, # difusao t citotóxica(procurar na literatura)
        "d_anti": d_anti, # difusao anticorpo(procurar na literatura)
        "lamb_f_m": lamb_f_m, # taxa de anticorpos consumidos durante o processo de opsonização pela micróglia
        "b_d": b_d, # taxa de ativacao de dc por odc destruidos(procurar na literatura)
        "r_dc": r_dc, # taxa de coleta de odc destruidos pelas DCs (procurar na literatura)
        "r_t": r_t, # agressividade de t citotoxica(procurar na literatura)

        "mu_dc": mu_dc, #Taxa de producao de células dendríticas (procurar na literatura)
        "gamma_D": gamma_D, #Taxa de migração de DC ativadas para o linfonodo (procurar na literatura)
        "gamma_F": gamma_F, #Taxa de migração de anticorpos para o tecido (procurar na literatura)
        "gamma_T": gamma_T, #Taxa de migração de T citotoxica para o tecido (procurar na literatura)

        "t_cito_media": 37,
        "dc_media": dc_media,
        "mic_media": mic_media,
        "odc_media": 400,


        "alpha_T_h": alpha_T_h,
        "alpha_T_c": alpha_T_c,
        "alpha_B": alpha_B,
        "b_T": b_T,
        "b_Tc": b_Tc,
        "b_rho": b_rho,
        "b_rho_b": b_rho_b,
        "rho_T": rho_T,
        "rho_Tc": rho_Tc,
        "rho_B": rho_B,
        "rho_F": rho_F,
        "estable_T_h": estable_T_h,
        "estable_B": estable_B,
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

    DL_vetor[0] = linfonodo_eqs[0]
    TL_c_vetor[0] = linfonodo_eqs[1]
    TL_h_vetor[0] = linfonodo_eqs[2]
    B_vetor[0] = linfonodo_eqs[3]
    FL_vetor[0] = linfonodo_eqs[4]

    for k in range(1,steps):
        dy = diferential(linfonodo_eqs, parameters)
        DL_atual = linfonodo_eqs[0] + h_t*dy[0]
        TL_c_atual = linfonodo_eqs[1] + h_t*dy[1]
        TL_h_atual = linfonodo_eqs[2] + h_t*dy[2]
        B_atual = linfonodo_eqs[3] + h_t*dy[3]
        FL_atual = linfonodo_eqs[4] + h_t*dy[4]
        
        for i in range(tam):
            for j in range(tam):
                oligo_destr = olide_anterior[i][j]
                microglia = mic_anterior[i][j]
                dc = dendritica_conv_anterior[i][j]
                da = dendritica_ativ_anterior[i][j]
                anticorpo = anticorpo_anterior[i][j]
                t_cito = t_cito_anterior[i][j]
                
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
                if microglia < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel microglia: " + str(microglia))
                if da < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel DA: " + str(da))
                if dc < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel dc: " + str(dc))
                if t_cito < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel t_cito: " + str(t_cito))
                if anticorpo < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel anticorpo: " + str(anticorpo))
                if oligo_destr < 0:
                    print("Tempo do Erro: " + str(k*h_t) + " - Variavel oligo_destr: " + str(oligo_destr))
                
        olide_anterior = np.copy(olide_atual)
        dendritica_conv_anterior = np.copy(dendritica_conv_atual)
        dendritica_ativ_anterior = np.copy(dendritica_ativ_atual)
        t_cito_anterior = np.copy(t_cito_atual)
        anticorpo_anterior = np.copy(anticorpo_atual)
        mic_anterior = np.copy(mic_atual)

        #calcula QoI
        aux_qoi = 0
        for i in range(tam):
            for j in range(tam):
                aux_qoi = aux_qoi + olide_atual[i][j]
        qoi[k] = aux_qoi
        #Atualização da concentração das populações que migram.
        #Valores das populaçoes que migram que estão em contato com os vasos sanguineos ou linfaticos
        DendriticasTecido = 0
        AnticorposTecido = 0
        TcitotoxicaTecido = 0

        for i in range(int(L/h_x)):
            for j in range(int(L/h_x)):
                if theta_LV[i][j] == 1:
                    if k%intervalo_figs ==0:
                        print("DA-ponto: " + str(dendritica_ativ_atual[i][j]))
                    DendriticasTecido += dendritica_ativ_atual[i][j]
                if theta_BV[i][j] == 1:
                    if k%intervalo_figs ==0:
                        print("AT-ponto: " + str(anticorpo_atual[i][j]))
                        print("TCD8-ponto: " + str(t_cito_atual[i][j]))
                    AnticorposTecido += anticorpo_atual[i][j]
                    TcitotoxicaTecido += t_cito_atual[i][j]

        parameters["TcitotoxicaTecido"] = TcitotoxicaTecido/V_BV
        parameters["DendriticasTecido"] = DendriticasTecido/V_LV
        parameters["AnticorposTecido"] = AnticorposTecido/V_BV
        
        linfonodo_eqs = [DL_atual, TL_c_atual, TL_h_atual, B_atual, FL_atual]
        DL_vetor[k] = DL_atual
        TL_c_vetor[k] = TL_c_atual
        TL_h_vetor[k] = TL_h_atual
        B_vetor[k] = B_atual
        FL_vetor[k] = FL_atual
    print("Terminei de rodar uma vez!")
    outputFile = open("returns.txt", "a")
    outputFile.write(str(qoi[-1]) + "\n")
    outputFile.close()
    return qoi[-1]
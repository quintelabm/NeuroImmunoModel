import numpy as np

def diferential(y, parameters):
    
    dy = np.zeros(5)
    DC = y[0]
    TC = y[1]
    TH = y[2]
    B = y[3]
    Igg = y[4]

    # Dendritic cells

    dy[0] = parameters["gamma_D"] * (parameters["DendriticasTecido"] - DC) * (parameters["V_LV"] / parameters["V_LN"])

    # Cytotoxic T cells

    ativacaoC = parameters["b_Tc"]*(parameters["rho_Tc"] * TC * DC - TC*DC)
    homeostaseC = parameters["alpha_T_c"] * (parameters["estable_T_c"] - TC)
    migracaoC = (parameters["gamma_T"] * (TC - parameters["TcitotoxicaTecido"])) * (parameters["V_BV"] / parameters["V_LN"])

    dy[1] = ativacaoC + homeostaseC - migracaoC

    # # Helper T cells
    
    ativacaoH = parameters["b_T"]*(parameters["rho_T"] * TH * DC - TH*DC)
    homeostaseH = parameters["alpha_T_h"] * (parameters["estable_T_h"] - TH)
    H_decaimento_ativacaoB = 0#(parameters["b_rho"] * TH * DC * B)
    dy[2] = ativacaoH - H_decaimento_ativacaoB + homeostaseH

    # # B cells
    ativacaoB = 0#(parameters["b_rho_b"] * ((parameters["rho_B"] * TH * DC) - (TH * DC * B)))
    homeostaseB = parameters["alpha_B"] * (parameters["estable_B"] - B)
    dy[3] = ativacaoB + homeostaseB
    
    # # Antibodies
    producaoF = parameters["rho_F"] * B
    migracaoF = ((parameters["gamma_F"] * (Igg - parameters["AnticorposTecido"])) * (parameters["V_BV"] / parameters["V_LN"]))
    dy[4] = producaoF - migracaoF

    return dy
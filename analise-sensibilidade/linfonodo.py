import numpy as np

def diferential(y, parameters):
    
    dy = np.zeros(5)
    
    # Dendritic cells
    dy[0] = parameters["gamma_D"] * (parameters["DendriticasTecido"] - y[0]) * (parameters["V_LV"] / parameters["V_LN"])
    # Cytotoxic T cells
    dy[1] = parameters["b_Tc"]*(parameters["rho_Tc"] * y[1] * y[0] - y[1]*y[0]) + parameters["alpha_T_c"] * (parameters["estable_T_c"] - y[1]) - (parameters["gamma_T"] * (y[1] - parameters["TcitotoxicaTecido"])) * (parameters["V_BV"] / parameters["V_LN"])
    # # Helper T cells
    dy[2] = parameters["b_T"]*(parameters["rho_T"] * y[2] * y[0] - y[2]*y[0]) - (parameters["b_rho"] * y[2] * y[0] * y[3]) + parameters["alpha_T_h"] * (parameters["estable_T_h"] - y[2])
    # # B cells
    dy[3] = (parameters["b_rho_b"] * ((parameters["rho_B"] * y[2] * y[0]) - (y[2] * y[0] * y[3]))) + parameters["alpha_B"] * (parameters["estable_B"] - y[3])
    # # Antibodies
    dy[4] = parameters["rho_F"] * y[3] - ((parameters["gamma_F"] * (y[4] - parameters["AnticorposTecido"])) * (parameters["V_BV"] / parameters["V_LN"]))

    return dy
import numpy as np
import Leachman

class ThermoEngine:
    def __init__(self, model_type="REAL", cp_ig=14.30, R_const=4.124, T0=298.15, P0=1.0):
        self.model = model_type.upper()
        self.cp = cp_ig
        self.R = R_const
        self.k = self.cp / (self.cp - self.R)
        self.T0 = T0
        self.P0 = P0

    def flash_TP(self, T_C, P_bar):
        if self.model == "IDEAL":
            T_K = T_C + 273.15
            h = self.cp * T_K
            s = self.cp * np.log(T_K / self.T0) - self.R * np.log(P_bar / self.P0)
            return {"T_C": T_C, "P_bar": P_bar, "H": h, "S": s}
        else:
            st = Leachman.calcular_estado("LEACHMAN", "H2", "TP", T_C, P_bar)
            return {"T_C": Leachman.K_para_C(st["T_K"]), "P_bar": P_bar, "H": st["H"], "S": st["S"]}

    def flash_PS(self, P_bar, S_val):
        if self.model == "IDEAL":
            T_K = self.T0 * np.exp((S_val + self.R * np.log(P_bar / self.P0)) / self.cp)
            T_C = T_K - 273.15
            h = self.cp * T_K
            return {"T_C": T_C, "P_bar": P_bar, "H": h, "S": S_val}
        else:
            st = Leachman.calcular_estado("LEACHMAN", "H2", "PS", P_bar, S_val)
            return {"T_C": Leachman.K_para_C(st["T_K"]), "P_bar": P_bar, "H": st["H"], "S": st["S"]}

    def flash_PH(self, P_bar, H_val):
        if self.model == "IDEAL":
            T_K = H_val / self.cp
            T_C = T_K - 273.15
            s = self.cp * np.log(T_K / self.T0) - self.R * np.log(P_bar / self.P0)
            return {"T_C": T_C, "P_bar": P_bar, "H": H_val, "S": s}
        else:
            st = Leachman.calcular_estado("LEACHMAN", "H2", "PH", P_bar, H_val)
            return {"T_C": Leachman.K_para_C(st["T_K"]), "P_bar": P_bar, "H": st["H"], "S": st["S"]}

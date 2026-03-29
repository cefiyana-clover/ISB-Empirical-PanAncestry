import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
import warnings
warnings.filterwarnings('ignore')

V_max_base, K_ATP, K_m, C_basal, k_leak = 5.0, 0.5, 1.0, 0.2, 0.05

def bioenergetic_system(state, t, PaCO2, R_stress, P_basal_genetic):
    A, G = state
    A, G = max(1e-5, A), max(1e-5, G)
    phi_vaso = np.exp(0.035 * (PaCO2 - 40))
    phi_bohr = 1 / (1 + 10**(0.1 * (40 - PaCO2)))
    P_ATP_current = P_basal_genetic * phi_vaso * phi_bohr
    V_GLT1 = V_max_base * (A / (K_ATP + A)) * (G / (K_m + G))
    dA_dt = P_ATP_current - C_basal - 4 * V_GLT1
    dG_dt = R_stress - V_GLT1 - k_leak * G
    return [dA_dt, dG_dt]

np.random.seed(42)
N = 10000

paco2_env = np.random.normal(39.2, 0.6, N) 

p_basal_eur, r_stress_eur = np.random.normal(8.0, 0.2, N), np.random.normal(0.6, 0.1, N)
p_basal_eas, r_stress_eas = np.random.normal(8.0, 0.2, N), np.random.normal(1.0, 0.15, N) 
p_basal_afr, r_stress_afr = np.random.normal(8.4, 0.2, N), np.random.normal(0.6, 0.1, N)

pop_paco2 = np.concatenate([paco2_env, paco2_env, paco2_env])
pop_pbasal = np.concatenate([p_basal_eur, p_basal_eas, p_basal_afr])
pop_rstress = np.concatenate([r_stress_eur, r_stress_eas, r_stress_afr])
pop_labels = np.array(['European (EUR)']*N + ['East Asian (EAS)']*N + ['African (AFR)']*N)

pop_paco2, pop_rstress = np.clip(pop_paco2, 20.0, 45.0), np.clip(pop_rstress, 0.1, 5.0)

final_atp, final_glu = np.zeros(30000), np.zeros(30000)
time_span, state0 = np.linspace(0, 100, 100), [3.0, 0.2]

for i in range(30000):
    sol = odeint(bioenergetic_system, state0, time_span, args=(pop_paco2[i], pop_rstress[i], pop_pbasal[i]))
    final_atp[i], final_glu[i] = sol[-1, 0], sol[-1, 1]

df_mega = pd.DataFrame({'Ancestry': pop_labels, 'ISB_Score': np.clip((3.0 - final_atp) / 3.0 * 5.0, 0.0, 5.0)})

plt.figure(figsize=(10, 6), dpi=300)
sns.violinplot(data=df_mega, x='Ancestry', y='ISB_Score', palette={'European (EUR)':'#1f77b4', 'East Asian (EAS)':'#d62728', 'African (AFR)':'#2ca02c'}, inner="quartile", alpha=0.8)
plt.axhline(2.5, color='black', linestyle='--', linewidth=2, label='Collapse Threshold (>2.5)')
plt.title('Genetic Vulnerability to Mitochondrial Shutdown (Borderline Stress)', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('Figure_Mega_PanAncestry_Borderline.png', dpi=300)

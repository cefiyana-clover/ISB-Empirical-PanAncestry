import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
import warnings
warnings.filterwarnings('ignore')

V_max_base, K_ATP, K_m, C_basal, k_leak, P_basal = 5.0, 0.5, 1.0, 0.2, 0.05, 8.0

def bioenergetic_system(state, t, PaCO2, R_stress):
    A, G = state
    A, G = max(1e-5, A), max(1e-5, G)
    phi_vaso = np.exp(0.035 * (PaCO2 - 40))
    phi_bohr = 1 / (1 + 10**(0.1 * (40 - PaCO2)))
    V_GLT1 = V_max_base * (A / (K_ATP + A)) * (G / (K_m + G))
    dA_dt = P_basal * phi_vaso * phi_bohr - C_basal - 4 * V_GLT1
    dG_dt = R_stress - V_GLT1 - k_leak * G
    return [dA_dt, dG_dt]

np.random.seed(42)
N_PER_GROUP = 500

paco2_pop = np.concatenate([np.random.normal(40.0, 2.0, N_PER_GROUP), np.random.normal(32.0, 3.0, N_PER_GROUP)])
r_stress_pop = np.concatenate([np.random.normal(0.8, 0.2, N_PER_GROUP), np.random.normal(2.8, 0.5, N_PER_GROUP)])
cohort_labels = np.concatenate([['Healthy Control']*N_PER_GROUP, ['Clinical (GERD+MDD)']*N_PER_GROUP])

paco2_pop, r_stress_pop = np.clip(paco2_pop, 20.0, 60.0), np.clip(r_stress_pop, 0.5, 4.0)

final_atp, final_glu = np.zeros(1000), np.zeros(1000)
time_span, state0 = np.linspace(0, 100, 500), [3.0, 0.2]

for i in range(1000):
    sol = odeint(bioenergetic_system, state0, time_span, args=(paco2_pop[i], r_stress_pop[i]))
    final_atp[i], final_glu[i] = sol[-1, 0], sol[-1, 1]

df = pd.DataFrame({'Cohort_Group': cohort_labels, 'PaCO2_Input': paco2_pop, 'Final_Glutamate_Load': final_glu, 'ISB_Score': np.clip((3.0 - final_atp) / 3.0 * 5.0, 0.0, 5.0)})

fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
sns.scatterplot(data=df, x='PaCO2_Input', y='Final_Glutamate_Load', hue='Cohort_Group', palette={'Healthy Control': '#2ca02c', 'Clinical (GERD+MDD)': '#d62728'}, alpha=0.7, ax=axes[0])
sns.kdeplot(data=df, x='ISB_Score', hue='Cohort_Group', fill=True, palette={'Healthy Control': '#2ca02c', 'Clinical (GERD+MDD)': '#d62728'}, alpha=0.5, ax=axes[1])
plt.tight_layout()
plt.savefig('Figure_ISB_Empirical_Bifurcation.png', dpi=300)

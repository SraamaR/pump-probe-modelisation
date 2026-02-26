import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Dossier des résultats [doit être créé au préalable !]
results_folder = "results/" #Dossier de sauvegarde des résultats, doit se terminer par "/"
data = []

pump_power = [1.2, 2.4, 4.4, 5.0, 6.1, 13.6, 39.6] #Puissance de la pompe en mW

colors = plt.cm.Dark2(np.linspace(0,1,len(pump_power)))

delta = 1.5
decalage = 2

#Lire les fichiers CSV
csv_0_filename = f'{results_folder}data/temporal_scan_2.00e+07.csv'
data.append(pd.read_csv(csv_0_filename, sep=';', decimal=','))

csv_1_filename = f'{results_folder}data/temporal_scan_2.90e+07.csv'
data.append(pd.read_csv(csv_1_filename, sep=';', decimal=','))

csv_2_filename = f'{results_folder}data/temporal_scan_4.00e+07.csv'
data.append(pd.read_csv(csv_2_filename, sep=';', decimal=','))

csv_3_filename = f'{results_folder}data/temporal_scan_4.20e+07.csv'
data.append(pd.read_csv(csv_3_filename, sep=';', decimal=','))

csv_4_filename = f'{results_folder}data/temporal_scan_4.70e+07.csv'
data.append(pd.read_csv(csv_4_filename, sep=';', decimal=','))

csv_5_filename = f'{results_folder}data/temporal_scan_7.00e+07.csv'
data.append(pd.read_csv(csv_5_filename, sep=';', decimal=','))

csv_6_filename = f'{results_folder}data/temporal_scan_1.10e+08.csv'
data.append(pd.read_csv(csv_6_filename, sep=';', decimal=','))

plt.figure(figsize=(3.36, 3.36), constrained_layout=True)

for i in range(0,7):
    plt.plot(data[i]['delay_ps']+delta*i+decalage, data[i]['R'], 
             linestyle='-', linewidth=2, label=f"{pump_power[i]} mW", color=colors[i])

plt.legend(fontsize=9, loc='best')
plt.ylabel('Reflectivity', fontsize=9)
plt.xlabel('Time (ps)', fontsize=9)
plt.grid(True, alpha=0.3)
plt.xlim(-0.25, 15)
plt.ylim(top=1.0)

plt.savefig(f'{results_folder}/comparaison/RoverTimex2.png', dpi=600)
print("\nGraphique sauvegardé sous 'reflectivity_vs_thickness.png'")

plt.show()
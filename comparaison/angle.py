import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#Dossier des résultats [doit être créé au préalable !]
results_folder = "results/" #Dossier de sauvegarde des résultats, doit se terminer par "/"

# Lire le fichier CSV
csv_simulation_vingt_filename = f'{results_folder}data/vingt_3_donnees_simulation_L_parallel_1.00e+05.csv'
data_simulation_vingt = pd.read_csv(csv_simulation_vingt_filename, sep=';', decimal=',')

csv_simulation_normale_filename = f'{results_folder}data/normale_3_donnees_simulation_L_parallel_1.00e+05.csv'
data_simulation_normale = pd.read_csv(csv_simulation_normale_filename, sep=';', decimal=',')


# Afficher les données pour vérification
print("Données de simulation vingt chargées :")
print(data_simulation_vingt.head())
print(f"\nNombre de lignes :  {len(data_simulation_vingt)}")

print("Données de simulation normale chargées :")
print(data_simulation_normale.head())
print(f"\nNombre de lignes :  {len(data_simulation_normale)}")



plt.figure(figsize=(3.36, 3.36), constrained_layout=True)
plt.plot(data_simulation_normale['Thickness (nm)'], data_simulation_normale['R_max'], 
         linestyle='-', linewidth=1,
         label='Simulation angle=0°')
plt.plot(data_simulation_vingt['Thickness (nm)'], data_simulation_vingt['R_max'], 
         linestyle='-', linewidth=1,
         label='Simulation angle=20°')

plt.xlabel('GaAs thickness (nm)', fontsize=9)
plt.ylabel('Maximum reflectivity', fontsize=9)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.2)
plt.xlim(250, 1000)
plt.ylim(0, 1.1)

ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_major_formatter(
    ticker.FuncFormatter(
        lambda x, pos: f"{x:g}" if pos % 2 == 0 else ""
    )
)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.grid(which='major', alpha=0.2, color='#666666', linestyle='-')
plt.grid(which='minor', alpha=0.1, color='#999999', linestyle='--')

# Sauvegarder le graphique
plt.savefig(f'{results_folder}/comparaison/comparaison_angle.png', dpi=300)
print("\nGraphique sauvegardé sous 'comparaison_angle.png'")

plt.plot()
plt.show()

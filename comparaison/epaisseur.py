import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#Dossier des résultats [doit être créé au préalable !]
results_folder = "results/" #Dossier de sauvegarde des résultats, doit se terminer par "/"

# Lire le fichier CSV
csv_experience_filename = 'comparaison/data_experimental/BEVEL2_summary.csv'
data_experience = pd.read_csv(csv_experience_filename, sep=';', decimal=',')

csv_simulation_filename = f'{results_folder}data/donnees_simulation_L_parallel_3.20e+07.csv'
data_simulation = pd.read_csv(csv_simulation_filename, sep=';', decimal=',')

csv_simulation_lineaire_filename = f'{results_folder}data/vingt_3_donnees_simulation_L_parallel_1.00e+05.csv'
data_simulation_lineaire = pd.read_csv(csv_simulation_lineaire_filename, sep=';', decimal=',')


# Afficher les données pour vérification
print("Données expérimentales chargées :")
print(data_experience.head())
print(f"\nNombre de lignes :  {len(data_experience)}")

print("Données de simulation chargées :")
print(data_simulation.head())
print(f"\nNombre de lignes :  {len(data_simulation)}")

print("Données de simulation linéaire chargées :")
print(data_simulation_lineaire.head())
print(f"\nNombre de lignes :  {len(data_simulation_lineaire)}")




#Thickness - R_max
plt.figure(figsize=(3.36, 3.36), constrained_layout=True)
plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_max'], 
         linestyle='-', linewidth=2,
         label='Simulation data')
plt.plot(data_experience['thickness_nm'], data_experience['Rmax'], 
         marker='.', linewidth=0, markersize=8,
         label='Experimental data')

xerr = data_experience['thickness_nm']*(np.abs(data_experience['thickness_nm']**2*5e-7-0.0005*data_experience['thickness_nm']+0.1151))  # Erreur horizontale basée sur une fonction quadratique
for i in range(len(xerr)):
    if xerr[i] < 0:
        xerr[i] = 0
plt.errorbar(
    data_experience['thickness_nm'],
    data_experience['Rmax'],
    xerr= xerr,      # erreur horizontale uniquement
    fmt='.',
    markersize=0,
    elinewidth=2,
    capsize=2,
    ecolor='darkorange'
)

plt.xlabel('GaAs thickness (nm)', fontsize=9)
plt.ylabel('Maximum reflectivity', fontsize=9)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.2)
plt.xlim(250, 700)
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
plt.savefig(f'{results_folder}/comparaison/reflectivity_vs_thickness.png', dpi=300)
print("\nGraphique sauvegardé sous 'reflectivity_vs_thickness.png'")






#Thickness - dR
plt.figure(figsize=(10, 6))

plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_max'] - data_simulation['R_linear (pump off)'], 
         linestyle='-', linewidth=2,
         label='$\Delta$R max simulation', color='blue')
plt.plot(data_experience['thickness_nm'], data_experience['dR'], 
         marker='o', linestyle='-', linewidth=2, markersize=6,
         label='$\Delta$R_max expérience', color='red')


plt.xlabel('Épaisseur GaAs (nm)', fontsize=12)
plt.ylabel('Réflectivité', fontsize=12)
plt.title('Réflectivité en fonction de l\'épaisseur', fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)
plt.xlim(150, 700)
plt.ylim(0, 1.1)

ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.grid(which='major', alpha=0.5, color='#666666', linestyle='-')
plt.grid(which='minor', alpha=0.2, color='#999999', linestyle='--')

# Ajuster la mise en page
plt.tight_layout()

# Sauvegarder le graphique
plt.savefig(f'{results_folder}comparaison/differential_reflectivity_vs_thickness.png', dpi=300, bbox_inches='tight')
print("\nGraphique sauvegardé sous 'differential_reflectivity_vs_thickness.png'")



#Thickness - R0
plt.figure(figsize=(3.36, 3.36), constrained_layout=True)
plt.plot(data_simulation_lineaire['Thickness (nm)'], data_simulation_lineaire['R_max'], 
         linestyle='-', linewidth=2,
         label='Simulation data')
plt.plot(data_experience['thickness_nm'], data_experience['Rmin'], 
         marker='.', linewidth=0, markersize=8,
         label='Experimental data')

plt.errorbar(
    data_experience['thickness_nm'],
    data_experience['Rmin'],
    xerr= xerr,      # erreur horizontale uniquement
    fmt='.',
    markersize=0,
    elinewidth=2,
    capsize=2,
    ecolor='darkorange'
)

plt.xlabel('GaAs thickness (nm)', fontsize=9)
plt.ylabel('Maximum reflectivity', fontsize=9)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.2)
plt.xlim(250, 700)
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

plt.savefig(f'{results_folder}comparaison/R0_vs_thickness.png', dpi=300)
print("\nGraphique sauvegardé sous 'R0_vs_thickness.png'")

#Thickness - R0 - POINT
plt.figure(figsize=(10, 6))

plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_linear (pump off)'], 
         marker='o', linestyle='-', linewidth=0, markersize=6,
         label='R0 simulation', color='blue')
plt.plot(data_experience['thickness_nm'], data_experience['Rmin'], 
         marker='o', linestyle='-', linewidth=0, markersize=6,
         label='R0 expérience', color='red')

plt.xlabel('Épaisseur GaAs (nm)', fontsize=12)
plt.ylabel('Réflectivité linéaire', fontsize=12)
plt.title('Réflectivité linéaire en fonction de l\'épaisseur', fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)
plt.xlim(150, 700)
plt.ylim(0, 1.1)

ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.grid(which='major', alpha=0.5, color='#666666', linestyle='-')
plt.grid(which='minor', alpha=0.2, color='#999999', linestyle='--')

# Ajuster la mise en page
plt.tight_layout()

# Sauvegarder le graphique
plt.savefig(f'{results_folder}comparaison/R0_vs_thickness_POINT.png', dpi=300, bbox_inches='tight')
print("\nGraphique sauvegardé sous 'R0_vs_thickness_POINT.png'")


#Thickness - dR/Rmax
plt.figure(figsize=(10, 6))

plt.plot(data_simulation['Thickness (nm)'], data_simulation['max(dR/R)'], 
         linestyle='-', linewidth=2,
         label='dR/R max simulation', color='blue')
plt.plot(data_experience['thickness_nm'], data_experience['dR/Rmax'], 
         marker='o', linestyle='-', linewidth=2, markersize=6,
         label='dR/R max expérience', color='red')

plt.xlabel('Épaisseur GaAs (nm)', fontsize=12)
plt.ylabel('Réflectivité', fontsize=12)
plt.title('Réflectivité en fonction de l\'épaisseur', fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)
plt.xlim(150, 700)
plt.ylim(0, 1.1)

ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.grid(which='major', alpha=0.5, color='#666666', linestyle='-')
plt.grid(which='minor', alpha=0.2, color='#999999', linestyle='--')

# Ajuster la mise en page
plt.tight_layout()

# Sauvegarder le graphique
plt.savefig(f'{results_folder}comparaison/dR_sur_R_vs_thickness.png', dpi=300, bbox_inches='tight')
print("\nGraphique sauvegardé sous 'dR_sur_R_vs_thickness.png'")



#Thickness - R0
plt.figure(figsize=(10, 6))

plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_linear (pump off)'], 
         linestyle='-', linewidth=2,
         label='R0 simulation', color='blue')
plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_max'], 
         linestyle='-', linewidth=2,
         label='Rmax simulation', color='red')
plt.plot(data_simulation['Thickness (nm)'], data_simulation['R_max'] - data_simulation['R_linear (pump off)'], 
         linestyle='-', linewidth=1,
         label='$\Delta$ R_max', color='green', alpha=0.5)

plt.xlabel('Épaisseur GaAs (nm)', fontsize=12)
plt.ylabel('Réflectivité linéaire', fontsize=12)
plt.title('Réflectivité linéaire en fonction de l\'épaisseur', fontsize=14)
plt.legend(fontsize=10, loc='best')
plt.grid(True, alpha=0.3)
plt.xlim(150, 700)
plt.ylim(0, 1.1)

ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.grid(which='major', alpha=0.5, color='#666666', linestyle='-')
plt.grid(which='minor', alpha=0.2, color='#999999', linestyle='--')

# Ajuster la mise en page
plt.tight_layout()

# Sauvegarder le graphique
plt.savefig(f'{results_folder}comparaison/R0RMaxDeltaR.png', dpi=300, bbox_inches='tight')
print("\nGraphique sauvegardé sous 'R0RMaxDeltaR.png'")

plt.show()
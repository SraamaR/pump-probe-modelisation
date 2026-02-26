import matplotlib.pyplot as plt
import numpy as np

# --------------------
# PARAMÈTRES UTILISATEURS A MODIFIER
# --------------------

#Renseigner ci-dessous les valeurs de la simulation à effectuer
#L'unité nécessaire pour le bon fonctionnement du code est indiquée entre crochet

#Paramètres du matériau
n0 = 3.7 #[adim]
wavelength = 780e-9 # [m]
gamma = 1/(50e-15) # [s^-1]
tau_c = 1e-12 # [s]
mr = 5.36e-32 #[kg] masse réduite effective
Egap = 2.272e-19 #[J]
d_cv = (1.12e-28)**2 #[C².m²] moment dipolaire, indiquer |d_cv|²

#Paramètre de simulation
L_objectif = 450e-9 #[m] épaisseur de semiconducteur à simuler
n_gold = 0.15 - 1j*4.74 #[adim] indice de réfraction du miroir d'or
couches_lineaires = [] #Liste des tuples (indice, epaisseur) des couches linéaires précédent le semiconducteur, laisser vide si seulement air

#Paramètres de discrétisation
NEt = 100 # nombre d’échantillons énergie 
Nt = 10000 # nombre de pas en temps
Nz = 11 # nombre de pas en espace
Etmax = 8.01e-21 #[J] borne supérieure de l'intégrale de polarisation

#Paramètres du pulse d'entrée
FWHM = 150e-15 #[s] durée à mi-hauteur du pulse d'entrée
tau0 = 1e-12 #[s] temps de pic du pulse d'entrée
amplitude = 1e8 #[V/m] amplitude maximale du champ électrique du pulse d'entrée

# --------------------
# CONSTANTES PHYSIQUES FONDAMENTALES
# --------------------

c = 3e8 #[m.s^-1]
h_bar = 1.05e-34 #[J.s]
m0 = 9.1e-31 #[kg] masse d'un électron
q = 1.6e-19 #[C] charge élémentaire
epsilon0 = 8.85e-12 
mu0 = 4*np.pi*1e-7

# -------------------
# TABLEAUX DE SORTIE 
# -------------------

#On liste ci-dessous les différents tableaux des variables à calculer.
#Ils peuvent être utilisés pour tracer les différentes courbes à la fin de la simulation.

#Tableaux des enveloppes du champ électrique
E_plus = np.zeros((Nz, Nt), dtype=complex) #Enveloppe du champ croissante
E_moins = np.zeros((Nz, Nt), dtype=complex) #Enveloppe du champ décroissante
E = np.zeros((Nz,Nt), dtype=complex) #Enveloppe croissante + décroissante du champ

#Tableaux des coefficients de Fresnel à l'entrée
t_in_array = np.zeros(Nt, dtype=complex)   # in : Air → GaAs
t_out_array = np.zeros(Nt, dtype=complex)  # out : GaAs → Air
r_in_array = np.zeros(Nt,dtype=complex)
r_out_array = np.zeros(Nt,dtype=complex)

#Tableaux de l'indice optique effectif et du coefficient d'absorption
n_eff_array = np.zeros(Nt, dtype=complex)
alpha_array = np.zeros(Nt) 

#Tableaux pour stocker l'évolution des populations
rho_e = np.zeros((NEt, Nz, Nt))
rho_h = np.zeros((NEt, Nz, Nt))

#Tableau pour stocker la valeur du produit de convolution
F_array = np.zeros((NEt, Nz, Nt), dtype=complex)
F_array_plus = np.zeros((NEt, Nz, Nt), dtype=complex)
F_array_moins = np.zeros((NEt, Nz, Nt), dtype=complex)

#Tableau pour stocker la valeur de la polarisation
polarisation= np.zeros((Nz,Nt),dtype=complex)
polarisation_plus= np.zeros((Nz,Nt),dtype=complex)
polarisation_moins= np.zeros((Nz,Nt),dtype=complex)

N_entree = np.zeros(Nt)

partie_reelle = np.zeros((Nt))
partie_imaginaire = np.zeros((Nt))

# --------------------
# VALEURS DERIVÉES
# --------------------
v_g = c/n0 #[m.s^-1] vitesse de groupe
omega0 = (2*np.pi*c/wavelength) #[rad.s^-1] pulsation de référence
k0 = n0*omega0/c

g = (d_cv*omega0)/(4*n0*c*epsilon0*h_bar) #[adim]
kappa = d_cv/(2*h_bar**2)  # force du couplage champ-populations utilisé dans l'évolution des populations
D0 = 1/(2*np.pi**2)*(2*mr/(h_bar)**2)**(3/2)

dt = L_objectif/(v_g*(Nz-1))  # pas en temps
dz = v_g*dt # pas en espace 

Et0 = h_bar*omega0-Egap

t = np.arange(Nt) * dt
z = np.arange(Nz) * dz 
Et = np.linspace(-Et0,Etmax*2,NEt)

a = 0.64/q #Paramètre de non-parabolicité
Dr = D0*np.sqrt((Et+Et0)*(1+a*(Et+Et0)))*(1+2*a*(Et+Et0)) #Densité d'états en 3D (il s'agit d'un tableau en Et)

#On précalcule le coefficient alpha utilisée dans le calcul du produit de convolution pour gagner du temps d'exécution.
alpha = np.exp(dt*(-1j*(Et)/h_bar-gamma))

sigma = FWHM / (2*np.sqrt(2*np.log(2)))

# --------------------
# DEFINITION DU PULSE D'ENTRÉE
# --------------------
base_pulse = np.exp(-((t-tau0)**2)/(2*sigma**2))
source_pulse = np.zeros((Nt), dtype=complex)
source_pulse = amplitude * base_pulse

# --------------------
# FONCTIONS 
# --------------------
def get_multistack_coeffs(layers_list):
    """
    Calcule les coefficients de réflexion et de transmission d'une structure stratifiée à partir de la matrice de transfert.
    La structure est définie par une liste de tuples (n, d) où n est l'indice de réfraction de la couche et d son épaisseur. La dernière interface est celle entre la dernière couche et un milieu effectif d'indice n_eff.
    Les coefficients sont calculés pour une incidence normale.

    Parameters
    ----------
    layers_list : list of tuples
        Liste des couches linéaires précédant le semiconducteur, sous la forme [(n1, d1), (n2, d2), ...].
    """
    M = np.eye(2, dtype=complex)
    n_curr = 1.0 #Indice de départ : Air
    
    for (n_next, d_next) in layers_list:
        #Première interface
        r = (n_curr-n_next)/(n_curr+n_next)
        t = 2*n_curr/(n_curr+n_next)
        #Mise à jour de la matrice M
        M_int = 1/t*np.array([[1,r],[r,1]], dtype=complex)
        M = np.dot(M,M_int)
        #Propagation dans la couche
        phi = k0*n_next*d_next
        M_prop = np.array([[np.exp(-1j*phi),0],[0,np.exp(1j*phi)]], dtype=complex)
        M = np.dot(M,M_prop)
        
        n_curr = n_next
        
    #Pour la dernière couche :
    r = (n_curr-n0)/(n_curr + n0)
    t = 2*n_curr / (n_curr+n0)
    M_int = 1/t*np.array([[1,r],[r,1]], dtype=complex)
    M = np.dot(M,M_int)
    
    #Extraction des coefficients à renvoyer
    r_in = M[1,0]/M[0,0]
    t_in = 1 / M[0,0]
    r_out = - M[0,1]/ M[0,0]
    t_out = (M[0,0]*M[1,1]-M[0,1]*M[1,0])/M[0,0]
    
    return r_in, t_in, r_out, t_out

def F_moins(m,n):
    F_array_moins[:, m,n] = alpha*F_array_moins[:,m,n-1] + 1/(gamma+(Et)/h_bar*1j)*(1-alpha)*E_moins[m,n-1]
    return

def F_plus(m,n):
    F_array_plus[:, m,n] = alpha*F_array_plus[:,m,n-1] + 1/(gamma+(Et)/h_bar*1j)*(1-alpha)*E_plus[m,n-1]
    return

def F(m, n):
    """
    Calcule la convolution récursive spectrale F(n) = α·F(n-1) + (1-α)/(γ + iE_t+E_t0/ℏ)·E(n-1).
    
    Parameters
    ----------
    m : int
        Index spatial.
    n : int
        Index temporel.
    """    
    F_array[:, m,n] = F_array_plus[:,m,n]*np.exp(-1j*k0*dz*m) + F_array_moins[:,m,n]*np.exp(1j*k0*dz*m)
    return
 
def f(m, n):
    """
    Intègre la polarisation P = ∫ (ρₑ + ρₕ - 1)·F(E,t)·D(E) dE.
    ATTENTION : ce n'est en réalité pas une polarisation car il manque le facteur multiplicatif en i|d_cv|**2/(2h_barGamma)
    Celui-ci est déjà comptabilisé dans le gain g0 mais pour calculer n_eff, il ne faut pas oublier de remettre ce facteur.
    
    Parameters
    ----------
    m : int
        Index spatial.
    n : int
        Index temporel.
    
    Returns
    -------
    complex
        Polarisation au point (m, n).
    """
    polarisation[m,n] = np.trapezoid((rho_e[:, m, n] + rho_h[:, m, n]-1) * F_array[:, m,n]*Dr,Et)
    return

def f_plus(m,n):
    polarisation_plus[m,n] = np.trapezoid((rho_e[:, m, n] + rho_h[:, m, n]-1) * F_array_plus[:, m,n]*Dr,Et)
    return polarisation_plus[m,n]

def f_moins(m,n):
    polarisation_moins[m,n] = np.trapezoid((rho_e[:, m, n] + rho_h[:, m, n]-1) * F_array_moins[:, m,n]*Dr,Et)
    return polarisation_moins[m,n]

def rho_calcul(n, m, E_value):
    """
    Résout dρ/dt = κ(1-ρₑ-ρₕ)|E·F| - ρ/τ_c avec clipping ρ ∈ [0, 0.5].
    
    Parameters
    ----------
    n : int
        Index temporel.
    m : int
        Index spatial.
    E_value : complex
        Champ électrique combiné.
    """
    
    interaction = np.real(np.conj(E_value) * F_array[:, m,n])
    drho =  kappa * (1-rho_e[:, m, n-1] - rho_h[:, m, n-1]) * interaction
    rho_e[:, m, n] = rho_e[:, m, n-1] + dt * (drho - 1/tau_c * rho_e[:, m, n-1])
    rho_h[:, m, n] = rho_h[:, m, n-1] + dt * (drho - 1/tau_c * rho_h[:, m, n-1])
    rho_e[:, m, n] = np.clip(rho_e[:, m, n], 0, 0.5)
    rho_h[:, m, n] = np.clip(rho_h[:, m, n], 0, 0.5)
    return


def propagation_croissante(n, m):
    """
    Propage l'onde progressive : E⁺(m+1, n) = E⁺(m, n-1) + dz·g·P+(m, n).
    
    Parameters
    ----------
    n : int
        Index temporel.
    m : int
        Index spatial.
    """
    
    pol_dynamique = f_plus(m, n)   
    E_plus[m+1, n] = E_plus[m, n-1] + dz * g * (pol_dynamique)  
    return 

def propagation_decroissante(n, m):
    """Propage l'onde progressive : E-(m-1, n) = E-(m, n-1) + dz·g·P-(m, n).
    
    Parameters
    ----------
    n : int
        Index temporel.
    m : int
        Index spatial.
    """
    
    pol_dynamique = f_moins(m, n)  
    E_moins[m-1, n] = E_moins[m, n-1] + dz * g * (pol_dynamique)  
    return

def calculate_total_neff(E_field, n, m):
    """
    Calcule n_eff en utilisant une stratégie Kramers-Kronig soustractive (ou scindée).
    - Partie Imaginaire : Calculée sur la population totale (saturation incluse).
    - Partie Réelle : Calculée UNIQUEMENT sur les porteurs excités (delta n), 
      le fond étant géré par n0 pour éviter l'erreur de cutoff.
    """
    
    # 1. Définition des variables spectrales
    denom = gamma**2 + (Et / h_bar)**2 
    
    # 2. Facteur de normalisation global
    prefactor = d_cv / (2 * h_bar * epsilon0)

    # ---------------------------------------------------------
    # PARTIE IMAGINAIRE (ABSORPTION / GAIN)
    # ---------------------------------------------------------
    pauli_imag = (rho_e[:, m, n] + rho_h[:, m, n] - 1)
    
    integrand_abs = pauli_imag * Dr * (gamma / denom)
    integral_abs = np.trapezoid(integrand_abs, Et)
    
    # Chi Imaginaire Total
    chi_imag = prefactor * integral_abs
    
    partie_imaginaire[n] = chi_imag 

    # ---------------------------------------------------------
    # PARTIE RÉELLE (INDICE DE RÉFRACTION)
    # ---------------------------------------------------------
    pauli_real = (rho_e[:, m, n] + rho_h[:, m, n]) # Pas de -1 !
    
    integrand_ref = pauli_real * Dr * ((Et / h_bar) / denom)
    integral_ref = np.trapezoid(integrand_ref, Et)
    
    # Chi Réel 
    chi_real = prefactor * integral_ref
    
    partie_reelle[n] = chi_real

    # 3. Reconstruction de la susceptibilité complexe hybride
    chi_resonant = chi_real + 1j * chi_imag

    # 4. Calcul de n_eff
    # n0 contient le background, chi_resonant contient le Delta_n (réel) et l'alpha total (imag)
    n_eff = np.sqrt(n0**2 + chi_resonant)
    
    # Sécurités numériques
    if np.imag(n_eff) < 0:
        n_eff = np.conj(n_eff)

    if np.isnan(n_eff) or np.isinf(n_eff):
        # Fallback : attention n_lin_complex doit être cohérent avec cette méthode
        # Idéalement juste n0 si pas de populations
        return n0 + 0j 
    
    if abs(n_eff.imag) > 1.0 or abs(n_eff.real) > 5.0 or abs(n_eff.real) < 1.0:
        return n0 + 0j

    return n_eff


def propagation():
    """Résout la propagation couplant Maxwell et Bloch optiques sur la grille (z,t).
    """
    
    #On initialise les différents coefficients à l'aide de l'indice linéaire complexe et de la TMM.
    r_eff_air_gaas, t_eff_air_gaas, r_eff_gaas_air, t_eff_gaas_air = get_multistack_coeffs(couches_lineaires)


    # On stocke les coefficients qui seront utilisés au pas de temps n=1
    r_out_array[0] = r_eff_gaas_air
    t_in_array[0] = t_eff_air_gaas
    r_in_array[0] = r_eff_air_gaas
    t_out_array[0] = t_eff_gaas_air
    
    for n in range(1,Nt-1):
        #Condition limite à l'interface Air/GaAs
        E_plus[0,n] = t_in_array[n-1] * source_pulse[n] + r_out_array[n-1] * E_moins[0,n-1]
        
        #Propagation croissante    
        for m in range(0,Nz-1):
            F_plus(m,n)
            F_moins(m,n)
            F(m,n)
            E_combined = E_plus[m,n-1]*np.exp(-1j*k0*dz*m) + E_moins[m,n-1]*np.exp(1j*k0*dz*m)
            rho_calcul(n, m, E_combined)
            propagation_croissante(n,m)
            
        m_back = Nz-1
            
        F_plus(m_back, n)
        F_moins(m_back, n)
        F(m_back, n)
        
        # 2. Calcul du champ combiné au fond (basé sur n-1 pour interaction causale)
        E_combined_back = E_plus[m_back, n-1]*np.exp(-1j*k0*dz*m_back) + E_moins[m_back, n-1]*np.exp(1j*k0*dz*m_back)
        
        # 3. Mise à jour des populations (rho) au fond pour l'instant n
        rho_calcul(n, m_back, E_combined_back)

        # Coefficient de réflexion GaAs → Or
        r_gold = (n0 - n_gold) / (n_gold + n0)
        
        L = (Nz-1)*dz
        dephasage = np.exp(-2j*k0*L)
        
        
        # Condition limite arrière réaliste
        E_moins[-1, n] = r_gold * E_plus[-1, n] * dephasage
        
        #Propagation décroissante
        for m in range(Nz-1,0,-1):
            f(m,n) #On met à jour la polarisation totale simplement pour nos analyses sur la SVEA
            propagation_decroissante(n,m) 
  
    
        E[0, n] = E_plus[0,n]+E_moins[0,n]
        
        F(0,n)
        
        n_eff = calculate_total_neff(E[0,n], n, 0)
            
        # Stockage de n_eff
        n_eff_array[n] = n_eff
        
        # Calcul du coefficient d'absorption alpha = 4*pi*k/lambda = 2*omega*Im(n)/c
        alpha_array[n] = 2 * omega0 * np.imag(n_eff) / c
    
    
        r_eff_air_gaas, t_eff_air_gaas, r_eff_gaas_air, t_eff_gaas_air = get_multistack_coeffs(couches_lineaires)
        
        t_in_array[n] = t_eff_air_gaas
        t_out_array[n] = t_eff_gaas_air 
        r_in_array[n] = r_eff_air_gaas
        r_out_array[n] = r_eff_gaas_air
        
        N_entree[n] = np.trapezoid(rho_e[:, 0, n] * Dr, Et)
        
        #Calcul de E total
        for m in range(1,Nz):
           E[m,n] = E_plus[m,n]*np.exp(-1j*k0*dz*m) + E_moins[m,n]*np.exp(1j*k0*dz*m)
    return

# --------------------
# RUN
# --------------------

propagation()


# --------------------
# PLOT
# --------------------
# Visualisation
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Propagation spatio-temporelle
im = axes[0,0].imshow(np.abs(E).T, aspect='auto', 
                      extent=[z[0]*1e9, z[-1]*1e9, t[0]*1e12, t[-1]*1e12],
                      origin='lower', cmap='plasma')
axes[0,0].set_xlabel('Position z (nm)')
axes[0,0].set_ylabel('Temps (ps)')
axes[0,0].set_title('|E(z,t)| - Propagation')
plt.colorbar(im, ax=axes[0,0])

# Évolution temporelle de l'onde
axes[0,1].plot(t*1e12, np.abs(E[0, :]), 'g-', label='Entrée', linewidth=2)
axes[0,1].plot(t*1e12, np.abs(source_pulse), 'y-', label='Pulse d\'entrée', linewidth=2)
axes[0,1].plot(t*1e12, np.abs(E[Nz//2, :]), 'b-', label='Milieu')
axes[0,1].plot(t*1e12, np.abs(E[-1, :]), 'r-', label='Sortie')
axes[0,1].set_xlabel('Temps (ps)')
axes[0,1].set_ylabel('|E| (V/m)')
axes[0,1].set_title('Évolution temporelle')
axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)

# Profil spatial
t_snap = Nt//3
axes[1,0].plot(z*1e9, np.abs(E[:, t_snap]), 'b-', linewidth=2)
axes[1,0].set_xlabel('Position z (nm)')
axes[1,0].set_ylabel('|E| (V/m)')
axes[1,0].set_title(f'Profil spatial à t={t[t_snap]*1e12:.1f} ps')
axes[1,0].grid(True, alpha=0.3)

# Populations
axes[1,1].plot(t*1e12, np.mean(rho_e[:, Nz//2, :], axis=0), 'r-', label='<ρₑ> milieu')
axes[1,1].plot(t*1e12, np.mean(rho_e[:, 0, :], axis=0), 'g-', label='<ρₑ> entrée')
axes[1,1].plot(t*1e12, np.mean(rho_e[:, -2, :], axis=0), 'b-', label='<ρₑ> fin')
axes[1,1].set_xlabel('Temps (ps)')
axes[1,1].set_ylabel('Population')
axes[1,1].set_title('Populations à plusieurs endroits')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.plot()

fig2, axes2 = plt.subplots(2, 1, figsize=(14, 10))

# Partie réelle de n_eff
axes2[0].plot(t[1:len(t)-1]*1e12, np.real(n_eff_array[1:len(t)-1]), 'b-', linewidth=2)
ax_alpha = axes2[0].twinx()
ax_alpha.plot(t[1:len(t)-1]*1e12, np.real(n_eff_array[1:len(t)-1]) - np.real(n_eff_array[1]) , 'r-', linewidth=2, label='α')
axes2[0].set_xlabel('Temps (ps)')
axes2[0].set_ylabel('Re(n_eff)')
ax_alpha.set_ylabel('Delta n')
ax_alpha.tick_params(axis='y')
axes2[0].set_title('Partie réelle de l\'indice effectif')
axes2[0].grid(True, alpha=0.3)

# Partie imaginaire de n_eff
axes2[1].plot(t[1:len(t)-1]*1e12, np.imag(n_eff_array[1:len(t)-1]), 'r-', linewidth=2)
ax_alpha = axes2[1].twinx()
ax_alpha.plot(t[1:len(t)-1]*1e12, alpha_array[1:len(t)-1] - alpha_array[1] , 'b-', linewidth=2, label='α')
axes2[1].set_xlabel('Temps (ps)')
axes2[1].set_ylabel('Im(n_eff)')
ax_alpha.set_ylabel('Delta alpha (m^{-1})')
axes2[1].set_title('Partie imaginaire de l\'indice effectif')
axes2[1].grid(True, alpha=0.3)


#Tracé pour vérifier la convergence de la partie imaginaire et réelle de l'intégrale de polarisation de n_eff
fig2, axes2 = plt.subplots(2, 1, figsize=(14, 10))
axes2[0].plot(t[1:len(t)-1]*1e12, partie_reelle[1:len(t)-1], 'b-', linewidth=2)
axes2[0].set_xlabel('Temps (ps)')
axes2[0].set_ylabel('Re(chi_eff)')
axes2[0].set_title('Partie réelle de chi_eff')
axes2[0].grid(True, alpha=0.3)

axes2[1].plot(t[1:len(t)-1]*1e12, partie_imaginaire[1:len(t)-1], 'b-', linewidth=2)
axes2[1].set_xlabel('Temps (ps)')
axes2[1].set_ylabel('Im(chi_eff)')
axes2[1].set_title('Partie imaginaire de chi_eff')
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.plot()

plt.show()


#############
# Test de convergence
#############
E_tmax_space = np.linspace(0,40*Etmax,100)
integrale_imaginaire = []
integrale_reelle = []

for i,cutoff in enumerate(E_tmax_space):
    print(f'Itération : {i+1}')
    Et = np.linspace(-Et0,cutoff,NEt)
    Dr = D0*np.sqrt((Et+Et0)*(1+a*(Et+Et0)))*(1+2*a*(Et+Et0))
    alpha = np.exp(dt*(-1j*(Et)/h_bar-gamma))
    propagation()
    integrale_imaginaire.append(np.max(np.abs(partie_imaginaire)))
    integrale_reelle.append(np.max(np.abs(partie_reelle)))
    

#Tracé pour vérifier la convergence de la partie imaginaire et réelle de l'intégrale de polarisation de n_eff
fig2, axes2 = plt.subplots(2, 1, figsize=(14, 10))
axes2[0].plot(E_tmax_space/q, integrale_reelle, 'b-', linewidth=2)
axes2[0].set_xlabel('Cutoff (eV)')
axes2[0].set_ylabel('Re(chi_eff)')
axes2[0].set_title('Partie réelle de chi_eff')
axes2[0].grid(True, alpha=0.3)

axes2[1].plot(E_tmax_space/q, integrale_imaginaire, 'b-', linewidth=2)
axes2[1].set_xlabel('Cutoff (eV)')
axes2[1].set_ylabel('Im(chi_eff)')
axes2[1].set_title('Partie imaginaire de chi_eff')
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.plot()
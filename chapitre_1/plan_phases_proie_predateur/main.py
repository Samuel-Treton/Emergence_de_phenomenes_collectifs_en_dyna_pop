#!/usr/bin/env python3
"""
Plan de phase Lotka–Volterra (proie–prédateur)

Modèle :
    dX/dt =  r1 * X  - alpha * X * Y
    dY/dt =  beta * X * Y - r2 * Y

où :
	X : Population de proies
	Y : Population de prédateurs
	r1 : Taux de croissance intrinsèque des proies
	r2 : Taux de mortalité intrinsèque des prédateurs
	alpha : Taux de prédation (efficacité des prédateurs)
	beta : Taux de conversion de la biomasse des proies en biomasse des prédateurs

Ce script :
	1) Importe les bibliothèques nécessaires
	2) Définit les paramètres du modèle
	3) Crée une grille pour le plan de phase
	4) Calcule le champ de vecteurs (dX/dt, dY/dt)
	5) Trace les isoclines pour dX/dt = 0 et dY/dt = 0
	6) Trouve les équilibres via la méthode de Newton
	7) Intègre une trajectoire du système à partir d'une condition initiale (X0, Y0)
	8) Affiche le plan de phase et l'évolution temporelle des populations

Dépendances :
	pip install numpy matplotlib scipy
    ou
	pip install -r requirements.txt
"""

# ======================================================================
# 1. Imports
# ======================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp

# ======================================================================
# 2. PARAMÈTRES
# ======================================================================
# Paramètres du modèle Lotka-Volterra
r1, r2 = 1.0, 1.0      # Taux de croissance intrinsèque des proies et de mortalité intrinsèque des prédateurs
alpha, beta = 0.5, 0.5 # Taux de prédation et de conversion de la biomasse

assert r1 > 0 and r2 > 0, "Les taux r1 et r2 doivent être strictement positifs."
assert alpha > 0 and beta > 0, "Les paramètres alpha et beta doivent être strictement positifs."

# Limites du plan de phase
XMIN, XMAX = -0.2, 5        # Limites du plan de phase en X
YMIN, YMAX = -0.2, 5        # Limites du plan de phase en Y
h = 0.05                  # Pas de discrétisation de la grille (en espace)

# Paramètres de visualisation
LW_PHASES      = 2.2      # Épaisseur des lignes de phase
LW_ISOCLINE    = 5.0      # Épaisseur des isoclines
EQ_SIZE        = 300.0    # Taille des marqueurs d'équilibre
PHASES_DENSITY = 1.5      # Densité des lignes de phase
ARROW_SIZE     = 2.2      # Taille des flèches

COLOR_EQ0   = '#ff5555ff' # (0,0)
COLOR_EQ1   = '#55caffff' # (r2/beta, r1/alpha)

# Paramètres de la trajectoire tracée
X0, Y0 = 1.0, 1.0         # Une valeur initiale pour simuler le système
dt = 0.01                 # Pas de temps pour l'intégration
TMAX = 25.0               # Temps maximum pour l'intégration
COLOR_TRAJECTORY = "#5571ffff" # Couleur de la trajectoire
LW_TRAJECTORY = 5.0      # Épaisseur de la trajectoire

# ======================================================================
# 3. GRILLE & CHAMP DE VECTEURS
# ======================================================================
X_vals = np.arange(XMIN, XMAX + h, h)
Y_vals = np.arange(YMIN, YMAX + h, h)
XX, YY = np.meshgrid(X_vals, Y_vals)

def f(X, Y): # d(X,Y)/dt = f(X,Y)
     X_prime = r1 * X - alpha * X * Y
     Y_prime = beta * X * Y - r2 * Y
     return X_prime, Y_prime

X_prime, Y_prime = f(XX, YY)

# ======================================================================
# 4. FIGURE + PHASES + MASQUE
# ======================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8)) # Création de la figure et des axes
# Plot de gauche : Plan de phase (ax1)
# Plot de droite : Évolution temporelle (ax2)

# Tracé des lignes de phase dans le plan de phase
PHASES = ax1.streamplot(XX, YY, X_prime, Y_prime,
                     color='black',
                     linewidth=LW_PHASES,
                     density=PHASES_DENSITY,
                     arrowsize=ARROW_SIZE)

# Ajout d'un masque pour les zones non physiques (X ≤ 0 ou Y ≤ 0)
# Création du masque pour les zones où les populations sont négatives ou nulles
mask_x_neg = XX <= 0
mask_y_neg = YY <= 0
mask_non_physical = mask_x_neg | mask_y_neg

# Application du masque gris transparent
ax1.contourf(
	XX, YY, 
	mask_non_physical.astype(int),  # Masque binaire : 1 pour zones non physiques, 0 sinon
	levels=[0.5, 1.5],              # Remplit uniquement là où mask == 1
	colors=['gray'],                # Couleur grise pour le masque
	alpha=0.7,                      # Transparence du masque
	zorder=-10                      # Derrière tous les autres éléments
)

# ======================================================================
# 5. ISOCLINES NUMÉRIQUES (NIVEAU 0) - PLOT DE GAUCHE
# ======================================================================
# Trace l'isocline pour dX/dt = 0 (courbe où la variation de la population de proies est nulle)
ISO_X = ax1.contour(
	XX, YY, X_prime,        # Grille et valeurs du champ pour dX/dt
	levels=[0],             # Niveau 0 : isocline
	colors=['red'],         # Couleur rouge pour cette isocline
	linestyles='-',         # Ligne pleine
	linewidths=LW_ISOCLINE, # Épaisseur définie plus haut
	zorder=0                # Derrière les autres éléments
)

# Trace l'isocline pour dY/dt = 0 (courbe où la variation de la population de prédateurs est nulle)
ISO_Y = ax1.contour(
	XX, YY, Y_prime,        # Grille et valeurs du champ pour dY/dt
	levels=[0],             # Niveau 0 : isocline
	colors=['green'],       # Couleur verte pour cette isocline
	linestyles='-',         # Ligne pleine
	linewidths=LW_ISOCLINE, # Épaisseur définie plus haut
	zorder=0                # Derrière les autres éléments
)

# Création de "proxy artists" pour ajouter les isoclines à la légende.
# Ces objets ne sont pas tracés sur la figure, mais servent uniquement à la légende.
proxy_dx = Line2D(
	[], [],                 # Pas de données réelles
	ls='-',                 # Ligne pleine
	lw=LW_ISOCLINE,         # Même épaisseur que l'isocline
	color='red',            # Même couleur que l'isocline dX/dt=0
	label=r'$\dot X = 0$'   # Légende en LaTeX
)
proxy_dy = Line2D(
	[], [],                 # Pas de données réelles
	ls='-',                 # Ligne pleine
	lw=LW_ISOCLINE,         # Même épaisseur que l'isocline
	color='green',          # Même couleur que l'isocline dY/dt=0
	label=r'$\dot Y = 0$'   # Légende en LaTeX
)

# ======================================================================
# 6. ÉQUILIBRES VIA NEWTON
# ======================================================================
def f_vec(z):
	# Besoin d'une fonction avec un seul argument pour root (scipy.optimize) (Newton)
	X, Y = z
	return np.array([r1*X - alpha*X*Y,
					 beta*X*Y - r2*Y])

# Recherche des équilibres via la méthode de Newton
initial_guess_trivial = [1e-8, 1e-8]  # Point de départ pour l'équilibre trivial
initial_guess_nontrivial = [2.0, 2.0]  # Point de départ pour l'équilibre non trivial

# Recherche numérique de l'équilibre trivial (X=0, Y=0) avec Newton
sol_trivial = root(f_vec, initial_guess_trivial)
# Recherche numérique de l'équilibre non trivial (X=r2/beta, Y=r1/alpha) avec Newton
sol_nontriv = root(f_vec, initial_guess_nontrivial)

# Vérification de la convergence des solutions
if not sol_trivial.success or not sol_nontriv.success:
    raise RuntimeError("Newton n'a pas convergé.")

# Extraction des coordonnées des équilibres
E_trivial_X, E_trivial_Y = sol_trivial.x
E_star_X, E_star_Y = sol_nontriv.x

# Calcul des équilibres théoriques
E_trivial_theorique = (0.0, 0.0)
E_star_theorique = (r2 / beta, r1 / alpha)

# Retour print des équilibres
print("Équilibres théoriques :")
print(f"  Trivial      : X = {E_trivial_theorique[0]:.6f}, Y = {E_trivial_theorique[1]:.6f}")
print(f"  Non trivial  : X = {E_star_theorique[0]:.6f}, Y = {E_star_theorique[1]:.6f}")
print()
print("Équilibres trouvés par Newton :")
print(f"  Trivial      : X = {E_trivial_X:.6f}, Y = {E_trivial_Y:.6f}")
print(f"  Non trivial  : X = {E_star_X:.6f}, Y = {E_star_Y:.6f}")

# Tracé des équilibres sur le plot de gauche
ax1.scatter(
	[E_trivial_X], [E_trivial_Y],                # Coordonnées de l'équilibre trivial (X=0, Y=0)
	s=EQ_SIZE,                                   # Taille du marqueur (définie plus haut)
	color=COLOR_EQ0,                             # Couleur du marqueur (équilibre trivial)
	edgecolors='black',                          # Couleur du contour du marqueur
	linewidths=2.0,                              # Épaisseur du contour
	label='Équilibre trivial (Newton)',          # Légende pour la figure
	zorder=10                                    # Priorité d'affichage (au-dessus des autres éléments)
)
ax1.scatter(
	[E_star_X], [E_star_Y],                      # Coordonnées de l'équilibre non trivial
	s=EQ_SIZE,                                   # Taille du marqueur (définie plus haut)
	color=COLOR_EQ1,                             # Couleur du marqueur (équilibre non trivial)
	edgecolors='black',                          # Couleur du contour du marqueur
	linewidths=2.0,                              # Épaisseur du contour
	label='Équilibre non trivial (Newton)',      # Légende pour la figure
	zorder=10                                    # Priorité d'affichage (au-dessus des autres éléments)
)

# ======================================================================
# 7. INTÉGRATION DE LA TRAJECTOIRE
# ======================================================================
def f_for_solve_ivp(t, Z):
	# solve_ivp prend en argument un champ de vecteurs de la forme
	# f_for_solve_ivp(t, Z)
    # Note : l'EDO est possiblement non autonome, d'où l'argument t
	X, Y = Z
	return f(X, Y)  # Retourne le champ de vecteurs (dX/dt, dY/dt)

# Intégration de la trajectoire avec une méthode d'ordre supérieur
t_range = (0, TMAX)
dt = 0.01
t_eval = np.arange(0, TMAX + dt, dt)
sol = solve_ivp(
	f_for_solve_ivp,      # Fonction du système différentiel
	t_range,              # Intervalle de temps d'intégration (début, fin)
	[X0, Y0],             # Conditions initiales : populations initiales de proies et prédateurs
	t_eval=t_eval,        # Points de temps où la solution est évaluée
	method='DOP853',      # DOP853 : Méthode de Dormand-Prince d'ordre 8 avec contrôle d'erreur adaptatif
	dense_output=True     # Permet d'obtenir une solution continue (interpolation)
)

# Ajout de la trajectoire dans le plan de phase (plot de gauche)
ax1.plot(
	sol.y[0], sol.y[1], 
	color=COLOR_TRAJECTORY,      # Couleur de la trajectoire dans le plan de phase
	linewidth=LW_TRAJECTORY,     # Épaisseur de la ligne de trajectoire
	label=f'Trajectoire depuis ({X0}, {Y0})',  # Légende pour la trajectoire
	zorder=1                     # Priorité d'affichage
)
ax1.scatter(
	[X0], [Y0], 
	s=200,                       # Taille du marqueur pour la condition initiale
	color=COLOR_TRAJECTORY,      # Même couleur que la trajectoire
	edgecolors='black',          # Contour noir du marqueur
	linewidths=2,                # Épaisseur du contour
	label=f'Une donnée initiale : ({X0}, {Y0})',  # Légende pour la condition initiale
	zorder=10                    # Priorité d'affichage
)

# Plot à droite : X(t) et Y(t) en fonction du temps
ax2.plot(
	sol.t, sol.y[0], 
	'r-',                        # Ligne rouge pour X(t)
	linewidth=3,                 # Épaisseur de la courbe
	label='X(t) - Population de proies'  # Légende
)
ax2.plot(
	sol.t, sol.y[1], 
	'g-',                        # Ligne verte pour Y(t)
	linewidth=3,                 # Épaisseur de la courbe
	label='Y(t) - Population de prédateurs'  # Légende
)

# ======================================================================
# 8. MISE EN FORME DES PLOTS
# ======================================================================
# Configuration du plot de gauche (plan de phase)
ax1.set_xlim(XMIN, XMAX)  # Limites de l'axe X (population de proies)
ax1.set_ylim(YMIN, YMAX)  # Limites de l'axe Y (population de prédateurs)
ax1.set_xlabel('Population de proies (X)', fontsize=16)  # Label de l'axe X
ax1.set_ylabel('Population de prédateurs (Y)', fontsize=16)  # Label de l'axe Y
ax1.set_title('Plan de phase', fontsize=18)  # Titre du graphique
ax1.tick_params(axis='both', which='major', labelsize=14)  # Taille des labels des axes
ax1.grid(True, linestyle=':', alpha=0.4)  # Grille en pointillés, légère

# Configuration du plot de droite (évolution temporelle)
ax2.set_xlim(0, TMAX)  # Limites de l'axe X (temps)
ax2.set_ylim(0, max(max(sol.y[0]), max(sol.y[1])) * 1.1)  # Limites de l'axe Y (populations)
ax2.set_xlabel('Temps (t)', fontsize=16)  # Label de l'axe X
ax2.set_ylabel('Taille des populations', fontsize=16)  # Label de l'axe Y
ax2.set_title(f'Évolution temporelle depuis ({X0}, {Y0})', fontsize=18)  # Titre du graphique
ax2.tick_params(axis='both', which='major', labelsize=14)  # Taille des labels des axes
ax2.grid(True, linestyle=':', alpha=0.4)  # Grille en pointillés, légère

# Ajoute les isoclines à la légende (proxies), puis affiche la légende complète
# Récupère les handles et labels existants de la légende du plot de gauche
handles1, labels1 = ax1.get_legend_handles_labels()
# Ajoute les proxies pour les isoclines (dX/dt=0 et dY/dt=0) à la légende
handles1.extend([proxy_dx, proxy_dy])
ax1.legend(handles=handles1, loc='upper right', frameon=True, facecolor='white', framealpha=1.0)  # Affiche la légende

# Légende pour le plot de droite
ax2.legend(loc='upper right', frameon=True, facecolor='white', framealpha=1.0)

# Affiche le plot
plt.tight_layout() # Ajuste la mise en page
plt.savefig("plan_phase_proie_predateur.svg", format="svg") # Sauvegarde en vectoriel (privilégier)
# plt.savefig("plan_phase_proie_predateur.png", format="png") # Sauvegarde en bitmap
plt.show()

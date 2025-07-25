#!/usr/bin/env python3
"""
Plan de phase pour le système proie-prédateur de Lotka–Volterra.

    dX/dt = r1 * X - alpha * X * Y
    dY/dt = beta * X * Y - r2 * Y

- Utilise scipy.integrate.odeint pour intégrer le système.
- Affiche :
    * le champ de vecteurs,
    * quelques trajectoires issues de conditions initiales variées,
    * les isoclines : dX/dt = 0 et dY/dt = 0,
    * les points d'équilibre (0,0) et (r2/beta, r1/alpha).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import root

# ========================================================================
# PARAMÈTRES DU MODÈLE
# ========================================================================
r1 = 1.0   	  # Taux de croissance intrinsèque des proies en absence de prédateurs
r2 = 1.0   	  # Taux de mortalité des prédateurs en absence de proies
alpha = 0.5   # Taux de prédation (impact des prédateurs sur les proies)
beta = 0.5    # Efficacité de conversion des proies en biomasse de prédateurs

assert r1 > 0 and r2 > 0, "Les taux de croissance et de mortalité doivent être positifs."
assert alpha > 0 and beta > 0, "Les taux de prédation et de conversion doivent être positifs."

dt = 1e-3  	  # Pas de temps pour l'intégration
TMAX = 1  # Temps maximum pour l'intégration

h = 1e-2      # Pas d'espace pour le champ de vecteurs

# Bornes du graphe
XMIN=0
XMAX=5
YMIN=0
YMAX=5

# ========================================================================
# FONCTION D'ÉVOLUTION DU SYSTÈME DE LOTKA-VOLTERRA
# ========================================================================
def f(X, Y):
	"""
	Équations différentielles du système Lotka–Volterra :
    renvoie le vecteur dérivé [dX/dt, dY/dt]
    en fonction des populations X et Y.
	"""
	global r1, r2, alpha, beta
	dXdt = r1 * X - alpha * X * Y
	dYdt = beta * X * Y - r2 * Y
	return np.array([dXdt, dYdt])

# ========================================================================
# GÉNÉRATION D'UNE GRILLE DE POINTS POUR LE CHAMP DE VECTEURS
# ========================================================================
X_vals = np.arange(XMIN, XMAX + h, h)
Y_vals = np.arange(YMIN, YMAX + h, h)
XX, YY = np.meshgrid(X_vals, Y_vals)

# Calcul du champ de vecteurs sur la grille (noter que tout est vectorisé)
dX = r1 * XX - alpha * XX * YY
dY = beta * XX * YY - r2 * YY

# --- 2.2. Tracé de base -------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))  # Taille raisonnable pour l'affichage

speed = np.sqrt(dX**2 + dY**2)  # Module du vecteur, pour moduler l'épaisseur
strm = ax.streamplot(XX, YY, dX, dY, color=speed, linewidth=1, density=1.2,
						arrowsize=1.2, cmap='viridis')
cbar = fig.colorbar(strm.lines, ax=ax, label='||f(X,Y)||')

# ========================================================================
# GÉNÉRATION DES ISOCLINES
# ========================================================================
# Recherche des isoclines en résolvant f(X, Y) = 0 sur la grille
# On cherche les points où dX/dt ≈ 0 et dY/dt ≈ 0

# Tolérance pour considérer une valeur comme "zéro"
tol = 1e-3

# Isocline dX/dt = 0
isocline_dX = np.abs(dX) < tol
ax.contour(XX, YY, isocline_dX, levels=[1], colors='tab:red', linestyles='--', linewidths=1.8, label='dX/dt = 0')

# Isocline dY/dt = 0
isocline_dY = np.abs(dY) < tol
ax.contour(XX, YY, isocline_dY, levels=[1], colors='tab:blue', linestyles='--', linewidths=1.8, label='dY/dt = 0')

# ========================================================================
# GÉNÉRATION DES ÉQUILIBRES
# ========================================================================
# Fonction pour scipy.optimize.root (Newton-Raphson)
def f_vec(z):
	X, Y = z
	return f(X, Y)

# Recherche des points d'équilibre avec Newton (méthode 'hybride' de scipy)
equilibres = []
# 1. Équilibre trivial
sol_trivial = root(f_vec, [0, 0], method='hybr')
if sol_trivial.success:
	X0, Y0 = sol_trivial.x
	ax.plot(X0, Y0, 'ko', ms=6, label='Équilibre trivial')
	equilibres.append((X0, Y0))

# 2. Équilibre non trivial (utilise une estimation initiale proche)
X_guess = r2 / beta
Y_guess = r1 / alpha
sol_nontrivial = root(f_vec, [X_guess, Y_guess], method='hybr')
if sol_nontrivial.success:
	X_star, Y_star = sol_nontrivial.x
	ax.plot(X_star, Y_star, 'ks', ms=6, label='Équilibre non trivial')
	equilibres.append((X_star, Y_star))

# ========================================================================
# MISE EN FORME ESTHÉTIQUE
# ========================================================================
ax.set_xlim(XMIN, XMAX)
ax.set_ylim(YMIN, YMAX)
ax.set_xlabel('Population de proies (X)')
ax.set_ylabel('Population de prédateurs (Y)')
ax.set_title('Plan de phase du système proie–prédateur (Lotka–Volterra)')
ax.grid(True, linestyle=':', alpha=0.5)

# plt.tight_layout()
plt.show()

# # --- 2.4. Trajectoires --------------------------------------------------
# # Temps pour l'intégration
# t = np.linspace(0, TMAX, int(TMAX / dt) + 1)

# # Génération automatique de conditions initiales (CI) :
# # On évite X=Y=0 (point trivial)
# rng = np.random.default_rng(seed=0)  # seed pour reproductibilité
# X0_list = rng.uniform(low=XMIN + 1e-3, high=XMAX, size=n_traj)
# Y0_list = rng.uniform(low=YMIN + 1e-3, high=YMAX, size=n_traj)

# for X0, Y0 in zip(X0_list, Y0_list):
# 	z0 = [X0, Y0]
# 	sol = odeint(lotka_volterra, z0, t, args=(r1, r2, alpha, beta))
# 	# On ne trace que la partie dans la fenêtre pour éviter du hors-champ disgracieux
# 	ax.plot(sol[:, 0], sol[:, 1], lw=1.2)

# # --- 2.5. Points d'équilibre -------------------------------------------
# # (0,0) et (r2/beta, r1/alpha) si alpha,beta ≠ 0
# ax.plot(0, 0, 'ko', ms=6, label='Équilibre trivial')
# if alpha != 0 and beta != 0:
# 	ax.plot(r2 / beta, r1 / alpha, 'ks', ms=6, label='Équilibre non trivial')

# # --- 2.6. Mise en forme esthétique -------------------------------------
# ax.set_xlim(XMIN, XMAX)
# ax.set_ylim(YMIN, YMAX)
# ax.set_xlabel('Population de proies (X)')
# ax.set_ylabel('Population de prédateurs (Y)')
# ax.set_title('Plan de phase du système proie–prédateur (Lotka–Volterra)')
# ax.grid(True, linestyle=':', alpha=0.5)

# # Légende : on ne veut pas d'entrées dupliquées
# handles, labels = ax.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# ax.legend(by_label.values(), by_label.keys(), loc='best')

# plt.tight_layout()
# return fig, ax

# # -------------------------------------------------------------------------
# # 3. Exemple d'utilisation (à supprimer si vous importez ce script ailleurs)
# # -------------------------------------------------------------------------
# if __name__ == "__main__":
# # Paramètres par défaut (modifiez à volonté)
# r1    = 1.0
# r2    = 1.0
# alpha = 0.5
# beta  = 0.5

# # Bornes du graphe
# XMIN, XMAX = 0, 5
# YMIN, YMAX = 0, 5

# # Appel à la fonction de tracé
# fig, ax = plot_phase_plane(r1, r2, alpha, beta,
# 							XMIN, XMAX, YMIN, YMAX,
# 							n_grid=25, n_traj=12, T=60, n_time=2000,
# 							stream=True)

# # Affichage
# plt.show()

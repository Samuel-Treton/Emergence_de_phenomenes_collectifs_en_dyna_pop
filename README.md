# Émergence de phénomènes collectifs en dynamique des populations

[![Licence](https://img.shields.io/badge/Licence-CC%20BY--SA%204.0-blue.svg)](https://creativecommons.org/licenses/by-sa/4.0/)
[![LaTeX](https://img.shields.io/badge/LaTeX-Document-brightgreen.svg)](M2_modelisation_poly_cours.pdf)

## 📚 À propos

Ce repository contient le matériel pédagogique du cours **"Émergence de phénomènes collectifs en dynamique des populations"** (XMS2MU100), dispensé dans le cadre du Master 2 MACS (Modélisation, Analyse Numérique et Calcul Scientifique) de l'Université de Nantes.

Le cours propose une introduction complète aux équations aux dérivées partielles de type réaction-diffusion, utilisées pour modéliser l'évolution spatio-temporelle de populations biologiques dans des contextes variés relevant de l'écologie ou de l'épidémiologie.

## 🎯 Objectifs pédagogiques

- Maîtriser les outils fondamentaux des équations différentielles ordinaires en dynamique des populations
- Comprendre la modélisation du mouvement de populations par diffusion
- Analyser les équations de réaction-diffusion et leurs applications biologiques
- Étudier en détail l'équation de Fisher-KPP et les phénomènes d'invasion

## 📖 Structure du cours

Le polycopié est organisé en quatre chapitres progressifs :

### **Chapitre 1 : Équations différentielles ordinaires (EDO)**
- Modélisation des variations de tailles de populations : $\dot{X} = f(X)$
- Outils d'analyse qualitative des solutions
- Exemples classiques en dynamique des populations

### **Chapitre 2 : Équations de diffusion (chaleur)**
- Modélisation du mouvement par dispersion diffusive : $\partial_t u = \Delta u$
- Interprétations individu-centrées de l'équation de la chaleur
- Propriétés essentielles en domaines bornés et non-bornés

### **Chapitre 3 : Équations de réaction-diffusion**
- Introduction aux équations : $\partial_t u = \Delta u + f(u)$
- Principe de comparaison et outils analytiques
- Applications à des exemples biologiques concrets

### **Chapitre 4 : Équation de Fisher-KPP**
- Étude approfondie de : $\partial_t u = \Delta u + u(1-u)$
- Analyse des solutions sur la droite réelle
- Ondes progressives et vitesse de propagation

## 📁 Contenu du repository

```
├── M2_modelisation_poly_cours.tex     # Document principal LaTeX
├── M2_modelisation_poly_cours.pdf     # Polycopié compilé
├── images/                            # Figures et diagrammes
│   ├── chap_1/                       # Illustrations chapitre 1
│   ├── annexe/                       # Figures annexes
│   └── DALL-E/                       # Images générées par IA
├── document_in_preparation/           # Matériel en développement
└── README.md                          # Ce fichier
```

## 🔧 Compilation du document

Le document LaTeX peut être compilé avec :

```bash
pdflatex M2_modelisation_poly_cours.tex
```

Ou utiliser le script fourni :
```bash
./LAUNCH.sh
```

## 👨‍🏫 Enseignant

**Samuel Tréton**  
Laboratoire de Mathématiques Jean Leray  
Nantes Université, France  
🌐 [Site web](https://www.samueltreton.fr/)

## 📚 Bibliographie recommandée

- **T. Giletti** - *EDP paraboliques pour la dynamique des populations* [(polycopié en ligne)](https://tgiletti.perso.math.cnrs.fr/pdf/M2MFA_EDP%20Dyn%20Pop.pdf)
- **L. Roques** - *Modèles de réaction-diffusion pour l'écologie spatiale*, Éditions Quae, 2013
- **R. A. Horn, C. R. Johnson** - *Matrix Analysis*, Cambridge University Press, 2012

## 🎓 Public cible

- Étudiants de Master 2 en mathématiques appliquées
- Chercheurs en modélisation mathématique en biologie
- Toute personne intéressée par les équations de réaction-diffusion

## 📊 Volume horaire

- **16h de cours magistraux**
- **12h de travaux dirigés**

## 🤝 Contributions

Les suggestions d'amélioration et corrections sont les bienvenues. N'hésitez pas à ouvrir une issue ou proposer une pull request.

## 📄 Licence

Ce matériel pédagogique est mis à disposition sous licence [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/).

---

*Dernière mise à jour : Juillet 2025*

C     ==================================================================
C     I     Approximation des Défauts Ponctuels Indépendants (ADPI)    I
C     I             pour la description thermodynamique                I
C     I                 d'un alliage/composé ordonné                   I
C     I             (composé intermétallique, minéral,...)             I
C     I        avec un, deux ou trois éléments intrinsèques            I
C     I                  et des nombres quelconques                    I
C     I          de sous-réseaux et d'éléments d'addition              I
C     ==================================================================
C     I                                                                I
C     I Grandeurs physiques calculées :                                I
C     I -----------------------------                                  I
C     I * quantités de défauts ponctuels                               I
C     I * quantités thermodynamiques                                   I
C     I (énergie libre, entropie, potentiels chimiques)                I
C     I en fonction des variables thermodynamiques intensives :        I
C     I composition = écart à la stoechiométrie, température           I
C     I (et éventuellement pression)                                   I
C     I                                                                I
C     I Les r1 premiers sous-réseaux contiennent l'espèce 1            I
C     I les r2 suivants l'espèce 2                                     I
C     I les r3 suivants l'espèce 3                                     I
C     I dans l'état fondamental sans défaut                            I
C     I                                                                I
C     I Aucun sous-réseau ne contient les éléments d'addition          I
C     I (i>3) à l'état fondamental                                     I
C     I                                                                I
C     I Les sous-réseaux d'indices > r1 + r2 + r3 sont interstitiels   I
C     I                                                                I
C     I On note p(r) le nombre de sites du sous-réseau r par maille    I
C     I                                                                I
C     I                 ============================                   I
C     I                 Trois modes de calcul ADPI :                   I
C     I                 ============================                   I
C     I                                                                I
C     I --------------------------------------------------             I
C     I Mode muVT --> ensemble grand canonique (mu(i),V,T)             I
C     I --------------------------------------------------             I
C     I                                                                I
C     I Relation de Gibbs-Duhem approchée pour la pression :           I
C     I à cause de cette relation P(mu(i),T) approchée,                I
C     I le balayage en potentiels chimiques résultant induit           I
C     I  une (légère) dérive en pression par rapport à la consigne.    I
C     I -> ajout d'une possibilité de boucle d'autocohérence sur mu(1) I
C     I pour éviter cette dérive en P.                                 I
C     I                                                                I
C     I En mode muVT, le programme ADPI effectue :                     I
C     I A) à partir de la relation de Gibbs-Duhem, le calcul approché  I
C     I    de la valeur de N(1)*mu(1) + N(2)*mu(2) + N(3)*mu(3)        I
C     I    correspondant à la pression prescrite                       I
C     I                                                                I
C     I    Attention : il s'agit de la pression externe                I
C     I                exercée sur le système par l'extérieur          I
C     I                                                                I
C     I B) un balayage en écarts de potentiels chimiques               I
C     I d_mu_i = mu(i_réf) - mu(i) (i = 1, 2 ou 3 <> i_réf)	       I
C     I             et en mu(i>3)                                      I
C     I								       I
C     I C) pour chaque valeur de d_mu_i et mu(i>3)                     I
C     I       l'écriture dans les fichiers ".adpi" des grandeurs :     I
C     I     (i) composition                                            I
C     I    (ii) quantités de défauts ponctuels x_déf, avec déf =       I
C     I  1(r) pour r1 < r <= r1 + r2 + r3                              I
C     I         ou r > r1 + r2 + r3 (1 interstitiel),                  I
C     I  2(r) pour r <= r1 ou r1 + r2 < r <= r1 + r2 + r3              I
C     I         ou r > r1 + r2 + r3 (2 interstitiel),                  I
C     I  3(r) pour r <= r1 + r2                                        I
C     I         ou r > r1 + r2 + r3 (3 interstitiel),                  I
C     I  L(r) pour r <= r1 + r2 + r3                                   I
C     I  i(r) pour tout r                                              I
C     I   (iii) énergie, volume, entropie de configuration             I
C     I         et enthalpie libre par atome ou par maille             I
C     I                                                                I
C     I Possibilité de fenêtres en composition pour l'écriture         I
C     I                                                                I
C     I Il a été ajouté ultérieurement                                 I
C     I la possibilité d'une boucle d'autocohérence sur mu(1)          I
C     I pour compenser la dérive en P                                  I
C     I induite par la relation approchée P(mu(i),T).                  I
C     I                                                                I
C     I -------------------------------------                          I
C     I Mode NPT=0 --> ensemble (N(i),P,T=0K)                          I
C     I -------------------------------------                          I
C     I                                                                I
C     I Minimisation de H (méthode du simplexe)                        I
C     I par rapport aux quantités de DP                                I
C     I sous les contraintes de quantités de matière constantes        I
C     I                                                                I
C     I Variantes disponibles pour le mode NPT=0 : P/p, x              I
C     I                                                                I
C     I --------------------------------                               I
C     I Mode NPT --> ensemble (N(i),P,T)                               I
C     I --------------------------------                               I
C     I                                                                I
C     I Résolution par la méthode de Newton-Raphson (NR)               I
C     I du système non linéaire d'inconnues ( M , x_d )                I
C     I correspondant à la minimisation de l'enthalpie libre           I
C     I                                                                I
C     I Variantes disponibles pour le mode NPT : P/p, T, x             I
C     I                                                                I
C     I Remarques :                                                    I
C     I un calcul NPT peut diverger si les paramètres sont mal choisis I
C     I => penser à jouer sur les paramètres suivants :                I
C     I * ordre (croissant / décroissant) du balayage en T ou x        I
C     I * finesse du balayage (nombre de points)                       I
C     I * valeurs initiales des inconnues                              I
C     I * fréquence de calcul de la matrice jacobienne                 I
C     I et également (plus rare)                                       I
C     I * nombre maximal d'itérations de NR                            I
C     I * précision pour l'arrêt de l'algorithme de NR                 I
C     I                                                                I
C     I ====================================================	       I
C     I          Possibilités de prise en compte :                     I
C     I * de DP complexes (approximative, en muVT seulement)           I
C     I * de DP chargés (en muVT"point" seulement)                     I
C     I ====================================================	       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 24/05/2022                              I
C     I Mardi des Rogations                                            I
C     ==================================================================
      PROGRAM ADPI
      USE CONSTANTES
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C############################
C############################
C#####Déclarations et formats
C############################
C############################
C###################################
C#####Tableaux de dimension variable
C###################################
C--------------------------------------------------------
C-----Commentaires écrits dans le fichier de compte-rendu
C--------------------------------------------------------
      CHARACTER * 100 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          COMMENTAIRE
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT
C-------------------------------------------------------------
C-----Coefficients des contraintes relatives
C-----aux compositions à écrire dans les fichiers de résultats
C-----et demi-largeurs de fenêtres correspondantes
C-----(ainsi qu'un tableau de valeurs initiales
C-----pour la commodité de lecture en raison des indices
C----- - absence de I_TYP_0)
C-------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          COEF_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_X_AT_INIT
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_X_AT
C---------------------------------------------------
C-----Matrices de décomposition LU pour le calcul
C-----des fractions atomiques issues des contraintes
C---------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          MAT_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          L_MAT_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          U_MAT_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          P_MAT_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          Q_MAT_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          MAT_INV_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          VECT_CTR
C-----------------------------------------------
C-----Fractions atomiques "centrales" et limites
C-----à écrire dans les fichiers de résultats
C-----------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_0_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_INF_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_SUP_CTR
C---------------------------------------------------
C-----Potentiels chimiques des éléments intrinsèques
C---------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_INTR
C-------------------------------------------------
C-----Potentiels chimiques des éléments d'addition
C-------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I
C--------------------------------------------
C-----Termes relatifs aux éléments d'addition
C--------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          ALPHA_I_R
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          SOMME_I
C------------------------------------------------------------------
C-----Valeurs initiales, nombres de pas et incréments
C-----des écarts de potentiels chimiques
C-----pour les éléments intrinsèques autres que celui de référence
C-----------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_POT_REF_INTR_INIT
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_D_POT_REF_INTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          PAS_D_POT_REF_INTR
C----------------------------------------------------
C-----Valeurs initiales, nombres de pas et incréments
C-----de potentiels chimiques des éléments d'addition
C----------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I_INIT
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_POT_I
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          PAS_POT_I
C---------------------------------------------------------
C-----Potentiels chimiques moyens des éléments d'addition,
C-----leurs carrés et les variances correspondantes
C-----pour les valeurs écrites
C---------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I_MOY
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I_2_MOY
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          VAR_POT_I
C--------------------------------------------------------------
C-----Produits partiels des nombres de pas
C-----pour les potentiels chimiques d'addition
C-----(utile pour la simulation de plusieurs boucles imbriquées
C-----en potentiels chimiques avec un seul indice)
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: 
     $                                          PROD_PART_N_POT_I
C------------------------------------------------
C-----Fractions d'antisites et lacunes (indice 0)
C-----dans les divers sous-réseaux
C------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_D_R
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-réseau
C-------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          P_R
C-----------------------------------------------------------------
C-----Nombre de sous-réseaux pour chaque espèce intrinsèque
C-----(utile pour le cas d'une extension à N espèces intrinsèques)
C-----------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_R_N
C-----------------------------------
C-----Nombres cumulés correspondants
C-----------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          RHO_N
C====================================================================
C====================================================================
C=====Paramètres GC pour des défauts simples neutres (INDIC_CHARGE=0)
C====================================================================
C====================================================================
C--------------------------------------------------------------
C-----Energies "brutes" des défauts sur les divers sous-réseaux
C--------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          E_B_D_R
C--------------------------------------------------------
C-----Energies GC des défauts sur les divers sous-réseaux
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          E_GC_D_R
C----------------------------------
C-----Enthalpies GC correspondantes
C----------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          H_GC_D_R
C-------------------------------------
C-----Mêmes quantités pour les volumes
C-------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          V_B_D_R
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          V_GC_D_R
C--------------------------------------------------------------------
C-----Enthalpies de formation des défauts sur les divers sous-réseaux
C--------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          H_FORM_D_R
C--------------------------------------
C-----Termes complémentaires d'entropie
C--------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          Z_TYP_R
C================================================
C================================================
C=====Paramètres pour des défauts simples chargés
C================================================
C================================================
C--------------------------------------------------
C-----Densité d'états électroniques :
C-----pour chaque point, valeur de l'énergie et DdE
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          TAB_DDE
C-------------------------------------------------------
C-----Tableau des nombres d'états de charges pour les DP
C-------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          NQ_D_R_Q
C---------------------------------------------------------------
C-----Charges (entières) des défauts sur les divers sous-réseaux
C---------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          Q_D_R_Q
C--------------------------------------------------------------
C-----Energies "brutes" des défauts sur les divers sous-réseaux
C--------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          E_B_D_R_Q
C--------------------------------------------------------
C-----Energies GC des défauts sur les divers sous-réseaux
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          E_GC_D_R_Q
C----------------------------------
C-----Enthalpies GC correspondantes
C----------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          H_GC_D_R_Q
C-------------------------------------
C-----Mêmes quantités pour les volumes
C-------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          V_B_D_R_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          V_GC_D_R_Q
C--------------------------------------------------------------------
C-----Enthalpies de formation des défauts sur les divers sous-réseaux
C--------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          H_FORM_D_R_Q
C===================================================================
C===================================================================
C=====Paramètres relatifs aux défauts complexes (pour ADPI sans chg)
C===================================================================
C===================================================================
C-----------------
C-----Multiplicité
C-----------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                MULTIPLICITE_COMPLEXE
C-----------------------------------------------------------
C-----Sous-réseau sur lequel est calculée cette multiplicité
C-----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                I_S_R_MULTIPLICITE_COMPLEXE
C--------------------------------
C-----Nombre de sites du complexe
C--------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NOMBRE_SITES_COMPLEXE
C--------------------------------------------------------------
C-----Numéros de sous-réseaux des sites occupés par le complexe
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_S_R_COMPLEXE
C-----------------------------------------------
C-----Types chimiques du complexe sur ces sites,
C-----dans le même ordre (0 = lacune)
C-----------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_TYPE_COMPLEXE
C-------------------------------------------------
C-----Energie et volume "bruts" de chaque complexe
C-------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  E_B_D_COMPLEXE 
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  V_B_D_COMPLEXE
C-------------------------------------------------------
C-----Energie, volume et enthalpie GC de chaque complexe
C-------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  E_GC_D_COMPLEXE                 
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  V_GC_D_COMPLEXE
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  H_GC_D_COMPLEXE
C--------------------------------------------------------
C-----Paramètres indicateurs de types et de sous-réseaux
C-----utiles au calcul des termes de potentiels chimiques  
C-----associés aux complexes
C--------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  ALPHA_TYPE_COMPLEXE
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  BETA_S_R_COMPLEXE
C--------------------------------------------------------------------
C-----Type chimique normal de chaque sous-réseau
C-----(0 pour sous-réseaux interstitiels)
C-----utile au calcul des termes de potentiels chimiques de complexes
C--------------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  I_TYPE_NORMAL_S_R
C----------------------------------------------------------
C-----Nombre de sites du sous-réseau r dans chaque complexe
C-----et nombre d'atomes de type i dans ce complexe
C----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  U_COMPLEXE_S_R
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  V_COMPLEXE_TYPE
C---------------------------
C-----Fractions de complexes
C---------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  X_D_COMPLEXE
C--------------------------------------------------
C-----Enthalpies de formation des défauts complexes
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          H_FORM_D_COMPLEXE
C-------------------------------------------
C-----Termes de somme sur les complexes
C-----relatifs aux éléments d'addition
C-----dans le calcul des fractions atomiques
C-------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          SOMME_COMPLEXE_V_TYPE_I
C===================================================================
C===================================================================
C=====Paramètres relatifs aux défauts complexes chargés (suffixe _Q)
C===================================================================
C===================================================================
C---------------------------------------
C-----Indicateur de DP complexes chargés
C---------------------------------------
       CHARACTER INDIC_COMPLEXES_Q
C--------------------------------------------------
C-----Nombre d'états de charge pour chaque complexe
C--------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NQ_COMPLEXE_Q
C--------------------------------------------------
C-----Nombre d'états de charge pour chaque complexe
C--------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                Q_COMPLEXE_Q
C-----------------
C-----Multiplicité
C-----------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                MULTIPLICITE_COMPLEXE_Q
C-----------------------------------------------------------
C-----Sous-réseau sur lequel est calculée cette multiplicité
C-----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                           I_S_R_MULTIPLICITE_COMPLEXE_Q
C--------------------------------
C-----Nombre de sites du complexe
C--------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NOMBRE_SITES_COMPLEXE_Q
C--------------------------------------------------------------
C-----Numéros de sous-réseaux des sites occupés par le complexe
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_S_R_COMPLEXE_Q
C------------------------------------------------------
C-----Types chimiques sur chacun des sites du complexe,
C-----dans le même ordre que les s-r (0 = lacune)
C------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_TYPE_COMPLEXE_Q
C-------------------------------------------------
C-----Energie et volume "bruts" de chaque complexe
C-----pour chaque état de charge
C-------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  E_B_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                              V_B_COMPLEXE_Q
C-------------------------------------------------------
C-----Energie, volume et enthalpie GC de chaque complexe
C-----pour chaque état de charge
C-------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  E_GC_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  V_GC_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  H_GC_COMPLEXE_Q
C----------------------------------------------------------
C-----Nombre de sites du sous-réseau r dans chaque complexe
C-----et nombre d'atomes de type i dans ce complexe
C----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  U_COMPLEXE_S_R_Q
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  V_COMPLEXE_TYPE_Q
C-------------------------------
C-----Fractions de complexes
C-----pour chaque état de charge
C-------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  X_COMPLEXE_Q
C--------------------------------------------------
C-----Enthalpies de formation des défauts complexes
C-----pour chaque état de charge
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                             H_FORM_COMPLEXE_Q
C============================================================
C============================================================
C=====Chaînes de caractères pour suffixes de noms de fichiers
C============================================================
C============================================================
C--------------------------------------------------------
C-----Chaîne de caractères pour les suffixes de DP
C-----(défaut ponctuel et sous-réseau) et leurs longueurs
C--------------------------------------------------------
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_R
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_TYP
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_R
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_TYP
C--------------------------------------------------------
C-----Chaîne de caractères pour les suffixes de complexes
C-----et leurs longueurs
C--------------------------------------------------------
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_COMPLEXE
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_COMPLEXE
C-----------------------------------------------------------------
C-----Tableau des noms de DP en fonction des types et sous-réseaux
C-----et longueurs correspondantes
C-----------------------------------------------------------------
      CHARACTER * 10 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          NOM_D_R_TYP
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          LONG_NOM_D_R_TYP
C---------------------------------------------------
C-----Noms des DP en fonction de leur indice unique,
C-----longueurs des noms correspondantes,
C-----sous-réseaux et types,
C-----énergies, volumes et enthalpies GC
C---------------------------------------------------
      CHARACTER * 10 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          NOM_D_IND
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                        LONG_NOM_D_IND
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                        I_TYP_R_D_IND
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                        E_GC_D_IND
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                        V_GC_D_IND
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                        H_GC_D_IND
C=============================================
C=============================================
C=====Tableaux relatifs à l'ADPI NPT à T = 0 K
C=============================================
C=============================================
C---------------------------------------------------------------
C-----Indice du DP en fonction de son sous-réseau et de son type
C---------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          IND_D_R_TYP 
C---------------------------------------------------------------------
C-----Tableau du simplexe contenant :
C-----(i) la fonction H à minimiser en première ligne
C-----(i) les contraintes dans les lignes 2 à N_TYP + 1
C-----(i) la fonction auxiliaire dans la ligne N_TYP + 2
C-----Ce tableau contient N_TYP_D_R + 2 colonnes :
C-----la première colonne contient les valeurs des contraintes,
C-----les autres colonnes contiennent les coefficients des contraintes
C---------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) :: TAB_SMPLX_D_R
C--------------------------------------
C-----Tableaux de résultats du simplexe
C--------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_1_SMPLX
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_2_SMPLX
C-----------------------------------------------
C-----En mode NPT=0 "balayage",
C-----tableau de résultats du simplexe précédent
C-----pour détection de limite de zone
C-----------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_2_SMPLX_PREC
C-------------------------------------------------------------
C-----Tableau des indices de la liste initiale I_2_SMPLX de DP
C-----classés par ordre croissant
C-------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: IND_INIT
C--------------------------------------------------------------
C-----Tableau des sous-réseaux et types des DP constitutionnels
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           I_TYP_R_DP_CONST
C------------------------------------------------------
C-----Matrice N_TYP x N_TYP entre les variables
C-----(N_maille, n_d_const) et ( N_I )
C-----utile au calcul des potentiels chimiques en NPT=0
C-----et inverse de cette matrice
C------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           MAT_CALC_POT
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           MAT_CALC_POT_INV
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           POT_CHIM_0K
C------------------------------------------------------------------
C-----Energies de référence des diverses espèces chimiques (eV/at.)
C------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           E_REF_TYP
C=====================================
C=====================================
C=====Tableaux relatifs à l'ADPI - NPT
C=====================================
C=====================================
C-----------------------------------------
C-----log_10 des fractions initiales de DP
C-----------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           LOG_X_D_R_INIT
C--------------------------------------------------------
C-----Vecteur, fonction vectorielle et matrice jacobienne
C-----du système non linéaire à résoudre en NPT
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           X_NPT
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           F_NPT
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           J_NPT
C----------------------------------------
C-----En NPT, vecteur initial à un indice
C----------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           X_NPT_INIT
C------------------------------------------------
C-----En NPT, ce même vecteur lu dans un fichier,
C-----pour chaque composition / température
C------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_NPT_INIT_FICH
C------------------------------------------------------------
C-----En NPT(Tx) : tableau pour la conservation des résultats 
C-----(aux diverses compo. = boucle interne)
C-----à la température courante en vue de l'initialisation
C-----à chaque compo. à la température suivante
C------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_NPT_COMPO_SVG
C###############################
C#####Tableaux de dimension fixe
C###############################
C-------------------------------------------------
C-----Chaîne de lecture dans le fichier de données
C-------------------------------------------------
      CHARACTER * 100 CHAINE_LECT
C---------------------------------------
C-----Base du nom des fichiers de sortie
C---------------------------------------
      CHARACTER * 150 FICH_FIN
C---------------------------------------------------------
C-----Indicateur de présence de sous-réseaux interstitiels
C---------------------------------------------------------
      CHARACTER * 1 INDIC_R_INTER
C-----------------------------------------------------------------
C-----Indicateur de prise en compte des interstitiels intrinsèques
C-----------------------------------------------------------------
      CHARACTER * 1 INDIC_INTER_INTR
C------------------------------------------------
C-----Indicateur de présence de défauts complexes
C------------------------------------------------
      CHARACTER * 1 INDIC_COMPLEXES
C---------------------------------------------------------
C-----Indicateur d'écriture des grandeurs thermodynamiques
C-----par atome (A/a) ou par maille (M/m)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_AT_MAILLE
C---------------------------------------------------------
C-----Indicateur d'écriture de l'énergie libre par atome :
C-----énergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_G
C------------------------------------------------------
C-----Indicateur du mode de calcul (muVT, NPT ou NPT=0)
C------------------------------------------------------
      CHARACTER * 10 INDIC_TYP_CALC
C-------------------------------------------------------------------
C-----Indicateur de type de calcul NPT(=0)
C----- * P/p : point (T,x) fixés (modes NPT et NPT=0)
C----- * T : balayage en température (mode NPT)
C----- * x : balayage en composition (modes NPT et NPT=0)
C----- * Tx et xT  : doubles balayages en température et composition
C-----  --> boucle externe sur T ou x pour Tx ou xT respectivement
C----- (mode NPT)
C-------------------------------------------------------------------
      CHARACTER * 2 INDIC_TYP_CALC_NPT
C---------------------------------------------------------------
C-----Produit courant des nombres de pas des éléments d'addition
C---------------------------------------------------------------
      INTEGER * 4 PROD_PART_COUR
C-----------------------------------------------------
C-----Caractères utiles à la constitution d'un format 
C-----pour l'écriture d'un nombre variable de colonnes
C-----(potentiel chimique de chaque espèce,
C-----fraction atomique de chaque espèce,
C-----fraction de DP sur le sous-réseau considéré
C-----et énergie de formation de ce DP)
C-----------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_X_DP
C------------------------------------------
C-----Format générique du titre associé
C-----et format particulier à chaque défaut
C------------------------------------------
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP_TYPE
C-----------------------------------------
C-----Idem pour le fichier d'énergie libre
C-----------------------------------------
      CHARACTER * 300 CAR_COL_VAR_E_L
      CHARACTER * 5000 CAR_TITRE_VAR_E_L
C-------------------------------------------------
C-----Idem pour le fichier de potentiels chimiques
C-------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_POT_CHIM
      CHARACTER * 5000 CAR_TITRE_VAR_POT_CHIM
C-------------------------------------------------
C-----Idem pour le fichier de pb de conv. de la BA
C-------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_PB_CONV_BA
      CHARACTER * 5000 CAR_TITRE_VAR_PB_CONV_BA
C--------------------------------------------------------------------
C-----Indicateur d'écriture des seuls points contenus dans la fenêtre
C--------------------------------------------------------------------
      CHARACTER * 1 INDIC_ECRIT_FENETRE
C-------------------------------------------------
C-----Ligne des fractions atomiques dans DATA.adpi
C-------------------------------------------------
      CHARACTER * 1000 LIGNE_X_AT
C--------------------------------------------------------
C-----Tableaux d'analyse de cette ligne chaîne par chaîne
C--------------------------------------------------------
      INTEGER * 4 I_DEBUT ( 1000 )
      INTEGER * 4 I_FIN ( 1000 )
      INTEGER * 4 LONG_CHAINE_X_AT ( 1000 )
      CHARACTER * 1000 CHAINE_X_AT ( 1000 )
C----------------------------------------------------
C-----Chaîne totale de suffixe de fractions atomiques
C----------------------------------------------------
      CHARACTER * 1000 CHAINE_TOTALE_X_AT
C--------------------------------
C-----Facteur thermodynamique k*T
C--------------------------------
      REAL * 8 K_T
C--------------------------------------
C-----Nombre d'atomes par maille (réel)
C--------------------------------------
      REAL * 8 N_AT_MAILLE
C---------------------------------
C-----Nombre d'atomes total (réel)
C-----(pour calcul NPT)
C---------------------------------
      REAL * 8 N_AT_TOT
C-------------------------------------
C-----Nombre de mailles initial (réel)
C-----(pour calcul NPT)
C-------------------------------------
      REAL * 8 N_MAILLES_INIT
C-------------------------------------------------
C-----Fréquence (nombre de pas) de calcul de J_NPT
C-------------------------------------------------
      INTEGER * 4 P_J_NPT
C----------------------------------------
C-----Fichier de valeurs initiales en NPT
C----------------------------------------
      CHARACTER * 200 FICH_VAL_INIT
C-----------------------------------------------------------
C-----Caractère indicateur d'écriture d'un fichier
C-----contenant les compositions où la minimisation a échoué
C-----(cas NPT=0)
C-----------------------------------------------------------
      CHARACTER * 1 CAR_FICH_550
C----------------------------------------------
C-----DP chargés : Répertoire et fichier de DdE
C----------------------------------------------
      CHARACTER * 200 REP_DDE
      CHARACTER * 100 FICH_DDE
C############
C#####Formats
C############
C--------------------
C-----Sauts de lignes
C--------------------
C      INCLUDE
C     $'format.inc'
C     ==================================================================
C     I								       I
C     I Formats de sauts de lignes				       I
C     I          						       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 21/06/2005                              I
C     I Saint Louis de Gonzague					       I
C     ==================================================================
C--------------------
C-----Sauts de lignes
C--------------------
    1 FORMAT ( / )
    2 FORMAT ( 2 / )
    3 FORMAT ( 3 / )
    4 FORMAT ( 4 / )
    5 FORMAT ( 5 / )
    6 FORMAT ( 6 / )
    7 FORMAT ( 7 / )
    8 FORMAT ( 8 / )
    9 FORMAT ( 9 / )
   10 FORMAT ( 10 / )
   11 FORMAT ( 11 / )
   12 FORMAT ( 12 / )
   13 FORMAT ( 13 / )
   14 FORMAT ( 14 / )
   15 FORMAT ( 15 / )
   16 FORMAT ( 16 / )
   17 FORMAT ( 17 / )
   18 FORMAT ( 18 / )
   19 FORMAT ( 19 / )
   20 FORMAT ( 20 / )
C---------------------
C-----Formats de réels
C---------------------
  406 FORMAT ( 6 ( 2X , G12.6 ) )
  407 FORMAT ( 7 ( 2X , G12.6 ) )
  408 FORMAT ( 8 ( 2X , G12.6 ) )
  409 FORMAT ( 9 ( 2X , G12.6 ) )
  433 FORMAT ( 3 ( 2X , G12.6 ) , 3 ( 5X , G12.6 ) )
  434 FORMAT ( 3 ( 2X , G12.6 ) , 4 ( 5X , G12.6 ) )
 1110 FORMAT ( 'x_at = ' , 100 ( 2X , G16.10 ) )
 1111 FORMAT ( 100 ( 2X , G16.10 ) )
 1112 FORMAT ( I7 , 100 ( 2X , G16.10 ) )
C--------------------------------
C-----Formats du fichier de liste
C--------------------------------
 1100 FORMAT ( 100 ( '=' ) )
 1200 FORMAT ( 100 ( '-' ) )
  501 FORMAT ( 2X , G12.6 )
  502 FORMAT ( 2X , I4 )
  503 FORMAT ( 3 ( 2X , I4 ) )
  504 FORMAT ( 4 ( 2X , I4 ) )
  505 FORMAT ( 3X , 'Valeur initiale (eV) = ' , G12.6 ,
     $         4X , 'Nombre de pas = ' , I9 ,
     $         4X , 'Incrément (eV) = ' , G12.6 )
  600 FORMAT ( "Nombre de types chimiques : " , 30X , I4 )
  601 FORMAT ( " Nombre d'éléments d'addition : " , 30X , I4 )
  602 FORMAT ( "   Nombre de sous-réseaux : " , 36X , I4 )
  603 FORMAT
     $ ( '      Nombre de sous-réseaux interstitiels :' , 26X , I4 )
  604 FORMAT
     $ ( '        Nombre de types de complexes pris en compte :' ,
     $   3X , I4 )
  705 FORMAT ( 15X , 100 ( 2X , I4 ) )
C-----Formats relatifs aux DP chargés
 2000 FORMAT
     $( 'Charges (unités |e-|) "n, A, p, D" pour le volume de SC =' ,
     $  F10.3 , ' A^3' )
 2001 FORMAT ( 5X , 'Itér.' , 4X , 'mu(e) (eV)' ,
     $ 6X , '-[Q(n)+Q(A)]' , 6X , '+Q(p)+Q(D)' )
 2002 FORMAT ( 5X , 'Itér.' , 4X , 'mu(e) (eV)' ,
     $ 12X , '-Q(n)', 12X , '-Q(A)' , 12X , '+Q(p)' , 12X , '+Q(D)' ,
     $ 12X , 'Q(totale)' , 6X , 'dQ(totale)/dmu_e' )
 3001 FORMAT ( 2X , I6 , 3 ( 2X , G16.7 ) )
 3002 FORMAT ( 2X , I6 , 7 ( 2X , G16.7 ) )
 3005 FORMAT ( '              T = ' , F10.2 , ' K' )
 3010 FORMAT ( 3 ( 5X , I3 ) , 4X , G16.7  )
 3020 FORMAT ( 2 ( 5X , I3 ) , 4X , G16.7  )
C####################################
C####################################
C#####Fin des déclarations et formats
C####################################
C####################################
C-------------------------------------------------------------------
C-----Ecriture à l'écran des références du programme avant exécution
C-------------------------------------------------------------------
        CALL REF_PROG
C       CALL REF_PROG_DIRECT
C&&&&& NOTE : la partie "lecture du fichier de données"
C&&&&& est commune aux cas sans et avec charges
C&&&&& => elle contient des tests sur la valeur de INDIC_CHARGE
C&&&&& A l'inverse, la suite du programme est scindée en deux parties,
C&&&&& la première pour DP non chargés, la seconde pour DP chargés,
C&&&&& et de tels tests y sont donc inutiles.
C##################################
C##################################
C#####Lecture du fichier de données
C##################################
C##################################
      OPEN ( 10 , FILE = 'DATA.adpi' )
C########################
C#####PARAMETRES GENERAUX
C########################
C---------------------------------
C-----Mode d'emploi et commentaire
C---------------------------------
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 11 ) .NE. 'Commentaire' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      READ ( 10 , '(A)' ) CHAINE_LECT
      CHAINE_LECT = ''
      N_LIGNES_COMMENTAIRE = - 1
      DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $          .NE.
     $'==================================================' )
        READ ( 10 , '(A)' ) CHAINE_LECT
        N_LIGNES_COMMENTAIRE = N_LIGNES_COMMENTAIRE + 1
      END DO
      REWIND ( 10 )
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 11 ) .NE. 'Commentaire' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      READ ( 10 , '(A)' ) CHAINE_LECT
      ALLOCATE ( COMMENTAIRE ( N_LIGNES_COMMENTAIRE ) )
      DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        READ ( 10 , '(A)' ) COMMENTAIRE ( I_LIGNE )
      END DO
C---------------------------------------
C-----Base du nom des fichiers de sortie
C---------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , '(A)' ) FICH_FIN
C--------------------------------------------
C-----Calcul de la vraie longueur de FICH_FIN
C--------------------------------------------
      LONG_FICH_FIN = INDEX ( FICH_FIN , ' ' ) - 1
C----------------------------------------------------
C-----Nombres de types chimiques,
C-----total et intrinsèque (hors éléments d'addition)
C----------------------------------------------------
      READ ( 10 , 6 )
      READ ( 10 , * ) N_TYP , N_TYP_INTR
C----------------------------------------------
C-----Interdiction de certains types de calculs
C----------------------------------------------
C     IF ( N_TYP .NE. N_TYP_INTR ) 
C    $ THEN
C       WRITE ( * , * ) '-----------------------------------'
C       WRITE ( * , * ) 'Type de calcul :'
C       WRITE ( * , * ) "éléments d'addition non disponibles"
C       WRITE ( * , * ) '-----------------------------------'
C       CALL INTERRUPTION
C     END IF
C     IF ( N_TYP .NE. 3 )
C    $ THEN
C       WRITE ( * , * ) '-----------------------------------'
C       WRITE ( * , * ) 'Type de calcul :'
C       WRITE ( * , * ) 'seul "3 types" possible'
C       WRITE ( * , * ) '-----------------------------------'
C       CALL INTERRUPTION
C     END IF
C----------------------------------------------------------------
C-----Nombre de sous-réseaux
C-----(entendu au sens d'ensemble d'atomes de même environnement)
C----------------------------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) N_R
C--------------------------------------------------------
C-----Présence de sous-réseaux interstitiels (O/o ou N/n)
C--------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_R_INTER
C----------------------------------------------------------------
C-----Prise en compte des interstitiels intrinsèques (O/o ou N/n)
C----------------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_INTER_INTR
      IF ( INDIC_R_INTER .EQ. 'N' .OR. INDIC_R_INTER .EQ. 'n' ) THEN
        IF ( INDIC_INTER_INTR .EQ. 'O'
     $  .OR. INDIC_INTER_INTR .EQ. 'o' )
     $  THEN
         WRITE ( * , * )
     $ '----------------------------------------------------'
         WRITE ( * , * )
     $ 'Indicateur de présence de sous-réseaux interstitiels = ' ,
     $   INDIC_R_INTER
         WRITE ( * , * )
     $ '=> ne pas sélectionner la prise en compte'
         WRITE ( * , * )
     $ 'des interstitiels intrinsèques'
         WRITE ( * , * )
     $ '----------------------------------------------------'
         CALL INTERRUPTION
        END IF
      END IF
C-----------------------------------------------------------------
C-----Nombres de sous-réseaux occupés par les espèces intrinsèques
C-----à l'état fondamental
C-----------------------------------------------------------------
      ALLOCATE ( N_R_N ( N_TYP_INTR ) )
        READ ( 10 , 4 )
        READ ( 10 , * ) 
     $ ( N_R_N ( I_TYP_INTR ) , I_TYP_INTR = 1 , N_TYP_INTR )
C-----------------------------------
C-----Calcul auxiliaire :
C-----nombres cumulés correspondants
C-----------------------------------
      ALLOCATE ( RHO_N ( 0 : N_TYP_INTR ) )
      RHO_N = 0
      RHO_N ( 0 ) = 0
      DO I_TYP = 1 , N_TYP_INTR
       DO J_TYP = 1 , I_TYP
        RHO_N ( I_TYP ) = RHO_N ( I_TYP )
     $                  + N_R_N ( J_TYP )
       END DO
      END DO
C---------------------------------------------------
C-----Calcul du nombre de sous-réseaux interstitiels
C---------------------------------------------------
      N_R_INTER = N_R
      DO I_TYP_INTR = 1 , N_TYP_INTR
        N_R_INTER = N_R_INTER - N_R_N ( I_TYP_INTR )
      END DO
      IF ( N_R_INTER .LT. 0 ) THEN
        WRITE ( * , * )
     $  '------------------------------------------------------------'
        WRITE ( * , * )
     $   'Le nombre total de sous-réseaux doit être au moins égal'
        WRITE ( * , * )
     $   'à la somme des nombres de sous-réseaux intrinsèques'
        WRITE ( * , * )
     $   '(la différence correspondant aux sous-réseaux interstitiels)'
        WRITE ( * , * )
     $   '------------------------------------------------------------'
        CALL INTERRUPTION
      END IF
      IF
     $     ( ( ( INDIC_R_INTER .EQ. 'N' .OR. INDIC_R_INTER .EQ. 'n' )
     $     .AND. N_R_INTER .GT. 0 )
     $  .OR. ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $     .AND. N_R_INTER .EQ. 0 ) )
     $ THEN
        WRITE ( * , * )
     $  '---------------------------------------------------------'
        WRITE ( * , * )
     $  'Incohérence entre le nombre de sous-réseaux interstitiels'
        WRITE ( * , * )
     $  ' = ' , N_R_INTER
        WRITE ( * , * )
     $  "et l'indicateur de présence de ces sous-réseaux"
        WRITE ( * , * )
     $  '---------------------------------------------------------'
        CALL INTERRUPTION
      END IF
C-----------------------------------------------------------------
C-----Etape auxiliaire :
C-----recopiage des nombres de sous-réseaux par espèce intrinsèque
C-----(la suite du programme ne gère pas encore N types intrinsèques)
C-----------------------------------------------------------------
      N_R_2 = 0
      N_R_3 = 0
      N_R_1 = N_R_N ( 1 )
      IF ( N_TYP_INTR .GT. 1 ) THEN
        N_R_2 = N_R_N ( 2 )
        IF ( N_TYP_INTR .GT. 2 ) THEN
          N_R_3 = N_R_N ( 3 )
        END IF
      END IF
C------------------------------------------------------------------
C-----Lecture du nombre de sites par maille pour chaque sous-réseau
C------------------------------------------------------------------
      ALLOCATE ( P_R ( N_R ) )
      READ ( 10 , 4 )
      READ ( 10 , * ) ( P_R ( I_R ) , I_R = 1 , N_R )
C-----------------------------------------------
C-----Présence de défauts complexes (O/o ou N/n)
C-----------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , '(A)' ) INDIC_COMPLEXES
C-----------------------------------------------------------
C-----Mise à zéro éventuelle du nombre de types de complexes
C-----------------------------------------------------------
       IF ( INDIC_COMPLEXES .EQ. 'N' .OR. INDIC_COMPLEXES .EQ. 'n' ) 
     $ THEN
        N_TYPES_COMPLEXES = 0
      END IF
C-----------------------------------------------------------
C-----Température prescrite (K) (non utilisée en mode NPT=0)
C-----------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) TEMPERATURE
C-----------------
C-----Facteur kB*T
C-----------------
        K_T = K_B * TEMPERATURE
C------------------------------
C-----Pression prescrite (kbar)
C------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) PRESSION
C------------------------------------------------------------------
C-----Energies de référence des diverses espèces chimiques (eV/at.)
C------------------------------------------------------------------
      ALLOCATE ( E_REF_TYP ( N_TYP ) )
      READ ( 10 , 4 )
      READ ( 10 , * ) ( E_REF_TYP ( I_TYP ) , I_TYP = 1 , N_TYP )
C----------------------------------------------------
C-----Choix d'écriture des grandeurs thermodynamiques 
C-----par atome (A/a) ou par maille (M/m)
C----------------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) INDIC_AT_MAILLE
      IF ( .NOT.
     $    ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm'
     $ .OR. INDIC_AT_MAILLE .EQ. 'A' .OR. INDIC_AT_MAILLE .EQ. 'a' ) )
     $ THEN
        WRITE ( * , * ) '---------------------------------------'
        WRITE ( * , * ) 'Ecriture des grandeurs thermodynamiques :'
        WRITE ( * , * ) 'choisir M/m (maille) ou A/a (atome)'
        WRITE ( * , * ) '---------------------------------------'
        CALL INTERRUPTION
      END IF
C---------------------------------------------------------
C-----Indicateur d'écriture de l'énergie libre par atome :
C-----énergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) INDIC_G
      IF ( .NOT.
     $    ( INDIC_G .EQ. 'T' .OR. INDIC_G .EQ. 't'
     $ .OR. INDIC_G .EQ. 'F' .OR. INDIC_G .EQ. 'f' ) )
     $ THEN
        WRITE ( * , * ) '---------------------------------------'
        WRITE ( * , * ) "Ecriture de l'énergie libre par atome  :"
        WRITE ( * , * ) 'choisir T/t (totale) ou F/f (formation)'
        WRITE ( * , * ) '---------------------------------------'
        CALL INTERRUPTION
      END IF
C----------------------------------------------------------------
C-----Le choix entre énergie libre totale et de formation
C-----n'est permis que si l'écriture par atome a été sélectionnée
C-----(sinon, l'option "totale" est seule permise)
C----------------------------------------------------------------
      IF ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm' )
     $ THEN
       IF ( INDIC_G .EQ. 'F' .OR. INDIC_G .EQ. 'f' ) THEN
        WRITE ( * , * ) 
     $ '------------------------------------------------------'
        WRITE ( * , * )
     $ 'Ecriture des grandeurs thermodynamiques par maille'
        WRITE ( * , * ) 
     $ "=> seule l'option d'énergie libre totale est autorisée"
        WRITE ( * , * ) 
     $ '------------------------------------------------------'
        CALL INTERRUPTION 
       END IF
      END IF
C----------------------------------------
C-----Type de calcul (muVT, NPT ou NPT=0)
C----------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_TYP_CALC
      LONG_INDIC_TYP_CALC = LEN_TRIM ( INDIC_TYP_CALC )
      IF ( .NOT.
     $    ( ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ .OR.
     $    ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' )  
     $ .OR.
     $    ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' ) )
     $   )
     $ THEN
        WRITE ( * , * ) '--------------------------'
        WRITE ( * , * ) 'Type de calcul :'
        WRITE ( * , * ) 'choisir muVT, NPT ou NPT=0'
        WRITE ( * , * ) '--------------------------'
        CALL INTERRUPTION
      END IF
C---------------------------------------------------   
C-----La prise en compte de complexes n'est possible
C-----que pour un calcul muVT
C---------------------------------------------------   
      IF ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )  
     $ .AND.
     $ ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) )
     $ THEN
        WRITE ( * , * )  
     $ '----------------------------------------------'
        WRITE ( * , * ) 
     $ "La prise en compte de complexes n'est possible"
        WRITE ( * , * )  
     $ 'que pour un calcul muVT'
        WRITE ( * , * )  
     $ '----------------------------------------------'
        CALL INTERRUPTION
      END IF 
C----------------------------------------------
C-----Interdiction de certains types de calculs
C-----(pour produire une version "bridée")
C----------------------------------------------
C     IF ( .NOT.
C    $ ( ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
C    $   )
C    $ THEN
C       WRITE ( * , * ) '--------------------'
C       WRITE ( * , * ) 'Type de calcul :'
C       WRITE ( * , * ) 'muVT seul disponible'
C       WRITE ( * , * ) '--------------------'
C       CALL INTERRUPTION
C     END IF
C-----------------------------------
C-----Indicateur (0/1) de DP chargés
C-----------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_CHARGE
      IF ( .NOT. ( INDIC_CHARGE .EQ. 0 .OR. INDIC_CHARGE .EQ. 1 ) ) THEN
        write (* , * ) "Vérifier INDIC_CHARGE = 0/1"
        stop
      END IF
C=======================================================
C=====Sous-section "paramètres généraux pour DP chargés"
C=====(lecture optionnelle si INDIC_CHARGE = 1)
C=======================================================
       IF ( INDIC_CHARGE .EQ. 1 ) THEN
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 13 : 57 )
     $   .NE. 'Sous-section générale relative aux DP chargés' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------------------
C-----Nombre maximum d'états de chargés considérés
C-----et type de calcul (1=balayage en mu_e, 2=neutralité électrique)
C--------------------------------------------------------------------
      READ ( 10 , 6 )
      READ ( 10 , * ) N_MAX_CHARGES , I_CALC_CHARGE
C-------------------------------------------------------
C-----Présence de défauts complexes chargés (O/o ou N/n)
C-------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , '(A)' ) INDIC_COMPLEXES_Q
C-----INDIC_COMPLEXES ne contrôle que les complexes non chargés
C-----et n'est pas utilisé en cas de DP chargés.
C-----Cependant, pour la clarté, on impose la cohérence
C-----entre INDIC_COMPLEXES et INDIC_COMPLEXES_Q
C-----(car INDIC_COMPLEXES est dans les "paramètres généraux"
C-----et est donc requis, comme INDIC_COMPLEXES_Q, en cas de DP chargés)
C----- INDIC_COMPLEXES_Q = O ==> INDIC_COMPLEXES = O
C----- INDIC_COMPLEXES_Q = N ==> INDIC_COMPLEXES = N
      IF ( INDIC_COMPLEXES_Q .EQ. 'o' .OR. INDIC_COMPLEXES_Q .EQ. 'O' )
     $ THEN
        IF ( INDIC_COMPLEXES .EQ. 'n' .OR. INDIC_COMPLEXES .EQ. 'N' )
     $  THEN
         write(*,*) "Présence de complexes chargés"
         write(*,*) "==> activer aussi INDIC_COMPLEXES (champ avant T)"
         WRITE(*,*)
         stop
        END IF
      END IF
      IF ( INDIC_COMPLEXES_Q .EQ. 'n' .OR. INDIC_COMPLEXES_Q .EQ. 'N' )
     $ THEN
        IF ( INDIC_COMPLEXES .EQ. 'o' .OR. INDIC_COMPLEXES .EQ. 'O' )
     $  THEN
         write(*,*) "Absence de complexes chargés"
         write(*,*)
     $   "==> désactiver aussi INDIC_COMPLEXES (champ avant T)"
         WRITE(*,*)
         stop
        END IF
      END IF
C-------------------------------------------------------------------
C-----Mise à zéro éventuelle du nombre de types de complexes chargés
C-------------------------------------------------------------------
       IF ( INDIC_COMPLEXES_Q .EQ. 'N' .OR. INDIC_COMPLEXES_Q .EQ. 'n' )
     $ THEN
        N_TYPES_COMPLEXES_Q = 0
      END IF
C-----------------------------------------------------
C-----Les DP chargés ne sont pris en compte qu'en muVT
C-----------------------------------------------------
      IF ( INDIC_CHARGE .EQ. 1 ) THEN
        IF
     $ ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
     $ THEN
       WRITE ( * , * ) '-------------------------'
       WRITE ( * , * ) 'Traitement des DP chargés :'
       WRITE ( * , * ) 'seulement en mode muVT'
       WRITE ( * , * ) '-------------------------'
       CALL INTERRUPTION
       END IF
      END IF
C-----------------------------------------------
C-----Energies du maximum de la bande de valence
C-----et du minimum de la bande de conduction
C-----------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) E_MAX_BV , E_MIN_BC
C--------------------------------------------------------
C-----Potentiel chimique électronique
C-----valeur initiale (mode 2) ou minimale (mode 1),
C-----puis incrément (mode 1) et valeur maximale (mode 1)
C--------------------------------------------------------
      READ ( 10 , 6 )
      READ ( 10 , * ) POT_CHIM_ELEC ,
     $                D_POT_CHIM_ELEC , POT_CHIM_ELEC_MAX
C---------------------------------
C-----Répertoire et fichier de DdE
C---------------------------------
      READ ( 10 , 4 )
      READ ( 10 , ' ( A ) ' ) REP_DDE
      LONG_REP_DDE = INDEX ( REP_DDE , ' ' ) - 1
      READ ( 10 , 4 )
      READ ( 10 , ' ( A ) ' ) FICH_DDE
      LONG_FICH_DDE = INDEX ( FICH_DDE , ' ' ) - 1
C-----------------------------------------------
C-----Précision d'arrêt et nombre maximum de pas
C-----pour l'algorithme NR "électronique"
C-----------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) EPS , N_MAX_PAS_NR_ELEC
C-----Fin de lecture optionnelle des paramètres généraux pour DP chargés
      END IF
C##########################################################
C#####Fin de la lecture de la section "PARAMETRES GENERAUX"
C##########################################################
C===================================================
C=====Ouverture préliminaire de tableaux nécessaires 
C=====à la lecture des sections suivantes
C===================================================
C---------------------------------------------------------------------
C-----Fractions atomiques dans l'alliage
C-----telles que fixées par les potentiels chimiques
C-----et fractions atomiques
C-----déduites de la valeur pour l'élément spécifié et des contraintes
C---------------------------------------------------------------------
      ALLOCATE ( X_AT ( N_TYP ) )
      ALLOCATE ( X_AT_0_CTR ( N_TYP ) )
C----------------------------------------------------------------
C-----Paramètres pertinents seulement si N_TYP > 2,
C-----et même seulement en muVT pour X_AT_INF_CTR et X_AT_SUP_CTR
C-----(tableaux ouverts dans tous les cas néanmoins)
C----------------------------------------------------------------
      ALLOCATE ( COEF_CTR ( N_TYP - 1 , N_TYP + 1 ) )
      ALLOCATE ( D_X_AT ( N_TYP ) )
      ALLOCATE ( D_X_AT_INIT ( N_TYP ) )
      ALLOCATE ( X_AT_INF_CTR ( N_TYP ) )
      ALLOCATE ( X_AT_SUP_CTR ( N_TYP ) )
C###########################################################
C#####PARAMETRES RELATIFS A UN CALCUL NPT / NPT=0
C#####Champs lus et utilisés seulement si calcul NPT / NPT=0
C###########################################################
C-----Test sur le type de calcul : autre que muVT ?
      IF
     $ ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
     $ THEN
C------------------------------------------------
C-----Recherche de l'en-tête de la grande section
C------------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 49 )
     $   .NE. 'PARAMETRES RELATIFS A UN CALCUL NPT ou NPT=0' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------------------------------------------------------------
C------------------------------------------------------------------
C-----Lecture de l'indicateur de type de calcul NPT(=0) :
C-----* P/p : point (T,x) fixés (modes NPT et NPT=0)
C-----* T : balayage en température (mode NPT)
C-----* x : balayage en composition (modes NPT et NPT=0)
C-----* Tx et xT  : doubles balayages en température et composition
C----- --> boucle externe sur T ou x pour Tx ou xT respectivement
C-----(mode NPT)
C------------------------------------------------------------------
C------------------------------------------------------------------
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 36 )
     $   .NE. 'Indicateur de type de calcul NPT(=0)' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 5 )
        READ ( 10 , '(A)' ) INDIC_TYP_CALC_NPT
        LONG_INDIC_TYP_CALC_NPT = LEN_TRIM ( INDIC_TYP_CALC_NPT )
C- - - - - - - - - - - - - - - - - - - - - - - 
C- - -Test sur la validité de cet indicateur lu
C- - - - - - - - - - - - - - - - - - - - - - - 
        IF ( .NOT.
     $    ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $ ) ) THEN
         WRITE ( * , * ) 
     $ '------------------------------------'
         WRITE ( * , * )
     $ 'Indicateur de type de calcul NPT(=0) :'
         WRITE ( * , * ) 
     $ 'P/p : point (T,x)'
         WRITE ( * , * )
     $  'T : balayage en température'
         WRITE ( * , * ) 
     $  'x : balayage en composition'
         WRITE ( * , * )
     $  'Tx ou xT : doubles balayages'
         WRITE ( * , * ) 
     $ '------------------------------------'
         CALL INTERRUPTION
        END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Fin du test sur la validité de cet indicateur lu
C- - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Avertissement concernant le mode NPT(Tx)
C     IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
C    $  .EQ. 'Tx' ) THEN
C       write(*,*) "************************************************"
C       write(*,*) "Vous avez sélectionné le mode NPT(Tx)."
C       write(*,*) "Cependant, à l'inverse de son homologue NPT(xT),"
C       write(*,*) "ce mode n'a pas été testé sur divers exemples."
C       write(*,*) "==> Il serait préférable de privilégier NPT(xT)."
C       write(*,*) "Si vous souhaitez poursuivre tout de même"
C       write(*,*) "en mode NPT(Tx), veuillez presser une touche."
C       write(*,*) "************************************************"
C       read(*,*)
C     END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Second test : en mode NPT=0, seules les options
C- - -"P/p" et "x" sont possibles.
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
       IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' )
     $ THEN
        IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $ ) THEN
         WRITE ( * , * )
         WRITE ( * , * )
     $ 'Mode NPT=0 ==> options possibles = "P/p" ou "x"'
         WRITE ( * , * )
         CALL INTERRUPTION
       END IF
      END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Troisième test : dans le cas d'un calcul NPT=0 avec N_TYP # 3
C- - -(section PARAMETRES GENERAUX), l'option "x" n'est pas possible
C- - -(seule l'option "P/p" est possible).
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
     $  .EQ. 'x' ) THEN
       IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' 
     $ .AND. N_TYP .NE. 3 ) THEN 
        WRITE ( * , * )
        WRITE ( * , * )
     $ "Système avec N_TYP # 3 types chimiques + mode NPT=0 :"
        WRITE ( * , * )
     $ 'option "x" impossible, choisir "P/p"'
        WRITE ( * , * )
        CALL INTERRUPTION
       END IF
      END IF
C----------------------------------------------------
C----------------------------------------------------
C-----Lecture optionnelle des champs relatifs aux cas
C-----impliquant un calcul en un point de composition,
C-----i.e. en NPT(T), NPT(p) et NPT=0(p)
C----- --> fractions atomiques de l'alliage
C-----dont la somme doit valoir 1
C-----------------------------------------------------
C-----------------------------------------------------
        IF ( 
     $                     INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $                .OR. INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $         .OR. (
     $                     INDIC_TYP_CALC
     $                   ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' 
     $               .AND. INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T' 
     $              ) 
     $      ) THEN
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. "Cas d'un calcul en un point de composition" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 2 )
        READ ( 10 , * )
     $ ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
C- - - - - - - - - - - - - - - - - - - - - - - -
C- - -Vérification de somme(x_at)=1
C- - -dans les cas où cette donnée est utilisée
C- - -i.e. pour NPT(T), NPT(p) et NPT=0(p)
C- - - - - - - - - - - - - - - - - - - - - - - -
         IF
     $ ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
     $ THEN
        IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $  .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $  .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $ ) THEN
        SOMME_VERIF = 0.D0
        DO I_TYP = 1 , N_TYP
          SOMME_VERIF = SOMME_VERIF + X_AT ( I_TYP )
        END DO
        IF ( DABS ( SOMME_VERIF - 1.D0 ) .GT. 1.D-12 ) THEN
          WRITE ( * , * ) 
     $ '----------------------------------------------'
          WRITE ( * , * )
     $ 'Calcul NPT(T) ou NPT(p/P) :'
          WRITE ( * , * )
     $ 'la somme des fractions atomiques doit valoir 1'
          WRITE ( * , * ) 
     $ '----------------------------------------------'
          CALL INTERRUPTION
        END IF
       END IF
      END IF
C- - - - - - - - - - - - - - - - - - - - - - - -
C- - -Fin de vérification de somme(x_at)=1
C- - -dans les cas où cette donnée est utilisée
C- - -i.e. pour NPT(T), NPT(p) et NPT=0(p)
C- - - - - - - - - - - - - - - - - - - - - - - -
C----------------------------------------------------------
C----------------------------------------------------------
C-----Fin de lecture optionnelle des champs relatifs au cas
C-----d'un calcul en un point de composition
C----------------------------------------------------------
C----------------------------------------------------------
      END IF
C-----------------------------------------------------
C-----------------------------------------------------
C-----Lecture optionnelle des champs pour le cas
C-----d'un balayage "r,theta" en composition
C-----autour de la stoechiométrie :
C-----nombre de pas angulaires entre 0 et 2 pi,
C-----"r" maximal (entre 0 et 1),
C-----nombre de pas entre 0 et r_max
C----- --> seulement en NPT=0, et pour un balayage "x"
C-----------------------------------------------------
C-----------------------------------------------------
      IF
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' )
     $ THEN
          IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
     $        .EQ. 'x' ) THEN
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. "Cas d'un balayage (r,theta) en composition" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 3 ) 
        READ ( 10 , * ) N_PAS_THETA , R_MAX , N_PAS_R
C------------------------------------------------------
C------------------------------------------------------
C-----Fin de lecture optionnelle des champs pour le cas
C-----d'un balayage "r,theta" en composition
C----- --> seulement en NPT=0, et pour un balayage "x"
C------------------------------------------------------
C------------------------------------------------------
       END IF
      END IF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C-----Lecture optionnelle des champs pour le cas
C-----d'un balayage en T ou en fraction atomique pour l'espèce choisie
C----- --> seulement en NPT, et pour un balayage "x", "T", "Tx" ou "xT"
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      IF
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
C- - - - - - - - - - - - - -
C- - -Cas d'un balayage "x"
C- - - - - - - - - - - - - -
       IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
     $ .EQ. 'x' ) THEN
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
          REWIND ( 10 )
          CHAINE_LECT = ''
          DO WHILE ( CHAINE_LECT ( 1 : 34 )
     $     .NE. "Nombre de pas en fraction atomique" )
            READ ( 10 , '(A)' ) CHAINE_LECT
          END DO
C------------
C-----Lecture
C------------
          READ ( 10 , 1 )
          READ ( 10 , * ) N_PAS_X_TYP_0
C- - - - - - - - - - - - - -
C- - -Cas d'un balayage "T"
C- - - - - - - - - - - - - -
        ELSE IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
     $  .EQ. 'T' ) THEN
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 18 )
     $   .NE. "Nombre de pas en T" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 1 )
        READ ( 10 , * ) N_PAS_T
C- - - - - - - - - - - - - - - - - -
C- - -Cas d'un balayage "Tx" ou "xT"
C- - - - - - - - - - - - - - - - - -
        ELSE IF ( INDIC_TYP_CALC_NPT
     $          ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $       .OR. INDIC_TYP_CALC_NPT
     $          ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' 
     $          )  THEN
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
          REWIND ( 10 )
          CHAINE_LECT = ''
          DO WHILE ( CHAINE_LECT ( 1 : 34 )
     $     .NE. "Nombre de pas en fraction atomique" )
            READ ( 10 , '(A)' ) CHAINE_LECT
          END DO
C------------
C-----Lecture
C------------
          READ ( 10 , 1 )
          READ ( 10 , * ) N_PAS_X_TYP_0
C----------------------------------------------
C-----Recherche de l'en-tête de la sous-section
C----------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 18 )
     $   .NE. "Nombre de pas en T" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 1 )
        READ ( 10 , * ) N_PAS_T
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C-----Fin de lecture optionnelle des champs pour le cas
C-----d'un balayage en T ou en fraction atomique pour l'espèce choisie
C----- --> seulement en NPT, et pour un balayage "x", "T", "Tx" ou "xT"
C----------------------------------------------------------------------
C----------------------------------------------------------------------
        END IF
       END IF
C-----Fin du test de type "autre que muVT"
      END IF
C##################################################
C#####Fin de la lecture de la section
C#####"PARAMETRES RELATIFS A UN CALCUL NPT / NPT=0"
C##################################################
C###########################################################
C#####PARAMETRES RELATIFS A UN CALCUL muVT ou NPT
C#####Champs lus et utilisés seulement si calcul muVT ou NPT
C###########################################################
C-----Test du type de calcul : NPT ou muVT
      IF ( .NOT.
     $    ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' ) )
     $ THEN
C============================================
C============================================
C=====Section lue seulement si DP non chargés
C============================================
C============================================
       IF ( INDIC_CHARGE .EQ. 0 ) THEN
C-------------------------------------------------------
C-----Recherche de l'en-tête de section
C-----(devenue inutile car en-têtes intérieurs utilisés) 
C-------------------------------------------------------
C       REWIND ( 10 )
C       CHAINE_LECT = ''
C       DO WHILE ( CHAINE_LECT ( 6 : 50 )
C    $   .NE. 'PARAMETRES RELATIFS A UN CALCUL muVT ou NPT' )
C         READ ( 10 , '(A)' ) CHAINE_LECT
C       END DO
C======================================================================
C=====1) Champs à inclure ssi N_TYP > 2, en NPT(x/Tx/xT) comme en muVT
C===== --> système de N_TYP - 2 contraintes
C=====y compris en muVT lorsque le filtre "fenêtres"
C=====n'est pas sélectionné 
C=====(cf. section suivante "PARAMETRES SPECIFIQUES A UN CALCUL muVT"),
C=====bien que ces champs soient alors inutilisés dans ce cas.
C======================================================================
         IF ( 
     $                       N_TYP .GT. 2 
     $   .AND. (
     $                       INDIC_TYP_CALC
     $                     ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT'
     $        .OR. (
     $                       INDIC_TYP_CALC
     $                     ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT'
     $          .AND. (
     $                       INDIC_TYP_CALC_NPT
     $                     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     ^                  .OR. INDIC_TYP_CALC_NPT
     $                     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $                  .OR. INDIC_TYP_CALC_NPT
     $                     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $                )
     $            )
     $        )
     $     ) THEN
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 44 )
     $   .NE. 'Coefficients (a,b) des N_TYP - 2 contraintes' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C---------------------------------------------------
C-----Coefficients de contraintes sur la composition
C-----sous la forme a(k,i)x(i)+b(k) = 0
C-----(seulement si N_TYP > 2)
C---------------------------------------------------
        READ ( 10 , 6 )
        DO I_TYP = 1 , N_TYP - 2
         READ ( 10 , * )
     $ ( COEF_CTR ( I_TYP , J_TYP ) , J_TYP = 1 , N_TYP + 1 )
        END DO
C----------------------------------------------------------------------
C-----La N_TYP - 1 ème contrainte est l'unité de la somme des fractions
C----------------------------------------------------------------------
       DO I_TYP = 1 , N_TYP
        COEF_CTR ( N_TYP - 1 , I_TYP ) = 1.D0
       END DO
       COEF_CTR ( N_TYP - 1 , N_TYP + 1 ) = 0.D0
C===================================================
C=====Fin des données 1) lues seulement si N_TYP > 2
C=====en NPT(x/Tx/xT) comme en muVT
C===================================================
      END IF
C======================================================================
C=====2) Champs à inclure (i) en NPT(x/Tx/xT),
C=====                    (ii) en muVT si N_TYP > 2,
C=====y compris en muVT lorsque le filtre "fenêtres"
C=====n'est pas sélectionné
C=====(cf. section suivante "PARAMETRES SPECIFIQUES A UN CALCUL muVT"),
C=====bien que ces champs soient alors inutilisés dans ce cas.
C======================================================================
         IF ( 
     $         ( 
     $                     INDIC_TYP_CALC
     $                   ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT'
     $        .AND. (
     $                     INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $                .OR. INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $                .OR. INDIC_TYP_CALC_NPT
     $                   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' 
     $              )
     $        )
     $   .OR. (
     $                    INDIC_TYP_CALC
     $                  ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT'
     $              .AND. N_TYP .GT. 2 
     $        )
     $     ) THEN
C-------------------------------------------------------
C-------------------------------------------------------
C-----Elément dont on étudie l'effet de l'enrichissement
C-----(toutes les autres fractions atomiques
C-----sont fixées dans des fenêtres à l'écriture)
C-----(toujours utilisé en NPT(x/Tx/xT),
C-----utilisé si N_TYP > 2 en muVT)
C-------------------------------------------------------
C-------------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $   .NE. "Elément dont on étudie l'effet de l'enrichissement" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 5 )
        READ ( 10 , * ) I_TYP_0
C--------------------------------------------------------------
C--------------------------------------------------------------
C-----Valeurs initiale et finale de fraction atomique écrite
C-----pour l'élément spécifié (dont on étudie l'enrichissement)
C-----(toujours utilisé en NPT(x/Tx/xT), 
C-----utilisé si N_TYP > 2 en muVT)
C--------------------------------------------------------------
C--------------------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $   .NE. "Valeurs extrémales du domaine de fraction atomique" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 5 )
        READ ( 10 , * ) X_TYP_0_INIT , X_TYP_0_FIN
C-----Test préliminaire
        IF ( X_TYP_0_INIT .LT. 0.D0
     $  .OR. X_TYP_0_FIN .LT. 0.D0 ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Valeurs initiale et finale de fraction atomique écrite'
          WRITE ( * , * )
          WRITE ( * , * )
     $ "pour l'élément spécifié (dont on étudie l'enrichissement) :"
          WRITE ( * , * )
     $ 'on doit avoir x_init > 0 et x_fin > 0'
          WRITE ( * , * )
     $ '---------------------------------------------------------'
        CALL INTERRUPTION
        END IF
C--------------------------------------------------------
C-----En plus de leur emploi en mode NPT(x/Tx/xT),
C-----X_TYP_0_INIT et X_TYP_0_FIN servent aussi à affecter
C-----X_AT_INF_CTR et X_AT_SUP_CTR, qui sont utilisés
C-----seulement en mode muVT si N_TYP > 2.
C--------------------------------------------------------
         X_AT_INF_CTR ( I_TYP_0 ) = X_TYP_0_INIT
         X_AT_SUP_CTR ( I_TYP_0 ) = X_TYP_0_FIN
C=======================================================
C=====Fin de la lecture optionnelle
C=====des données 2) toujours utilisées en NPT(x/Tx/xT),
C=====ou utilisées si N_TYP > 2 en muVT
C=======================================================
        END IF
C===================================
C===================================
C=====Fin du test "DP non chargés ?"
C===================================
C===================================
       END IF
C-----Fin du test sur le type de calcul (NPT ou muVT)
      END IF
C##################################################
C#####Fin de la lecture de la section
C#####"PARAMETRES RELATIFS A UN CALCUL muVT ou NPT"
C##################################################
C####################################################
C#####PARAMETRES SPECIFIQUES A UN CALCUL muVT
C#####Champs lus et utilisés seulement si calcul muVT
C####################################################
C-----Test sur le type de calcul : muVT ?
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ THEN  
C-------------------------------------------------------
C-----Recherche de l'en-tête de section
C-----(devenue inutile car en-têtes intérieurs utilisés) 
C-------------------------------------------------------
C       REWIND ( 10 )
C       CHAINE_LECT = ''
C       DO WHILE ( CHAINE_LECT ( 6 : 44 )
C    $   .NE. 'PARAMETRES SPECIFIQUES A UN CALCUL muVT' )
C         READ ( 10 , '(A)' ) CHAINE_LECT
C       END DO
C============================================
C============================================
C=====Premier cas de lecture : DP non chargés
C===== --> balayage en mu_i
C============================================
C============================================
       IF ( INDIC_CHARGE .EQ. 0 ) THEN
C--------------------------------------------------------
C--------------------------------------------------------
C-----Type chimique intrinsèque de référence
C-----pour les potentiels chimiques,
C-----précision et nombre maximum d'itérations
C-----pour l'arrêt de la boucle d'autocohérence sur POT_1
C--------------------------------------------------------
C--------------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. '(i) Type chimique intrinsèque de référence' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
       READ ( 10 , 3 )
       READ ( 10 , * ) I_TYP_REF_MU , PRECISION_MU_1 , N_ITER_MAX_MU_1
C-----------------
C-----Vérification
C-----------------
      IF ( N_ITER_MAX_MU_1 .LT. 1 ) THEN
       write(*,*)
     $ " *** PROGRAMME INTERROMPU ***"
       write(*,*)
     $ "Autocohérence sur mu1 -> vérifier que N_ITER_MAX_MU_1 > 0"
       write(*,*)
     $ "(N_ITER_MAX_MU_1 = 1 ==> pas d'autocohérence)"
        stop
      END IF
C----------------------------------------
C-----La suite du programme ne fonctionne
C-----que si le type de référence est 1
C----------------------------------------
      if ( I_TYP_REF_MU .ne. 1 ) then
        write ( * , * ) 'La présente version ne fonctionne'
        write ( * , * ) 'que pour type_réf.(mu) = 1'
        call interruption
      end if
C========================================
C=====Données lues seulement si N_TYP > 2
C========================================
       IF ( N_TYP .GT. 2 ) THEN
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C-----Indicateur d'écriture des seuls points contenus dans la fenêtre
C-----(indicateur lu seulement pour plus de 2 types)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 47 )
     $   .NE. "Indicateur d'écriture des seuls points contenus" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 1 )
        READ ( 10 , '(A)' ) INDIC_ECRIT_FENETRE
        IF ( INDIC_ECRIT_FENETRE .EQ. 'O'
     $  .OR. INDIC_ECRIT_FENETRE .EQ. 'o' ) THEN 
          I_ECRIT_FENETRE = 1
        ELSE IF ( INDIC_ECRIT_FENETRE .EQ. 'N'
     $       .OR. INDIC_ECRIT_FENETRE .EQ. 'n' ) THEN
          I_ECRIT_FENETRE = 0
        ELSE
          WRITE ( * , * ) 
     $          '----------------------------------------------'
          WRITE ( * , * ) 
     $          "Indicateur d'écriture des points  : O/o ou N/n" 
          WRITE ( * , * ) 
     $          '----------------------------------------------'
          CALL INTERRUPTION
        END IF
C------------------------------------------------
C------------------------------------------------
C-----Demi-largeurs des fenêtres en composition
C-----pour les éléments autres que celui spécifié
C-----(seulement si N_TYP > 2)
C------------------------------------------------
C------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. "Demi-largeurs des fenêtres en composition" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 3 )
        READ ( 10 , * )
     $ ( D_X_AT_INIT ( I_TYP ) , I_TYP = 1 , N_TYP - 1 ) 
        I = 0
        DO I_TYP = 1 , I_TYP_0 - 1
          I = I + 1
          D_X_AT ( I_TYP ) = D_X_AT_INIT ( I )
        END DO
        DO I_TYP = I_TYP_0 + 1 , N_TYP
          I = I + 1
          D_X_AT ( I_TYP ) = D_X_AT_INIT ( I )
        END DO
C----------------------------------------------------------
C-----Dans le cas où N_TYP = 2, rien n'est lu
C-----mais il faut tout de même initialiser I_ECRIT_FENETRE
C----------------------------------------------------------
      ELSE
        I_ECRIT_FENETRE = 0
C=============================================================
C=====Fin de la lecture des paramètres pertinents si N_TYP > 2
C=============================================================
      END IF
C-----------------------------------------------------------
C-----Recherche de l'en-tête de la sous-section muVT
C-----sur les balayages en mu(intrinsèques) et mu(additions)
C-----(permet de sauter la section optionnelle précédente)
C-----------------------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 15 : 54 )
     $ .NE. 'Sous-section muVT relative aux balayages' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
        READ ( 10 , '(A)' ) CHAINE_LECT
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C-----Propriétés des séries d'écarts de potentiels chimiques
C-----pour les éléments intrinsèques autres que celui de référence
C-----------------------------------------------------------------
C-----------------------------------------------------------------
      ALLOCATE ( D_POT_REF_INTR_INIT ( N_TYP_INTR ) )
      ALLOCATE ( N_D_POT_REF_INTR ( N_TYP_INTR ) )
      ALLOCATE ( PAS_D_POT_REF_INTR ( N_TYP_INTR ) )
C--------------------------------------------------
C-----Initialisations (à 1 pour les nombres de pas)
C--------------------------------------------------
      D_POT_REF_INTR_INIT = 0.D0
      N_D_POT_REF_INTR = 1
      PAS_D_POT_REF_INTR = 0.D0
C------------
C-----Lecture
C------------
      DO I_TYP = 1 , N_TYP_INTR
       IF ( I_TYP .NE. I_TYP_REF_MU ) THEN 
        READ ( 10 , 7 )
        READ ( 10 , * ) D_POT_REF_INTR_INIT ( I_TYP ) ,
     $                  N_D_POT_REF_INTR ( I_TYP ) ,
     $                  PAS_D_POT_REF_INTR ( I_TYP )
        END IF
       END DO
C-------------
C-----Ecriture
C-------------
C     DO I_TYP = 1 , N_TYP_INTR
C       write ( * , * ) D_POT_REF_INTR_INIT ( I_TYP ) ,
C    $                  N_D_POT_REF_INTR ( I_TYP ) ,
C    $                  PAS_D_POT_REF_INTR ( I_TYP )
C      END DO
C-------------------------------------------------------------------
C-----Connexion avec la suite du programme, dont la présente version
C-----ne fonctionne pas pour plus de 3 types intrinsèques
C-----(même problème que plus haut pour N_R_1, N_R_2 et N_R_3)
C-------------------------------------------------------------------
      IF ( N_TYP_INTR .GE. 2 ) THEN
        D_POT_INIT_2_1 =  D_POT_REF_INTR_INIT ( 2 )
        N_D_POT_2_1 = N_D_POT_REF_INTR ( 2 )
        PAS_D_POT_2_1 = PAS_D_POT_REF_INTR ( 2 )
        IF ( N_TYP_INTR .GT. 2 ) THEN
         D_POT_INIT_3_1 =  D_POT_REF_INTR_INIT ( 3 )
         N_D_POT_3_1 = N_D_POT_REF_INTR ( 3 )
         PAS_D_POT_3_1 = PAS_D_POT_REF_INTR ( 3 )
        ELSE
         N_D_POT_3_1 = 1
        END IF
      ELSE
        N_D_POT_2_1 = 1
        N_D_POT_3_1 = 1
      END IF
C---------------------------------------------------
C---------------------------------------------------
C-----Propriétés des séries des mu_i (additions)
C-----La boucle de N_TYP_INTR + 1 à N_TYP ci-dessous
C-----permet de ne rien lire
C-----en cas d'alliage sans addition
C---------------------------------------------------
C---------------------------------------------------
      ALLOCATE ( POT_I_INIT ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( N_POT_I ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( PAS_POT_I ( N_TYP_INTR + 1 : N_TYP ) ) 
      DO I_TYP = N_TYP_INTR + 1 , N_TYP 
        READ ( 10 , 7 )
        READ ( 10 , * ) POT_I_INIT ( I_TYP ) ,
     $                  N_POT_I ( I_TYP ) ,
     $                  PAS_POT_I ( I_TYP )
       END DO
C=========================================
C=========================================
C=====Deuxième cas de lecture : DP chargés
C===== --> pas de balayage en mu_i
C=========================================
C=========================================
      ELSE
C---------------------------------------------
C-----Recherche du début de cette sous-section
C---------------------------------------------
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 51 )
     $   .NE. 'Cas "muVT+charges" : spécification directe des mu_i' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
        DO I = 1 , 3
         READ ( 10 , * )
        END DO
C----------------------------------------------------------------
C-----mu(2) (eV) et éventuellement mu(3) (si ternaire intrinsèque)
C----------------------------------------------------------------
       IF ( N_TYP_INTR .EQ. 2 ) THEN
         READ ( 10 , * ) POT_2
       ELSE
         READ ( 10 , * ) POT_2 , POT_3
       END IF
       IF ( N_TYP_INTR .LT. N_TYP ) THEN
        READ ( 10 , * )
     $ ( POT_I_INIT ( I ) , I = N_TYP_INTR + 1 , N_TYP )
       END IF
C===============================================================
C===============================================================
C=====Fin de test pour "muVT+balayage" (en l'absence de charges)
C=====ou "muVT+point" (avec DP chargés)
C===============================================================
C===============================================================
       END IF
C-----Fin du test de type "muVT"
      END IF
C##############################################
C#####Fin de la lecture de la section
C#####"PARAMETRES SPECIFIQUES A UN CALCUL muVT"
C##############################################
C###################################################
C#####PARAMETRES SPECIFIQUES A UN CALCUL NPT
C#####Champs lus et utilisés seulement si calcul NPT
C###################################################
C-----Test "mode NPT ?" pour la lecture de cette grande section
      IF
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
C---------------------------------------------
C-----Recherche de l'en-tête de grande section
C---------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 43 )
     $   .NE. 'PARAMETRES SPECIFIQUES A UN CALCUL NPT' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-------------------------------------------------
C-------------------------------------------------
C-----Valeurs extrémales du domaine de température
C----- --> champs utilisés en NPT(T/Tx/xT)
C-------------------------------------------------
C-------------------------------------------------
        IF (          INDIC_TYP_CALC_NPT
     $              ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. "T"
     $           .OR. INDIC_TYP_CALC_NPT 
     $              ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. "Tx"
     $           .OR. INDIC_TYP_CALC_NPT
     $              ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. "xT"
     $     ) THEN
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 44 )
     $   .NE. "Valeurs extrémales du domaine de température" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 1 )
        READ ( 10 , * ) T_INIT , T_FIN
C-----Test préliminaire 1
        IF ( T_INIT .LT. 0.D0
     $  .OR. T_FIN .LT. 0.D0 ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Domaine de température en NPT(T/Tx/xT) :'
          WRITE ( * , * )
     $ ' choisir T_min > 0 et T_max > 0'
          WRITE ( * , * )
     $ '---------------------------------------------------------'
        CALL INTERRUPTION
        END IF
C-----Test préliminaire 2
        IF ( T_INIT .LT. T_FIN ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Domaine de température en NPT(T) :'
          WRITE ( * , * )
     $ ' veuillez choisir T_init > T_fin'
          WRITE ( * , * )
     $ ' (convergence très sensible aux valeurs initiales sinon).'
          WRITE ( * , * )
     $ '---------------------------------------------------------'
        CALL INTERRUPTION
        END IF
C---------------------------
C-----Fin du test "NPT(T) ?"
C---------------------------
       END IF
C==========================================
C==========================================
C=====Mode NPT : paramètres de l'algorithme
C=====toujours lus en NPT (p, x ou T)
C==========================================
C==========================================
C-------------------------------------------
C-----Recherche de l'en-tête de sous-section
C-------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 17 : 56 )
     $   .NE. "Sous-section NPT relative à l'algorithme" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-----Lecture séquentielle des champs (pas de recherche d'en-tête)
C-----dans cette sous-section
C--------------------------------------------------------------------
C-----Nombre maximal d'itérations pour l'algorithme de Newton-Raphson
C--------------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) N_ITER_MAX_NPT
C---------------------------------------------------------------------
C-----Précision requise pour l'arrêt de l'algorithme de Newton-Raphson
C-----(maximum des valeurs absolues
C-----des composantes de l'écart entre deux pas)
C---------------------------------------------------------------------
        READ ( 10 , 5 )
        READ ( 10 , * ) PRECISION_NPT
C-----------------------------------------------------------------
C-----Fréquence (nombre de pas) de calcul de la matrice jacobienne
C-----------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) P_J_NPT
C----------------------------------------
C-----Valeur du paramètre alpha dans NRCG
C----------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) ALPHA_NRCG_NPT
C----------------------------------------------------------------------
C-----NRCG :
C----- 1) indicateur à l'étape préliminaire pour "élimination de x_1<0"
C-----= 1 si réduction du pas suivant les seules composantes x_1(i)<0,
C-----= 2 si réduction de toutes les composantes du pas
C----- 2) valeur du coefficient de réduction
C----------------------------------------------------------------------
        READ ( 10 , 7 )
        READ ( 10 , * ) INDIC_TYPE_REDUC_NRCG_NPT , COEF_REDUC_NRCG_NPT
C--------------------------------------------------------------------
C-----Valeur minimale "lambda_min" pour la réduction du pas dans NRCG
C--------------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) VALEUR_LAMBDA_MIN_NRCG_NPT
C=================================
C=================================
C=====Mode NPT : valeurs initiales
C=================================
C=================================
C---------------------------------------------------------
C-----La variable "nombre de mailles" est initialisée à 1,
C-----car sa valeur ne joue pas sur la résolution
C-----(inconnues = variables intensives).
C---------------------------------------------------------
         N_MAILLES_INIT = 1.D0
C-----------------------------------------------------
C-----Ouverture du tableau des valeurs initiales de DP
C-----------------------------------------------------
      ALLOCATE ( LOG_X_D_R_INIT ( 0 : N_TYP , N_R ) )
C-------------------------------------------
C-----Recherche de l'en-tête de sous-section
C-------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 12 : 58 )
     $   .NE. 'Sous-section NPT relative aux valeurs initiales' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-----Lecture séquentielle (pas de recherche d'en-tête)
C-----pour le champ suivant, qui est lu dans tous les cas NPT,
C-----i.e. NPT(p), NPT(x/Tx/xT) et NPT(T)
        READ ( 10 , 4 )
        READ ( 10 , * ) INDIC_LECT_VAL_INIT_NPT
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Test de compatibilité : la lecture des valeurs initiales
C- - -dans un fichier n'est possible que pour NPT(x/Tx/xT)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        IF (
     $                   INDIC_LECT_VAL_INIT_NPT .EQ. 1
     $             .AND. INDIC_TYP_CALC_NPT
     $                 ( 1 : LONG_INDIC_TYP_CALC_NPT ) .NE. "x"
     $             .AND. INDIC_TYP_CALC_NPT
     $                 ( 1 : LONG_INDIC_TYP_CALC_NPT ) .NE. "Tx"
     $             .AND. INDIC_TYP_CALC_NPT
     $                 ( 1 : LONG_INDIC_TYP_CALC_NPT ) .NE. "xT"
     $    ) THEN
         write ( * , * ) "Lecture des valeurs initiales dans un fichier"
         write ( * , * ) "possible seulement en mode NPT(x/Tx/xT)"
         stop
        END IF
C----------------------------------------------------------------------
C-----Cas "NPT(xT) et indic. lect. val. init. fich. = 0" :
C-----indicateur de réinitialisation des variables NPT
C-----par les valeurs "DATA.adpi" à chaque début d'itération externe,
C-----i.e. entre les boucles externe (x) et interne (T)
C-----Valeur 1 : les valeurs "DATA.adpi" sont utilisées
C-----à chaque incrément de x, avant le début de la boucle en T,
C-----à l'intérieur de laquelle les val. init. sont affectées
C-----proche en proche "PeP" de (x_n,T_p) à (x_n,T_p+1)
C-----Valeur 0 : les valeurs "DATA.adpi" sont utilisées seulement
C-----avant le début des boucles, puis l'affectation se fait toujours
C-----en mode "PeP", même au changement de composition (x_n,T_P) ->
C----- (x_n+1,T_1) (N et P étant les nombres de pas en x et T (resp.).
C-----Remarque : ce champ n'est lu (et ne doit donc être inclus en mode
C-----NPT(xT)), que si n'est pas sélectionnée ci-dessus l'option
C-----"lect. val. init. dans fich.", cette option prenant alors
C-----le pas sur les valeurs "DATA.adpi", et étant utilisée
C-----à chaque changement de x, avant le début de la boucle sur T
C-----(c'est-à-dire avec la même fréquence que les valeurs "DATA.adpi"
C-----lorsque le champ présent vaut 1).
C----------------------------------------------------------------------
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. "xT" ) THEN
C-----L'indicateur est préalablement initialisé, pour le cas où
C-----"lect. val. init. fich." i.e. INDIC_LECT_VAL_INIT_NPT = 1
C-----(puisque dans ce cas, la lecture ci-dessous n'a pas lieu).
C-----(Cette initialisation n'est peut-être pas nécessaire,
C-----mais elle est tout de même faite "par précaution").
         IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 0 ) THEN
C-----Recherche de l'en-tête
         REWIND ( 10 )
         CHAINE_LECT = ''
         DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $    .NE. 'Cas "NPT(xT) et indic. lect. val. init. fich. = 0"' )
           READ ( 10 , '(A)' ) CHAINE_LECT
         END DO
C-----Lecture
          READ ( 10 , 16 )
          READ ( 10 , * ) INDIC_LECT_VAL_INIT_XT
C-----Si "lect. val. init. fich." i.e. INDIC_LECT_VAL_INIT_NPT = 1,
C-----l'indicateur est tout de même initialisé,
C-----puisque dans ce cas, la lecture ci-dessus n'a pas lieu.
C-----(Cette initialisation n'est peut-être pas nécessaire,
C-----mais elle est tout de même faite "par précaution").
         ELSE
           INDIC_LECT_VAL_INIT_XT = 0
         END IF
        END IF
C---------------------------------------------------------------
C---------------------------------------------------------------
C-----1) Cas de lecture dans un fichier de valeurs initiales :
C-----cette lecture est effectuée seulement en mode NPT(x/Tx/xT),
C-----en raison du test précédent (sinon, le programme s'arrête)
C---------------------------------------------------------------
C---------------------------------------------------------------
        IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
C- - -D'après la remarque ci-dessus, on se trouve forcément ici
C- - -dans l'un des cas NPT(x/Tx/xT).
C-----Recherche de l'en-tête
         REWIND ( 10 )
         CHAINE_LECT = ''
         DO WHILE ( CHAINE_LECT ( 1 : 35 )
     $    .NE. 'Nom du fichier de valeurs initiales' )
           READ ( 10 , '(A)' ) CHAINE_LECT
         END DO
C-----Lecture séquentielle (pas de recherche d'en-tête intermédiaire)
C-----du nom de ce fichier et de son nombre de lignes
         READ ( 10 , '(A)' ) CHAINE_LECT
         READ ( 10 , '(A)' ) FICH_VAL_INIT
          LONG_FICH_VAL_INIT
     $  = INDEX ( FICH_VAL_INIT , ' ' ) - 1
         READ ( 10 , 4 )
         READ ( 10 , * ) N_LIGNES_FICH_VAL_INIT
C- - -Test sur la cohérence entre la longueur du fichier
C- - -et le nombre de pas en composition :
C- - - le nombre de pas et le nombre de lignes
C- - -doivent être multiples ou sous-multiples l'un de l'autre.
C- - -D'après la remarque ci-dessus, on se trouve forcément ici
C- - -dans l'un des cas NPT(x/Tx/xT)
C- - - ==> c'est N_PAS_X_TYP_0 lu qui est utilisé
          IF ( MOD ( N_LIGNES_FICH_VAL_INIT ,
     $               N_PAS_X_TYP_0 ) .NE. 0
     $   .AND. MOD ( N_PAS_X_TYP_0 ,
     $               N_LIGNES_FICH_VAL_INIT ) .NE. 0 ) THEN
         write(*,*)
     $   "Calcul NPT avec valeurs initiales lues dans un fichier :"
         write(*,*)
     $  "   -> veuillez choisir N_pas et N_lignes tels que"
         write(*,*)
     $  "      mod(N_pas;N_lignes) ou mod(N_lignes;N_pas) = 0"
         write(*,*)
         stop
         END IF
C-------------------------------------------------------------
C-------------------------------------------------------------
C-----2) Cas de lecture des valeurs initiales dans DATA.adpi :
C-----a lieu aussi bien en NPT(p) que NPT(x) ou NPT(T)
C-------------------------------------------------------------
C-------------------------------------------------------------
      ELSE
C-----Recherche de l'en-tête
         REWIND ( 10 )
         CHAINE_LECT = ''
         DO WHILE ( CHAINE_LECT ( 1 : 30 )
     $     .NE. 'log_10 des quantités initiales' )
          READ ( 10 , '(A)' ) CHAINE_LECT
         END DO
C--------------------------------------------------------------
C-----"log_10 des quantités initiales", partie 1/3 :
C-----lecture des quantités initiales d'antisites et de lacunes
C--------------------------------------------------------------
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Lecture "synthétique" (adaptée à N types intrinsèques)
C- - -des valeurs initiales pour les antisites
C- - -des diverses espèces intrinsèques sur les divers sous-réseaux
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( LOG_X_D_R_INIT ( J_TYP , I_R ) ,
     $     I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP ) )
C-----Test pour éviter des valeurs initiales aberrantes
C-----(induit erreur dans l'algorithme de résolution)
           do I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP )
            if ( LOG_X_D_R_INIT ( J_TYP , I_R ) .GE. 0.D0 ) then
             write(*,*)
             write(*,*)
     $      "Valeurs initiales en NPT : choisir Log(x_DP) < 0"
             write(*,*)
             stop
            end if
           end do
          END IF
        END DO
       END DO
C- - - - - - - - - - - - - - - - - - - - - - -
C- - -Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C- - - - - - - - - - - - - - - - - - - - - - -
      READ ( 10 , 2 )
      READ ( 10 , * )
     $ ( LOG_X_D_R_INIT ( 0 , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
         do I_R = 1 , RHO_N ( N_TYP_INTR )
          if ( LOG_X_D_R_INIT ( 0 , I_R ) .GE. 0.D0 ) then
            write(*,*)
            write(*,*)
     $     "Valeurs initiales en NPT : choisir Log(x_DP) < 0"
            write(*,*)
            stop
           end if
          end do
C-------------------------------------------------
C-----"log_10 des quantités initiales", partie 2/3 :
C-----lecture optionnelle des  quantités initiales
C-----relatives aux interstitiels intrinsèques
C-------------------------------------------------
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( LOG_X_D_R_INIT ( I_TYP_INTR , I_R ) ,
     $   I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
          do I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R 
           if ( LOG_X_D_R_INIT ( I_TYP_INTR , I_R ) .GE. 0.D0 ) then
             write(*,*)
             write(*,*)
     $      "Valeurs initiales en NPT : choisir Log(x_DP) < 0"
             write(*,*)
             stop
            end if
           end do
       END DO
      END IF
C----------------------------------------------------------------------
C-----"log_10 des quantités initiales", partie 3/3 :
C-----lecture optionnelle des quantités initiales d'éléments d'addition
C-----(substitutionnels et interstitiels)
C----------------------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         READ ( 10 , * )
     $   ( LOG_X_D_R_INIT ( I_TYP , I_R ) ,
     $     I_R = 1 , RHO_N ( N_TYP_INTR ) )
          do I_R = 1 , RHO_N ( N_TYP_INTR )
           if ( LOG_X_D_R_INIT ( I_TYP , I_R ) .GE. 0.D0 ) then
             write(*,*)
             write(*,*)
     $      "Valeurs initiales en NPT : choisir Log(x_DP) < 0"
             write(*,*)
             stop
            end if
           end do
C- - - - - - - - - - - - - - - - - - - - - - - -
C- - -Cas de présence d'interstitiels d'addition
C- - - - - - - - - - - - - - - - - - - - - - - -
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( LOG_X_D_R_INIT ( I_TYP , I_R ) ,
     $     I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
          do I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R 
           if ( LOG_X_D_R_INIT ( I_TYP , I_R ) .GE. 0.D0 ) then
             write(*,*)
             write(*,*)
     $      "Valeurs initiales en NPT : choisir Log(x_DP) < 0"
             write(*,*)
             stop
            end if
           end do
         END IF
      END DO
C-----------------------------------------
C-----Lecture du nombre de mailles initial
C-----------------------------------------
C     CHAINE_LECT = ''
C     DO WHILE ( CHAINE_LECT ( 1 : 25 )
C    $     .NE. 'Nombre de mailles initial' )
C       READ ( 10 , '(A)' ) CHAINE_LECT
C     END DO
C     READ ( 10 , '(A)' ) CHAINE_LECT
C     READ ( 10 , * ) N_MAILLES_INIT
C--------------------------------------------
C--------------------------------------------
C-----Fin du test sur INDIC_LECT_VAL_INIT_NPT
C--------------------------------------------
C--------------------------------------------
       END IF
C=====================================
C=====================================
C=====Fin de la lecture optionnelle
C=====des valeurs initiales (mode NPT)
C=====================================
C=====================================
C-----"Fin du test "mode NPT ?" pour la lecture de cette grande section
      END IF
C#############################################
C#####Fin de la lecture de la section
C#####"PARAMETRES SPECIFIQUES A UN CALCUL NPT"
C#############################################
C#######################################################
C#####PARAMETRES RELATIFS A LA SUPERCELLULE DE REFERENCE
C#######################################################
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 47 )
     $ .NE. 'PARAMETRES DE LA SUPERCELLULE DE REFERENCE' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------
C-----Energie de référence de la cellule sans défaut (eV)
C--------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) E_REF
C------------------------------------------------------------
C-----Volume de référence de la cellule sans défaut (A * * 3)
C------------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) V_REF
C---------------------------------------------------
C-----Nombre de mailles contenues dans cette cellule
C---------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) N_MAILLE_REF
C#########################################################
C#####Fin de lecture de la section
C#####"PARAMETRES RELATIFS A LA SUPERCELLULE DE REFERENCE"
C#########################################################
C######################################################
C#####     Lecture des paramètres de SC des DP
C======================================================
C#####CAS DE DP NON CHARGES : DP simples (et complexes)
C======================================================
C######################################################
C-----Test "DP non chargés ?"
      IF ( INDIC_CHARGE .EQ. 0 ) THEN
C###################################################
C#####PARAMETRES RELATIFS AUX DP SIMPLES NON CHARGES
C###################################################
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 59 )
     $ .NE. 'PARAMETRES DES SUPERCELLULES DE DP SIMPLES NON CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-------------------------------------------------------
C-----Ouverture des tableaux
C-----d'énergies "brutes" des cellules avec défauts (eV)
C-----et de même taille que la cellule de référence
C-----et énergies GC correspondantes
C-------------------------------------------------------
      ALLOCATE ( E_B_D_R ( 0 : N_TYP , N_R ) )
      ALLOCATE ( E_GC_D_R ( 0 : N_TYP , N_R ) )
C-------------------------------------
C-----Mêmes quantités pour les volumes
C-------------------------------------
      ALLOCATE ( V_B_D_R ( 0 : N_TYP , N_R ) )
      ALLOCATE ( V_GC_D_R ( 0 : N_TYP , N_R ) )
C----------------------------------------------------------
C-----Enthalpies GC des défauts sur les divers sous-réseaux
C----------------------------------------------------------
      ALLOCATE ( H_GC_D_R ( 0 : N_TYP , N_R ) )
C---------------------------------------------------------------
C-----Enthalpies de formation des DP sur les divers sous-réseaux
C---------------------------------------------------------------
      ALLOCATE ( H_FORM_D_R ( 0 : N_TYP , N_R ) )
C--------------------------------------
C-----Termes complémentaires d'entropie
C--------------------------------------
      ALLOCATE ( Z_TYP_R ( N_R ) )
C-----------------------------------------------------------
C-----Initialisation des énergies et volumes de supercellule
C-----------------------------------------------------------
      E_B_D_R = 0.D0
      V_B_D_R = 0.D0
C===================================================
C=====Lecture des énergies des cellules avec défauts
C===================================================
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 17 ) .NE. 'Energies "brutes"' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      DO I = 1 , 2
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
C------------------------------------------------------------
C------------------------------------------------------------
C-----Partie 1/3 : 
C-----lecture des énergies "brutes" d'antisites et de lacunes
C------------------------------------------------------------
C------------------------------------------------------------
C---------------------------------------------------------------
C-----Lecture "synthétique" (adaptée à N types intrinsèques)
C-----des énergies d'antisites des diverses espèces intrinsèques
C----- sur les divers sous-réseaux
C---------------------------------------------------------------
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( E_B_D_R ( J_TYP , I_R ) ,
     $     I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP ) )
          END IF
        END DO
       END DO
C-------------
C-----Ecriture
C-------------
C      DO I_TYP = 1 , N_TYP_INTR
C        DO J_TYP = 1 , N_TYP_INTR
C         IF ( J_TYP .NE. I_TYP ) THEN
C           write ( * , * )
C    $   ( E_B_D_R ( J_TYP , I_R ) ,
C    $     I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP ) )
C         END IF
C       END DO
C      END DO
C---------------------------------------------
C-----Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C---------------------------------------------
      READ ( 10 , 2 )
      READ ( 10 , * )
     $ ( E_B_D_R ( 0 , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C------------------------------------------------------------------
C------------------------------------------------------------------
C-----Partie 2/3 :
C-----lecture optionnelle des énergies d'interstitiels intrinsèques
C------------------------------------------------------------------
C------------------------------------------------------------------
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR 
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( E_B_D_R ( I_TYP_INTR , I_R ) ,
     $   I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
       END DO
      END IF
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----Partie 3/3 :
C-----lecture optionnelle des énergies d'éléments d'addition
C-----(substitutionnels et interstitiels)
C-----------------------------------------------------------
C-----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         READ ( 10 , * )
     $   ( E_B_D_R ( I_TYP , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C-----------------------------------------------
C-----Cas de présence d'interstitiels d'addition
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( E_B_D_R ( I_TYP , I_R ) ,
     $     I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
        END IF
      END DO
C==================================================
C=====Lecture des volumes des cellules avec défauts
C==================================================
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 15 ) .NE. 'Volumes "bruts"' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      DO I = 1 , 2
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
C----------------------------------------------------------
C----------------------------------------------------------
C-----Partie 1/3 :
C-----lecture des volumes "bruts" d'antisites et de lacunes
C----------------------------------------------------------
C----------------------------------------------------------
C--------------------------------------------------------------
C-----Lecture "synthétique" (adaptée à N types intrinsèques)
C-----des volumes d'antisites des diverses espèces intrinsèques
C----- sur les divers sous-réseaux
C--------------------------------------------------------------
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( V_B_D_R ( J_TYP , I_R ) ,
     $     I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP ) )
          END IF
        END DO
       END DO
C---------------------------------------------
C-----Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C---------------------------------------------
      READ ( 10 , 2 )
      READ ( 10 , * )
     $ ( V_B_D_R ( 0 , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C-----------------------------------------------------------------
C-----------------------------------------------------------------
C-----Partie 2/3 :
C-----lecture optionnelle des volumes d'interstitiels intrinsèques
C-----------------------------------------------------------------
C-----------------------------------------------------------------
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( V_B_D_R ( I_TYP_INTR , I_R ) ,
     $   I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
       END DO
      END IF
C----------------------------------------------------------
C----------------------------------------------------------
C-----Partie 3/3 :
C-----Lecture optionnelle des volumes d'éléments d'addition
C-----(substitutionnels et interstitiels)
C----------------------------------------------------------
C----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         READ ( 10 , * )
     $   ( V_B_D_R ( I_TYP , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C-----------------------------------------------
C-----Cas de présence d'interstitiels d'addition
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( V_B_D_R ( I_TYP , I_R ) ,
     $     I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
         END IF
      END DO
C#####################################################
C#####Fin de lecture optionnelle de la section
C#####"PARAMETRES RELATIFS AUX DP SIMPLES NON CHARGES"
C#####################################################
C###########################################################
C#####PARAMETRES RELATIFS AUX DP COMPLEXES NON CHARGES
C#####(lecture optionnelle ou simple allocation pour G_ADPI)
C###########################################################
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' )
     $ THEN
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 46 )
     $ .NE. 'PARAMETRES RELATIFS AUX DEFAUTS COMPLEXES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------------------------------------
C-----Nombre maximum de sites des complexes
C------------------------------------------
       READ ( 10 , 4 )
       READ ( 10 , * ) N_SITES_COMPLEXES_MAX
C--------------------------------------------------------
C-----Nombre de types de défauts complexes pris en compte
C--------------------------------------------------------
       READ ( 10 , 4 )
       READ ( 10 , * ) N_TYPES_COMPLEXES
C--------------------------------------------------
C-----Ouverture des tableaux relatifs aux complexes
C-----dont les valeurs vont être lues
C--------------------------------------------------
       ALLOCATE ( MULTIPLICITE_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( I_S_R_MULTIPLICITE_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( NOMBRE_SITES_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( I_S_R_COMPLEXE
     $          ( N_TYPES_COMPLEXES , N_SITES_COMPLEXES_MAX ) )
       ALLOCATE ( I_TYPE_COMPLEXE
     $          ( N_TYPES_COMPLEXES , N_SITES_COMPLEXES_MAX ) )
       ALLOCATE ( E_B_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( V_B_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
C-----------------------------------
C-----Positionnement dans le fichier
C-----------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 35 )
     $ .NE. 'Caractéristiques des complexes' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------
C-----Boucle sur les types de complexes 
C--------------------------------------
       DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
        READ ( 10 , 5 )
        READ ( 10 , * ) I_COMPLEXE_VERIF ,
     $                  MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ,
     $                  I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ,
     $                  NOMBRE_SITES_COMPLEXE ( I_COMPLEXE )
        IF ( I_COMPLEXE_VERIF .NE. I_COMPLEXE ) THEN
         WRITE ( * , * ) '--------------------------------'
         WRITE ( * , * ) 'Numéroter les complexes de 1 à N'
         WRITE ( * , * ) '--------------------------------'
         CALL INTERRUPTION
        END IF
        READ ( 10 , 2 )
        READ ( 10 , * )
     $ ( I_S_R_COMPLEXE ( I_COMPLEXE , I_SITE ) , 
     $   I_SITE = 1 , NOMBRE_SITES_COMPLEXE ( I_COMPLEXE ) )
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( I_TYPE_COMPLEXE ( I_COMPLEXE , I_SITE ) ,
     $   I_SITE = 1 , NOMBRE_SITES_COMPLEXE ( I_COMPLEXE ) )           
         READ ( 10 , 4 )
         READ ( 10 , * )
     $   E_B_D_COMPLEXE ( I_COMPLEXE ) ,
     $   V_B_D_COMPLEXE ( I_COMPLEXE )
       END DO
C--------------------------------------------------------------
C-----Fin de la lecture optionnelle des paramètres de complexes
C--------------------------------------------------------------
      ELSE
C----------------------------------------------------------------
C-----Sinon, allocation de tableaux passés en paramètres à G_ADPI
C----------------------------------------------------------------
        N_TYPES_COMPLEXES = 0
        ALLOCATE ( MULTIPLICITE_COMPLEXE ( N_TYPES_COMPLEXES ) )
        ALLOCATE ( I_S_R_MULTIPLICITE_COMPLEXE ( N_TYPES_COMPLEXES ) )
       END IF
C#######################################################
C#####Fin de lecture optionnelle de la section
C#####"PARAMETRES RELATIFS AUX DP COMPLEXES NON CHARGES"
C#######################################################
C##################################################
C#####     Lecture des paramètres de SC des DP
C==================================================
C#####CAS DE DP CHARGES : DP simples (et complexes)
C==================================================
C##################################################
      ELSE
C=========================================
C=========================================
C=====Première partie : DP simples chargés
C=========================================
C=========================================
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 47 )
     $   .NE. 'PARAMETRES RELATIFS AUX DP SIMPLES CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------------
C-----Recherche de la section des paramètres (GC) de DP chargés
C--------------------------------------------------------------
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 21 )
     $  .NE. 'Pour chaque DP chargé' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      DO I = 1 , 3
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO    
C===========================================
C=====Ouverture des tableaux pour DP chargés
C===========================================
C-------------------------------------------------------
C-----Tableau des nombres d'états de charges pour les DP
C-------------------------------------------------------
      ALLOCATE ( NQ_D_R_Q ( 0 : N_TYP , N_R ) )
C-------------------
C-----Charges des DP
C-------------------
      ALLOCATE ( Q_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C--------------------------------------------------
C-----Energies "brutes" des cellules DP chargés
C-----et de même taille que la cellule de référence
C-----et énergies GC correspondantes
C--------------------------------------------------
      ALLOCATE ( E_B_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
      ALLOCATE ( E_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C-------------------------------------
C-----Mêmes quantités pour les volumes
C-------------------------------------
      ALLOCATE ( V_B_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
      ALLOCATE ( V_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C----------------------------------------------------------
C-----Enthalpies GC des défauts sur les divers sous-réseaux
C----------------------------------------------------------
      ALLOCATE ( H_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C---------------------------------------------------------------
C-----Enthalpies de formation des DP sur les divers sous-réseaux
C---------------------------------------------------------------
      ALLOCATE ( H_FORM_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C-------------------------------------------------------
C-----Initialisation à 0 des paramètres de DP chargés
C-----(inutile pour le calcul, mais permet d'identifier
C-----les triplets (i_q,i_typ,i_r) non affectés à un DP)
C-------------------------------------------------------
      NQ_D_R_Q = 0
      E_B_D_R_Q = 0.D0
      V_B_D_R_Q = 0.D0
      Q_B_D_R_Q = 0.D0
C=====================================================
C=====Lecture des paramètres des cellules avec défauts
C=====================================================
C------------------------------------------------------------
C------------------------------------------------------------
C-----Partie 1/3 : 
C-----lecture des paramètres d'antisites et de lacunes
C------------------------------------------------------------
C------------------------------------------------------------
C---------------------------------------------------------------
C-----Lecture "synthétique" (adaptée à N types intrinsèques)
C-----des énergies d'antisites des diverses espèces intrinsèques
C----- sur les divers sous-réseaux
C---------------------------------------------------------------
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           READ ( 10 , 2 )
           DO I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP )
            READ ( 10 , * ) NQ_D_R_Q ( J_TYP , I_R )
C            write ( * , * ) J_TYP , I_R , NQ_D_R_Q ( J_TYP , I_R )
            DO I_Q = 1 , NQ_D_R_Q ( J_TYP , I_R )
           READ ( 10 , * )
     $   Q_D_R_Q ( I_Q , J_TYP , I_R ) ,
     $   E_B_D_R_Q ( I_Q , J_TYP , I_R ) ,
     $   V_B_D_R_Q ( I_Q , J_TYP , I_R )
          END DO
         END DO
           END IF
          END DO
      END DO
C---------------------------------------------
C-----Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C---------------------------------------------
      READ ( 10 , 2 )
      DO I_R = 1 , RHO_N ( N_TYP_INTR )
         READ ( 10 , * ) NQ_D_R_Q ( 0 , I_R )
         DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R )
          READ ( 10 , * )
     $   Q_D_R_Q ( I_Q , 0 , I_R ) ,
     $   E_B_D_R_Q ( I_Q , 0 , I_R ) ,
     $   V_B_D_R_Q ( I_Q , 0 , I_R )
         END DO
      END DO
C------------------------------------------------------------------
C------------------------------------------------------------------
C-----Partie 2/3 :
C-----lecture optionnelle des énergies d'interstitiels intrinsèques
C------------------------------------------------------------------
C------------------------------------------------------------------
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR
         READ ( 10 , 2 )
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R
          READ ( 10 , * ) NQ_D_R_Q ( I_TYP_INTR , I_R )
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP_INTR , I_R )
           READ ( 10 , * )
     $   Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) ,
     $   E_B_D_R_Q ( I_Q , I_TYP_INTR , I_R ) ,
     $   V_B_D_R_Q ( I_Q , I_TYP_INTR , I_R )
         END DO
        END DO
       END DO
      END IF
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----Partie 3/3 :
C-----lecture optionnelle des énergies d'éléments d'addition
C-----(substitutionnels et interstitiels)
C-----------------------------------------------------------
C-----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         DO I_R = 1 , RHO_N ( N_TYP_INTR )
          READ ( 10 , * ) NQ_D_R_Q ( I_TYP , I_R )
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           READ ( 10 , * )
     $     Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $     E_B_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $     V_B_D_R_Q ( I_Q , I_TYP , I_R )
          END DO
         END DO
C-----------------------------------------------
C-----Cas de présence d'interstitiels d'addition
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
           READ ( 10 , 2 )
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R
          READ ( 10 , * ) NQ_D_R_Q ( I_TYP , I_R )
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           READ ( 10 , * )
     $     Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $     E_B_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $     V_B_D_R_Q ( I_Q , I_TYP , I_R )
           END DO
          END DO
        END IF
C---------------------------------------------
C---------------------------------------------
C-----Fin de boucle sur les types (partie 3/3)
C---------------------------------------------
C---------------------------------------------
      END DO
C===========================================
C===========================================
C=====Deuxième partie : DP complexes chargés
C=====(optionnelle)
C===========================================
C===========================================
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 54 )
     $ .NE. 'PARAMETRES RELATIFS AUX DEFAUTS COMPLEXES CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------
C-----Nombre maximum de sites des complexes chargés
C--------------------------------------------------
       READ ( 10 , 5 )
       READ ( 10 , * ) N_MAX_SITES_COMPLEXES_Q
C----------------------------------------------------------------
C-----Nombre de types de défauts complexes chargés pris en compte
C----------------------------------------------------------------
       READ ( 10 , 4 )
       READ ( 10 , * ) N_TYPES_COMPLEXES_Q
C----------------------------------------------------------
C-----Ouverture des tableaux relatifs aux complexes chargés
C-----dont les valeurs vont être lues
C----------------------------------------------------------
       ALLOCATE ( NQ_COMPLEXE_Q ( N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( Q_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                            N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( MULTIPLICITE_COMPLEXE_Q ( N_TYPES_COMPLEXES_Q ) )
       ALLOCATE
     $ ( I_S_R_MULTIPLICITE_COMPLEXE_Q ( N_TYPES_COMPLEXES_Q ) )        
       ALLOCATE ( NOMBRE_SITES_COMPLEXE_Q ( N_TYPES_COMPLEXES_Q ) )     
       ALLOCATE ( I_S_R_COMPLEXE_Q
     $          ( N_TYPES_COMPLEXES_Q , N_MAX_SITES_COMPLEXES_Q ) )
       ALLOCATE ( I_TYPE_COMPLEXE_Q
     $          ( N_TYPES_COMPLEXES_Q , N_MAX_SITES_COMPLEXES_Q ) )
       ALLOCATE ( E_B_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( V_B_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( E_GC_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( V_GC_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( H_GC_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( H_FORM_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                               N_TYPES_COMPLEXES_Q ) )
       ALLOCATE ( X_COMPLEXE_Q ( N_MAX_CHARGES ,
     $                           N_TYPES_COMPLEXES_Q ) )
C-----------------------------------
C-----Positionnement dans le fichier
C-----------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 17 : 54 )
     $ .NE. 'Caractéristiques des complexes chargés' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C----------------------------------------------
C-----Boucle sur les types de complexes chargés
C----------------------------------------------
       DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
         READ ( 10 , 5 )
         READ ( 10 , * ) I_COMPLEXE_Q_VERIF ,
     $                   MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) ,
     $           I_S_R_MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) ,       
     $                   NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q )
         write ( * , * ) I_COMPLEXE_Q_VERIF ,
     $                   MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) ,
     $           I_S_R_MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) ,       
     $                   NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q )
         IF ( I_COMPLEXE_Q_VERIF .NE. I_COMPLEXE_Q ) THEN
          WRITE ( * , * ) '----------------------------------------'
          WRITE ( * , * ) 'Numéroter les complexes chargés de 1 à N'
          WRITE ( * , * ) '----------------------------------------'
          CALL INTERRUPTION
         END IF
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( I_S_R_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE ) ,
     $   I_SITE = 1 , NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q ) )  
         READ ( 10 , 2 )
         READ ( 10 , * )
     $ ( I_TYPE_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE ) ,
     $   I_SITE = 1 , NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q ) )   
         READ ( 10 , 4 )
C-----Lecture des énergies et volumes "bruts" pour chaque état de charge
        READ ( 10 , * ) NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
        DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
         READ ( 10 , * )
     $   Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) ,
     $   E_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) ,
     $   V_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
        END DO
C-----Fin de boucle sur complexes chargés
       END DO   
C----------------------------------------------------------------------
C-----Fin de la lecture optionnelle des paramètres de complexes chargés
C----------------------------------------------------------------------
       END IF
C-----Fin du test "DP non chargés ?"
C-----(conditionne le choix des blocs de données lus)
      END IF
C##################################################
C##### Fin de lecture des paramètres de SC des DP
C==================================================
C#####     CAS DE DP NON CHARGES OU CHARGES
C==================================================
C#####        DP simples (et complexes)
C##################################################
C############################################
C############################################
C#####Fin de la lecture du fichier de données
C############################################
C############################################
C##############################################################
C&&&&& A partir d'ici, le programme est constitué  de la mise 
C&&&&& bout à bout des deux cas (i) sans et (ii) avec charges,
C&&&&& avec simple "branchement" à l'un ou l'autre cas.
C&&&&& Bien que cela corresponde à des doublons d'instructions
C&&&&& (e.g pour ALPHA_TYPE_COMPLEXE ou I_TYPE_NORMAL_S_R),
C&&&&& cela permet aussi les doublons d'allocations
C&&&&& (ce qui serait source d'erreur si ces allocations
C&&&&& étaient réellement effectuées deux fois).
C##############################################################
C&&&&& Une amélioration ultérieure pourrait consister
C&&&&& à "factoriser"' ces doublons.
C##############################################################
C&&&&& Pour l'instant, seule la lecture du fichier de données
C&&&&& est commune aux deux cas.
C##############################################################
      IF ( INDIC_CHARGE .EQ. 0 ) THEN
       WRITE ( * , * )
       WRITE ( * , * )
     $ '                         ##########################'
       WRITE ( * , * )
     $ '                         ##### DP NON CHARGES #####'
       WRITE ( * , * )
     $ '                         ##########################'
       WRITE ( * , * )
      ELSE
       WRITE ( * , * )
       WRITE ( * , * )
     $ '                         ######################'
       WRITE ( * , * )
     $ '                         ##### DP CHARGES #####'
       WRITE ( * , * )
     $ '                         ######################'
       WRITE ( * , * )
C-----Branchement à la section "DP chargés"
      GOTO 3333
      END IF
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&& Cas "DP NON CHARGES"
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&& 
C------------------------
C-----Ecritures à l'écran
C------------------------
      WRITE ( * , 1100 )
      WRITE ( * , 600 ) N_TYP
      WRITE ( * , 601 ) N_TYP - N_TYP_INTR
      WRITE ( * , 602 ) N_R
      WRITE ( * , 603 ) N_R_INTER
      WRITE ( * , * )
     $ 'Nombre de sites par maille pour chaque sous-réseau : '
      WRITE ( * , * )
     $ ( P_R ( I_R ) , I_R = 1 , N_R )
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' )
     $ THEN
       WRITE ( * , * )
     $ '--------------------------------------------------------------'
        WRITE ( * , * )
     $ 'Option "défauts complexes" sélectionnée'
       WRITE ( * , * )
     $ '- - - - - - - - - - - - - - - - - - - -'
       WRITE ( * , * )
     $ '   * La présence de DP complexes peut induire'
       WRITE ( * , * )
     $ '     certaines fractions atomiques irréalistes (> 1).'
       WRITE ( * , * )
     $ '     Ceci provient de taux de DP complexes > 1'
      WRITE ( * , * )
     $ "     pour certains potentiels chimiques, mais n'empêche pas"
      WRITE ( * , * )
     $ "     l'utilisation pratique de l'ADPI+complexes"
      WRITE ( * , * )
     $ "     qui reste valable autour de la stoechiométrie."
       WRITE ( * , * )
     $ '   * Un comportement irréaliste (xDPcmplx, Gat(x))'
       WRITE ( * , * )
     $ "     peut aussi survenir pour des xi compatibles avec l'ADPI,"
       WRITE ( * , * )
     $ "     e.g. pour xB = 0 dans FeAl-B avec complexes incluant B." 
       WRITE ( * , * )
     $ "     Ces valeurs xB = 0 sont invalides,"
       WRITE ( * , * )
     $ "     car elles proviennent d'un forçage à 0 par le programme"
       WRITE ( * , * )
     $ "     de fractions atomiques trouvées initialement < 0,"
       WRITE ( * , * )
     $ "     à cause du terme X (HDR, éq. 2.133) de dénominateur de xi."
       WRITE ( * , * )
     $ "     Là encore, l'utilisation des résultats reste possible,"
       WRITE ( * , * )
     $ "     moyennant l'élimination de ces points xB = 0 artificiels"
       WRITE ( * , * )
     $ "     e.g. via un filtre xB > 10e-10 au lieu de 0."
        WRITE ( * , * )
     $ '- - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - -'
        WRITE ( * , 604 ) N_TYPES_COMPLEXES
        WRITE ( * , * )
     $ '- - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - -'
       ELSE
        WRITE ( * , * )
     $ '- - - - - - - - - - - - - - - - - - - -'
        WRITE ( * , * )
     $ 'Pas de défauts complexes pris en compte'
        WRITE ( * , 1100 )
       END IF
C-------------------------------
C-----Format variable d'écriture
C-------------------------------
         CALL
     $   FORMAT_VARIABLE_ADPI
     $ ( N_TYP , INDIC_AT_MAILLE , INDIC_G ,
     $   CAR_COL_VAR_X_DP , CAR_TITRE_VAR_X_DP ,
     $   CAR_COL_VAR_E_L , CAR_TITRE_VAR_E_L ,
     $   CAR_COL_VAR_POT_CHIM , CAR_TITRE_VAR_POT_CHIM ,
     $   CAR_COL_VAR_PB_CONV_BA , CAR_TITRE_VAR_PB_CONV_BA )
C     write(*,*) "CAR_TITRE_VAR_POT_CHIM =",
C    $           CAR_TITRE_VAR_POT_CHIM
C     write(*,*) "CAR_COL_VAR_POT_CHIM =",
C    $           CAR_COL_VAR_POT_CHIM
C---------------------------------------------------
C-----Potentiels chimiques des éléments intrinsèques
C---------------------------------------------------
      ALLOCATE ( POT_INTR ( N_TYP_INTR ) )
C-------------------------------------------------
C-----Potentiels chimiques des éléments d'addition
C-------------------------------------------------
      ALLOCATE ( POT_I ( N_TYP_INTR + 1 : N_TYP ) )
C--------------------------------------------
C-----Termes relatifs aux éléments d'addition
C--------------------------------------------
      ALLOCATE ( ALPHA_I_R ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( SOMME_I ( N_TYP_INTR + 1 : N_TYP ) )
C---------------------------------------------------------------------
C-----Calcul des nombres d'atomes de chaque espèce par maille
C-----(utiles au calcul de la somme pondérée des potentiels chimiques)
C---------------------------------------------------------------------
      N_1_MAILLE = 0
      N_2_MAILLE = 0
      N_3_MAILLE = 0
      DO I_R = 1 , N_R_1
       N_1_MAILLE = N_1_MAILLE + P_R ( I_R )
      END DO
      DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
       N_2_MAILLE = N_2_MAILLE + P_R ( I_R )
      END DO
      DO I_R = 1 + N_R_2 + N_R_1 , N_R_3 + N_R_2 + N_R_1
       N_3_MAILLE = N_3_MAILLE + P_R ( I_R )
      END DO
C------------------------------
C-----Fractions correspondantes
C------------------------------
      X_1_MAILLE = DFLOAT ( N_1_MAILLE )
     $           / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
      X_2_MAILLE = DFLOAT ( N_2_MAILLE )
     $           / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
      X_3_MAILLE = DFLOAT ( N_3_MAILLE )
     $           / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
C------------------------------------------------------------------
C-----Calcul de l'énergie et du volume de référence par maille (eV)
C------------------------------------------------------------------
      E_REF_MAILLE = E_REF / DFLOAT ( N_MAILLE_REF )
      V_REF_MAILLE = V_REF / DFLOAT ( N_MAILLE_REF )
C---------------------------------------------------------------------
C-----Calcul des enthalpies totale et par maille de la cellule sans DP
C---------------------------------------------------------------------
      H_REF = E_REF - V_REF * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
      H_REF_MAILLE = H_REF / DFLOAT ( N_MAILLE_REF )
C	write ( * , * ) V_REF_MAILLE
C------------------------------------------------------------------
C-----Calcul approché de la somme pondérée des potentiels chimiques
C-----pour une maille à partir de la relation de Gibbs-Duhem
C------------------------------------------------------------------
      S_POT = V_REF_MAILLE * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
     $      + E_REF_MAILLE
C==================================================
C==================================================
C=====Calcul optionnel de grandeurs préliminaires
C=====relatives aux complexes éventuels NON CHARGES
C==================================================
C==================================================
C----------------------------------------------
C----------------------------------------------
C-----Test sur la prise en compte des complexes
C----------------------------------------------
C----------------------------------------------
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' )
     $ THEN
C--------------------------------------------------------------------
C-----Calcul du nombre de sites du sous-réseau r dans chaque complexe
C-----et du nombre d'atomes de type i dans ce complexe
C--------------------------------------------------------------------
        ALLOCATE ( U_COMPLEXE_S_R ( N_TYPES_COMPLEXES , N_R ) )
        ALLOCATE ( V_COMPLEXE_TYPE ( N_TYPES_COMPLEXES , N_TYP ) )
        U_COMPLEXE_S_R = 0
        V_COMPLEXE_TYPE = 0
C-----------------------------
C-----Boucle sur les complexes
C-----------------------------
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
C-------------------------------------
C-----Boucle sur les sites du complexe
C-------------------------------------
          DO I_SITE = 1 , NOMBRE_SITES_COMPLEXE ( I_COMPLEXE )
           DO I_R = 1 , N_R
            IF ( I_S_R_COMPLEXE ( I_COMPLEXE , I_SITE ) .EQ. I_R )
     $      THEN
            U_COMPLEXE_S_R ( I_COMPLEXE , I_R )
     $    = U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) + 1
            END IF
           END DO
           DO I_TYP = 1 , N_TYP
            IF ( I_TYPE_COMPLEXE ( I_COMPLEXE , I_SITE ) .EQ. I_TYP )
     $      THEN
            V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP )
     $    = V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP ) + 1
            END IF
           END DO
C-----------------------------------------------
C-----Fin de la boucle sur les sites du complexe
C-----------------------------------------------
          END DO
C---------------------------------------
C-----Fin de la boucle sur les complexes
C---------------------------------------
        END DO
C--------------------------------------------------
C-----Ecriture de U_COMPLEXE_S_R et V_COMPLEXE_TYPE
C--------------------------------------------------
        write ( * , * )
     $'       Paramètres u_d(s-r) et v_d(type) des complexes :'
        do i_complexe = 1 , n_types_complexes
         write ( * , * ) '               Complexe ' , i_complexe
         write ( * , * ) '    U_COMPLEXE_S_R' ,
     $                  ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) ,
     $                    i_r = 1 , n_r )
         write ( * , * ) '   V_COMPLEXE_TYPE' ,
     $                  ( V_COMPLEXE_TYPE ( I_COMPLEXE , I_typ ) ,
     $                    i_typ = 1 , n_typ )
        end do
        WRITE ( * , 1100 )
C------------------------------------
C------------------------------------
C-----Fin du test sur INDIC_COMPLEXES
C------------------------------------
C------------------------------------
       END IF
C=======================================================
C=======================================================
C=====Fin du calcul optionnel de grandeurs préliminaires
C=====relatives aux complexes éventuels NON CHARGES
C=======================================================
C=======================================================
C######################################
C#####Ecriture dans le fichier de liste
C#####des conditions du calcul
C######################################
      OPEN ( 100 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                           // '_adpi_'
     $                           // INDIC_TYP_CALC
     $                            ( 1 : LONG_INDIC_TYP_CALC )
     $                           // '_liste' )
C------------------------ 
C-----Paramètres généraux
C------------------------ 
      WRITE ( 100 , 1100 )
      DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        WRITE ( 100 , * ) COMMENTAIRE ( I_LIGNE )
      END DO
      WRITE ( 100 , 1100 )
      WRITE ( 100 , * ) 'Nombre de types : '
      WRITE ( 100 , 502 ) N_TYP
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) "Nombre d'éléments d'addition : "
      WRITE ( 100 , 502 ) N_TYP - N_TYP_INTR
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) 'Nombre de sous-réseaux :'
      WRITE ( 100 , 502 ) N_R
      WRITE ( 100 , 1200 )
      IF ( N_R_3 .GT. 0 ) THEN
        WRITE ( 100 , * )
     $ 'Nombres de sous-réseaux 1, 2, 3 et interstitiels :'
        WRITE ( 100 , 504 ) N_R_1 , N_R_2 , N_R_3 , N_R_INTER
        WRITE ( 100 , 1200 )
      ELSE
        WRITE ( 100 , * )
     $ 'Nombres de sous-réseaux 1, 2 et interstitiels :'
        WRITE ( 100 , 503 ) N_R_1 , N_R_2 , N_R_INTER
        WRITE ( 100 , 1200 )
      END IF
      WRITE ( 100 , * )
     $   'Nombre de sites par maille pour chacun des ' , N_R ,
     $   ' sous-réseaux : '
      WRITE ( 100 , * )
     $ ( P_R ( I_R ) , I_R = 1 , N_R )
      WRITE ( 100 , 1200 )
C-------------------------------
C-----Formule du composé parfait
C-------------------------------
      IF ( N_TYP_INTR .GE. 2 ) THEN
 2212 FORMAT ( 1X , 'A(' , F4.2 , ')B(' , F4.2 , ')' )
 2213 FORMAT ( 1X , 'A(' , F4.2 , ')B(' , F4.2 , ')C(' , F4.2 , ')' )
        WRITE ( 100 , * )
     $     'Formule du composé stoechiométrique : '
        IF ( N_TYP_INTR .EQ. 2 ) THEN
          WRITE ( 100 , 2212 ) X_1_MAILLE , X_2_MAILLE
        ELSE
          WRITE ( 100 , 2213 ) X_1_MAILLE , X_2_MAILLE , X_3_MAILLE
        END IF
        WRITE ( 100 , 1200 )
      END IF
      WRITE ( 100 , * ) 'Pression (kbar) :'
      WRITE ( 100 , 501 ) PRESSION
      WRITE ( 100 , 1200 )
C---------------------------------------
C-----Ecriture optionnelle des complexes
C---------------------------------------
      IF ( N_TYPES_COMPLEXES .GT. 0 ) THEN
      WRITE ( 100 , * )
     $   'Défauts complexes :'
      WRITE ( 100 , * )
     $   '-----------------'
      WRITE ( 100 , * )
     $   '     Nombre de complexes pris en compte :'
      WRITE ( 100 , * ) N_TYPES_COMPLEXES
      WRITE ( 100 , * ) '         Numéro :'
      WRITE ( 100 , 705 )
     $ ( I_COMPLEXE ,
     $   I_COMPLEXE = 1 , N_TYPES_COMPLEXES )
      WRITE ( 100 , * ) '         Multiplicité :'
      WRITE ( 100 , 705 )
     $ ( MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ,
     $   I_COMPLEXE = 1 , N_TYPES_COMPLEXES )
      WRITE ( 100 , * ) '         Nombre de sites :'
      WRITE ( 100 , 705 )
     $ ( NOMBRE_SITES_COMPLEXE ( I_COMPLEXE ) ,
     $   I_COMPLEXE = 1 , N_TYPES_COMPLEXES )
        WRITE ( 100 , * )
     $ '            ----------------------------'
        WRITE ( 100 , * )
     $ '            Informations complémentaires :'
        WRITE ( 100 , * )
     $ '            ----------------------------'
       DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
        WRITE ( 100 , * )
     $ '            ------------------'
        WRITE ( 100 , * )
     $ '            Numéro du complexe :'
        WRITE ( 100 , * )
     $ '            ------------------'
        WRITE ( 100 , * ) I_COMPLEXE
        WRITE ( 100 , * ) 
     $ '            Nombre de sites par sous-réseau :'
        WRITE ( 100 , * )
     $ ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) , I_R = 1 , N_R ) 
        WRITE ( 100 , * )
     $ "            Nombre d'atomes par type chimique :"
        WRITE ( 100 , * )
     $ ( V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP ) , I_TYP = 1 , N_TYP )
       END DO
       ELSE
        WRITE ( 100 , * )
     $ 'Pas de défauts complexes pris en compte'
       END IF
      WRITE ( 100 , 1200 )
C---------------------------------------------
C-----Fin d'écriture optionnelle des complexes
C---------------------------------------------
C-------------
C-----Cas muVT
C-------------
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ THEN
      WRITE ( * , * )
     $ '                                   ##### CALCUL muVT #####'
      WRITE ( * , 1100 )
      WRITE ( 100 , 1 )
      WRITE ( 100 , 1100 )
      WRITE ( 100 , * )
     $ '                                   ##### CALCUL muVT #####'
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) 'Température (K) : '
      WRITE ( 100 , 501 ) TEMPERATURE
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * )
     $ "Somme pondérée des potentiels chimiques (eV) (fonction de P)"
      WRITE ( 100 , * )
     $ '(sur les éléments intrinsèques) :'
      WRITE ( 100 , 501 ) S_POT
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * )
     $  'Caractéristiques des séries de potentiels chimiques :'
      WRITE ( 100 , * )
     $  '---------------------------------------------------'
      IF ( N_TYP_INTR .GE. 2 ) THEN
         WRITE ( 100 , * ) 
     $   'mu(1) - mu(2) :'
         WRITE ( 100 , * )
     $   '--------------'
         WRITE ( 100 , 505 )
     $   D_POT_INIT_2_1 , N_D_POT_2_1 , PAS_D_POT_2_1
         IF ( N_TYP_INTR .GT. 2 ) THEN
           WRITE ( 100 , * )
     $     'mu(1) - mu(3) :'
           WRITE ( 100 , * )
     $     '--------------'
           WRITE ( 100 , 505 )
     $     D_POT_INIT_3_1 , N_D_POT_3_1 , PAS_D_POT_3_1
          END IF
        END IF
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
        WRITE ( 100 , * )
     $ '-------------------'
        WRITE ( 100 , * ) 
     $  "Elément d'addition " , I_TYP , ':'
        WRITE ( 100 , * )
     $ '-------------------'
      WRITE ( 100 , 505 )
     $  POT_I_INIT ( I_TYP ) , N_POT_I ( I_TYP ) , PAS_POT_I ( I_TYP )
      END DO
        WRITE ( 100 , 1200 )
C--------------
C-----Cas NPT=0
C--------------
      ELSE IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC )
     $     .EQ. 'NPT=0' )
     $ THEN
       WRITE ( * , * ) '##### CALCUL NPT=0 #####'
       WRITE ( * , 1100 )
       WRITE ( 100 , 1 )
       WRITE ( 100 , 1100 )
       WRITE ( 100 , * ) '##### CALCUL NPT=0 #####'
       WRITE ( 100 , 1200 )
  500 FORMAT ( 100 ( 2X , F8.4 ) )
          IF ( INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P' 
     $    .OR. INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p' ) THEN
           WRITE ( 100 , * ) 'Fractions atomiques :'
           WRITE ( 100 , 500 ) ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
           WRITE ( 100 , 1200 )
          END IF
      END IF
C------------
C-----Cas NPT
C------------
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
      WRITE ( * , * ) '##### CALCUL NPT #####'
      WRITE ( * , 1100 )
      WRITE ( 100 , 1 )
      WRITE ( 100 , 1100 )
      WRITE ( 100 , * ) '##### CALCUL NPT #####'
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) 'Température (K) : '
      WRITE ( 100 , 501 ) TEMPERATURE
      WRITE ( 100 , 1200 )
      END IF
C-----------------------------------------------------------
C-----Suffixes de sous-réseau et de type pour les DP simples
C-----------------------------------------------------------
      ALLOCATE ( W_R ( N_R ) )
      ALLOCATE ( L_W_R ( N_R ) )
      ALLOCATE ( W_TYP ( 0 : N_TYP ) )
      ALLOCATE ( L_W_TYP ( 0 : N_TYP ) )
      DO I_R = 1 , N_R
        WRITE ( UNIT = W_R ( I_R ) , FMT = '(I4)' ) I_R
        DO I = LEN ( W_R ( I_R ) ) , 1 , - 1
         IF ( ICHAR ( W_R ( I_R ) ( I : I ) ) .EQ. 32 ) THEN
           L_W_R ( I_R ) = I + 1
           GOTO 1263
         END IF
        END DO
 1263   CONTINUE
      END DO
      W_TYP ( 0 ) = 'L'
      L_W_TYP ( 0 ) = 1
      DO I_TYP = 1 , N_TYP
        WRITE ( UNIT = W_TYP ( I_TYP ) , FMT = '(I4)' ) I_TYP
        DO I = LEN ( W_TYP ( I_TYP ) ) , 1 , - 1
         IF ( ICHAR ( W_TYP ( I_TYP ) ( I : I ) ) .EQ. 32 ) THEN
           L_W_TYP ( I_TYP ) = I + 1
           GOTO 2263
         END IF
        END DO
 2263   CONTINUE
      END DO
C--------------------------------------------------------
C-----Suffixe de numéro pour les DP complexes (optionnel)
C--------------------------------------------------------
      IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) 
     $     THEN
      ALLOCATE ( L_W_COMPLEXE ( N_TYPES_COMPLEXES ) )
      ALLOCATE ( W_COMPLEXE ( 1 : N_TYPES_COMPLEXES ) )
       DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
        WRITE ( UNIT = W_COMPLEXE ( I_COMPLEXE ) , FMT = '(I4)' )
     $ I_COMPLEXE
        DO I = LEN ( W_COMPLEXE ( I_COMPLEXE ) ) , 1 , - 1
         IF ( ICHAR ( W_COMPLEXE ( I_COMPLEXE ) ( I : I ) ) .EQ. 32 )
     $   THEN
           L_W_COMPLEXE ( I_COMPLEXE ) = I + 1
           GOTO 3263
         END IF
        END DO
 3263   CONTINUE
       END DO
      END IF
C================================================================
C=====Matrices relatives à l'inversion LU du sytème d'équations
C=====(de dimension N_TYP - 1 ) donnant les fractions atomiques
C=====des éléments autres que I_TYP_0 en fonction des contraintes
C=====(en nombre N_TYP - 2) et de la fraction atomique de I_TYP_0
C=====pour le point courant en potentiels chimiques
C================================================================
C----------------------------------------
C-----Dimension de ce système d'équations
C----------------------------------------
      N_SYS_CTR = N_TYP - 1
      ALLOCATE ( MAT_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( L_MAT_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( U_MAT_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( P_MAT_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( Q_MAT_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( MAT_INV_CTR ( N_SYS_CTR , N_SYS_CTR ) )
      ALLOCATE ( VECT_CTR ( N_SYS_CTR ) ) 
C======================================================================
C=====Ouverture des fichiers de résultats par sous-réseau et par DP
C=====ainsi que des fichiers d'énergie libre et de potentiels chimiques
C=====(sauf pour le cas NPT=0)
C======================================================================
      IF ( .NOT.
     $   ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' ) )
     $ THEN
C--------------------------------
C-----Boucle sur les sous-réseaux
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------------------------
C-----Boucle sur les types (plus 0 = lacune)
C-------------------------------------------
      DO I_TYP = 0 , N_TYP
C----------------------
C-----Numéro de fichier
C----------------------
        I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'écriture de l'interstitiel
C-----pour le sous-réseau et l'espèce donnés
C--------------------------------------------
       INDIC_ECRIT_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR ) THEN
          INDIC_ECRIT_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o') )
     $  THEN
         INDIC_ECRIT_INTER = 1
        END IF
       END IF
C---------------------------------------------
C-----Ouverture du fichier s'il existe
C-----(toutes les combinaisons ( I_R , I_TYP )
C-----ne correspondent pas à un DP)
C---------------------------------------------
        IF ( ( I_R .LE. N_R_1
     $   .AND. I_TYP .NE. 1 )
     $  .OR. ( I_R .GT. N_R_1
     $   .AND. I_R .LE. N_R_1 + N_R_2
     $   .AND. I_TYP .NE. 2 ) 
     $  .OR. ( I_R .GT. N_R_1 + N_R_2
     $   .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3
     $   .AND. I_TYP .NE. 3 )
     $  .OR. ( I_R .GT. N_R_1 + N_R_2 + N_R_3  
     $   .AND. INDIC_ECRIT_INTER .EQ. 1 ) ) THEN
         IF ( I_TYP .EQ. 0 ) THEN
          OPEN ( I_FICH , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                      // '_adpi_'
     $                      // INDIC_TYP_CALC
     $                       ( 1 : LONG_INDIC_TYP_CALC )
     $                      // '_'
     $                      //  W_TYP ( I_TYP )
     $                           ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) ) 
     $                      // '_r_'
     $                      //  W_R ( I_R )
     $                   ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) ) )
          ELSE
          OPEN ( I_FICH , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                      // '_adpi_'
     $                      // INDIC_TYP_CALC
     $                       ( 1 : LONG_INDIC_TYP_CALC )
     $                      // '_'
     $                      //  W_TYP ( I_TYP ) 
     $              ( L_W_TYP ( I_TYP ) : LEN ( W_TYP ( I_TYP ) ) )
     $                      // '_r_'
     $                      //  W_R ( I_R )
     $              ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) ) )
          END IF
          WRITE ( I_FICH , * ) '-----------'
          WRITE ( I_FICH , * ) 'Sous-réseau = ' , I_R 
          WRITE ( I_FICH , * ) '---------------'
          WRITE ( I_FICH , * ) 'Défaut ponctuel = ' ,
     $                       W_TYP ( I_TYP )
     $                     ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) )
          WRITE ( I_FICH , * ) '---------------'
C-------------------------------------------
C-----Titre spécifique suivant le type de DP
C-------------------------------------------
          J = 1
          DO WHILE ( CAR_TITRE_VAR_X_DP ( J : J + 1 ) .NE. 'DP' )
           J = J + 1
          END DO
          J = J + 2
          CAR_TITRE_VAR_X_DP_TYPE = CAR_TITRE_VAR_X_DP
          IF ( I_TYP .NE. 0 ) THEN
          L = LEN ( W_TYP ( I_TYP ) ) - L_W_TYP ( I_TYP ) + 1 
     $      + LEN ( W_R ( I_R ) ) - L_W_R ( I_R ) + 1
     $      + 2
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '='
     $ // W_TYP ( I_TYP )
     $ ( L_W_TYP ( I_TYP ) : LEN ( W_TYP ( I_TYP ) ) )
     $ // '_'
     $ //  W_R ( I_R ) ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
           J = LEN_TRIM ( CAR_TITRE_VAR_X_DP_TYPE ) + 1
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '='
     $ // W_TYP ( I_TYP ) 
     $ ( L_W_TYP ( I_TYP ) : LEN ( W_TYP ( I_TYP ) ) )
     $ // '_'
     $ //  W_R ( I_R ) ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
C---------------------
C-----Cas d'une lacune
C---------------------
          ELSE
          L = LEN_TRIM ( W_TYP ( I_TYP ) ) 
     $      + LEN ( W_R ( I_R ) ) - L_W_R ( I_R ) + 1
     $      + 2
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '='
     $ // W_TYP ( I_TYP ) ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) )
     $ // '_'
     $ //  W_R ( I_R ) ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
           J = LEN_TRIM ( CAR_TITRE_VAR_X_DP_TYPE ) + 1
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '='
     $ // W_TYP ( I_TYP ) ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) )
     $ // '_'
     $ //  W_R ( I_R ) ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
          END IF
          WRITE ( I_FICH , * ) CAR_TITRE_VAR_X_DP_TYPE
     $          ( 1 : LEN_TRIM ( CAR_TITRE_VAR_X_DP_TYPE ) )
C----------------------------------
C-----Fin du test d'existence du DP
C----------------------------------
        END IF
C-----------------------------------
C-----Fin de la boucle sur les types
C-----------------------------------
       END DO
C------------------------------------------
C-----Fin de la boucle sur les sous-réseaux
C------------------------------------------
      END DO
C----------------------------------------------------
C-----Ouverture optionnelle des fichiers de complexes
C----------------------------------------------------
         IF
     $ ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) 
     $   THEN
C-----------------------------
C-----Boucle sur les complexes
C-----------------------------
      DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
          I_FICH = 1000000 + I_COMPLEXE
          OPEN
     $ ( I_FICH , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $           // '_adpi_'
     $           // INDIC_TYP_CALC
     $            ( 1 : LONG_INDIC_TYP_CALC )
     $           // '_'
     $           // 'cmplx_'
     $           //  W_COMPLEXE ( I_COMPLEXE )
     $             ( L_W_COMPLEXE ( I_COMPLEXE )
     $             : LEN ( W_COMPLEXE ( I_COMPLEXE ) ) ) )
C----------------------
C-----Ecriture du titre
C----------------------
          WRITE ( I_FICH , * ) '---------------'
          WRITE ( I_FICH , * ) 'Défaut complexe = ' ,
     $                       W_COMPLEXE ( I_COMPLEXE )
     $           ( 1 : LEN_TRIM ( W_COMPLEXE ( I_COMPLEXE ) ) )
          WRITE ( I_FICH , * ) '---------------'
C------------------------------------
C-----Ecriture des titres de colonnes
C------------------------------------
          J = 1
          DO WHILE ( CAR_TITRE_VAR_X_DP ( J : J + 1 ) .NE. 'DP' )
           J = J + 1
          END DO
          J = J + 2
C	  write ( * , * ) 'Cmplx = ' , I_COMPLEXE
C         write ( * , * ) 'J = ' , J
          CAR_TITRE_VAR_X_DP_TYPE = CAR_TITRE_VAR_X_DP
          L = LEN ( W_COMPLEXE ( I_COMPLEXE ) )
     $      - L_W_COMPLEXE ( I_COMPLEXE ) + 1
     $      + 2
C         write ( * , * ) 'LEN ( W_COMPLEXE ) = ' ,
C    $    LEN ( W_COMPLEXE ( I_COMPLEXE ) )
C        write ( * , * ) 'LEN_W_COMPLEXE = ' ,
C    $    L_W_COMPLEXE ( I_COMPLEXE ) 
C         write ( * , * ) 'L = ' , L
          CAR_TITRE_VAR_X_DP_TYPE = CAR_TITRE_VAR_X_DP
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '=c_'
     $ // W_COMPLEXE ( I_COMPLEXE )
     $ ( L_W_COMPLEXE ( I_COMPLEXE )
     $ : LEN ( W_COMPLEXE ( I_COMPLEXE ) ) )
           J = LEN_TRIM ( CAR_TITRE_VAR_X_DP_TYPE ) + 1
          CAR_TITRE_VAR_X_DP_TYPE ( J : J + L )
     $  = '=c_'
     $ // W_COMPLEXE ( I_COMPLEXE )
     $ ( L_W_COMPLEXE ( I_COMPLEXE )
     $ : LEN ( W_COMPLEXE ( I_COMPLEXE ) ) )
         WRITE ( I_FICH , * ) CAR_TITRE_VAR_X_DP_TYPE
     $     ( 1 : LEN_TRIM ( CAR_TITRE_VAR_X_DP_TYPE ) )
C------------------------------------
C-----Fin de boucle sur les complexes
C------------------------------------
       END DO
C----------------------------------------------------------
C-----Fin d'ouverture optionnelle des fichiers de complexes
C----------------------------------------------------------
      END IF
C-----------------------------------------
C-----Ouverture du fichier d'énergie libre
C-----------------------------------------
      OPEN ( 110 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                 // '_adpi_'
     $                 // INDIC_TYP_CALC
     $                  ( 1 : LONG_INDIC_TYP_CALC )
     $                 // '_'
     $                 // 'e_l' )
      WRITE ( 110 , * ) CAR_TITRE_VAR_E_L
     $ ( 1 : LEN_TRIM ( CAR_TITRE_VAR_ E_L ) )
C-------------------------------------------------
C-----Ouverture du fichier de potentiels chimiques
C-------------------------------------------------
      OPEN ( 120 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                 // '_adpi_'
     $                 // INDIC_TYP_CALC
     $                  ( 1 : LONG_INDIC_TYP_CALC )
     $                 // '_'
     $                 // 'pot_chim' )
      WRITE ( 120 , * ) CAR_TITRE_VAR_POT_CHIM
     $ ( 1 : LEN_TRIM ( CAR_TITRE_VAR_ POT_CHIM ) )
C====================================================
C=====Fin des ouvertures de fichiers pour muVT et NPT
C====================================================
      END IF
C	write ( * , * ) 'Fin ouverture fichiers'
C	stop
C===================================================
C=====Repérage des défauts et de leurs quantités GC
C=====par un indice unique (utile pour NPT=0 et NPT)
C===================================================
C-----------------------------------------------
C-----Indice de chaque DP
C-----en fonction de son sous-réseau de son type
C-----------------------------------------------
C-------------------------------------------------
C-----Cet indice est mis à zéro également
C-----pour les DP "réels" mais "non utiles",
C-----i.e. tels que "E_DP >= 0".
C-----Le nombre de DP manipulés N_TYP_D_R est aussi
C-----réduit d'une unité à chaque fois
C-----(ce nombre fixe la dimension du système NPT :
C-----N_NPT = N_TYP_D_R + 1).
C-------------------------------------------------
      ALLOCATE ( IND_D_R_TYP ( 0 : N_TYP , N_R ) )
         CALL
     $   INDICE_D_R_TYP
     $ ( N_TYP , N_TYP_INTR ,
     $   N_R , N_R_1 , N_R_2 , N_R_3 ,
     $   INDIC_R_INTER , INDIC_INTER_INTR ,
     $   E_B_D_R ,
     $   N_TYP_D_R , IND_D_R_TYP )
C-----Ecriture
C     WRITE ( 100 , * ) 'Nombre de types de DP =' , N_TYP_D_R
C	write ( * , * ) 'Fin INDICE_D_R_TYP'
      DO I_R = 1 , N_R
       DO I_TYP = 0 , N_TYP
        write ( * , * ) "I_R=",I_R," I_TYP=",I_TYP," INDICE=",
     $  IND_D_R_TYP ( I_TYP , I_R ) 
       END DO
      END DO
C=====================================================
C=====Noms de DP en fonction des types et sous-réseaux
C=====et longueurs correspondantes
C=====================================================
      ALLOCATE ( NOM_D_R_TYP ( 0 : N_TYP , N_R ) )
      ALLOCATE ( LONG_NOM_D_R_TYP ( 0 : N_TYP , N_R ) )
C--------------------------------------
C-----Boucles de type et de sous-réseau
C--------------------------------------
      DO I_TYP = 0 , N_TYP
       DO I_R = 1 , N_R
C--------------------------------------------
C-----Indicateur d'écriture de l'interstitiel
C-----pour le sous-réseau et l'espèce donnés
C--------------------------------------------
       INDIC_ECRIT_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR ) THEN
          INDIC_ECRIT_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o') )
     $  THEN
          INDIC_ECRIT_INTER = 1
        END IF
       END IF
C---------------------------
C-----Test d'existence du DP
C---------------------------
        IF ( ( I_R .LE. N_R_1
     $   .AND. I_TYP .NE. 1 )
     $  .OR. ( I_R .GT. N_R_1
     $   .AND. I_R .LE. N_R_1 + N_R_2
     $   .AND. I_TYP .NE. 2 )
     $  .OR. ( I_R .GT. N_R_1 + N_R_2
     $   .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3
     $   .AND. I_TYP .NE. 3 )
     $  .OR. ( I_R .GT. N_R_1 + N_R_2 + N_R_3
     $   .AND. INDIC_ECRIT_INTER .EQ. 1 ) ) THEN
         IF ( I_TYP .EQ. 0 ) THEN
          NOM_D_R_TYP ( I_TYP , I_R )
     $                      = 'N_'
     $                      //  W_TYP ( I_TYP )
     $                           ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) )
     $                      // '_r_'
     $                      //  W_R ( I_R )
     $                   ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
          ELSE
          NOM_D_R_TYP ( I_TYP , I_R )
     $                      = 'N_'
     $                      //  W_TYP ( I_TYP )
     $            ( L_W_TYP ( I_TYP ) : LEN ( W_TYP ( I_TYP ) ) )
     $                      // '_r_'
     $                      //  W_R ( I_R )
     $           ( L_W_R ( I_R ) : LEN ( W_R ( I_R ) ) )
          END IF
          LONG_NOM_D_R_TYP ( I_TYP , I_R )
     $  = INDEX ( NOM_D_R_TYP ( I_TYP , I_R ) , ' ' ) - 1
        END IF
C----------------------------------------------
C-----Fin des boucles de type et de sous-réseau
C----------------------------------------------
       END DO
      END DO
C	write ( * , * ) 'Fin des noms de DP'
C-------------
C-----Ecriture
C-------------
  423 FORMAT
     $ ( ' s-r = ' , I2 , ' - type = ' , I2 ,
     $   ' - ind_DP = ' , I2 ,  ' - symb_DP = "' , A , '"' )
        WRITE ( 100 , 1200 )
        DO I_R = 1 , N_R
         DO I_TYP = 0 , N_TYP
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
           WRITE ( 100 , 423 ) I_R ,
     $                         I_TYP ,
     $                         IND_D_R_TYP ( I_TYP , I_R ) ,
     $                         NOM_D_R_TYP ( I_TYP , I_R )
     $         ( 1 : LONG_NOM_D_R_TYP ( I_TYP , I_R ) )
          END IF
         END DO
        END DO
        WRITE ( 100 , 1200 )
C##################################################
C#####Calcul des énergies, volumes et enthalpies GC
C#####des DP simples (et éventuellement complexes)
C#####               DP NON CHARGES
C##################################################
      E_GC_D_R = 0.D0
      V_GC_D_R = 0.D0
      H_GC_D_R = 0.D0
      DO I_TYP = 0 , N_TYP
        DO I_R = 1 , N_R
         IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
           E_GC_D_R ( I_TYP , I_R )
     $   = E_B_D_R ( I_TYP , I_R ) - E_REF
           V_GC_D_R ( I_TYP , I_R )
     $   = V_B_D_R ( I_TYP , I_R ) - V_REF
           H_GC_D_R ( I_TYP , I_R )
     $   = E_GC_D_R ( I_TYP , I_R )
     $   + V_GC_D_R ( I_TYP , I_R )
     $   * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
         ELSE
          H_GC_D_R ( I_TYP , I_R ) = 100.D0
         END IF
        END DO
      END DO
C	write ( * , * ) 'Fin du calcul des quantités GC'
C----------------------------------------------------------------
C-----Calcul optionnel des quantités GC des complexes non chargés
C----------------------------------------------------------------
       ALLOCATE ( E_GC_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( V_GC_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
       ALLOCATE ( H_GC_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' )
     $ THEN
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
         E_GC_D_COMPLEXE ( I_COMPLEXE )
     $ = E_B_D_COMPLEXE ( I_COMPLEXE ) - E_REF
         V_GC_D_COMPLEXE ( I_COMPLEXE )
     $ = V_B_D_COMPLEXE ( I_COMPLEXE ) - V_REF
         H_GC_D_COMPLEXE ( I_COMPLEXE )
     $ = E_GC_D_COMPLEXE ( I_COMPLEXE )
     $ + V_GC_D_COMPLEXE ( I_COMPLEXE )
     $ * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
        END DO
C	write ( * , * ) 'Fin du calcul des quantités GC(cmplx)'
       END IF
C==================================================
C=====Noms de DP en fonction de leur indice unique,
C=====longueurs correspondantes,
C=====sous-réseaux et types,
C=====et énergies, volumes et enthalpies GC
C==================================================
      ALLOCATE ( NOM_D_IND ( 0 : N_TYP_D_R ) )
      ALLOCATE ( LONG_NOM_D_IND ( 0 : N_TYP_D_R ) )
      ALLOCATE ( I_TYP_R_D_IND ( 0 : N_TYP_D_R , 2 ) )
      ALLOCATE ( E_GC_D_IND ( 0 : N_TYP_D_R ) )
      ALLOCATE ( V_GC_D_IND ( 0 : N_TYP_D_R ) )
      ALLOCATE ( H_GC_D_IND ( 0 : N_TYP_D_R ) )
      DO I_R = 1 , N_R
        DO I_TYP = 0 , N_TYP
C------------------------------------------------
C-----Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C-----ne correspondent à aucun DP ( I_TYP , I_R )
C------------------------------------------------
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN 
C	write ( * , * ) I_R , I_TYP , IND_D_R_TYP ( I_TYP , I_R )
            NOM_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) ) 
     $    = NOM_D_R_TYP ( I_TYP , I_R )
            LONG_NOM_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) )
     $    = LONG_NOM_D_R_TYP ( I_TYP , I_R )
            I_TYP_R_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) , 1 )
     $    = I_TYP
            I_TYP_R_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) , 2 )
     $    = I_R
            E_GC_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) )
     $    = E_GC_D_R ( I_TYP , I_R )
            V_GC_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) )
     $    = V_GC_D_R ( I_TYP , I_R )
            H_GC_D_IND ( IND_D_R_TYP ( I_TYP , I_R ) )
     $    = H_GC_D_R ( I_TYP , I_R )
          END IF
        END DO
      END DO
      NOM_D_IND ( 0 ) = 'N_mailles'
      LONG_NOM_D_IND ( 0 ) = LEN_TRIM ( NOM_D_IND ( 0 ) ) 
C     do i_d_p = 1 , N_TYP_D_R
C	write ( * , * ) ( I_TYP_R_D_IND ( I_D_P , I ) , I = 1 , 2 ) 
C     end do
C=========================================================
C=====Fin du repérage des défauts et de leurs quantités GC
C=====par un indice unique (utile pour NPT=0 et NPT)
C=========================================================
C--------------------------------------------------------
C-----Paramètres indicateurs de types et de sous-réseaux
C-----utiles au calcul des termes de potentiels chimiques
C-----relatifs aux complexes éventuels
C--------------------------------------------------------
      ALLOCATE ( ALPHA_TYPE_COMPLEXE ( 0 : N_TYP ) )
      ALLOCATE ( BETA_S_R_COMPLEXE ( 1 : N_R ) ) 
      ALPHA_TYPE_COMPLEXE = 1
      ALPHA_TYPE_COMPLEXE ( 0 ) = 0
      BETA_S_R_COMPLEXE = 1
      DO I_R = N_R_1 + N_R_2 + N_R_3 + 1 , N_R
        BETA_S_R_COMPLEXE ( I_R ) = 0
      END DO
C	write ( * , * ) 'alpha(t) : ' ,
C    $ ( ALPHA_TYPE_COMPLEXE ( I_TYP ) , I_TYP = 0 , N_TYP )
C       write ( * , * ) 'beta(r) : ' ,
C    $ ( BETA_S_R_COMPLEXE ( I_R ) , I_R = 1 , N_R )
C--------------------------------------------------------------------
C-----Type chimique normal de chaque sous-réseau
C-----(0 pour sous-réseaux interstitiels)
C-----utile au calcul des termes de potentiels chimiques de complexes
C--------------------------------------------------------------------
      ALLOCATE ( I_TYPE_NORMAL_S_R ( N_R ) )
      I_TYPE_NORMAL_S_R = 0
      DO I_R = 1 , N_R_1
        I_TYPE_NORMAL_S_R ( I_R ) = 1
      END DO
      DO I_R = N_R_1 + 1 , N_R_1 + N_R_2
        I_TYPE_NORMAL_S_R ( I_R ) = 2
      END DO
      DO I_R = N_R_1 + N_R_2 + 1 , N_R_1 + N_R_2 + N_R_3
        I_TYPE_NORMAL_S_R ( I_R ) = 3
      END DO
C-----------------------------
C-----Ecriture de vérification
C-----------------------------
C       write ( * , * ) 't_0(r) : ' ,
C    $ ( I_TYPE_NORMAL_S_R ( I_R ) , I_R = 1 , N_R )
C-----------------------------------------------
C-----Ouverture des tableaux de fractions de DP
C-----et initialisation à 0.
C-----En mode NPT, ces valeurs seront inchangées
C-----pour les DP "non utiles", i.e. ceux pour
C-----lesquels est spécifié "E_DP >= 0"
C-----(en particulier, certains interstitiels)
C-----------------------------------------------
      ALLOCATE ( X_D_R ( 0 : N_TYP , N_R ) )
      X_D_R = 0.D0
C--------------------------------------------------
C-----Ouverture des tableaux de complexes éventuels
C--------------------------------------------------
C     write ( * , * ) 'N_TYPES_COMPLEXES = ' , N_TYPES_COMPLEXES
      ALLOCATE ( X_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
      ALLOCATE ( H_FORM_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
      ALLOCATE ( SOMME_COMPLEXE_V_TYPE_I ( N_TYP_INTR + 1 : N_TYP ) )
      X_D_COMPLEXE = 0.D0
C#################################################
C#################################################
C#####Traitement du cas muVT (DP et énergie libre)
C#################################################
C#################################################
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ THEN
C	write ( *  , * ) 'Début traitement muVT'
C----------------------------------------------------
C-----Nombre total de pas sur les éléments d'addition
C-----et produits partiels
C----------------------------------------------------
      IF ( N_TYP .GT. N_TYP_INTR ) THEN
       ALLOCATE ( PROD_PART_N_POT_I ( N_TYP_INTR + 1 : N_TYP ) )
       PROD_PART_N_POT_I ( N_TYP_INTR + 1 ) = N_POT_I ( N_TYP_INTR + 1 )
       DO I_TYP = N_TYP_INTR + 2 , N_TYP
           PROD_PART_N_POT_I ( I_TYP ) 
     $   = PROD_PART_N_POT_I ( I_TYP - 1 ) * N_POT_I ( I_TYP )
       END DO
       N_TOT_POT_I = PROD_PART_N_POT_I ( N_TYP )
      END IF
      IF ( N_TYP .EQ. N_TYP_INTR ) N_TOT_POT_I = 1
C       write ( * , * ) 'N_TOT_POT_I = ' , N_TOT_POT_I
C-----------------------------------------------------------------------
C-----Initialisation du nombre de points en composition écrits
C-----(contenus dans les fenêtres choisies),
C-----des valeurs moyennes de mu(1) - mu(2) et mu(1) - mu(3)
C-----ainsi que de leurs carrés
C-----et idem pour les éléments d'addition
C-----(calcul de l'écart-type) pour ces points
C-----(utile pour affiner le balayage en mu(1) - mu(2) et mu(1) - mu(3))
C-----------------------------------------------------------------------
      ALLOCATE ( POT_I_MOY ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( POT_I_2_MOY ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( VAR_POT_I ( N_TYP_INTR + 1 : N_TYP ) )
      N_POINTS = 0
      D_POT_MOY_2_1 = 0.D0
      D_POT_2_MOY_2_1 = 0.D0
      D_POT_MOY_3_1 = 0.D0
      D_POT_2_MOY_3_1 = 0.D0
      POT_I_MOY = 0.D0
      POT_I_2_MOY = 0.D0
C	write ( * , * ) 'Début des boucles muVT'
C-----------------------------------------------------------
C-----Ecriture concernant l'autocohérence éventuelle sur mu1
C-----------------------------------------------------------
      IF ( N_ITER_MAX_MU_1 .EQ. 1 ) THEN
       write(*,*)
       write(*,*)
     $ "*** Vous avez choisi N_ITER_MAX_MU_1 = 1"
       write(*,*) " ==> pas d'autocohérence sur mu1 ***"
       write(*,*)
      ELSE
       write(*,*)
       write(*,*)
     $ "*** Vous avez choisi N_ITER_MAX_MU_1 > 1"
       write(*,*) " ==> boucle d'autocohérence (BA) sur mu1 ***"
       write(*,*) "  NB : la mise à jour de POT_1 est faite"
       write(*,*) "       en fin d'itération de BA,"
       write(*,*) "       et l'itér. ITER_MU_1 = 1 correspond"
       write(*,*) "       aux propriétés avant BA"
       write(*,*) "       (e.g. pour x_at et G_AT écrits à l'écran)."
       write(*,*) "  NB2 : les points du balayage en pot. chimiques"
       write(*,*) "        pour lesquels la BA a rencontré"
       write(*,*) "        une difficulté de convergence sont écrits"
       write(*,*) "        eux aussi dans les fichiers de sortie."
       write(*,*) "        Ils sont également indiqués dans le fichier"
       write(*,*) '        "pb_conv_BA", en distinguant'
       write(*,*)
     $ " * nombre maxi. d'itérations atteint (INDIC_PB_CONV_BA=1)"
       write(*,*)
     $ ' * comportement divergent "NaN" (INDIC_PB_CONV_BA=2)'
       write(*,*)
      OPEN ( 450 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                 // '_adpi_'
     $                 // INDIC_TYP_CALC
     $                  ( 1 : LONG_INDIC_TYP_CALC )
     $                 // '_'
     $                 // 'pb_conv_BA' )
      WRITE ( 450 , * ) CAR_TITRE_VAR_PB_CONV_BA
     $ ( 1 : LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA ) )
      END IF
C-------------------------------------------------------------
C-----Initialisation du nombre de points
C-----pour lesquels la BA ne s'est pas déroulée convenablement
C-------------------------------------------------------------
      NB_POINTS_PB_CONV_BA = 0
C###################################################
C#####Début des boucles sur les potentiels chimiques
C###################################################
C---------------------------------------
C-----Début du balayage en mu(1) - mu(2)
C---------------------------------------
        I_POUR_CENT = 0
        DO K_D_POT_2_1 = 1 , N_D_POT_2_1
C----------------------
C-----Balayage linéaire
C----------------------
          D_POT_2_1 = D_POT_INIT_2_1
     $              + PAS_D_POT_2_1 * DFLOAT ( K_D_POT_2_1 )
C---------------------------
C-----Balayage logarithmique
C---------------------------
C         D_POT_2_1 = D_POT_INIT_2_1
C    $          + PAS_D_POT_2_1 * K_T * DLOG ( DFLOAT ( K_D_POT_2_1 ) )
C---------------------------------------
C-----Début du balayage en mu(1) - mu(3)
C---------------------------------------
        DO K_D_POT_3_1 = 1 , N_D_POT_3_1
C----------------------
C-----Balayage linéaire
C----------------------
          D_POT_3_1 = D_POT_INIT_3_1
     $              + PAS_D_POT_3_1 * DFLOAT ( K_D_POT_3_1 )
C-----------------------------------------------------------
C-----Calcul des potentiels chimiques des éléments 1, 2 et 3
C-----Pour 1, qui est l'élément de référence, il s'agit de
C-----l'estimation initiale, avant boucle d'autocohérence.
C-----------------------------------------------------------
             POT_1_INIT = ( S_POT
     $                  + DFLOAT ( N_2_MAILLE ) * D_POT_2_1
     $                  + DFLOAT ( N_3_MAILLE ) * D_POT_3_1 )
     $             / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
             POT_2 = POT_1_INIT - D_POT_2_1
             POT_3 = POT_1_INIT - D_POT_3_1
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C-----Début optionnel du balayage en mu ( i > N_TYP_INTR ) (additions)
C---------------------------------------------------------------------
C---------------------------------------------------------------------
        DO K_TOT_POT_I = 1 , N_TOT_POT_I
C----------------------------------------------------------------
C-----Indices des éléments d'addition à partir de l'indice unique
C-----et potentiels chimiques correspondants
C-----(simulation de plusieurs boucles de potentiels chimiques
C-----imbriquées à l'aide d'un seul indice K_TOT_POT_I)
C----------------------------------------------------------------
         J_COUR = K_TOT_POT_I - 1
C-------------------------------
C-----Eléments >= N_TYP_INTR + 2
C-------------------------------
         K_COUR_POT_I = MOD ( K_TOT_POT_I - 1 , N_TOT_POT_I )
         DO I_TYP = N_TYP , N_TYP_INTR + 2 , - 1
            PROD_PART_COUR = PROD_PART_N_POT_I ( I_TYP - 1 )
            K_COUR_POT_I
     $   =  J_COUR / PROD_PART_COUR + 1
            J_COUR = MOD ( J_COUR , PROD_PART_COUR ) 
C----------------------
C-----Balayage linéaire
C----------------------
            POT_I ( I_TYP )
     $    = POT_I_INIT ( I_TYP )
     $    + PAS_POT_I ( I_TYP )
     $    * DFLOAT ( K_COUR_POT_I )
C---------------------------
C-----Balayage logarithmique
C---------------------------
C           POT_I ( I_TYP )
C    $    = POT_I_INIT ( I_TYP )
C    $    + PAS_POT_I ( I_TYP ) * K_T
C    $    * DLOG ( DFLOAT ( K_COUR_POT_I ) )
         END DO
C---------------------------
C-----Elément N_TYP_INTR + 1
C---------------------------
       IF ( N_TYP .GT. N_TYP_INTR ) THEN
            K_COUR_POT_I = J_COUR + 1
C----------------------
C-----Balayage linéaire
C----------------------
            POT_I ( N_TYP_INTR + 1 )
     $    = POT_I_INIT ( N_TYP_INTR + 1 )
     $    + PAS_POT_I ( N_TYP_INTR + 1 )
     $    * DFLOAT ( K_COUR_POT_I )
C---------------------------
C-----Balayage logarithmique
C---------------------------
C           POT_I ( N_TYP_INTR + 1 )
C    $    = POT_I_INIT ( N_TYP_INTR + 1 )
C    $    + PAS_POT_I ( N_TYP_INTR + 1 ) * K_T
C    $    * DLOG ( DFLOAT ( K_COUR_POT_I ) )
       END IF
C--------------------------------------------------------
C-----Précision et nombre maximum d'itérations
C-----pour l'arrêt de la boucle d'autocohérence sur POT_1
C-----(lus dans le fichier DATA.adpi)
C--------------------------------------------------------
C      PRECISION_MU_1 = 1.D-8
C      N_ITER_MAX_MU_1 = 100
C-----Ligne ci-dessous uniquement pour test de la BA :
C-----lorsque celle-ci est divergente, cela est-il dû
C-----à la valeur initiale de POT_1, attribuée systématiquement
C-----via l'estimation avec S_POT (cf. supra) ?
C-----d'après les premiers tests sur Al_cfc(B,Ti),
C-----système pour lequel S_POT => POT_1_INIT -3,7477 eV,
C-----la modification "manuelle" de POT_1_INIT comme ci-dessous
C-----(POT_1_INIT entre  -3,60 et -3,85 eV) ne change pas
C-----le comportement qui reste divergent.
C      POT_1_INIT = -3.60
C-------------------------------------------------------------
C-----Initialisation du potentiel chimique de référence POT_1,
C-----du nombre d'itérations,
C-----de la valeur sauvegardée (itér. BA précédente) de POT_1,
C-----et de l'indicateur de "problème de convergence de la BA"
C-----avant la boucle d'autocohérence (BA)
C-------------------------------------------------------------
       ITER_MU_1 = 0
       POT_1 = POT_1_INIT
       POT_1_PREC = 1.D100
       INDIC_PB_CONV_BA = 0
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&Boucle d'autocohérence sur POT_1
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
       DO WHILE ( DABS ( POT_1 - POT_1_PREC ) .GT. PRECISION_MU_1
     $ .AND. ITER_MU_1 .LT. N_ITER_MAX_MU_1 )
       ITER_MU_1 = ITER_MU_1 + 1
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        IF ( ITER_MU_1 .EQ. N_ITER_MAX_MU_1 ) THEN
          INDIC_PB_CONV_BA = 1
        END IF
       END IF
C-----Ecriture
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        IF ( ITER_MU_1 .EQ. 1 ) THEN
         write(*,*) '----------'
         write ( * , * )
     $ "### Avant boucle d'autocohérence sur mu1 : " ,
     $         ' POT_1(initial) = ' , POT_1
         write(*,*) '----------'
        END IF
        write(*,*)
     $ '                ***** ITER_MU_1 = ' , ITER_MU_1 , ' *****'
       END IF
C--------------------------------------------------------
C--------------------------------------------------------
C-----Termes cumulés sur les sous-réseaux
C-----pour le calcul de la fraction atomique de l'alliage
C--------------------------------------------------------
C--------------------------------------------------------
C-------------------------------------------------
C-----Termes au numérateur des fractions atomiques
C-------------------------------------------------
          SOMME_1_TYP_2 = 0.D0
          SOMME_2_TYP_2 = 0.D0
          SOMME_1_TYP_3 = 0.D0
          SOMME_2_TYP_3 = 0.D0
          SOMME_I = 0.D0
C---------------------------------------------------
C-----Termes au dénominateur des fractions atomiques
C---------------------------------------------------
          SOMME_0 = 0.D0
          SOMME_MAILLE = 0.D0
C==================================
C==================================
C=====Boucle sur les sous-réseaux 1
C==================================
C==================================
          DO I_R = 1 , N_R_1
C	     write ( * , * ) 'I_R = ' , I_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------
C-----Lacunes sur les sous-réseaux 1
C-----------------------------------
           ALPHA_0_R
     $   = DEXP ( - ( H_GC_D_R ( 0 , I_R )
     $              + POT_1 ) / K_T )
C----------------------------------------------------------
C-----Antisites éventuels 2(r), 3(r) sur les sous-réseaux 1
C----------------------------------------------------------
       ALPHA_2_R = 0.D0
       ALPHA_3_R = 0.D0
       IF ( N_TYP_INTR .GT. 1 ) THEN
           ALPHA_2_R
     $   = DEXP ( - ( H_GC_D_R ( 2 , I_R )
     $              + POT_1 - POT_2 ) / K_T )
        IF ( N_TYP_INTR .GT. 2  ) THEN 
           ALPHA_3_R
     $   = DEXP ( - ( H_GC_D_R ( 3 , I_R )
     $              + POT_1 - POT_3 ) / K_T )
        END IF
       END IF
C-----------------
C-----Dénominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_2_R + ALPHA_3_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
           ALPHA_I_R = 0.D0
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
               ALPHA_I_R ( I_TYP )
     $       = DEXP ( - ( H_GC_D_R ( I_TYP , I_R )
     $                  + POT_1 - POT_I ( I_TYP ) ) / K_T )
            DENOMINATEUR_X_D_R = DENOMINATEUR_X_D_R
     $                         + ALPHA_I_R ( I_TYP )
           END DO
C=======================================
C=====Fractions de défauts ponctuels
C=====et énergies de formation de ces DP
C=======================================
C------------
C-----Lacunes
C------------
           X_D_R ( 0 , I_R ) = ALPHA_0_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 0 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 0 , I_R ) = 1.D-100
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_1
C-----Rappel : les termes "SOMME" sont également calculés
C-----car utiles plus bas au calcul des fractions atomiques em muVT
           SOMME_0
     $   = SOMME_0
     $   + DFLOAT ( P_R ( I_R ) )
     $   * ( 1.D0 - X_D_R ( 0 , I_R ) )
C--------------------------
C-----Antisites 2 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 1 ) THEN
           X_D_R ( 2 , I_R ) = ALPHA_2_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 2 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 2 , I_R ) = 1.D-100
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) + POT_1 - POT_2
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 2 , I_R )
       END IF
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           X_D_R ( 3 , I_R ) = ALPHA_3_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 3 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 3 , I_R ) = 1.D-100
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) + POT_1 - POT_3
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 3 , I_R )
       END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             X_D_R ( I_TYP , I_R ) = ALPHA_I_R ( I_TYP )
     $                             / DENOMINATEUR_X_D_R
             IF ( X_D_R ( I_TYP , I_R ) .LE. 0.D0 )
     $          X_D_R ( I_TYP , I_R ) = 1.D-100
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_1 - POT_I ( I_TYP )
             SOMME_I ( I_TYP )
     $     = SOMME_I ( I_TYP )
     $     + DFLOAT ( P_R ( I_R ) ) * X_D_R ( I_TYP , I_R )
           END DO
        END DO
C==================================
C==================================
C=====Boucle sur les sous-réseaux 2
C=====(lorsqu'ils existent)
C==================================
C==================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
C	     write ( * , * ) 'I_R = ' , I_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------------------------
C-----Lacunes et antisites 1(r) sur les sous-réseaux 2
C-----------------------------------------------------
           ALPHA_0_R
     $   = DEXP ( - ( H_GC_D_R ( 0 , I_R )
     $              + POT_2 ) / K_T )
           ALPHA_1_R
     $   = DEXP ( - ( H_GC_D_R ( 1 , I_R )
     $              - POT_1 + POT_2 ) / K_T )
C----------------------------------------------------
C-----Antisites éventuels 3(r) sur les sous-réseaux 2
C----------------------------------------------------
       ALPHA_3_R = 0.D0
       IF ( N_TYP_INTR .GT. 2 ) THEN
           ALPHA_3_R
     $   = DEXP ( - ( H_GC_D_R ( 3 , I_R )
     $              - POT_3 + POT_2 ) / K_T )
        END IF
C-----------------
C-----Dénominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_3_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
           ALPHA_I_R = 0.D0
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             ALPHA_I_R ( I_TYP )
     $     = DEXP ( - ( H_GC_D_R ( I_TYP , I_R )
     $                + POT_2 - POT_I ( I_TYP ) ) / K_T )
             DENOMINATEUR_X_D_R = DENOMINATEUR_X_D_R
     $                          + ALPHA_I_R ( I_TYP )
           END DO
C=======================================
C=====Fractions de défauts ponctuels
C=====et énergies de formation de ces DP
C=======================================
C------------
C-----Lacunes
C------------
           X_D_R ( 0 , I_R ) = ALPHA_0_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 0 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 0 , I_R ) = 1.D-100
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_2
           SOMME_0
     $   = SOMME_0
     $   + DFLOAT ( P_R ( I_R ) )
     $   * ( 1.D0 - X_D_R ( 0 , I_R ) )
C----------------
C-----Antisites 1
C----------------
           X_D_R ( 1 , I_R ) = ALPHA_1_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 1 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 1 , I_R ) = 1.D-100
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_1 + POT_2
           SOMME_1_TYP_2_R
     $   = 1.D0 - X_D_R ( 1 , I_R ) - X_D_R ( 0 , I_R )
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           X_D_R ( 3 , I_R ) = ALPHA_3_R / DENOMINATEUR_X_D_R
C          write ( * , * ) I_R , X_D_R ( 3 , I_R )
           IF ( X_D_R ( 3 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 3 , I_R ) = 1.D-100
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) - POT_3 + POT_2
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 3 , I_R )
           SOMME_1_TYP_2_R
     $   = SOMME_1_TYP_2_R - X_D_R ( 3 , I_R )
      END IF
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
            X_D_R ( I_TYP , I_R ) = ALPHA_I_R ( I_TYP )
     $                            / DENOMINATEUR_X_D_R
            IF ( X_D_R ( I_TYP , I_R ) .LE. 0.D0 )
     $           X_D_R ( I_TYP , I_R ) = 1.D-100
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_2 - POT_I ( I_TYP )
            SOMME_1_TYP_2_R
     $    = SOMME_1_TYP_2_R
     $    - X_D_R ( I_TYP , I_R )
            SOMME_I ( I_TYP )
     $    = SOMME_I ( I_TYP )
     $    + DFLOAT ( P_R ( I_R ) ) * X_D_R ( I_TYP , I_R )
         END DO
         SOMME_1_TYP_2 = SOMME_1_TYP_2
     $                 + SOMME_1_TYP_2_R * DFLOAT ( P_R ( I_R ) )
        END DO
C==================================
C==================================
C=====Boucle sur les sous-réseaux 3
C=====(lorsqu'ils existent)
C==================================
C==================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------------------------------
C-----Antisites 1(r), 2(r) et lacunes sur les sous-réseaux 3
C-----------------------------------------------------------
           ALPHA_1_R
     $   = DEXP ( - ( H_GC_D_R ( 1 , I_R )
     $              - POT_1 + POT_3 ) / K_T )
           ALPHA_2_R
     $   = DEXP ( - ( H_GC_D_R ( 2 , I_R )
     $              - POT_2 + POT_3 ) / K_T )
           ALPHA_0_R
     $   = DEXP ( - ( H_GC_D_R ( 0 , I_R )
     $              + POT_3 ) / K_T )
C-----------------
C-----Dénominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_2_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
           ALPHA_I_R = 0.D0
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             ALPHA_I_R ( I_TYP )
     $     = DEXP ( - ( H_GC_D_R ( I_TYP , I_R )
     $                + POT_3 - POT_I ( I_TYP ) ) / K_T )
             DENOMINATEUR_X_D_R = DENOMINATEUR_X_D_R
     $                          + ALPHA_I_R ( I_TYP )
           END DO
C=======================================
C=====Fractions de défauts ponctuels
C=====et énergies de formation de ces DP
C=======================================
C------------
C-----Lacunes
C------------
           X_D_R ( 0 , I_R ) = ALPHA_0_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 0 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 0 , I_R ) = 1.D-100
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_3
           SOMME_0
     $   = SOMME_0
     $   + DFLOAT ( P_R ( I_R ) )
     $   * ( 1.D0 - X_D_R ( 0 , I_R ) )
C----------------
C-----Antisites 1
C----------------
           X_D_R ( 1 , I_R ) = ALPHA_1_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 1 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 1 , I_R ) = 1.D-100
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_1 + POT_3
           SOMME_1_TYP_3_R
     $   = 1.D0 - X_D_R ( 1 , I_R ) - X_D_R ( 0 , I_R )
C----------------
C-----Antisites 2
C----------------
           X_D_R ( 2 , I_R ) = ALPHA_2_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 2 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 2 , I_R ) = 1.D-100
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) - POT_2 + POT_3
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 2 , I_R )
           SOMME_1_TYP_3_R
     $   = SOMME_1_TYP_3_R - X_D_R ( 2 , I_R )
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR 
C-------------------------------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
            X_D_R ( I_TYP , I_R ) = ALPHA_I_R ( I_TYP )
     $                            / DENOMINATEUR_X_D_R
            IF ( X_D_R ( I_TYP , I_R ) .LE. 0.D0 )
     $           X_D_R ( I_TYP , I_R ) = 1.D-100
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_3 - POT_I ( I_TYP )
            SOMME_1_TYP_3_R
     $    = SOMME_1_TYP_3_R
     $    - X_D_R ( I_TYP , I_R )
            SOMME_I ( I_TYP )
     $    = SOMME_I ( I_TYP )
     $    + DFLOAT ( P_R ( I_R ) ) * X_D_R ( I_TYP , I_R )
         END DO
C-------------------------------------
C-----Terme de somme 1 pour l'espèce 3
C-------------------------------------
         SOMME_1_TYP_3 = SOMME_1_TYP_3
     $                 + SOMME_1_TYP_3_R * DFLOAT ( P_R ( I_R ) )
        END DO
C==========================================================
C==========================================================
C=====Boucle optionnelle sur les sous-réseaux interstitiels
C==========================================================
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------
C-----Défauts 1(r)
C-----------------
           ALPHA_1_R
     $   = DEXP ( - ( H_GC_D_R ( 1 , I_R ) - POT_1 ) / K_T )
C-----------------------------------
C-----Défauts éventuels 2(r) et 3(r)
C-----------------------------------
           ALPHA_2_R = 0.D0
           ALPHA_3_R = 0.D0
           IF ( N_TYP_INTR .GT. 1 ) THEN
            ALPHA_2_R
     $    = DEXP ( - ( H_GC_D_R ( 2 , I_R ) - POT_2 ) / K_T )
             IF ( N_TYP_INTR .GT. 2 ) THEN
               ALPHA_3_R
     $       = DEXP ( - ( H_GC_D_R ( 3 , I_R ) - POT_3 ) / K_T )
             END IF
           END IF
C        write(*,*)'H_GC_D_R ( 1 , I_R )=',H_GC_D_R ( 1 , I_R )
C        write(*,*)'H_GC_D_R ( 2 , I_R )=',H_GC_D_R ( 2 , I_R )
C        write(*,*)'H_GC_D_R ( 3 , I_R )=',H_GC_D_R ( 3 , I_R )
C        write(*,*)'ALPHA_1_R=',ALPHA_1_R
C        write(*,*)'ALPHA_2_R=',ALPHA_2_R
C        write(*,*)'ALPHA_3_R=',ALPHA_3_R
C-----------------
C-----Dénominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_2_R + ALPHA_3_R
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
           ALPHA_I_R = 0.D0
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             ALPHA_I_R ( I_TYP )
     $     = DEXP ( - ( H_GC_D_R ( I_TYP , I_R )
     $                - POT_I ( I_TYP ) ) / K_T )
             DENOMINATEUR_X_D_R = DENOMINATEUR_X_D_R
     $                          + ALPHA_I_R ( I_TYP )
           END DO
C=======================================
C=====Fractions de défauts ponctuels
C=====et énergies de formation de ces DP
C=======================================
C---------------------------------------------
C-----Eléments 1, 2 et 3 (ce dernier éventuel)
C---------------------------------------------
           X_D_R ( 1 , I_R ) = ALPHA_1_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 1 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 1 , I_R ) = 1.D-100
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_1
           X_D_R ( 2 , I_R ) = ALPHA_2_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 2 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 2 , I_R ) = 1.D-100
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) - POT_2
           IF ( N_TYP_INTR .GT. 2 ) THEN
           X_D_R ( 3 , I_R ) = ALPHA_3_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 3 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 3 , I_R ) = 1.D-100
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) - POT_3
           END IF
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 2 , I_R )
           IF ( N_TYP_INTR .GT. 2 ) THEN
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_D_R ( 3 , I_R )
           END IF
           SOMME_MAILLE_R
     $   = X_D_R ( 1 , I_R ) + X_D_R ( 2 , I_R )
           IF ( N_TYP_INTR .GT. 2 ) THEN
           SOMME_MAILLE_R
     $   = SOMME_MAILLE_R + X_D_R ( 3 , I_R )  
           END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             X_D_R ( I_TYP , I_R ) = ALPHA_I_R ( I_TYP )
     $                             / DENOMINATEUR_X_D_R
             IF ( X_D_R ( I_TYP , I_R ) .LE. 0.D0 )
     $            X_D_R ( I_TYP , I_R ) = 1.D-100
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) - POT_I ( I_TYP )
             SOMME_I ( I_TYP )
     $     = SOMME_I ( I_TYP )
     $     + DFLOAT ( P_R ( I_R ) ) * X_D_R ( I_TYP , I_R )
             SOMME_MAILLE_R
     $     = SOMME_MAILLE_R + X_D_R ( I_TYP , I_R )
C------------------------------------------------
C-----Fin de la boucle sur les types > N_TYP_INTR
C------------------------------------------------
          END DO
          SOMME_MAILLE = SOMME_MAILLE
     $                 + SOMME_MAILLE_R * DFLOAT ( P_R ( I_R ) )
C========================================================
C========================================================
C-----Fin de la boucle sur les sous-réseaux interstitiels
C========================================================
C========================================================
        END DO
C----------------------------------------------------------
C-----Fin du test d'existence de sous-réseaux interstitiels
C----------------------------------------------------------
       END IF
C---------------------------------------------
C-----Initialisation des termes de sommes
C-----utiles au calcul des fractions atomiques
C-----(toujours effectuée, sinon erreur
C-----dans le calcul des fractions atomiques)
C---------------------------------------------
       SOMME_COMPLEXE_U_TYPE_2 = 0.D0
       SOMME_COMPLEXE_V_TYPE_2 = 0.D0
       SOMME_COMPLEXE_U_TYPE_3 = 0.D0
       SOMME_COMPLEXE_V_TYPE_3 = 0.D0
       SOMME_COMPLEXE_V_TYPE_I = 0.D0
       SOMME_COMPLEXE_U_TOTALE = 0.D0
       SOMME_COMPLEXE_V_TOTALE = 0.D0
C==============================================================
C==============================================================
C=====Calcul optionnel des quantités de défauts complexes
C=====(rapportées au sous-réseau spécifié pour la multiplicité)
C=====ainsi que des enthalpies de formation
C=====et des termes de sommes pour les fractions atomiques
C==============================================================
C==============================================================
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) 
     $ THEN
C=============================
C=====Boucle sur les complexes
C=============================
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
          DELTA_MU_COMPLEXE = 0.D0
C-------------------------------------
C-------------------------------------
C-----Boucle sur les sites du complexe
C-------------------------------------
C-------------------------------------
          DO I_SITE = 1 , NOMBRE_SITES_COMPLEXE ( I_COMPLEXE )
                I_S_R_COUR = I_S_R_COMPLEXE ( I_COMPLEXE , I_SITE )
                I_TYPE_COUR = I_TYPE_COMPLEXE ( I_COMPLEXE , I_SITE )
                I_TYPE_NORMAL_COUR = I_TYPE_NORMAL_S_R ( I_S_R_COUR )
C------------------------------------------------------------------
C-----Calcul du potentiel chimique de l'élément sur le site courant
C-----et du potentiel chimique "normal" sur le site courant
C------------------------------------------------------------------
C--------------------------------
C-----Cas 1 : un type intrinsèque
C--------------------------------
           IF ( N_TYP_INTR .EQ. 1 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
                POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
                POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
                POT_CHIM_NORMAL_SITE = 0.D0
              END IF
C------------------------------------
C-----Cas 2 : deux types intrinsèques
C------------------------------------
           ELSE IF ( N_TYP_INTR .EQ. 2 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .EQ. 2 ) THEN
                POT_CHIM_SITE = POT_2
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
               POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 2 ) THEN
               POT_CHIM_NORMAL_SITE = POT_2
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
               POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
               POT_CHIM_NORMAL_SITE = 0.D0
              END IF
C-------------------------------------
C-----Cas 3 : trois types intrinsèques
C-------------------------------------
           ELSE IF ( N_TYP_INTR .EQ. 3 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .EQ. 2 ) THEN
                POT_CHIM_SITE = POT_2
              ELSE IF ( I_TYPE_COUR .EQ. 3 ) THEN
                POT_CHIM_SITE = POT_3
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
                POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 2 ) THEN
                POT_CHIM_NORMAL_SITE = POT_2
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 3 ) THEN
                POT_CHIM_NORMAL_SITE = POT_3
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
                POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
               POT_CHIM_NORMAL_SITE = 0.D0
              END IF
           END IF
C-------------------------------------------------------------------
C-----Calcul du terme de potentiel chimique pour le complexe courant
C-------------------------------------------------------------------
           DELTA_MU_COMPLEXE
     $   = DELTA_MU_COMPLEXE
     $   + DFLOAT ( BETA_S_R_COMPLEXE ( I_S_R_COUR ) )
     $   * POT_CHIM_NORMAL_SITE  
     $   - DFLOAT ( ALPHA_TYPE_COMPLEXE ( I_TYPE_COUR ) )
     $   * POT_CHIM_SITE
C-----------------------------
C-----Ecriture de vérification
C-----------------------------
C      if ( I_COMPLEXE .ge. 6 ) then
C         write ( * , * ) 'Complexe ' , I_COMPLEXE
C	  write ( * , * ) 'site : ' , I_SITE
C	  write ( * , * ) 'sous-réseau : ' , I_S_R_COUR
C	  write ( * , * ) 'type : ' , I_TYPE_COUR
C         write ( * , * ) 'beta(r) : ' ,
C    $           BETA_S_R_COMPLEXE ( I_S_R_COUR )
C         write ( * , * ) 'mu_normal : ' , POT_CHIM_NORMAL_SITE
C         write ( * , * ) 'alpha(t) : ' ,
C    $		 ALPHA_TYPE_COMPLEXE ( I_TYPE_COUR )
C	   write ( * , * ) 'mu : ' , POT_CHIM_SITE
C       end if
C-----------------------------------------------
C-----------------------------------------------
C-----Fin de la boucle sur les sites du complexe
C-----------------------------------------------
C-----------------------------------------------
          END DO
C-----------------------------
C-----Ecriture de vérification
C-----------------------------
C      if ( I_COMPLEXE .ge. 6 ) then
C         write ( * , * ) 'Complexe ' , I_COMPLEXE
C	  write ( * , * ) 'POT_1 = ' , POT_1
C	  write ( * , * ) 'POT_2 = ' , POT_2
C	  write ( * , * ) 'POT_I = ' , 
C    $ ( POT_I ( I_TYP ) , I_TYP = N_TYP_INTR + 1 , N_TYP )
C	  write ( * , * ) 'delta_mu = ' , DELTA_MU_COMPLEXE
C	 write ( * , * ) '----------------------------------------'
C	end if
C--------------------------
C-----Facteur intermédiaire
C--------------------------
          FACT = H_GC_D_COMPLEXE ( I_COMPLEXE ) + DELTA_MU_COMPLEXE
C---------------------------------------
C-----Enthalpie de formation du complexe
C---------------------------------------
          H_FORM_D_COMPLEXE ( I_COMPLEXE ) = FACT
          FACT = FACT / K_T
C-------------------------
C-----Fraction de complexe
C-------------------------
          X_D_COMPLEXE ( I_COMPLEXE )
     $  = DFLOAT ( MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) )
     $  * DEXP ( - FACT )   
C-------------------------------------------------------
C-----DP complexes :
C-----termes de sommes relatifs à l'espèce 2 intrinsèque
C-----pour le calcul des fractions atomiques
C-------------------------------------------------------
       IF ( N_TYP_INTR .GT. 1 ) THEN
         SOMME_INTER_U = 0.D0
         DO I_R = N_R_1 + 1 , N_R_1 + N_R_2
          SOMME_INTER_U
     $  = SOMME_INTER_U 
     $  + DFLOAT ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
         END DO
         SOMME_COMPLEXE_U_TYPE_2
     $ = SOMME_COMPLEXE_U_TYPE_2
     $ + SOMME_INTER_U
     $ * X_D_COMPLEXE ( I_COMPLEXE )
         SOMME_COMPLEXE_V_TYPE_2
     $ = SOMME_COMPLEXE_V_TYPE_2
     $  + DFLOAT ( V_COMPLEXE_TYPE ( I_COMPLEXE , 2 ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
     $ * X_D_COMPLEXE ( I_COMPLEXE )
       END IF
C-------------------------------------------------------
C-----DP complexes :
C-----termes de sommes relatifs à l'espèce 3 intrinsèque
C-----pour le calcul des fractions atomiques
C-------------------------------------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
         SOMME_INTER_U = 0.D0
         DO I_R = N_R_1 + N_R_2 + 1 , N_R_1 + N_R_2 + N_R_3
          SOMME_INTER_U
     $  = SOMME_INTER_U
     $  + DFLOAT ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
         END DO
         SOMME_COMPLEXE_U_TYPE_3
     $ = SOMME_COMPLEXE_U_TYPE_3
     $ + SOMME_INTER_U
     $ * X_D_COMPLEXE ( I_COMPLEXE )
         SOMME_COMPLEXE_V_TYPE_3
     $ = SOMME_COMPLEXE_V_TYPE_3
     $  + DFLOAT ( V_COMPLEXE_TYPE ( I_COMPLEXE , 3 ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
     $ * X_D_COMPLEXE ( I_COMPLEXE )
       END IF
C------------------------------------------------------
C-----DP complexes :
C-----termes de sommes relatifs aux éléments d'addition
C-----pour le calcul des fractions atomiques
C------------------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         SOMME_COMPLEXE_V_TYPE_I ( I_TYP )
     $ = SOMME_COMPLEXE_V_TYPE_I ( I_TYP )
     $  + DFLOAT ( V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
     $ * X_D_COMPLEXE ( I_COMPLEXE )
       END DO
C---------------------------------------------------------------
C-----DP complexes :
C-----termes de sommes relatifs aux quantités de matière totales
C---------------------------------------------------------------
         SOMME_INTER_U = 0.D0
         DO I_R = 1 , N_R_1 + N_R_2 + N_R_3
          SOMME_INTER_U
     $  = SOMME_INTER_U
     $  + DFLOAT ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
         END DO
         SOMME_COMPLEXE_U_TOTALE
     $ = SOMME_COMPLEXE_U_TOTALE
     $ + SOMME_INTER_U
     $ * X_D_COMPLEXE ( I_COMPLEXE )
         SOMME_INTER_V = 0.D0
         DO I_TYP = 1 , N_TYP
          SOMME_INTER_V
     $  = SOMME_INTER_V
     $  + DFLOAT ( V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP ) )
         END DO
         SOMME_COMPLEXE_V_TOTALE
     $ = SOMME_COMPLEXE_V_TOTALE
     $  + SOMME_INTER_V
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE ( I_COMPLEXE ) ) )
     $ * X_D_COMPLEXE ( I_COMPLEXE )
C=======================================
C=====Fin de la boucle sur les complexes
C=======================================
        END DO
C=======================================================
C=======================================================
C=====Fin du calcul optionnel des fractions de complexes
C=====et de leurs enthalpies de formation
C=======================================================
C=======================================================
       END IF
C========================
C========================
C=====Fractions atomiques
C========================
C========================
      IF ( N_TYP_INTR .GT. 1 ) THEN
C	 write ( * , * ) 'Fractions atomiques'
C---------------------------
C-----Fraction atomique de 2
C---------------------------
C      write ( * , * ) SOMME_2_TYP_2
C      write ( * , * ) SOMME_1_TYP_2
C      write ( * , * ) SOMME_COMPLEXE_U_TYPE_2
C      write ( * , * ) SOMME_COMPLEXE_V_TYPE_2
C      write ( * , * ) SOMME_0
C      write ( * , * ) SOMME_MAILLE
C      write ( * , * ) SOMME_COMPLEXE_U_TOTALE
C      write ( * , * ) SOMME_COMPLEXE_V_TOTALE
           X_AT ( 2 )
     $ = ( SOMME_2_TYP_2 + SOMME_1_TYP_2
     $   - SOMME_COMPLEXE_U_TYPE_2 + SOMME_COMPLEXE_V_TYPE_2 )
     $ / ( SOMME_0 + SOMME_MAILLE 
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
C          write ( * , * ) 'X_AT ( 2 ) = ' , X_AT ( 2 )
         IF ( X_AT ( 2 ) .LT. 0.D0 ) X_AT ( 2 ) = 0.D0  
C          write ( * , * ) 'X_AT ( 2 ) = ' , X_AT ( 2 )
           IF ( N_TYP_INTR .GT. 2 ) THEN
C---------------------------
C-----Fraction atomique de 3
C---------------------------
C      write ( * , * ) SOMME_2_TYP_3
C      write ( * , * ) SOMME_1_TYP_3
C      write ( * , * ) SOMME_COMPLEXE_U_TYPE_3
C      write ( * , * ) SOMME_COMPLEXE_V_TYPE_3
C      write ( * , * ) SOMME_0
C      write ( * , * ) SOMME_MAILLE
C      write ( * , * ) SOMME_COMPLEXE_U_TOTALE
C      write ( * , * ) SOMME_COMPLEXE_V_TOTALE
           X_AT ( 3 )
     $ = ( SOMME_2_TYP_3 + SOMME_1_TYP_3
     $   - SOMME_COMPLEXE_U_TYPE_3 + SOMME_COMPLEXE_V_TYPE_3 )
     $ / ( SOMME_0 + SOMME_MAILLE
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
C          write ( * , * ) 'X_AT ( 3 ) = ' , X_AT ( 3 )
           IF ( X_AT ( 3 ) .LT. 0.D0 ) X_AT ( 3 ) = 0.D0
C          write ( * , * ) 'X_AT ( 3 ) = ' , X_AT ( 3 )
          END IF
      END IF
C------------------------------------------------
C-----Fractions atomiques des éléments d'addition  
C------------------------------------------------
       DO I_TYP = N_TYP_INTR + 1 , N_TYP
C	     write ( * , * ) 'SOMME_I = ' , SOMME_I ( I_TYP )
C	     write ( * , * ) 'SOMME_COMPLEXE_V_TYPE_I = ' ,
C     $  SOMME_COMPLEXE_V_TYPE_I ( I_TYP )
C	     write ( * , * ) 'SOMME_0 = ' , SOMME_0
C           write ( * , * ) 'SOMME_MAILLE = ' , SOMME_MAILLE
C	     write ( * , * ) 'SOMME_COMPLEXE_U_TOTALE = ' ,
C     $   SOMME_COMPLEXE_U_TOTALE
C	     write ( * , * ) 'SOMME_COMPLEXE_V_TOTALE = ' ,
C     $   SOMME_COMPLEXE_V_TOTALE
           X_AT ( I_TYP )
     $ = ( SOMME_I ( I_TYP ) + SOMME_COMPLEXE_V_TYPE_I ( I_TYP ) )
     $ / ( SOMME_0 + SOMME_MAILLE 
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
C	  write ( * , * ) I_TYP , X_AT ( I_TYP )
        IF ( X_AT ( I_TYP ) .LT. 0.D0 ) X_AT ( I_TYP ) = 0.D0
       END DO
C---------------------------
C-----Fraction atomique de 1
C---------------------------
       X_AT ( 1 ) = 1.D0
       DO I_TYP = 2 , N_TYP
        X_AT ( 1 ) = X_AT ( 1 ) - X_AT ( I_TYP )
       END DO
       IF ( X_AT ( 1 ) .LT. 0.D0 ) X_AT ( 1 ) = 0.D0
C	write ( * , * ) ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
C==============================================================
C==============================================================
C=====Calcul des quantités thermodynamiques dans l'ADPI           
C=====(énergie et volume par maille, entropie de configuration, 
C=====énergie libre par maille, quantités par atome)         
C==============================================================
C==============================================================
         CALL
     $   G_ADPI
     $ ( N_TYP , N_TYP_INTR ,
     $   N_R_1 , N_R_2 , N_R_3 , N_R ,
     $   P_R ,
     $   N_1_MAILLE , N_2_MAILLE , N_3_MAILLE ,
     $   E_REF_MAILLE , V_REF_MAILLE ,
     $   E_REF_TYP ,
     $   X_AT ,
     $   X_D_R ,
     $   E_GC_D_R , V_GC_D_R ,
     $   INDIC_COMPLEXES , N_TYPES_COMPLEXES ,
     $   MULTIPLICITE_COMPLEXE ,
     $   I_S_R_MULTIPLICITE_COMPLEXE ,
     $   E_GC_D_COMPLEXE , V_GC_D_COMPLEXE ,
     $   X_D_COMPLEXE ,
     $   TEMPERATURE , PRESSION ,
     $   N_AT_MAILLE ,
     $   E_MAILLE , V_MAILLE ,
     $   S_CONF_MAILLE , G_MAILLE ,
     $   E_AT , V_AT ,
     $   S_CONF_AT , G_AT ,
     $   G_AT_FORM )
C-----------------------------------------------------------
C-----Avant mise à jour de POT_1,
C-----calcul de Somme(x_i*mu_i) (pour comparaison avec G_AT)
C-----------------------------------------------------------
       SOMME_X_I_MU_I = POT_1 * X_AT ( 1 )
         IF ( N_TYP_INTR .GT. 1 ) THEN
          SOMME_X_I_MU_I = SOMME_X_I_MU_I + X_AT ( 2 ) * POT_2
          IF ( N_TYP_INTR .GT. 2 ) THEN
           SOMME_X_I_MU_I = SOMME_X_I_MU_I + X_AT ( 3 ) * POT_3
          END IF
         END IF
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
           SOMME_X_I_MU_I = SOMME_X_I_MU_I
     $                    + X_AT ( I_TYP ) * POT_I ( I_TYP )
         END DO
C----------------------------------------------------------------
C-----Sauvegarde de POT_1, faite dans tous les cas (BA ou non),
C-----quoique utile seulement si BA (i.e. si N_ITER_MAX_MU_1 > 1)
C-----(=> elle pourrait n'être faite que si N_ITER_MAX_MU_1 > 1)
C----------------------------------------------------------------
       POT_1_PREC = POT_1
C--------------------------------------------------------
C-----Si la BA est activée (i.e. si N_ITER_MAX_MU_1 > 1),
C-----mise à jour de POT_1, grâce à g_at=Somme(x_i*mu_i)
C--------------------------------------------------------
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
         POT_1 = G_AT
         IF ( N_TYP_INTR .GT. 1 ) THEN
          POT_1  = POT_1 - X_AT ( 2 ) * POT_2
          IF ( N_TYP_INTR .GT. 2 ) THEN
           POT_1  = POT_1 - X_AT ( 3 ) * POT_3
          END IF
         END IF
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
           POT_1 = POT_1 - X_AT ( I_TYP ) * POT_I ( I_TYP )
         END DO
         POT_1 = POT_1 / X_AT ( 1 )
       END IF
C-------------
C-----Ecriture
C-------------
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
       write(*,*) 'ITER_MU_1 = ' , ITER_MU_1 ,
     $ " SOMME_X_I_MU_I = " , SOMME_X_I_MU_I ,
     $ " avant mise à jour de POT_1"
       write ( * , * ) 'ITER_MU_1 = ' , ITER_MU_1 , ' x_at : ' ,
     $ ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
       write(*,*) 'ITER_MU_1 = ' , ITER_MU_1 , ' G_AT = ' , G_AT
       write(*,*) 'ITER_MU_1 = ' , ITER_MU_1 ,
     $           ' POT_1 = ' , POT_1 , 'après mise à jour (BA) via G_AT'
       write(*,*) '----------'
       END IF
C----------------------------------------------------
C-----Test de comportement divergent ("NaN") de la BA
C----------------------------------------------------
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        IF ( X_AT ( 1 ) .NE. X_AT ( 1 )
     $  .OR. G_AT .NE. G_AT
     $  .OR. POT_1 .NE. POT_1 )
     $  THEN
          INDIC_PB_CONV_BA = 2
        END IF
       END IF
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&Fin de la boucle d'autocohérence sur POT_1
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      END DO
      IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
       write(*,*)
     $ "###### Fin de la boucle d'autocohérence sur mu_1 après " ,
     $ ITER_MU_1 , " itérations ######"
       write(*,*)
      END IF
C====================================================================
C====================================================================
C=====A partir des contraintes sur la composition
C=====et de la fraction atomique calculée précédemment
C=====pour l'élément spécifié I_TYP_0,
C=====inversion LU du système d'équations de contraintes
C=====pour obtenir les valeurs des autres fractions atomiques
C=====autour desquelles doivent être centrées les fenêtres d'écriture
C=====(toujours utilisé en NPT, seulement si N_TYP > 2 en muVT)
C====================================================================
C====================================================================
         IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT'
     $ .OR. ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT'
     $  .AND. N_TYP .GT. 2 ) ) THEN
C----------------------------------------------------
C----------------------------------------------------
C-----Boucle sur les lignes du système de contraintes
C----------------------------------------------------
C----------------------------------------------------
      DO I_TYP = 1 , N_SYS_CTR
C---------------------------------------------------------
C-----Boucle sur les colonnes du système de contraintes
C-----en sautant la colonne I_TYP_0
C-----(élément dont on étudie l'effet de l'enrichissement)
C---------------------------------------------------------
         I_CTR = 0
         DO J_TYP = 1 , I_TYP_0 - 1
           I_CTR = I_CTR + 1
           MAT_CTR ( I_TYP , I_CTR )
     $   = COEF_CTR ( I_TYP , J_TYP )
         END DO
         DO J_TYP = I_TYP_0 + 1 , N_TYP
           I_CTR = I_CTR + 1
           MAT_CTR ( I_TYP , I_CTR )
     $   = COEF_CTR ( I_TYP , J_TYP )
         END DO
         IF ( I_TYP .LT. N_TYP - 1 ) THEN
           VECT_CTR ( I_TYP ) = - COEF_CTR ( I_TYP , N_TYP + 1 )
     $                          - COEF_CTR ( I_TYP , I_TYP_0 )
     $                          * X_AT ( I_TYP_0 )
         ELSE
           VECT_CTR ( I_TYP ) = 1.D0 - X_AT ( I_TYP_0 )
         END IF
C--------------------------------------------------------------
C--------------------------------------------------------------
C-----Fin de la boucle sur les lignes du système de contraintes
C--------------------------------------------------------------
C--------------------------------------------------------------
      END DO
C     write ( * , * ) '----------'
C     write ( * , * ) 'I_TYP_0 = ' , I_TYP_0
C     write ( * , * ) '----------'
C     write ( * , * ) '-------------------'
C     write ( * , * ) 'X_AT ( I_TYP_0 ) = ' , X_AT ( I_TYP_0 )
C     write ( * , * ) '-------------------'
C     write ( * , * ) '-------'
C     write ( * , * ) 'MAT_CTR'
C     write ( * , * ) '-------'
C     DO I = 1 , N_SYS_CTR
C	WRITE ( * , * ) ( MAT_CTR ( I , J ) , J = 1 , N_SYS_CTR )
C     END DO	
C     write ( * , * ) '--------'
C     write ( * , * ) 'VECT_CTR'
C     write ( * , * ) '--------'
C     WRITE ( * , * ) ( VECT_CTR ( I ) , I = 1 , N_SYS_CTR )	
C	stop
C-------------------------------------------
C-------------------------------------------
C-----Inversion LU du système de contraintes
C-------------------------------------------
C-------------------------------------------
         CALL
     $   FACT_LU
     $ ( N_SYS_CTR ,
     $   MAT_CTR ,
     $   P_MAT_CTR , Q_MAT_CTR ,
     $   L_MAT_CTR , U_MAT_CTR )
         CALL
     $   INVERSE_LU
     $ ( N_SYS_CTR ,
     $   L_MAT_CTR , U_MAT_CTR ,
     $   P_MAT_CTR , Q_MAT_CTR ,
     $   MAT_INV_CTR )
      X_AT_0_CTR = 0.D0
      I_CTR = 0
      DO I_TYP = 1 , I_TYP_0 - 1
       I_CTR = I_CTR + 1
       DO J_TYP = 1 , N_TYP - 1 
         X_AT_0_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP )
     $ + MAT_INV_CTR ( I_CTR , J_TYP ) * VECT_CTR ( J_TYP ) 
        END DO
      END DO  
      DO I_TYP = I_TYP_0 + 1 , N_TYP
       I_CTR = I_CTR + 1
       DO J_TYP = 1 , N_TYP - 1
         X_AT_0_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP )
     $ + MAT_INV_CTR ( I_CTR , J_TYP ) * VECT_CTR ( J_TYP ) 
        END DO
      END DO
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----Fractions atomiques limites
C-----(pour l'écriture dans les fichiers de résultats)
C-----déduites de la valeur courante pour l'élément spécifié
C-----et des contraintes sur les fractions atomiques
C-----------------------------------------------------------
C-----------------------------------------------------------
      DO I_TYP = 1 , I_TYP_0 - 1
         X_AT_INF_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP ) - D_X_AT ( I_TYP )
         X_AT_SUP_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP ) + D_X_AT ( I_TYP )
       END DO
      DO I_TYP = I_TYP_0 + 1 , N_TYP
         X_AT_INF_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP ) - D_X_AT ( I_TYP )
         X_AT_SUP_CTR ( I_TYP )
     $ = X_AT_0_CTR ( I_TYP ) + D_X_AT ( I_TYP )
       END DO
C     write ( * , * ) '----------'
C     write ( * , * ) 'I_TYP_0 = ' , I_TYP_0
C     write ( * , * ) '----------'
C     write ( * , * ) '-------------------'
C     write ( * , * ) 'X_AT ( I_TYP_0 ) = ' , X_AT ( I_TYP_0 )
C     write ( * , * ) '-------------------'
C      DO I_TYP = 1 , I_TYP_0 - 1
C	write ( * , * ) 'Type = ' , I_TYP
C       write ( * , * ) 'X_AT_0_CTR = ' , X_AT_0_CTR ( I_TYP )
C       write ( * , * ) 'X_AT_INF_CTR = ' , X_AT_INF_CTR ( I_TYP )
C       write ( * , * ) 'X_AT_SUP_CTR = ' , X_AT_SUP_CTR ( I_TYP )
C      END DO
C      DO I_TYP = I_TYP_0 + 1 , N_TYP
C       write ( * , * ) 'Type = ' , I_TYP 
C       write ( * , * ) 'X_AT_0_CTR = ' , X_AT_0_CTR ( I_TYP )
C       write ( * , * ) 'X_AT_INF_CTR = ' , X_AT_INF_CTR ( I_TYP )
C       write ( * , * ) 'X_AT_SUP_CTR = ' , X_AT_SUP_CTR ( I_TYP )
C      END DO
C	stop
C--------------------------------------------------------
C-----Fin du test d'utilisation du système de contraintes
C--------------------------------------------------------
       END IF
C===============================================================
C===============================================================
C=====Fin de l'utilisation optionnelle du système de contraintes
C===============================================================
C===============================================================
C==================================================================
C==================================================================
C=====Ecriture des quantités de défauts ponctuels
C=====dans les fichiers de résultats
C=====par sous-réseau r et par défaut ponctuel :
C=====(i)  si r <= N_R_1,
C=====                   DP = L(r), 2(r), 3(r)  ou i(r) (i > 3) ;
C=====(ii) si N_R_1 < r <= N_R_1 + N_R_2,
C=====                   DP = L(r), 1(r), 3(r), i(r) (i > 3) ;
C=====(iii) si N_R_1 + N_R_2 < r <= N_R_1 + N_R_2 + N_R_3,
C=====                   DP = L(r), 1(r), 2(r), i(r) (i > 3) ;
C=====(iv) si N_R_1 + N_R_2 + N_R_3 < r,
C=====                   DP = i(r) (i > 3)
C==================================================================
C==================================================================
C----------------------------------------------------
C-----Seules sont écrites les valeurs de compositions
C-----comprises dans les fenêtres spécifiées
C-----si N_TYP > 2
C----- -> indicateur de ce que le point de compo.  est
C----- dans la fenêtre (toujours vrai si 2 types) 
C----------------------------------------------------
        INDIC_FENETRE = 1
        IF ( N_TYP .GT. 2 ) THEN
         DO I_TYP = 1 , N_TYP
          IF ( X_AT ( I_TYP ) .GT. X_AT_SUP_CTR ( I_TYP )
     $    .OR. X_AT ( I_TYP ) .LT. X_AT_INF_CTR ( I_TYP ) )
     $     INDIC_FENETRE = 0
         END DO
        END IF
C-------------------------
C-----Ecriture optionnelle
C-------------------------
         IF ( ( I_ECRIT_FENETRE .EQ. 1 .AND. INDIC_FENETRE .EQ. 1 )
     $     .OR. I_ECRIT_FENETRE .EQ. 0 ) THEN
          D_POT_MOY_2_1 = D_POT_MOY_2_1 + D_POT_2_1
          D_POT_2_MOY_2_1 = D_POT_2_MOY_2_1 + D_POT_2_1 * * 2
          D_POT_MOY_3_1 = D_POT_MOY_3_1 + D_POT_3_1
          D_POT_2_MOY_3_1 = D_POT_2_MOY_3_1 + D_POT_3_1 * * 2
          DO J_TYP = N_TYP_INTR + 1 , N_TYP
           POT_I_MOY ( J_TYP ) = POT_I_MOY ( J_TYP ) + POT_I ( J_TYP )
           POT_I_2_MOY ( J_TYP ) = POT_I_2_MOY ( J_TYP )
     $                           + POT_I ( J_TYP ) * * 2
          END DO
          N_POINTS = N_POINTS + 1
        END IF
C--------------------------------
C--------------------------------
C-----Boucle sur les sous-réseaux
C--------------------------------
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------
C-----Boucle sur les types
C-------------------------
        DO I_TYP = 0 , N_TYP
C----------------------
C-----Numéro du fichier
C----------------------
          I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'écriture de l'interstitiel
C-----pour le sous-réseau et l'espèce donnés
C--------------------------------------------
       INDIC_ECRIT_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR ) THEN
          INDIC_ECRIT_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o') )
     $  THEN
          INDIC_ECRIT_INTER = 1
        END IF
       END IF
C--------------------------------------------
C-----Toutes les combinaisons ( I_R , I_TYP )
C-----ne correspondent pas à un DP
C--------------------------------------------
          IF ( ( I_R .LE. N_R_1
     $     .AND. I_TYP .NE. 1 )
     $    .OR. ( I_R .GT. N_R_1
     $     .AND. I_R .LE. N_R_1 + N_R_2
     $     .AND. I_TYP .NE. 2 )
     $    .OR. ( I_R .GT. N_R_1 + N_R_2
     $     .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3
     $     .AND. I_TYP .NE. 3 )
     $    .OR. ( I_R .GT. N_R_1 + N_R_2 + N_R_3
     $     .AND. INDIC_ECRIT_INTER .EQ. 1 ) ) THEN
C---------------------------------------------------
C-----Potentiels chimiques des éléments intrinsèques
C---------------------------------------------------
        POT_INTR ( 1 ) = POT_1
         IF ( N _TYP_INTR .GT. 1 ) THEN
          POT_INTR ( 2 ) = POT_2
          IF ( N _TYP_INTR .GT. 2 ) THEN
           POT_INTR ( 3 ) = POT_3
          END IF
        END IF
C------------------------------------------------------------------
C-----Cas binaire (tous types confondus : intrinsèques + additions)
C------------------------------------------------------------------
        IF ( N_TYP .EQ. 2 ) THEN
          WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_R ( I_TYP , I_R ) , H_FORM_D_R ( I_TYP , I_R )
        ELSE
C------------------------------------------------------
C-----Cas ternaire et plus 
C-----(tous types confondus : intrinsèques + additions)
C----- -> fenêtres de composition éventuelles
C------------------------------------------------------
         IF ( ( I_ECRIT_FENETRE .EQ. 1 .AND. INDIC_FENETRE .EQ. 1 ) 
     $     .OR. I_ECRIT_FENETRE .EQ. 0 ) THEN
             WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_R ( I_TYP , I_R ) , H_FORM_D_R ( I_TYP , I_R )
           END IF
         END IF
C--------------------------------------
C-----Fin du test d'existence de défaut
C--------------------------------------
        END IF
C----------------------------------------------
C----------------------------------------------
C-----Fin des boucles sur types et sous-réseaux
C----------------------------------------------
C----------------------------------------------
       END DO
      END DO
C===================================================
C===================================================
C=====Ecriture éventuelle des quantités de complexes
C===================================================
C===================================================
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) 
     $ THEN
C-----------------------------
C-----Boucle sur les complexes
C-----------------------------
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
C----------------------
C-----Numéro du fichier
C----------------------
          I_FICH = 1000000 + I_COMPLEXE
C----------------
C-----Cas binaire
C----------------
        IF ( N_TYP .EQ. 2 ) THEN
          WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_COMPLEXE ( I_COMPLEXE ) ,
     $       H_FORM_D_COMPLEXE ( I_COMPLEXE )
        ELSE
C---------------------------------------------------
C-----Cas ternaire et plus (fenêtres de composition)
C---------------------------------------------------
         IF ( ( I_ECRIT_FENETRE .EQ. 1 .AND. INDIC_FENETRE .EQ. 1 )
     $     .OR. I_ECRIT_FENETRE .EQ. 0 ) THEN
             WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_COMPLEXE ( I_COMPLEXE ) , 
     $       H_FORM_D_COMPLEXE ( I_COMPLEXE )
           END IF
         END IF
C---------------------------------------
C-----Fin de la boucle sur les complexes
C---------------------------------------
        END DO
C-----------------------------------------
C-----Fin du test d'existence de complexes
C-----------------------------------------
      END IF
C==================================================================
C==================================================================
C=====Ecriture de l'énergie, du volume et de l'énergie libre,
C=====en distinguant le cas binaire (fenêtres de composition sinon)
C==================================================================
C==================================================================
       IF (
     $      (
     $        ( ( I_ECRIT_FENETRE .EQ. 1 .AND. INDIC_FENETRE .EQ. 1 )
     $       .OR. I_ECRIT_FENETRE .EQ. 0 )
     $                                   .AND. ( N_TYP .NE. 2 )
     $      )
     $  .OR. ( N_TYP .EQ. 2 )
     $    ) 
     $ THEN
C-------------------------------------------------
C-----Cas d'écriture de l'énergie libre par maille
C-------------------------------------------------
            IF 
     $    ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm' )
     $      THEN
             WRITE ( 110 , CAR_COL_VAR_E_L )
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       E_MAILLE , V_MAILLE , S_CONF_MAILLE , G_MAILLE
C------------------------------------------------
C-----Cas d'écriture de l'énergie libre par atome
C------------------------------------------------
            ELSE
C- - - - - - - - - - - - - - - - - - - - - - -
C- - -Dans le cas G/atome, deux possibilités :
C- - -G totale ou G de formation
C- - - - - - - - - - - - - - - - - - - - - - -
             IF ( INDIC_G .EQ. 'T' .OR. INDIC_G .EQ. 't' )
     $       THEN
              WRITE ( 110 , CAR_COL_VAR_E_L )
     $      ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $        TEMPERATURE ,
     $        E_AT , V_AT , S_CONF_AT , G_AT
             ELSE
              WRITE ( 110 , CAR_COL_VAR_E_L )
     $      ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $        TEMPERATURE ,
     $        E_AT , V_AT , S_CONF_AT , G_AT_FORM
C- - - - - - - - - - - - - - - - - -
C- - -Fin du test "totale/formation"
C- - - - - - - - - - - - - - - - - -
             END IF
C-------------------------------
C-----Fin du test "atome/maille"
C-------------------------------
            END IF
C----------------------------------------
C-----Fin du test sur "N_TYP et fenêtres"
C----------------------------------------
       END IF
C==================================================================
C==================================================================
C=====Ecriture des potentiels chimiques
C=====en distinguant le cas binaire (fenêtres de composition sinon)
C==================================================================
C==================================================================
C----------------
C-----Cas binaire
C----------------
        IF ( N_TYP .EQ. 2 ) THEN
          WRITE ( 120 , CAR_COL_VAR_POT_CHIM )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE
        ELSE
C---------------------------------------------------
C-----Cas ternaire et plus (fenêtres de composition)
C---------------------------------------------------
         IF ( ( I_ECRIT_FENETRE .EQ. 1 .AND. INDIC_FENETRE .EQ. 1 )
     $     .OR. I_ECRIT_FENETRE .EQ. 0 ) THEN
          WRITE ( 120 , CAR_COL_VAR_POT_CHIM )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE
           END IF
         END IF
C===================================================================
C=====En cas de boucle d'autocohérence,
C=====écriture dans le fichier "pb_conv_BA" des points pour lesquels
C=====la BA ne s'est pas effectuée convenablement
C===================================================================
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        IF ( INDIC_PB_CONV_BA .NE. 0 ) THEN
          NB_POINTS_PB_CONV_BA = NB_POINTS_PB_CONV_BA + 1
          WRITE ( 450 , CAR_COL_VAR_PB_CONV_BA )
     $       POT_1_INIT ,
     $     ( POT_INTR ( J_TYP ) , J_TYP = 2 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $       INDIC_PB_CONV_BA
        END IF
       END IF
C-----------------------------------
C-----Pourcentage du calcul effectué
C-----et écriture correspondante
C-----si un seul type intrinsèque
C-----------------------------------
        IF ( N_TYP_INTR .EQ. 1 ) THEN
          I_POUR_CENT_NOUV = 100 * K_TOT_POT_I / N_TOT_POT_I
          I_ECRIT_POUR_CENT = 0
          IF ( I_POUR_CENT_NOUV .NE. I_POUR_CENT ) THEN
            I_POUR_CENT = I_POUR_CENT_NOUV
            I_ECRIT_POUR_CENT = 1
          END IF
          IF ( I_ECRIT_POUR_CENT .EQ. 1 ) THEN
           WRITE ( * , * )
     $      '           # # # # # ' , I_POUR_CENT ,
     $ ' % du calcul muVT effectués        # # # # #'
          END IF
        END IF
C#################################################
C#####Fin des boucles sur les potentiels chimiques
C#################################################
C=================================
C=====Fin de la boucle sur mu(i>3)
C=================================
        END DO
C=======================================
C=====Fin de la boucle sur mu(1) - mu(3)
C=======================================
       END DO
C----------------------------------------
C-----Pourcentage du calcul effectué
C-----et écriture correspondante
C-----si au moins deux types intrinsèques
C----------------------------------------
        IF ( N_TYP_INTR .GT. 1 ) THEN
          I_POUR_CENT_NOUV = 100 * K_D_POT_2_1 / N_D_POT_2_1
          I_ECRIT_POUR_CENT = 0
          IF ( I_POUR_CENT_NOUV .NE. I_POUR_CENT ) THEN
            I_POUR_CENT = I_POUR_CENT_NOUV
            I_ECRIT_POUR_CENT = 1
          END IF
          IF ( I_ECRIT_POUR_CENT .EQ. 1 ) THEN
           WRITE ( * , * )
     $      '           # # # # # ' , I_POUR_CENT ,
     $ ' % du calcul muVT effectués        # # # # #'
          END IF
        END IF
C=======================================
C=====Fin de la boucle sur mu(1) - mu(2)
C=======================================
      END DO
C-------------------------------------------------------------------
C-----En cas de boucle d'autocohérence, écriture du nombre de points
C-----pours lesquels la BA ne s'est pas déroulée convenablement
C-------------------------------------------------------------------
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        write(*,*)
     $ "   ----------------------------------------------------------"
        write(*,*)
     $ "      Boucle d'autocohérence sur mu_1 : nombre de points"
        write(*,*)
     $ "     pour lesquels la BA a rencontré un pb. = " ,
     $  NB_POINTS_PB_CONV_BA 
        write(*,*)
     $ '   Ces points sont répertoriés dans le fichier ""pb_conv_BA".'
        write(*,*)
     $ "   ----------------------------------------------------------"
       END IF
C-----------------------------------------------------------------------
C-----Valeurs moyennes de mu(1) - mu(2), mu(1) - mu(3)
C-----et de leurs carrés (calcul de l'écart-type) pour ces points
C-----(utile pour affiner le balayage en mu(1) - mu(2) et mu(1) - mu(3))
C-----------------------------------------------------------------------
        D_POT_MOY_2_1 = D_POT_MOY_2_1 / DFLOAT ( N_POINTS )  
        D_POT_2_MOY_2_1 = D_POT_2_MOY_2_1 / DFLOAT ( N_POINTS )
        D_POT_MOY_3_1 = D_POT_MOY_3_1 / DFLOAT ( N_POINTS )
        D_POT_2_MOY_3_1 = D_POT_2_MOY_3_1 / DFLOAT ( N_POINTS )
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
          POT_I_MOY ( I_TYP ) = POT_I_MOY ( I_TYP )
     $                        / DFLOAT ( N_POINTS )
          POT_I_2_MOY ( I_TYP ) = POT_I_2_MOY ( I_TYP ) 
     $                          / DFLOAT ( N_POINTS )
        END DO
        VAR_D_POT_2_1 = D_POT_2_MOY_2_1 - D_POT_MOY_2_1 * * 2
        VAR_D_POT_3_1 = D_POT_2_MOY_3_1 - D_POT_MOY_3_1 * * 2
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
          VAR_POT_I ( I_TYP )
     $  = POT_I_2_MOY ( I_TYP ) - POT_I_MOY ( I_TYP ) * * 2
        END DO
C-------------
C-----Ecriture
C-------------
       IF ( N_POINTS .GT. 0 ) THEN
        WRITE ( 100 , * ) 'Nombre de points écrits :'
        WRITE ( 100 , * ) N_POINTS
C--------------------------------------
C-----Ecriture dans le fichier de liste
C--------------------------------------
        IF ( N_TYP_INTR .GE. 2 ) THEN
          WRITE ( 100 , 1200 )
          WRITE ( 100 , * )
     $   'Pour les ' , N_POINTS , 'points écrits :'
          WRITE ( 100 , * )
     $   'valeur moyenne de mu(1) - mu(2) = ' , D_POT_MOY_2_1 , 'eV'
          WRITE ( 100 , * )
     $   'écart-type = ' , DSQRT ( VAR_D_POT_2_1 )
         IF ( N_TYP_INTR .GT. 2 ) THEN
          WRITE ( 100 , * )
     $   'valeur moyenne de mu(1) - mu(3) = ' , D_POT_MOY_3_1 , 'eV'
          WRITE ( 100 , * )
     $   'écart-type = ' , DSQRT ( VAR_D_POT_3_1 )
         END IF
        END IF
C-----------------------
C-----Ecriture à l'écran
C-----------------------
        IF ( N_TYP_INTR .GE. 2 ) THEN
          WRITE ( * , 1100 ) 
          WRITE ( * , * )
     $   'Pour les ' , N_POINTS , 'points écrits :'
          WRITE ( * , * )
     $   'valeur moyenne de mu(1) - mu(2) = ' , D_POT_MOY_2_1 , 'eV'
          WRITE ( * , * )
     $   'écart-type = ' , DSQRT ( VAR_D_POT_2_1 )
         IF ( N_TYP_INTR .GT. 2 ) THEN
          WRITE ( * , * )
     $   'valeur moyenne de mu(1) - mu(3) = ' , D_POT_MOY_3_1 , 'eV'
          WRITE ( * , * )
     $   'écart-type = ' , DSQRT ( VAR_D_POT_3_1 )
         END IF
        END IF
C-------------------------------------------
C-----Ecritures pour les éléments d'addition
C-------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
C	 IF ( I_TYP .NE. I_TYP_0 ) THEN
          WRITE ( * , * )
     $   'Pour les ' , N_POINTS , 'points écrits :'
          WRITE ( 100 , * )
     $   "valeur moyenne de mu pour l'élément d'addition " , I_TYP ,
     $   '= ' , POT_I_MOY ( I_TYP )
        WRITE ( 100 , * )
     $   'écart-type = ' , DSQRT ( VAR_POT_I ( I_TYP ) )
          WRITE ( * , * )
     $   "valeur moyenne de mu pour l'élément d'addition " , I_TYP ,
     $   '= ' , POT_I_MOY ( I_TYP )
        WRITE ( * , * )
     $   'écart-type = ' , DSQRT ( VAR_POT_I ( I_TYP ) )
C   	 END IF
        END DO
        WRITE ( 100 , * )
     $   'Utiliser ces valeurs pour affiner les séries éventuelles'
        WRITE ( 100 , * )
     $   'mu(1) - mu(2), mu(1) - mu(3) et mu(additions)'
        WRITE ( 100 , 1100 )
        WRITE ( * , * )
     $   'Utiliser ces valeurs pour affiner les séries éventuelles'
        WRITE ( * , * )
     $   'mu(1) - mu(2), mu(1) - mu(3) et mu(additions)'
        WRITE ( * , 1100 )
C--------------------------------------
C-----Cas où aucun point n'a été trouvé
C--------------------------------------
      ELSE
        WRITE ( * , * )
     $   '------------------------------------------------------------'
        WRITE ( * , * )
     $   "Aucun point n'a été trouvé. Vous pouvez :"
        WRITE ( * , * )
     $   '* modifier les propriétés des séries de potentiels chimiques'
        WRITE ( * , * )
     $   '* choisir des fenêtres plus larges'
        WRITE ( * , * )
     $   "   pour l'affichage des fractions atomiques"
        WRITE ( * , * )
     $   '------------------------------------------------------------'
      END IF
C###########################
C###########################
C#####Fin du traitement muVT
C###########################
C###########################
      END IF
C#################################
C#################################
C#####Traitement du cas NPT=0 (DP)
C#################################
C#################################
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC )
     $     .EQ. 'NPT=0' ) THEN
C	write ( * , * ) 'H_REF_MAILLE = ' , H_REF_MAILLE
C       write ( * , * ) 'H_GC_D_R = ' , H_GC_D_R
C---------------------------------------------------
C-----Ouverture de tableaux de résultats du simplexe
C---------------------------------------------------
      ALLOCATE ( I_1_SMPLX ( N_TYP_D_R + 1 ) )
      ALLOCATE ( I_2_SMPLX ( N_TYP ) )
C--------------------------------------------------------------------
C-----Tableau des indices initiaux des DP classés par ordre croissant
C--------------------------------------------------------------------
      ALLOCATE ( IND_INIT ( N_TYP ) )
C----------------------------------------
C-----Simplexe précédent (balayage NPT=0)
C----------------------------------------
      ALLOCATE ( I_2_SMPLX_PREC ( N_TYP ) )
C--------------------------------------------------
C-----Types et sous-réseaux des DP constitutionnels
C--------------------------------------------------
      ALLOCATE ( I_TYP_R_DP_CONST ( N_TYP - 1 , 2 ) )
C------------------------------------------------------
C-----Matrice N_TYP x N_TYP entre les variables
C-----(N_maille, n_d_const) et ( N_I )
C-----utile au calcul des potentiels chimiques en NPT=0
C-----et inverse de cette matrice
C------------------------------------------------------
      ALLOCATE ( MAT_CALC_POT ( N_TYP , N_TYP ) )
      ALLOCATE ( MAT_CALC_POT_INV ( N_TYP , N_TYP ) )
C-----------------------------------------------
C-----Tableau des potentiels chimiques à T = 0 K
C-----------------------------------------------
      ALLOCATE ( POT_CHIM_0K ( N_TYP ) )
C===============================================================
C=====Séquence 1 d'opérations relatives au cas de calcul "point"
C===============================================================
       IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $    ) THEN 
C-------------------------------------------------------
C-----Relecture dans DATA.adpi des fractions atomiques
C-----sous forme d'une chaîne de caractères
C-----(pour utilisation comme suffixe de nom de fichier)
C-------------------------------------------------------
C---------------------------
C-----Recherche de l'en-tête
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 49 )
     $   .NE. 'PARAMETRES RELATIFS A UN CALCUL NPT ou NPT=0' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 32 )
     $   .NE. "fractions atomiques de l'alliage" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------------
C-----Fractions atomiques de l'alliage (la somme doit valoir 1)
C--------------------------------------------------------------
        READ ( 10 , '(A)' ) CHAINE_LECT
        READ ( 10 , '(A)' ) LIGNE_X_AT
C----------------------------------------
C-----Analyse des fractions atomiques
C-----sous forme de chaînes de caractères
C----------------------------------------
         CALL
     $   ANALYSE_LIGNE
     $ ( LIGNE_X_AT ,
     $   N_CHAINES_X_AT ,
     $   I_DEBUT , I_FIN ,
     $   LONG_CHAINE_X_AT , CHAINE_X_AT )
C        write ( * , * ) 'La ligne ' ,
C    $   ligne_x_at ( 1 : len_trim (ligne_x_at ) )
C        write ( * , * ) 'contient ' , N_CHAINES_x_at , 'chaînes'
C        do i = 1 , N_CHAINES_X_AT
C         write ( * , * ) 'chaîne, début, fin, longueur, contenu' ,
C    $                   i, I_DEBUT ( i ) , I_FIN ( i ) ,
C    $                   LONG_CHAINE_X_AT ( I ) ,
C    $            CHAINE_X_AT ( I ) ( 1 : LONG_CHAINE_X_AT ( I ) )
C        end do
C-----------------
C-----Vérification
C-----------------
      IF ( N_CHAINES_X_AT .NE. N_TYP ) THEN
        WRITE ( * , * ) '------------------------------------------'
        WRITE ( * , * ) 'Nombre de types chimiques = ' , N_TYP
        WRITE ( * , * ) '=> nombre de fractions atomiques incorrect'
        WRITE ( * , * ) "(vérifier l'absence de tabulation)"
        WRITE ( * , * ) '------------------------------------------'
        CALL INTERRUPTION
      END IF
C-----------------------------------------------------
C-----Constitution d'une chaîne unique de suffixe X_AT
C-----------------------------------------------------
      CHAINE_TOTALE_X_AT = '' 
      DO I_TYP = 1 , N_TYP
          CHAINE_TOTALE_X_AT 
     $  = CHAINE_TOTALE_X_AT
     $ ( 1 : LEN_TRIM ( CHAINE_TOTALE_X_AT ) )
     $ // 'x_' // W_TYP ( I_TYP ) 
     $          ( L_W_TYP ( I_TYP ) : LEN ( W_TYP ( I_TYP ) ) )
     $ // '=' //  CHAINE_X_AT ( I_TYP )
     $          ( 1 : LONG_CHAINE_X_AT ( I_TYP ) )
        IF ( I_TYP .LT. N_TYP ) THEN
          CHAINE_TOTALE_X_AT
     $  = CHAINE_TOTALE_X_AT ( 1 : LEN_TRIM ( CHAINE_TOTALE_X_AT ) )
     $ // '_' 
        END IF
      END DO
      LONG_CHAINE_TOTALE_X_AT = LEN_TRIM ( CHAINE_TOTALE_X_AT )
C     write ( * , * ) CHAINE_TOTALE_X_AT ( 1 : LONG_CHAINE_TOTALE_X_AT )
C-------------
C-----Ecriture
C-------------
          OPEN ( 500 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                      // '_adpi_' 
     $                      // INDIC_TYP_CALC
     $                       ( 1 : LONG_INDIC_TYP_CALC )
     $                      // '_'
     $                      // CHAINE_TOTALE_X_AT 
     $                     ( 1 : LONG_CHAINE_TOTALE_X_AT ) )
       WRITE ( 500 , 1100 )
       DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        WRITE ( 500 , * ) COMMENTAIRE ( I_LIGNE )
       END DO
       WRITE ( 500 , 1100 )
       WRITE ( 500 , * ) '##### CALCUL NPT=0 #####'
       WRITE ( 500 , 1200 )
       WRITE ( 500 , * ) 'Fractions atomiques :'
       WRITE ( 500 , 500 ) ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
       WRITE ( 500 , 1200 )
C-----------------------------------------------------
C-----Nombres de pas en r et theta = 1 en mode "point"
C-----------------------------------------------------
        N_PAS_R = 1
        N_PAS_THETA = 0
C======================================================================
C=====Fin de séquence 1 d'opérations relatives au cas de calcul "point"
C======================================================================
       ELSE
C====================================================================
C=====Séquence 1 d'opérations relatives au cas de calcul "balayage x"
C====================================================================
          OPEN ( 500 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                      // '_adpi_' 
     $                      // INDIC_TYP_CALC
     $                       ( 1 : LONG_INDIC_TYP_CALC ) )
       WRITE ( 500 , 1100 )
       DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        WRITE ( 500 , * ) COMMENTAIRE ( I_LIGNE )
       END DO
       WRITE ( 500 , 1100 )
       WRITE ( 500 , * ) '##### CALCUL NPT=0 #####'
       WRITE ( 500 , 1200 )
  335 FORMAT ( 1X ,
     $ 'Identification des DP (rotation dans le sens trigonométrique)' )
  435  FORMAT ( 1X ,
     $ 24 ( '-' ) ,  'Limites de zones de DP' , 24 ( '-' ) ,
     $ '|' , 6 ( '-' ) ,
     $ 'DP et pot. chim. juste "après" ces limites' , 6 ( '-' ) )
  436  FORMAT ( 8X , 'x(ron)' , 8X , 'y(ron)' ,
     $          7X , 'x_at(1)' , 
     $          7X , 'x_at(2)' ,
     $          7X , 'x_at(3)' , 1X , '|' , 1X ,
     $          10X , 'DP' ,
     $          9X , 'mu(1)(eV)' ,
     $          3X , 'mu(2)(eV)' ,
     $          3X , 'mu(3)(eV)' )
       WRITE ( 500 , 435 )
       WRITE ( 500 , 436 )
C-------------------------------------------------------------
C-----Ouverture du fichier de "minimisation impossible pour H"
C-----(optionnel)
C-------------------------------------------------------------
       WRITE ( * , * ) 
     $ '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
       WRITE ( * , * ) 'Calcul NPT=0 :'
       WRITE ( * , * )
     $ "Désirez-vous un fichier annexe contenant les compositions"
       WRITE ( * , * ) 
     $ "pour lesquelles la minimisation a échoué"
       WRITE ( * , * ) 
     $ "(o pour oui, autre caractère pour non) ?"
       WRITE ( * , * ) 
     $ '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
       READ ( * , * ) CAR_FICH_550
       IF ( CAR_FICH_550 .EQ. 'o' ) THEN
          OPEN ( 550 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                      // '_adpi_'
     $                      // INDIC_TYP_CALC
     $                       ( 1 : LONG_INDIC_TYP_CALC )
     $                      // '_'
     ^                      // 'pb_min_H' )
       WRITE ( 550 , 1100 )
       DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        WRITE ( 550 , * ) COMMENTAIRE ( I_LIGNE )
       END DO
       WRITE ( 550 , 1100 )
       WRITE ( 550 , * ) '##### CALCUL NPT=0 #####'
       WRITE ( 550 , 1200 )
       WRITE ( 550 , * ) 'Fractions atomiques avec pb de min. H:'
       WRITE ( 550 , 1200 )
       END IF
C------------------------------------------------
C-----Ecritures dans le fichier de liste
C-----(un seul cercle pour identification des DP)
C------------------------------------------------
       WRITE ( 100 , 1100 )
       WRITE ( 100 , 335 )
       WRITE ( 100 , 435 )
       WRITE ( 100 , 436 )
       END IF
C-------------------------------------
C-----Ouverture du tableau du simplexe
C-------------------------------------
      ALLOCATE ( TAB_SMPLX_D_R ( N_TYP + 2 , N_TYP_D_R + 2 ) )
C--------------------------------------------------------------
C-----Indicateur d'existence d'au moins un point de composition
C-----pour lequel (i) soit la fonction est non bornée,
C-----(ii) soit il n'y a pas de solution 
C--------------------------------------------------------------
      I_ERREUR_SIMPLEXE = 0
C#############################################
C#############################################
C#####Boucle sur les pas de balayage "r,theta"
C#####(un seul si mode "point")
C#############################################
C#############################################
      I_POUR_CENT = 0
C----------------- 
C----------------- 
C-----Boucle sur r
C----------------- 
C----------------- 
      DO I_PAS_R = 1 , N_PAS_R
          I_POUR_CENT_NOUV = 100 * I_PAS_R / N_PAS_R
          I_ECRIT_POUR_CENT = 0
          IF ( I_POUR_CENT_NOUV .NE. I_POUR_CENT ) THEN
            I_POUR_CENT = I_POUR_CENT_NOUV
            I_ECRIT_POUR_CENT = 1
          END IF
          IF ( I_ECRIT_POUR_CENT .EQ. 1 ) THEN
           WRITE ( * , * )
     $      I_POUR_CENT , ' % du calcul NPT=0 effectués'
          END IF
C-------------------------------------------------------
C-------------------------------------------------------
C-----Boucle sur theta (jusqu'à N_PAS_THETA + 1
C-----pour pouvoir comparer les points 1 et N_PAS_THETA)
C-------------------------------------------------------
C-------------------------------------------------------
        DO I_PAS_THETA = 1 , N_PAS_THETA + 1
C-------------------------
C-----Composition courante
C-------------------------
          X_COORD_1_STOECHIO
     $  = 1.D0
     $  - ( DFLOAT ( N_1_MAILLE ) - DFLOAT ( N_2_MAILLE ) )
     $  / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
          X_COORD_1_STOECHIO
     $  = X_COORD_1_STOECHIO / 2.D0
          X_COORD_2_STOECHIO
     $  = 1.D0
     $  - ( DFLOAT ( N_1_MAILLE ) + DFLOAT ( N_2_MAILLE ) )
     $  / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
          X_COORD_2_STOECHIO
     $  = X_COORD_2_STOECHIO * DSQRT ( 3.D0 ) / 2.D0
         THETA_COUR = DFLOAT ( I_PAS_THETA )
     $              / DFLOAT ( N_PAS_THETA )
     $              * 2.D0 * PI
         X_COORD_1 = X_COORD_1_STOECHIO
     $             + DFLOAT ( I_PAS_R ) / DFLOAT ( N_PAS_R )
     $             * R_MAX
     $             * DCOS ( THETA_COUR )
         X_COORD_2 = X_COORD_2_STOECHIO
     $             + DFLOAT ( I_PAS_R ) / DFLOAT ( N_PAS_R )
     $             * R_MAX
     $             * DSIN ( THETA_COUR )
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' ) THEN
         X_AT ( 1 ) = 1.D0 - ( X_COORD_1 + X_COORD_2 / DSQRT ( 3.D0 ) )
         X_AT ( 2 ) = X_COORD_1 - X_COORD_2 / DSQRT ( 3.D0 )
         X_AT ( 3 ) = 1.D0 - X_AT ( 1 ) - X_AT ( 2 )
C--------------------------------------------------
C-----Vérification de la validité de la composition
C-----ainsi obtenue par balayage
C--------------------------------------------------
         I_COMPO = 1
         DO I_TYP = 1 , N_TYP
           IF ( X_AT ( I_TYP ) .GT. 1.D0
     $     .OR. X_AT ( I_TYP ) .LT. 0.D0 ) THEN
            I_COMPO = 0
           END IF
         END DO
       END IF
C#####################################################################
C#####Remplissage du tableau du simplexe contenant :
C#####(i) la fonction à minimiser en première ligne
C#####(i) les contraintes dans les lignes 2 à N_TYP + 1
C#####(i) la fonction auxiliaire dans la ligne N_TYP + 2
C#####Ce tableau contient N_TYP_D_R + 2 colonnes :
C#####la première colonne contient les valeurs des contraintes,
C#####les autres colonnes contiennent les coefficients des contraintes
C#####################################################################
C------------------------------------------
C-----Initialisation du tableau du simplexe
C------------------------------------------
      TAB_SMPLX_D_R = 0.D0
      TAB_SMPLX_D_R ( 1 , 2 ) = H_REF_MAILLE
C=================================
C=====Boucles sur les sous-réseaux
C=================================
      DO I_R = 1 , N_R
C---------------------------------------------------------
C-----Partie du tableau correspondant au nombre de mailles
C-----(deuxième colonne)
C---------------------------------------------------------
C--------------------------
C-----Nombre total d'atomes
C--------------------------
       IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           TAB_SMPLX_D_R ( 2 , 2 )
     $   = TAB_SMPLX_D_R ( 2 , 2 )
     $   - DFLOAT ( P_R ( I_R ) )
       END IF
C-----------------------
C-----Elément 2 éventuel
C-----------------------
       IF ( N_TYP_INTR .GE. 2 ) THEN
        IF ( I_R .GT. N_R_1 .AND. I_R .LE. N_R_1 + N_R_2 ) THEN
           TAB_SMPLX_D_R ( 3 , 2 )
     $   = TAB_SMPLX_D_R ( 3 , 2 )
     $   + DFLOAT ( P_R ( I_R ) ) * ( X_AT ( 2 ) - 1.D0 )
        ELSE IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           TAB_SMPLX_D_R ( 3 , 2 )
     $   = TAB_SMPLX_D_R ( 3 , 2 )
     $   + DFLOAT ( P_R ( I_R ) ) * X_AT ( 2 )
        END IF
       END IF
C-----------------------
C-----Elément 3 éventuel
C-----------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        IF ( I_R .GT. N_R_1 + N_R_2
     $ .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           TAB_SMPLX_D_R ( 4 , 2 )
     $   = TAB_SMPLX_D_R ( 4 , 2 )
     $   + DFLOAT ( P_R ( I_R ) ) * ( X_AT ( 3 ) - 1.D0 )
       ELSE IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           TAB_SMPLX_D_R ( 4 , 2 )
     $   = TAB_SMPLX_D_R ( 4 , 2 )
     $   + DFLOAT ( P_R ( I_R ) ) * X_AT ( 3 )
        END IF
       END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
       DO J_TYP = N_TYP_INTR + 1 , N_TYP
         IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           TAB_SMPLX_D_R ( J_TYP + 1 , 2 )
     $   = TAB_SMPLX_D_R ( J_TYP + 1 , 2 )
     $   + DFLOAT ( P_R ( I_R ) ) * X_AT ( J_TYP )
         END IF
       END DO             
C==========================
C=====Boucles sur les types
C==========================
        DO I_TYP = 0 , N_TYP
C----------------------------------------------------
C-----Indicateur de prise en compte de l'interstitiel
C-----pour le sous-réseau et le type courant
C----------------------------------------------------
       INDIC_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR ) THEN
          INDIC_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o') )
     $  THEN
          INDIC_INTER = 1
          END IF
         END IF
        IF ( ( I_R .LE. N_R_1
     $   .AND. I_TYP .NE. 1 )
     $  .OR. ( I_R .GT. N_R_1
     $   .AND. I_R .LE. N_R_1 + N_R_2
     $   .AND. I_TYP .NE. 2 )
     $  .OR. ( I_R .GT. N_R_1 + N_R_2
     $   .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3
     $   .AND. I_TYP .NE. 3 )
     $  .OR. ( I_R .GT. N_R_1 + N_R_2 + N_R_3
     $   .AND. INDIC_INTER .EQ. 1 ) ) THEN
C==================================================
C=====Remplissage de la première ligne (fonction H)
C==================================================
         TAB_SMPLX_D_R ( 1 , IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $ = H_GC_D_R ( I_TYP , I_R )
C=========================================
C=====Remplissage des lignes 2 à N_TYP + 1
C=====(contraintes sous forme "f = 0")
C=========================================
C-----------------------------------------
C-----Ligne de la quantité totale d'atomes
C-----------------------------------------
       TAB_SMPLX_D_R ( 2 , 1 ) = 1.D0
       IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( 2 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = 1.D0
           END IF
       ELSE 
           IF ( I_TYP .NE. 0 ) THEN
                    TAB_SMPLX_D_R ( 2 , 
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - 1.D0
           END IF
       END IF
C---------------------------------------
C-----Ligne éventuelle de l'élément 2
C-----(si au moins 2 types intrinsèques)
C---------------------------------------
       IF ( N_TYP_INTR .GE. 2 ) THEN 
        IF ( I_R .GT. N_R_1 .AND. I_R .LE. N_R_1 + N_R_2 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = 1.D0 - X_AT ( 2 )
           ELSE IF ( I_TYP .NE. 2 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = 1.D0
           END IF
       ELSE IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - X_AT ( 2 )
          ELSE IF ( I_TYP .EQ. 2 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - 1.D0
           END IF
       ELSE
           IF ( I_TYP .EQ. 2 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( 2 ) - 1.D0 
           ELSE IF ( I_TYP .NE. 0 ) THEN
                    TAB_SMPLX_D_R ( 3 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( 2 )
           END IF
        END IF
       END IF
C----------------------------------
C-----Ligne de l'élément 3 éventuel
C-----(si 3 types intrinsèques)
C----------------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        IF ( I_R .GT. N_R_1 + N_R_2
     $ .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = 1.D0 - X_AT ( 3 )
           ELSE IF ( I_TYP .NE. 3 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = 1.D0
           END IF
       ELSE IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - X_AT ( 3 )
          ELSE IF ( I_TYP .EQ. 3 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - 1.D0
           END IF
       ELSE
           IF ( I_TYP .EQ. 3 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( 3 ) - 1.D0
           ELSE IF ( I_TYP .NE. 0 ) THEN
                    TAB_SMPLX_D_R ( 4 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( 3 )
           END IF
        END IF
       END IF
C---------------------------------------------
C-----Lignes des éléments d'addition éventuels
C---------------------------------------------
       DO J_TYP = N_TYP_INTR + 1 , N_TYP
         IF ( I_R .LE. N_R_1 + N_R_2 + N_R_3 ) THEN
           IF ( I_TYP .EQ. 0 ) THEN
                    TAB_SMPLX_D_R ( J_TYP + 1 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - X_AT ( J_TYP )
           ELSE IF ( I_TYP .EQ. J_TYP ) THEN
                    TAB_SMPLX_D_R ( J_TYP + 1 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = - 1.D0
           END IF
       ELSE
           IF ( I_TYP .EQ. J_TYP ) THEN
                    TAB_SMPLX_D_R ( J_TYP + 1 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( J_TYP ) - 1.D0 
           ELSE IF ( I_TYP .NE. 0 ) THEN
                    TAB_SMPLX_D_R ( J_TYP + 1 ,
     $                              IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $            = X_AT ( J_TYP )
           END IF
          END IF
C-------------------------------------------------
C-----Fin de la boucle sur les éléments d'addition
C-------------------------------------------------
         END DO
C---------------------------------- 
C-----Fin du test d'existence de DP
C---------------------------------- 
         END IF
C----------------------------------------------------------
C-----Fin des boucles de sous-réseaux et de types chimiques
C----------------------------------------------------------
        END DO
      END DO
C==============================================
C=====Fin du remplissage du tableau du simplexe
C==============================================
C-----------------------------------------------------------
C-----La procédure du simplexe recherchant le maximum
C-----de la fonction lue en première ligne de TAB_SMPLX_D_R,
C-----on transforme le problème du minimum d'enthalpie
C-----en maximum (- enthalpie)
C-----------------------------------------------------------
      DO I_VAR = 1 , N_TYP_D_R + 2
         TAB_SMPLX_D_R ( 1 , I_VAR ) =  - TAB_SMPLX_D_R ( 1 , I_VAR )
      END DO
C--------------------------------------------------------------
C-----Procédure du simplexe :
C-----il y a zéro contraintes du type M1 ( < ) ou M2 ( > )
C-----et N_TYP contraintes du type M3 ( = )
C-----Les arguments de la procédures sont, dans l'ordre :
C----- * le tableau du simplexe ;
C----- * les nombres N_E de contraintes et N_VAR de variables ;
C----- * les nombres de lignes et de colonnes du tableau
C-----   (respectivement N_E + 2 , N_VAR + 1)
C----- * les nombres de contraintes du type (<),  (>) et (=)
C----- * l'indicateur d'absence de solution (IC)
C-----   et les deux tableaux de résultats
C--------------------------------------------------------------
         CALL
     $   SIMPLX
     $ ( TAB_SMPLX_D_R ,
     $   N_TYP , N_TYP_D_R + 1 ,
     $   N_TYP + 2 , N_TYP_D_R + 2 ,
     $   0 , 0 , N_TYP ,
     $   IC , I_1_SMPLX , I_2_SMPLX )
C-----------------------------------------------------------------------
C-----Une erreur de minimisation (fonction non bornée, pas de solution)
C-----est retenue pour écriture d'un avertissement en fin de calcul.
C-----Remarque : seuls les points de composition corrects sont concernés
C-----------------------------------------------------------------------
         IF ( IC .NE. 0 .AND. I_COMPO .EQ. 1 ) THEN
           I_ERREUR_SIMPLEXE = 1
         END IF
C-----------------------------------------------------------------
C-----En mode "balayage x", comparaison du point avec le précédent
C-----(comparaison sur theta pour un même r)
C-----en excluant le premier point en theta,
C-----car il correspond à r différent
C-----Mise à jour de I_2_SMPLX_PREC
C-----------------------------------------------------------------
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' ) THEN
        INDIC_LIMITE_ZONE = 0
         IF ( I_PAS_THETA .GT. 1 ) THEN
         DO I_TYP = 1 , N_TYP
          INDIC_MEME_DP = 0
          DO J_TYP = 1 , N_TYP
           IF ( I_2_SMPLX_PREC ( I_TYP ) .EQ. I_2_SMPLX ( J_TYP )  )
     $     THEN
            INDIC_MEME_DP = 1
           END IF 
          END DO
          IF ( INDIC_MEME_DP .EQ. 0 ) THEN
            INDIC_LIMITE_ZONE = 1
          END IF
         END DO
        END IF
        IF ( I_COMPO .EQ. 0 ) THEN
          INDIC_LIMITE_ZONE = 0
        END IF
        I_2_SMPLX_PREC = I_2_SMPLX
       END IF
C===============================================
C=====Si l'on se trouve en limite de zone,
C=====calcul des potentiels chimiques de la zone
C===============================================
        IF ( INDIC_LIMITE_ZONE .EQ. 1 ) THEN
C-------------------------------------------------------------------
C-----Classement par ordre croissant des indices des défauts trouvés
C-------------------------------------------------------------------
         CALL
     $   TRI_CROISSANT_ENTIER
     $ ( N_TYP , I_2_SMPLX , IND_INIT )
C----------------------------------------------------------------
C-----Recherche des types et sous-réseaux des DP constitutionnels
C----------------------------------------------------------------
        I_DP_CONST = 0
        DO I_TYP = 2 , N_TYP
          I_DP_CONST = I_DP_CONST + 1
          I_TYP_R_DP_CONST ( I_DP_CONST , 1 )
     $  = I_TYP_R_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 , 1 )
          I_TYP_R_DP_CONST ( I_DP_CONST , 2 )
     $  = I_TYP_R_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 , 2 )
        END DO
C-----------------------------------------------------------------
C-----Constitution de la matrice N_TYP x N_TYP entre les variables
C-----(N_maille, n_d_const) et ( N_I )
C-----utile au calcul des potentiels chimiques en NPT=0
C-----et inverse de cette matrice
C-----------------------------------------------------------------
        MAT_CALC_POT = 0.D0
C=================================
C=====Première colonne (N_mailles)
C=================================
        MAT_CALC_POT ( 1 , 1 ) = DFLOAT ( N_1_MAILLE )
        IF ( N_TYP_INTR .GT. 1 ) THEN
          MAT_CALC_POT ( 2 , 1 ) = DFLOAT ( N_2_MAILLE )
        END IF
        IF ( N_TYP_INTR .GT. 2 ) THEN
          MAT_CALC_POT ( 3 , 1 ) = DFLOAT ( N_3_MAILLE )
        END IF
        DO I_TYP = N_TYP_INTR + 1 , N_TYP 
          MAT_CALC_POT ( I_TYP , 1 ) = 0.D0
        END DO 
C=============================
C=====Deuxième colonne (n_d_1)
C=============================
C----------------------
C-----Elément ( 1 , 2 )
C----------------------
        IF ( I_TYP_R_DP_CONST ( 1 , 1 ) .EQ. 1 ) THEN
         IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .GT. N_R_1 ) THEN
           MAT_CALC_POT ( 1 , 2 ) = 1.D0
         END IF 
        ELSE
         IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .LE. N_R_1 ) THEN
           MAT_CALC_POT ( 1 , 2 ) = - 1.D0
         END IF
        END IF
C---------------------------------------------------------
C-----Elément ( 2 , 2 ) pour au moins 2 types intrinsèques
C---------------------------------------------------------
        IF ( N_TYP_INTR .GT. 1 ) THEN
          IF ( I_TYP_R_DP_CONST ( 1 , 1 ) .EQ. 2 ) THEN
            IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .LE. N_R_1 
     $      .OR. I_TYP_R_DP_CONST ( 1 , 2 ) .GT. N_R_1 + N_R_2 ) THEN
               MAT_CALC_POT ( 2 , 2 ) = 1.D0
            END IF
          ELSE
            IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .GT. N_R_1 
     $      .AND. I_TYP_R_DP_CONST ( 1 , 2 ) .LE. N_R_1 + N_R_2 ) THEN
               MAT_CALC_POT ( 2 , 2 ) = - 1.D0
            END IF
          END IF
        END IF
C------------------------------------------------
C-----Elément ( 3 , 2 ) pour 3 types intrinsèques
C------------------------------------------------
        IF ( N_TYP_INTR .GT. 2 ) THEN
          IF ( I_TYP_R_DP_CONST ( 1 , 1 ) .EQ. 3 ) THEN
            IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .LE. N_R_1 + N_R_2
     $      .OR. I_TYP_R_DP_CONST ( 1 , 2 ) .GT. N_R_1 + N_R_2 + N_R_3 )
     $      THEN
               MAT_CALC_POT ( 3 , 2 ) = 1.D0
            END IF
          ELSE
            IF ( I_TYP_R_DP_CONST ( 1 , 2 ) .GT. N_R_1 + N_R_2
     $     .AND. I_TYP_R_DP_CONST ( 1 , 2 ) .LE. N_R_1 + N_R_2 + N_R_3 ) 
     $      THEN
             MAT_CALC_POT ( 3 , 2 ) = - 1.D0
            END IF
          END IF
        END IF
C--------------------------------------------
C-----Eléments ( 2 , 2 ) et/ou ( 3 , 2 )
C-----pour moins de 2 ou 3 types intrinsèques  
C--------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
          IF ( I_TYP_R_DP_CONST ( 1 , 1 ) .EQ. I_TYP ) THEN
              MAT_CALC_POT ( I_TYP , 2 ) = 1.D0
          END IF
        END DO
C==============================
C=====Troisième colonne (n_d_2)
C==============================
C----------------------
C-----Elément ( 1 , 3 )
C----------------------
        IF ( I_TYP_R_DP_CONST ( 2 , 1 ) .EQ. 1 ) THEN
         IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .GT. N_R_1 ) THEN
           MAT_CALC_POT ( 1 , 3 ) = 1.D0
         END IF
        ELSE
         IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .LE. N_R_1 ) THEN
           MAT_CALC_POT ( 1 , 3 ) = - 1.D0
         END IF
        END IF
C---------------------------------------------------------
C-----Elément ( 2 , 3 ) pour au moins 2 types intrinsèques
C---------------------------------------------------------
        IF ( N_TYP_INTR .GT. 1 ) THEN
          IF ( I_TYP_R_DP_CONST ( 2 , 1 ) .EQ. 2 ) THEN
           IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .LE. N_R_1
     $      .OR. I_TYP_R_DP_CONST ( 2 , 2 ) .GT. N_R_1 + N_R_2 ) THEN
             MAT_CALC_POT ( 2 , 3 ) = 1.D0
           END IF
          ELSE
            IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .GT. N_R_1
     $      .AND. I_TYP_R_DP_CONST ( 2 , 2 ) .LE. N_R_1 + N_R_2 ) THEN
             MAT_CALC_POT ( 2 , 3 ) = - 1.D0
            END IF
          END IF
        END IF
C------------------------------------------------
C-----Elément ( 3 , 3 ) pour 3 types intrinsèques
C------------------------------------------------
        IF ( N_TYP_INTR .GT. 2 ) THEN
          IF ( I_TYP_R_DP_CONST ( 2 , 1 ) .EQ. 3 ) THEN
            IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .LE. N_R_1 + N_R_2
     $      .OR. I_TYP_R_DP_CONST ( 2 , 2 ) .GT. N_R_1 + N_R_2 + N_R_3 )
     $      THEN
             MAT_CALC_POT ( 3 , 3 ) = 1.D0
            END IF
          ELSE
            IF ( I_TYP_R_DP_CONST ( 2 , 2 ) .GT. N_R_1 + N_R_2
     $     .AND. I_TYP_R_DP_CONST ( 2 , 2 ) .LE. N_R_1 + N_R_2 + N_R_3 ) 
     $      THEN
             MAT_CALC_POT ( 3 , 3 ) = - 1.D0
            END IF
          END IF
        END IF
C--------------------------------------------
C-----Eléments ( 2 , 3 ) et/ou ( 3 , 3 )
C-----pour moins de 2 ou 3 types intrinsèques
C--------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
          IF ( I_TYP_R_DP_CONST ( 2 , 1 ) .EQ. I_TYP ) THEN
              MAT_CALC_POT ( I_TYP , 3 ) = 1.D0
          END IF
        END DO
  122  FORMAT ( 100 ( 2X , F12.5 ) )
C       write ( * , 1200 )
C       do i = 1 , n_typ
C	  write ( * , 122 ) ( MAT_CALC_POT ( i , j ) , j = 1 , n_typ ) 
C	end do
C----------------------------
C-----Inversion de la matrice
C----------------------------
         CALL
     $   INV_MAT_3_3
     $ ( MAT_CALC_POT , MAT_CALC_POT_INV )
C----------------------------------------------
C-----Calcul des potentiels chimiques à T = 0 K
C----------------------------------------------
       POT_CHIM_0K = 0.D0
       DO I_TYP = 1 , N_TYP
         POT_CHIM_0K ( I_TYP ) 
     $ = POT_CHIM_0K ( I_TYP )
     $ + MAT_CALC_POT_INV ( 1 , I_TYP )
     $ * H_REF_MAILLE
     $ + MAT_CALC_POT_INV ( 2 , I_TYP )
     $ * H_GC_D_R ( I_TYP_R_DP_CONST ( 1 , 1 ) ,
     $              I_TYP_R_DP_CONST ( 1 , 2 ) )
     $ + MAT_CALC_POT_INV ( 3 , I_TYP )
     $ * H_GC_D_R ( I_TYP_R_DP_CONST ( 2 , 1 ) ,
     $              I_TYP_R_DP_CONST ( 2 , 2 ) )
       END DO
C     write ( * , 122 ) ( POT_CHIM_0K ( i_typ ) , i_typ = 1 , n_typ )
C     write ( * , 1100 )
C======================================================
C=====Fin de calcul des potentiels chimiques de la zone
C======================================================
       END IF
C----------------------------------
C----------------------------------
C-----Ecriture pour le mode "point"
C----------------------------------
C----------------------------------
       IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $    ) THEN
  542 FORMAT ( 10X , 'N_mailles' , 3X , 100 ( A ) )
  543 FORMAT ( 100 ( F9.3 , 1X ) )
      WRITE ( 500 , * ) 'Tableau du simplexe initial :'
      WRITE ( 500 , 1200 )
      WRITE ( 500 , 542 )
     $ ( NOM_D_IND ( I ) , I = 1 , N_TYP_D_R )
      DO I_TYP = 1 , N_TYP + 2
          WRITE ( 500 , 543 )
     $  ( TAB_SMPLX_D_R ( I_TYP , I_D_R ) ,
     $    I_D_R = 1 , N_TYP_D_R + 2 )
      END DO
      WRITE ( 500 , 1200 )
      WRITE ( 500 , * ) 'Recherche des DP constitutionnels :'
      WRITE ( 500 , 1200 )
      IF ( IC .EQ. 1 ) THEN
       WRITE ( 500 , * ) "La fonction n'est pas bornée"
       WRITE ( 500 , 1200 )
      ELSE IF ( IC .EQ. - 1 ) THEN
       WRITE ( 500 , * ) 'Aucune solution avec ces contraintes'
       WRITE ( 500 , 1200 )
      ELSE IF ( IC .EQ. 0 ) THEN
       WRITE ( 500 , * ) 'Une solution trouvée :'
       WRITE ( 500 , 1200 ) 
       WRITE ( 500 , * ) 'Numéros des variables non nulles :'
C-------------------------------------------------------------
C-----Il y a autant de variables non nulles que de contraintes
C-----(et donc que de types chimiques)
C-------------------------------------------------------------
  700 FORMAT ( 10X , '* ' , A )
       WRITE ( 500 , 1200 )
       WRITE ( 500 , * ) ( I_2_SMPLX ( I_TYP ) , I_TYP = 1 , N_TYP )
       WRITE ( 500 , 1200 )
       WRITE ( 500 , '(A)' ) '  Leurs types :'
       DO I_TYP = 1 , N_TYP
         IF ( I_2_SMPLX ( I_TYP ) .EQ. 1 ) THEN
           WRITE ( 500 , 700 ) NOM_D_IND ( I_2_SMPLX ( I_TYP ) - 1 )
         END IF
       END DO
       DO I_TYP = 1 , N_TYP
         IF ( I_2_SMPLX ( I_TYP ) .NE. 1 ) THEN
           WRITE ( 500 , 700 ) 
     $      NOM_D_IND ( I_2_SMPLX ( I_TYP ) - 1 )
     $    ( 1 : LONG_NOM_D_IND ( I_2_SMPLX ( I_TYP ) - 1 ) )
         END IF 
       END DO
       WRITE ( 500 , 1200 ) 
      END IF
C----------------------------------------
C----------------------------------------
C-----Ecriture pour le mode "balayage"
C-----(seulement si composition correcte)
C----------------------------------------
C----------------------------------------
       ELSE
  437  FORMAT ( 5 ( 2X , F12.7 ) , 1X , '|' , 1X ,
     $          2 ( 2X , A ) ,  3 ( 2X , F10.5 ) )
       IF ( I_COMPO .EQ. 1 ) THEN
        IF ( IC .EQ. 0 ) THEN
         IF ( INDIC_LIMITE_ZONE .EQ. 1 ) THEN
          WRITE ( 500 , 437 ) X_COORD_1 , X_COORD_2 , 
     $                     ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP ) ,
     $    ( NOM_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 )
     $ ( 1 : LONG_NOM_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 ) ) ,
     $      I_TYP = 2 , N_TYP ) ,
     $    ( POT_CHIM_0K ( I_TYP ) , I_TYP = 1 , N_TYP )
C       do i_typ = 1 , n_typ - 1
C        write ( 500 , * ) 
C    $ ( I_TYP_R_DP_CONST ( I_TYP , I ) , I = 1 , 2 )
C       end do
C       write ( 500 , 1200 )
C       do i = 1 , n_typ
C         write ( 500 , 122 ) ( MAT_CALC_POT ( i , j ) , j = 1 , n_typ )
C       end do
         END IF
C----------------------------------------------------------------
C-----Cas d'écriture du point de composition à  "erreur simplexe"
C-----dans un fichier annexe
C----------------------------------------------------------------
        ELSE
         IF ( CAR_FICH_550 .EQ. 'o' ) THEN
          WRITE ( 550 , 500 ) ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
         END IF
        END IF
C-------------------------------
C-----Fin du test de composition
C-------------------------------
       END IF
C---------------------------------
C---------------------------------
C-----Fin du test "point/balayage"
C---------------------------------
C---------------------------------
       END IF
C----------------------------------------------------
C-----Au dernier pas en "r",
C-----écriture des résultats dans le fichier de liste
C----------------------------------------------------
       IF ( I_PAS_R .EQ. N_PAS_R ) THEN
        IF ( IC .EQ. 0 ) THEN
         IF ( INDIC_LIMITE_ZONE .EQ. 1 ) THEN
         WRITE ( 100 , 437 ) X_COORD_1 , X_COORD_2 ,
     $                     ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP ) ,
     $    ( NOM_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 )
     $ ( 1 : LONG_NOM_D_IND ( I_2_SMPLX ( IND_INIT ( I_TYP ) ) - 1 ) ) ,
     $      I_TYP = 2 , N_TYP ) ,
     $    ( POT_CHIM_0K ( I_TYP ) , I_TYP = 1 , N_TYP )
         END IF
        END IF
       END IF
C##################################################
C#####Fin de boucle sur les pas de balayage "theta"
C##################################################
        END DO
C##############################################
C#####Fin de boucle sur les pas de balayage "r"
C##############################################
       END DO
C-------------------------------------------------------------        
C-----Ecriture du compte-rendu d'erreur éventuelle du simplexe
C-------------------------------------------------------------        
       IF ( I_ERREUR_SIMPLEXE .NE. 0 ) THEN
         WRITE ( 100 , 1200 )
         WRITE ( 100 , * )
     $ 'ATTENTION : lors de ce calcul, pour certaines compositions :'
         WRITE ( 100 , * )
     $ 'la minimisation de H était impossible :'
         WRITE ( 100 , * )
     $ "(i) soit la fonction n'est pas bornée"
         WRITE ( 100 , * )
     $ '(ii) soit aucune solution avec ces contraintes'
         WRITE ( 100 , 1200 )
       END IF
C############################
C############################
C#####Fin du traitement NPT=0
C############################
C############################
      END IF
C################################################
C################################################
C#####Traitement du cas NPT (DP et énergie libre)
C################################################
C################################################
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
C-----------------------------------
C-----Nombre d'inconnues du problème
C-----------------------------------
       N_NPT = N_TYP_D_R + 1
       write ( * , * )
     $ "Dimension du système d'équations NPT = " , N_NPT
       write ( * , * ) 
     $ '------------------------------------'
C--------------------------------------------
C-----Vecteur, fonction et matrice jacobienne
C-----du problème non linéaire
C--------------------------------------------
       ALLOCATE ( X_NPT ( 1 : N_NPT ) )
       ALLOCATE ( F_NPT ( 1 : N_NPT ) ) 
       ALLOCATE ( J_NPT ( 1 : N_NPT , 1 : N_NPT ) ) 
C-------------------------------------------------
C-----Sauvegarde du vecteur initial à un indice
C-----(seulement pour comparaison au vecteur final
C----- --> critère de "convergence" NRCG)
C-------------------------------------------------
       ALLOCATE ( X_NPT_INIT ( 1 : N_NPT ) )
C- - - - - - - - - - - - - - - - - - - - - - - -
C- - -Le vecteur initial est conservé
C- - -(utile seulement pour test de convergence)
C- - - - - - - - - - - - - - - - - - - - - - - -
      DO I_R = 1 , N_R
       DO I_TYP = 0 , N_TYP
        IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
           X_NPT_INIT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $   = LOG_X_D_R_INIT ( I_TYP , I_R )
        END IF
       END DO
      END DO     
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Variable "nombre de mailles" (format = réel)
C- - -Rappel : N_MAILLES_INIT, qui ne joue pas sur les résultats,
C- - -(xDP = variables intensives) a été initialisé à 1.D0
C- - -au début du traitement NPT
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      X_NPT_INIT ( N_NPT ) = N_MAILLES_INIT
C---------------------------------------------------------------------
C-----Nombre initial d'atomes par maille (entier),
C-----dont est déduit ci-dessous N_AT_TOT (fonction de N_MAILLES)
C-----N_AT_TOT = variable passée en argument de la procédure SYS_NR
C-----Rappel : N_AT_MAILLE_INIT est une constante,
C-----comme N_1_MAILLE, N_2_MAILLE, N_3_MAILLE calculés dès la lecture
C-----du fichier de données, et relatifs au système sans DP.
C---------------------------------------------------------------------
       N_AT_MAILLE_INIT = N_1_MAILLE
     $                  + N_2_MAILLE
     $                  + N_3_MAILLE
C---------------------------------------------------------
C-----Nombre total d'atomes (format = réel)
C-----= variable passée en argument de la procédure SYS_NR
C---------------------------------------------------------
      N_AT_TOT = N_MAILLES_INIT * DFLOAT ( N_AT_MAILLE_INIT )
C-----------------------------------------------------------------------
C-----Affectation des paramètres lus 
C-----(nombres de pas et valeurs extrémales pour les boucles sur T et x)
C-----aux variables des boucles interne et externe.
C-----Dans le cas "x" ou "T", c'est la boucle interne qui est parcourue.
C-----------------------------------------------------------------------
      IF ( INDIC_TYP_CALC_NPT
     $   ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' ) THEN
        N_PAS_BOUCLE_EXTERNE = 1
        N_PAS_BOUCLE_INTERNE = N_PAS_X_TYP_0
        PAS_BOUCLE_EXTERNE = 0.D0
        PAS_BOUCLE_INTERNE = ( X_TYP_0_FIN - X_TYP_0_INIT ) 
     $                     / DFLOAT ( N_PAS_X_TYP_0 )
       ELSE IF ( INDIC_TYP_CALC_NPT
     $         ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T' ) THEN
        N_PAS_BOUCLE_EXTERNE = 1
        N_PAS_BOUCLE_INTERNE = N_PAS_T
        PAS_BOUCLE_EXTERNE = 0.D0
        PAS_BOUCLE_INTERNE = ( T_FIN - T_INIT )
     $                     / DFLOAT ( N_PAS_T )
       ELSE IF ( INDIC_TYP_CALC_NPT
     $         ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
        N_PAS_BOUCLE_EXTERNE = N_PAS_T
        N_PAS_BOUCLE_INTERNE = N_PAS_X_TYP_0
        PAS_BOUCLE_EXTERNE = ( T_FIN - T_INIT )
     $                     / DFLOAT ( N_PAS_T )
        PAS_BOUCLE_INTERNE = ( X_TYP_0_FIN - X_TYP_0_INIT )
     $                     / DFLOAT ( N_PAS_X_TYP_0 )
       ELSE IF ( INDIC_TYP_CALC_NPT
     $         ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
        N_PAS_BOUCLE_EXTERNE = N_PAS_X_TYP_0
        N_PAS_BOUCLE_INTERNE = N_PAS_T
        PAS_BOUCLE_EXTERNE = ( X_TYP_0_FIN - X_TYP_0_INIT )
     $                     / DFLOAT ( N_PAS_X_TYP_0 )
        PAS_BOUCLE_INTERNE = ( T_FIN - T_INIT )
     $                     / DFLOAT ( N_PAS_T )
       ELSE IF ( INDIC_TYP_CALC_NPT
     $         ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p' 
     $      .OR. INDIC_TYP_CALC_NPT
     $         ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P' ) THEN
        N_PAS_BOUCLE_EXTERNE = 1
        N_PAS_BOUCLE_INTERNE = 1
        PAS_BOUCLE_EXTERNE = 0.D0
        PAS_BOUCLE_INTERNE = 0.D0
       END IF
C--------------------------------------------------------------------
C-----Ouverture des tableaux de conservation des résultats
C-----d'une itération sur l'autre
C-----(pour affectation des valeurs initiales).
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Ce tableau n'est utile qu'en NPT(Tx).
C-----La boucle interne (premier indice du tableau) contient alors
C-----les val. init. aux diverses compo. à la température précédente.
C--------------------------------------------------------------------
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
         ALLOCATE ( X_NPT_COMPO_SVG ( N_PAS_BOUCLE_INTERNE , N_NPT ) )
       END IF
C=====================================================================
C=====Initialisation préalable du vecteur d'inconnues (x_d ; M)  
C---------------------------------------------------------------------
C=====En NPT(T/x/Tx/xT), l'initialisation est d'abord faite ici
C=====(i.e. avant la double boucle sur T et x)
C=====à l'aide des valeurs lues dans DATA.adpi.
C=====* En NPT(T) : cette initialisation concerne la T initiale.
C=====  Ensuite, initialisation "de proche en proche"
C===== (--> résultat du point T(N) = amorce du point T(N+1))
C=====* En NPT(x/Tx/xT) : cette initialisation de proche en proche
C=====  pourra être supplantée par une autre initialisation
C=====  effectuée  à l'intérieur des boucles sur T et/ou x :
C=====    soit à l'aide des valeurs de DATA.adpi,
C=====    soit via un fichier de valeurs initiales (f(x)) (cf. infra).
C=====   [f(x) signifie que ces val. init. dépendent de la compo.]
C=====================================================================
C----------------------------------------------------------------
C-----Variable "nombre de mailles" (format = réel)
C-----Rappel : N_MAILLES_INIT, qui ne joue pas sur les résultats,
C-----(xDP = variables intensives) a été initialisé à 1.D0
C-----au début du traitement NPT
C----------------------------------------------------------------
      X_NPT ( N_NPT ) = N_MAILLES_INIT
C-----------------------------
C-----Variables "log_10(x_DP)"
C-----------------------------
      DO I_R = 1 , N_R
        DO I_TYP = 0 , N_TYP
C- - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C- - -ne correspondent à aucun DP ( I_TYP , I_R )
C- - -ou à un DP non utile (E_DP > 0 ou = 0)
C- - -et ne sont donc pas initialisées.
C- - - - - - - - - - - - - - - - - - - - - - - - 
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
                X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $        = LOG_X_D_R_INIT ( I_TYP , I_R )
          END IF
        END DO
      END DO
C===========================================================
C=====Fin de l'initialisation préalable du vecteur (x_d ; M)  
C===========================================================
C---------------------------------------------------------------------
C-----Fichier où seront écrites pour chaque composition
C-----les valeurs convergées des variables pour le calcul NPT en cours
C----- --> fichier ouvert dans tous les cas NPT(P/p/T/x/Tx/xT)
C---------------------------------------------------------------------
          OPEN ( 210 , FILE = FICH_FIN ( 1 : LONG_FICH_FIN )
     $                     // '_adpi_'
     $                     // INDIC_TYP_CALC
     $                      ( 1 : LONG_INDIC_TYP_CALC )
     $                     // '_'
     $                     // 'varNPT' )
C-----------------------------------------------------------------------
C-----Seulement dans les cas NPT(x/Tx/xT)
C-----et si l'option "fich. val. init." est choisie :
C-----ouverture et lecture du fichier de valeurs initiales
C-----(dont le contenu est placé dans X_NPT_INIT_FICH),
C-----puis ouverture du fichier où sera écrit le "suivi"
C-----de l'utilisation de ce tableau (le numéro de ligne pour chaque x).
C-----------------------------------------------------------------------
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' 
     $  ) THEN
        IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
         ALLOCATE ( X_NPT_INIT_FICH ( N_LIGNES_FICH_VAL_INIT , N_NPT ) )
         OPEN ( 310 , FILE = FICH_VAL_INIT ( 1 : LONG_FICH_VAL_INIT ) )
         DO I_VAL_INIT = 1 , N_LIGNES_FICH_VAL_INIT
         READ ( 310 , * ) ITER_VERIF ,
     $   ( X_NPT_INIT_FICH ( I_VAL_INIT , I_NPT ) , I_NPT = 1 , N_NPT ) 
         END DO
         OPEN ( 410 , FILE = "fichier_NPT_suivi_lect_val_init" )
        ENDIF
       END IF
C########################################################
C#####Double boucle sur la température
C#####et sur la fraction atomique pour l'élément spécifié
C########################################################
      I_POUR_CENT = 0
      I_COMPTEUR_PAS = 0
C# # # # # # # # # #
C# # #Boucle externe
C# # # # # # # # # #
      DO I_PAS_BOUCLE_EXTERNE = 1 , N_PAS_BOUCLE_EXTERNE
C       write(*,*) "I_PAS_BOUCLE_EXTERNE=",I_PAS_BOUCLE_EXTERNE
C--------------------------------------------------------------
C-----En NPT(xT) seul : à chaque début de composition,
C-----réinitialisation optionnelle avec les valeurs de DATA.adpi
C-----(optionnelle car seulement si INDIC_LECT_VAL_INIT_XT = 1)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Si INDIC_LECT_VAL_INIT_NPT = 1,
C-----cette réinitialisation sera supplantée plus loin
C-----à l'aide des valeurs lues dans le fichier de val. init. :
C-----à l'intérieur des boucles, juste avant la résolution,
C-----et toujours en début de nouvelle compo. c'est-à-dire
C-----pour la première température (I_PAS_BOUCLE_INTERNE = 1).
C--------------------------------------------------------------
          IF ( INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
           IF ( INDIC_LECT_VAL_INIT_XT .EQ. 1 ) THEN
            X_NPT ( N_NPT ) = N_MAILLES_INIT
            DO I_R = 1 , N_R
             DO I_TYP = 0 , N_TYP
               IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
                     X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $             = LOG_X_D_R_INIT ( I_TYP , I_R )
               END IF
             END DO
            END DO
           END IF
          END IF
C# # # # # # # # # #
C# # #Boucle interne
C# # # # # # # # # #
       DO I_PAS_BOUCLE_INTERNE = 1 , N_PAS_BOUCLE_INTERNE
C       write(*,*) "I_PAS_BOUCLE_INTERNE=",I_PAS_BOUCLE_INTERNE
        I_COMPTEUR_PAS = I_COMPTEUR_PAS + 1
C===================================================
C===================================================
C=====Opérations préliminaires à la résolution
C=====et à faire pour chaque pas en composition et T
C===================================================
C===================================================
C=======================================================================
C=====1) Affectation de la composition courante :
C=====* cas NPT(x/Tx/xT) : calcul des fractions atomiques courantes X_AT
C=====à partir de la valeur courante de X_AT ( I_TYP_0 )
C=====et des contraintes en composition
C=====* cas NPT(p/P/T) : la composition X_AT n'est pas modifiée
C=====et garde la valeur lue dans DATA.adpi
C=======================================================================
C= = = = = = = = = = = = 
C= = =* Cas NPT(x/Tx/xT)
C= = = = = = = = = = = = 
       IF ( INDIC_TYP_CALC_NPT 
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' 
     $ .OR. INDIC_TYP_CALC_NPT 
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' 
     $ .OR. INDIC_TYP_CALC_NPT 
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $    ) THEN
        X_AT = 0.D0
C-----Cas x/Tx : le balayage en x est traité par la boucle interne
          IF ( INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $    .OR. INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
          X_AT_0_COUR = X_TYP_0_INIT
     $                + DFLOAT ( I_PAS_BOUCLE_INTERNE )
     $                * PAS_BOUCLE_INTERNE
C-----Cas xT : le balayage en x est traité par la boucle externe
          ELSE
     $    IF ( INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
          X_AT_0_COUR = X_TYP_0_INIT
     $                + DFLOAT ( I_PAS_BOUCLE_EXTERNE )
     $                * PAS_BOUCLE_EXTERNE
          END IF
          X_AT ( I_TYP_0 ) = X_AT_0_COUR
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =En mode NPT(x/Tx/xT), dans les cas ternaires ou plus,
C= = =à partir des contraintes sur la composition
C= = =et de la fraction atomique courante
C= = =pour l'élément spécifié I_TYP_0,
C= = =inversion LU du système d'équations de contraintes
C= = =pour obtenir les valeurs des autres fractions atomiques
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      IF ( N_TYP .GT. 2 ) THEN
C----------------------------------------------------
C-----Boucle sur les lignes du système de contraintes
C----------------------------------------------------
       DO I_TYP = 1 , N_SYS_CTR
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Boucle sur les colonnes du système de contraintes
C- - -en sautant la colonne I_TYP_0
C- - -(élément dont on étudie l'effet de l'enrichissement)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         I_CTR = 0
         DO J_TYP = 1 , I_TYP_0 - 1
           I_CTR = I_CTR + 1
           MAT_CTR ( I_TYP , I_CTR )
     $   = COEF_CTR ( I_TYP , J_TYP )
         END DO
         DO J_TYP = I_TYP_0 + 1 , N_TYP
           I_CTR = I_CTR + 1
           MAT_CTR ( I_TYP , I_CTR )
     $   = COEF_CTR ( I_TYP , J_TYP )
         END DO
         IF ( I_TYP .LT. N_TYP - 1 ) THEN
           VECT_CTR ( I_TYP ) = - COEF_CTR ( I_TYP , N_TYP + 1 )
     $                          - COEF_CTR ( I_TYP , I_TYP_0 )
     $                          * X_AT ( I_TYP_0 )
         ELSE
           VECT_CTR ( I_TYP ) = 1.D0 - X_AT ( I_TYP_0 )
         END IF
C-----------------------------------------------------------
C-----Fin de boucle sur les lignes du système de contraintes
C-----------------------------------------------------------
      END DO
C---------------------------------------------------------
C-----Inversion LU de la matrice du système de contraintes
C---------------------------------------------------------
         CALL
     $   FACT_LU
     $ ( N_SYS_CTR ,
     $   MAT_CTR ,
     $   P_MAT_CTR , Q_MAT_CTR ,
     $   L_MAT_CTR , U_MAT_CTR )
         CALL
     $   INVERSE_LU
     $ ( N_SYS_CTR ,
     $   L_MAT_CTR , U_MAT_CTR ,
     $   P_MAT_CTR , Q_MAT_CTR ,
     $   MAT_INV_CTR )
C------------------------------------------------------
C-----Fractions atomiques correspondant aux contraintes
C-----et à la valeur courante de X_AT ( I_TYP_0 )
C------------------------------------------------------
      I_CTR = 0
      DO I_TYP = 1 , I_TYP_0 - 1
       I_CTR = I_CTR + 1
       DO J_TYP = 1 , N_TYP - 1 
         X_AT ( I_TYP )
     $ = X_AT ( I_TYP )
     $ + MAT_INV_CTR ( I_CTR , J_TYP ) * VECT_CTR ( J_TYP ) 
        END DO
      END DO  
      DO I_TYP = I_TYP_0 + 1 , N_TYP
       I_CTR = I_CTR + 1
       DO J_TYP = 1 , N_TYP - 1
         X_AT ( I_TYP )
     $ = X_AT ( I_TYP )
     $ + MAT_INV_CTR ( I_CTR , J_TYP ) * VECT_CTR ( J_TYP ) 
        END DO
       END DO
C     write ( * , * ) '-------'
C     write ( * , 1110 ) x_at
C     write ( * , * ) '-------'
C------------------------------------
C-----Cas de l'alliage binaire
C-----(pas de contrainte dans ce cas)
C------------------------------------
       ELSE
        X_AT ( 3 - I_TYP_0 ) = 1.D0 - X_AT ( I_TYP_0 )
C---------------------------------------------
C-----Fin du test sur binaire/ternaire ou plus
C---------------------------------------------
       END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = =Fin de l'obtention des autres fractions atomiques
C= = =pour les cas NPT(x/Tx/xT)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = = = = = = = = = =
C= = =* Cas NPT(p/P/T)
C= = = = = = = = = = =
       ELSE
     $ IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P' ) THEN
C-----Dans ce cas, rien n'est fait :
C-----la composition garde la valeur lue dans DATA.adpi
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = =Fin du test "cas NPT(x/Tx/xT) ou NPT(p/P/T) ?"
C= = =pour l'affectation de la composition courante
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
      END IF
C=======================================================
C=====Fin de l'affectation 1) de la composition courante
C=======================================================
C=======================================================================
C=====2) Affectation de la température courante :
C=====* cas NPT(T/Tx/xT) : calcul de la valeur courante
C=====* cas NPT(p/P/x) : on garde la valeur constante lue dans DATA.adpi
C=======================================================================
C= = = = = = = = = = = = 
C= = =* Cas NPT(T/Tx/xT)
C= = = = = = = = = = = = 
C-----Cas où le balayage en T est traité par la boucle interne
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $  .OR. INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
        TEMPERATURE = T_INIT
     $              + DFLOAT ( I_PAS_BOUCLE_INTERNE )
     $              * PAS_BOUCLE_INTERNE
         K_T = K_B * TEMPERATURE
C-----Cas où le balayage en T est traité par la boucle externe
        ELSE
     $  IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
        TEMPERATURE = T_INIT
     $              + DFLOAT ( I_PAS_BOUCLE_EXTERNE )
     $              * PAS_BOUCLE_EXTERNE
         K_T = K_B * TEMPERATURE
C= = = = = = = = = = = = = = = = = = = = =
C= = =Fin de l'obtention de la température
C= = =pour les cas NPT(T/Tx/xT)
C= = = = = = = = = = = = = = = = = = = = =
C= = = = = = = = = = = 
C= = =* Cas NPT(p/P/x)
C= = = = = = = = = = = 
       ELSE
     $ IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' 
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p' 
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P' ) THEN
C-----Dans ce cas, rien n'est fait :
C-----la température garde la valeur lue dans DATA.adpi
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = =Fin du test "cas NPT(T/Tx/xT) ou NPT(p/P/x) ?"
C= = =pour l'affectation de la température courante
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
      END IF
C=======================================================
C=====Fin de l'affectation 2) de la température courante
C=======================================================
C=====================================================================
C=====Si l'option INDIC_LECT_VAL_INIT_NPT = 1 est choisie,
C=====et donc dans les cas NPT(x/Tx/xT) seulement
C=====[puisque le programme s'arrête si
C=====INDIC_LECT_VAL_INIT_NPT = 1 et NPT(p/P/T)],
C=====utilisation du fichier de valeurs initiales
C===== --> prend le pas sur l'initialisation par DATA.adpi
C=====réalisée plus haut en dehors des boucles sur T et/ou x.
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
C=====Cette opération a lieu (ssi INDIC_LECT_VAL_INIT_NPT = 1) :
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(x) à tous les pas en compo. (avec val. init. f(x))
C=====   [f(x) signifie que ces val. init. dépendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(Tx) : à la première température seulement
C=====               (--> premier pas de boucle externe)
C=====                pour toutes les compositions
C=====               (--> tous les pas de boucle interne)
C=====               (avec la bonne ligne du fichier).
C=====             Aux températures suivantes T(i),
C=====             les valeurs initiales sont celles (f(x)) à T(i-1)
C=====             conservées dans le tableau X_NPT_COMPO_SVG.
C=====   [f(x) signifie que ces val. init. dépendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(xT) : à chaque composition
C=====               (--> chaque pas de boucle externe)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C=====Rappel : si INDIC_LECT_VAL_INIT_NPT = 0 :
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(T/x), l'initialisation se fait à l'aide de DATA.adpi
C=====à la première T ou x, puis de proche en proche aux autres T ou x.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(Tx), l'initialisation se fait
C=====(i) à la première T : par DATA.adpi (quelle que soit la compo.)
C=====(ii) aux autres T(i) : à l'aide des val. init. (=f(x))
C=====     conservées dans le tableau X_NPT_COMPO_SVG.
C=====(i.e. le schéma d'initialisation "Tx" est similaire dans les cas
C=====INDIC_LECT_VAL_INIT_NPT = 0 et INDIC_LECT_VAL_INIT_NPT = 1,
C=====en remplaçant, à la première, T "DATA.adpi" par "fich. val. init")
C=====   [f(x) signifie que ces val. init. dépendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(xT), l'initialisation est effectuée
C=====(i) à l'aide des valeurs de DATA.adpi au début de chaque compo.
C=====(ii) de proche en proche à chaque température dans cette compo.
C=====================================================================
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = =Premier cas : NPT(x) -> la boucle sur x est la boucle interne
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
       IF ( 
     $      INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $    ) THEN
C----------------------------------------------
C-----Si INDIC_LECT_VAL_INIT_NPT = 1 et NPT(x),
C-----affectation de l'indice de ligne I_COUR
C-----dans le fichier de valeurs initiales
C----------------------------------------------
        IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Instructions permettant un changement de pas en composition
C- - -entre le fichier de val. init. et le calcul en cours
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Rappel : le nombre de pas et le nombre de lignes
C- - -doivent être multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arrête (test à la lecture 
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autorisés sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilisée successivement
C- - -     pour N_FRAC points en composition consécutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 100
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> seule une ligne toutes les N_FRAC est lue dans le fichier. 
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Ainsi, un cas tel que N_PAS_X_TYP_0 = 300
C- - -et N_LIGNES_FICH_VAL_INIT = 200 est interdit.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( N_PAS_BOUCLE_INTERNE .GE. N_LIGNES_FICH_VAL_INIT ) THEN
               N_FRAC = N_PAS_BOUCLE_INTERNE / N_LIGNES_FICH_VAL_INIT 
               I_COUR = ( I_PAS_BOUCLE_INTERNE - 1 ) / N_FRAC + 1
          ELSE
               N_FRAC = N_LIGNES_FICH_VAL_INIT / N_PAS_BOUCLE_INTERNE
               I_COUR = 1 + ( I_PAS_BOUCLE_INTERNE - 1 ) * N_FRAC
          END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Test supplémentaire de précaution pour éviter
C- - -d'éventuels dépassements dans le numéro de ligne
C- - -(même si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui précèdent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C----------------------------------------------------------------------
C-----Pour NPT(x) et si INDIC_LECT_VAL_INIT_NPT = 1, initialisation
C-----à l'aide des valeurs lues dans le fichier de valeurs initiales
C-----(prend alors le pas sur l'initialisation faite plus haut,
C-----en dehors de la boucle sur T ou x, à l'aide des valeurs initiales
C-----lues dans DATA.adpi)
C----------------------------------------------------------------------
          DO I_NPT = 1 , N_NPT
             X_NPT ( I_NPT ) = X_NPT_INIT_FICH ( I_COUR , I_NPT )
          END DO
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fin du double test pour vérifier que l'on est en mode NPT(x)
C= = =et que INDIC_LECT_VAL_INIT_NPT = 1
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        END IF
      END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Deuxième cas : NPT(Tx) -> la boucle sur x est la boucle interne
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Le fichier de val. init. n'est utilisé qu'à la première T.
C= = =Aux T suivantes, on utilise un tableau conservant les résultats
C= = =pour les diverses compositions à la T précédente.
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $    ) THEN
        IF ( I_PAS_BOUCLE_EXTERNE .EQ. 1 ) THEN
C----------------------------------------------------------------------
C-----Si NPT(Tx) et s'il s'agit du premier pas en T, réinitialisation :
C-----* si INDIC_LECT_VAL_INIT_NPT = 1 :
C-----par utilisation du fichier de valeurs initiales (= f(x))
C-----via l'indice de ligne I_COUR fonction de I_PAS_BOUCLE_INTERNE,
C-----* si INDIC_LECT_VAL_INIT_NPT = 0 :
C-----par utilisation de DATA.adpi (-> val. init. indépendantes de x).
C----------------------------------------------------------------------
C- - - - - - - - - - - - - - - - - - - - - - - - -
C- - - - - - - - - - - - - - - - - - - - - - - - -
C- - -1) Cas "affectation par fichier val. init."
C- - - - - - - - - - - - - - - - - - - - - - - - -
C- - - - - - - - - - - - - - - - - - - - - - - - -
         IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Instructions permettant un changement de pas en composition
C- - -entre le fichier de val. init. et le calcul en cours
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Rappel : le nombre de pas et le nombre de lignes
C- - -doivent être multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arrête (test à la lecture
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autorisés sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilisée successivement
C- - -     pour N_FRAC points en composition consécutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 100
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> seule une ligne toutes les N_FRAC est lue dans le fichier.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Ainsi, un cas tel que N_PAS_X_TYP_0 = 300
C- - -et N_LIGNES_FICH_VAL_INIT = 200 est interdit.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( N_PAS_BOUCLE_INTERNE .GE. N_LIGNES_FICH_VAL_INIT ) THEN
               N_FRAC = N_PAS_BOUCLE_INTERNE / N_LIGNES_FICH_VAL_INIT
               I_COUR = ( I_PAS_BOUCLE_INTERNE - 1 ) / N_FRAC + 1
          ELSE
               N_FRAC = N_LIGNES_FICH_VAL_INIT / N_PAS_BOUCLE_INTERNE
               I_COUR = 1 + ( I_PAS_BOUCLE_INTERNE - 1 ) * N_FRAC
          END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Test supplémentaire de précaution pour éviter
C- - -d'éventuels dépassements dans le numéro de ligne
C- - -(même si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui précèdent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x/xT/Tx) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x/xT/Tx) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -NPT(Tx) : au premier pas en T, si INDIC_LECT_VAL_INIT_NPT = 1,
C- - -initialisation (f(x)) à l'aide du fichier de valeurs initiales
C- - -via I_COUR déduit de I_PAS_BOUCLE_INTERNE
C- - - --> initialisation (f(x)) qui a lieu à toutes les compo.
C- - -   [f(x) signifie que ces val. init. dépendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO I_NPT = 1 , N_NPT
             X_NPT ( I_NPT ) = X_NPT_INIT_FICH ( I_COUR , I_NPT )
          END DO
C- - - - - - - - - - - - - - - - - - - - 
C- - - - - - - - - - - - - - - - - - - - 
C- - -2) Cas "affectation par DATA.adpi"
C- - - - - - - - - - - - - - - - - - - - 
C- - - - - - - - - - - - - - - - - - - - 
         ELSE
C- - - - - - - - - - - - - - - - -
C- - -Variable "nombre de mailles"
C- - - - - - - - - - - - - - - - -
           X_NPT ( N_NPT ) = N_MAILLES_INIT
C- - - - - - - - - - - - - - -
C- - -Variables "log_10(x_DP)"
C- - - - - - - - - - - - - - -
           DO I_R = 1 , N_R
            DO I_TYP = 0 , N_TYP
C- - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C- - -ne correspondent à aucun DP ( I_TYP , I_R )
C- - -ou à un DP non utile (E_DP > 0 ou = 0)
C- - -et ne sont donc pas initialisées.
C- - - - - - - - - - - - - - - - - - - - - - - - 
              IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
                    X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $            = LOG_X_D_R_INIT ( I_TYP , I_R )
              END IF
            END DO
           END DO
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Fin du test "INDIC_LECT_VAL_INIT_NPT = 0 ou 1"
C- - -(pour l'initialisation "Tx" au premier pas en T)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
         END IF
C--------------------------------------------------------
C-----En NPT(Tx), s'il ne s'agit pas du premier pas en T,
C-----utilisation des valeurs de T précédentes,
C-----conservées dans  X_NPT_COMPO_SVG
C--------------------------------------------------------
        ELSE
         DO I = 1 , N_NPT
          X_NPT ( I ) = X_NPT_COMPO_SVG ( I_PAS_BOUCLE_INTERNE , I )
         END DO
C-------------------------------------
C-----Fin du test "premier pas en T ?"
C-------------------------------------
        END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fin du test pour vérifier que l'on est en mode NPT(Tx)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Troisième cas : NPT(xT) -> la boucle sur x est la boucle externe
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Le fichier de val. init. est utilisé à chaque compo.
C= = =(i.e. à chaque pas de boucle externe),
C= = =avant le démarrage de la boucle en T
C= = =(i.e. au premier pas de la boucle interne)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = = Rque : cette réinitialisation "fich. val. init" propre à "xT"
C= = =n'étant effectuée qu'au premier pas de boucle interne,
C= = =elle pourrait donc être extraite de la boucle interne
C= = =et regroupée avec l'initialisation "xT" par DATA.adpi
C= = =située plus haut entre les deux boucles.
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $    ) THEN
C----------------------------------------------------------
C-----Si NPT(xT) et INDIC_LECT_VAL_INIT_NPT = 1,
C-----et s'il s'agit du premier pas en T (boucle interne) :
C-----affectation de l'indice de ligne I_COUR
C-----dans le fichier de valeurs initiales
C----------------------------------------------------------
        IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
         IF ( I_PAS_BOUCLE_INTERNE .EQ. 1 ) THEN
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Instructions permettant un changement de pas en composition
C- - -entre le fichier de val. init. et le calcul en cours
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Rappel : le nombre de pas et le nombre de lignes
C- - -doivent être multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arrête (test à la lecture
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autorisés sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilisée successivement
C- - -     pour N_FRAC points en composition consécutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou égal à)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 100
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> seule une ligne toutes les N_FRAC est lue dans le fichier.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Ainsi, un cas tel que N_PAS_X_TYP_0 = 300
C- - - et N_LIGNES_FICH_VAL_INIT = 200 est interdit.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( N_PAS_BOUCLE_EXTERNE .GE. N_LIGNES_FICH_VAL_INIT ) THEN
               N_FRAC = N_PAS_BOUCLE_EXTERNE / N_LIGNES_FICH_VAL_INIT
               I_COUR = ( I_PAS_BOUCLE_EXTERNE - 1 ) / N_FRAC + 1
          ELSE
               N_FRAC = N_LIGNES_FICH_VAL_INIT / N_PAS_BOUCLE_EXTERNE
               I_COUR = 1 + ( I_PAS_BOUCLE_EXTERNE - 1 ) * N_FRAC
          END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Test supplémentaire de précaution pour éviter
C- - -d'éventuels dépassements dans le numéro de ligne
C- - -(même si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui précèdent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x/xT/Tx) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x/xT/Tx) - itér. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C----------------------------------------------------------------------
C-----Pour NPT(xT) et si INDIC_LECT_VAL_INIT_NPT = 1, initialisation
C-----à l'aide des valeurs lues dans le fichier de valeurs initiales
C-----(prend alors le pas sur l'initialisation faite plus haut,
C-----en dehors de la boucle sur T ou x, à l'aide des valeurs initiales
C-----lues dans DATA.adpi)
C----------------------------------------------------------------------
          DO I_NPT = 1 , N_NPT
             X_NPT ( I_NPT ) = X_NPT_INIT_FICH ( I_COUR , I_NPT )
          END DO
C- - -Fin du test pour vérifier que c'est le premier pas en T
         END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fin du double test pour vérifier que l'on est en mode NPT(xT)
C= = =et que INDIC_LECT_VAL_INIT_NPT = 1
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        END IF
      END IF
C============================================================
C=====Fin d'utilisation optionnelle des val. init. du fichier
C=====(les 3 cas NPT(x/Tx/xT) et INDIC_LECT_VAL_INIT_NPT = 1)
C============================================================
C=====================================================
C=====================================================
C=====Fin des opérations préliminaires à la résolution
C=====================================================
C=====================================================
C======================================
C======================================
C=====Résolution du système d'équations
C=====à chaque pas de T et/ou de x
C======================================
C======================================
C	write ( * , * ) '###################'
C	write ( * , * ) 'Début résolution NR'	
C	write ( * , * ) '###################'
         CALL
     $   SYS_NR
     $ (
C---------------------------------------
C-----Paramètres dont dépend la fonction
C---------------------------------------
     $   N_TYP , N_TYP_INTR ,
     $   N_R_1 , N_R_2 , N_R_3 , N_R ,
     $   P_R ,
     $   IND_D_R_TYP ,
     $   H_GC_D_R ,
     $   K_T ,
     $   X_AT , N_AT_TOT ,
     $   H_REF_MAILLE ,
C----------------------------------
C-----Paramètres généraux de sys_NR
C----------------------------------
     $   N_NPT , X_NPT , F_NPT , J_NPT , P_J_NPT ,
     $   ALPHA_NRCG_NPT ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT ,
     $   INDIC_TYPE_REDUC_NRCG_NPT , COEF_REDUC_NRCG_NPT ,
     $   N_ITER_MAX_NPT , PRECISION_NPT )
C	write ( * , * ) '#################'
C	write ( * , * ) 'Fin résolution NR'	
C	write ( * , * ) '#################'
C================================================
C================================================
C=====Fin de la résolution du système d'équations
C================================================
C================================================
C--------------------------------------------------------
C-----Vérification de la convergence de l'algorithme NRCG
C---------------------------------------------------------
      ECART_FIN = 0.D0
      DO I = 1 , N_NPT
        IF ( DABS ( X_NPT ( I ) - X_NPT_INIT ( I ) )
     $  .GT. ECART_FIN )
     $  ECART_FIN = DABS ( X_NPT ( I ) - X_NPT_INIT ( I ) )
      END DO
      IF ( ECART_FIN .EQ. 0.D0 ) THEN
        WRITE ( * , * )
     $ '-------------------------------------------------'
        WRITE ( * , * )
     $ 'Algorithme NRCG :'
        WRITE ( * , * )
     $ 'problème de convergence'
        WRITE ( * , * )
     $ '(le nouveau point est identique au point initial)'
        WRITE ( * , * )
     $ '==> modifier le point initial'
        WRITE ( * , * )
     $ '-------------------------------------------------'
        CALL INTERRUPTION
      END IF
C------------------------------------------------------------------
C-----Ecriture de la solution à l'écran et dans le fichier "varNPT"
C----- --> une ligne = le vecteur convergé à chaque pas T ou x
C------------------------------------------------------------------
C       write ( * , 1200 )
C       write ( * , * ) 'Solution :'
C       write ( * , 1200 )
C       write ( * , * )
C    $ 'Quantités de DP (log_10) - Nombre de mailles'
C	write ( * , 1111 )
C    $ ( X_NPT ( I ) , I = 1 , N_NPT )
C       write ( * , 1200 )
C-------------------------------------------------------------------
C-----Ecriture dans le fichier "varNPT" :
C-----* en modes NPT(Tx/xT) -> pour les divers points de composition
C-----à la température visée (la dernière du balayage)
C-----* en mode NPT(x) -> pour les divers points de composition
C-----à la température de consigne
C-------------------------------------------------------------------
C-----Rque : la variable I_LIGNE_FICH_VAR_NPT ne renvoie à rien !
C-----Ce n'est pas gênant, car cet indice n'est pas utilisé
C-----à la lecture d'un fichier de val. init. du type "varNPT".
        INDIC_ECRIT_FICH_VAR_NPT = 0
        IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x' ) THEN
         INDIC_ECRIT_FICH_VAR_NPT = 1
        END IF
        IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
         IF ( I_PAS_BOUCLE_INTERNE .EQ. N_PAS_BOUCLE_INTERNE ) THEN
         INDIC_ECRIT_FICH_VAR_NPT = 1
         END IF
        END IF
        IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
         IF ( I_PAS_BOUCLE_EXTERNE .EQ. N_PAS_BOUCLE_EXTERNE ) THEN
         INDIC_ECRIT_FICH_VAR_NPT = 1
         END IF
        END IF
        IF ( INDIC_ECRIT_FICH_VAR_NPT .EQ. 1) THEN
           WRITE ( 210 , 1112 ) I_LIGNE_FICH_VAR_NPT ,
     $                        ( X_NPT ( I ) , I = 1 , N_NPT )
        END IF
C==============================================================
C==============================================================
C=====Calcul des quantités thermodynamiques dans l'ADPI
C=====(énergie et volume par maille, entropie de configuration,
C=====énergie libre par maille, quantités par atome)
C==============================================================
C==============================================================
C-------------------------------------------------
C-----Passage des variables du système d'équations
C-----aux quantités de DP
C-------------------------------------------------
      DO I_R = 1 , N_R
        DO I_TYP = 0 , N_TYP
C------------------------------------------------
C-----Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C-----ne correspondent à aucun DP ( I_TYP , I_R )
C-----ou à un DP non "utile" ("E_DP >= 0")
C------------------------------------------------
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
            X_D_R ( I_TYP , I_R )
     $    = DEXP ( X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $           * DLOG ( 10.D0 ) )
          END IF
        END DO
      END DO
C-----------------------------------------------------------------
C-----Procédure de calcul des quantités thermodynamiques de l'ADPI
C-----------------------------------------------------------------
         CALL
     $   G_ADPI
     $ ( N_TYP , N_TYP_INTR ,
     $   N_R_1 , N_R_2 , N_R_3 , N_R ,
     $   P_R ,
     $   N_1_MAILLE , N_2_MAILLE , N_3_MAILLE ,
     $   E_REF_MAILLE , V_REF_MAILLE ,
     $   E_REF_TYP ,
     $   X_AT ,
     $   X_D_R ,
     $   E_GC_D_R , V_GC_D_R ,
     $   INDIC_COMPLEXES , N_TYPES_COMPLEXES ,
     $   MULTIPLICITE_COMPLEXE ,
     $   I_S_R_MULTIPLICITE_COMPLEXE ,
     $   E_GC_D_COMPLEXE , V_GC_D_COMPLEXE ,
     $   X_D_COMPLEXE ,
     $   TEMPERATURE , PRESSION ,
     $   N_AT_MAILLE ,
     $   E_MAILLE , V_MAILLE ,
     $   S_CONF_MAILLE , G_MAILLE ,
     $   E_AT , V_AT ,
     $   S_CONF_AT , G_AT ,
     $   G_AT_FORM )
C        write ( * , * ) '     Fin de G_ADPI'
C=============================================================
C=============================================================
C=====Fin du calcul des quantités thermodynamiques dans l'ADPI
C=============================================================
C=============================================================
C================================================================
C================================================================
C=====Potentiels chimiques (= multiplicateurs de Lagrange en NPT)
C================================================================
C================================================================
C=====Calcul préliminaire des termes complémentaires d'entropie
C=====(même calcul que dans F_NR_AVEC_DERIV)
      Z_TYP_R = 1.D0
C-------------------
C-----Sous-réseaux 1
C-------------------
      DO I_R = 1 , N_R_1
         DO I_TYP = 0 , N_TYP
	   IF ( I_TYP .NE. 1 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
	   END IF
         END DO
      END DO
C-----------------------------
C-----Sous-réseaux 2 éventuels
C-----------------------------
      DO I_R = N_R_1 + 1 , N_R_1 + N_R_2
         DO I_TYP = 0 , N_TYP
           IF ( I_TYP .NE. 2 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C-----------------------------
C-----Sous-réseaux 3 éventuels
C-----------------------------
      DO I_R = N_R_1 + N_R_2 + 1 , N_R_1 + N_R_2 + N_R_3
         DO I_TYP = 0 , N_TYP
           IF ( I_TYP .NE. 3 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C-----------------------------------------
C-----Sous-réseaux interstitiels éventuels
C-----------------------------------------
      DO I_R = N_R_1 + N_R_2 + N_R_3 + 1 , N_R
         DO I_TYP = 0 , N_TYP
           IF ( I_TYP .NE. 0 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C=====Potentiels chimiques des éléments intrinsèques
         POT_INTR ( 1 ) 
     $ = - H_GC_D_R ( 0 , N_R_1 )
     $   - K_T
     $   * DLOG ( X_D_R ( 0 , N_R_1 )
     $          / Z_TYP_R ( N_R_1 ) )
      IF ( N_TYP_INTR .GT. 1 ) THEN
         POT_INTR ( 2 ) 
     $ = - H_GC_D_R ( 0 , N_R_1 + N_R_2 )
     $   - K_T
     $   * DLOG ( X_D_R ( 0 , N_R_1 + N_R_2 )
     $          / Z_TYP_R ( N_R_1 + N_R_2 ) )
      END IF
      IF ( N_TYP_INTR .GT. 2 ) THEN
         POT_INTR ( 3 ) 
     $ = - H_GC_D_R ( 0 , N_R_1 + N_R_2 + N_R_3 )
     $   - K_T
     $   * DLOG ( X_D_R ( 0 , N_R_1 + N_R_2 + N_R_3 )
     $          / Z_TYP_R ( N_R_1 + N_R_2 + N_R_3 ) )
      END IF
C=====Potentiels chimiques des éléments d'addition
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
           POT_I ( I_TYP )
     $   = H_GC_D_R ( I_TYP , N_R_1 )
     $   + K_T
     $   * DLOG ( X_D_R ( I_TYP , N_R_1 )
     $          / Z_TYP_R ( N_R_1 ) )
     $   + POT_INTR ( 1 )
      END DO
C===========================================
C===========================================
C=====Fin du calcul des potentiels chimiques
C===========================================
C===========================================
C==============================================================
C==============================================================
C=====Enthalpies de formation des DP (fonctions des pot. chim.)
C==============================================================
C==============================================================
C==================================
C=====Boucle sur les sous-réseaux 1
C==================================
          DO I_R = 1 , N_R_1
C------------
C-----Lacunes
C------------
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_INTR ( 1 )
C--------------------------
C-----Antisites 2 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 1 ) THEN
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) + POT_INTR ( 1 ) - POT_INTR ( 2 )
       END IF
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) + POT_INTR ( 1 ) - POT_INTR ( 3 )
       END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_INTR ( 1 ) - POT_I ( I_TYP )
           END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 1
C=========================================
        END DO
C==================================
C=====Boucle sur les sous-réseaux 2
C=====(lorsqu'ils existent)
C==================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_INTR ( 2 )
C----------------
C-----Antisites 1
C----------------
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_INTR ( 1 ) + POT_INTR ( 2 )
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) - POT_INTR ( 3 ) + POT_INTR ( 2 )
      END IF
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR
C-------------------------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_INTR ( 2 ) - POT_I ( I_TYP )
         END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 2
C=========================================
        END DO
C==================================
C=====Boucle sur les sous-réseaux 3
C=====(lorsqu'ils existent)
C==================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_INTR ( 3 )
C----------------
C-----Antisites 1
C----------------
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_INTR ( 1 ) + POT_INTR ( 3 )
C----------------
C-----Antisites 2
C----------------
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) - POT_INTR ( 2 ) + POT_INTR ( 3 )
C-------------------------------------------------
C-----Traitement éventuel des espèces > N_TYP_INTR 
C-------------------------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
           H_FORM_D_R ( I_TYP , I_R )
     $   = H_GC_D_R ( I_TYP , I_R ) +  POT_INTR ( 3 ) - POT_I ( I_TYP )
         END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 3
C=========================================
        END DO
C==========================================================
C=====Boucle optionnelle sur les sous-réseaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C-----------------------------------------------------
C-----Eléments 1, 2 et 3 (ces deux derniers éventuels)
C-----------------------------------------------------
           H_FORM_D_R ( 1 , I_R )
     $   = H_GC_D_R ( 1 , I_R ) - POT_INTR ( 1 )
           IF ( N_TYP_INTR .GT. 1 ) THEN
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) -  POT_INTR ( 2 )
            END IF
           IF ( N_TYP_INTR .GT. 2 ) THEN
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) - POT_INTR ( 3 )
           END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) - POT_I ( I_TYP )
          END DO
C========================================================
C=====Fin de la boucle sur les sous-réseaux interstitiels
C========================================================
        END DO
C==========================================================
C=====Fin du test d'existence de sous-réseaux interstitiels
C==========================================================
       END IF
C==============================================
C==============================================
C=====Fin du calcul des enthalpies de formation
C==============================================
C==============================================
C============================================
C============================================
C=====Ecriture dans les fichiers de résultats
C============================================
C============================================
C-----------------------------------------
C-----Indicateur d'écriture conditionnelle
C-----suivant le cas NPT(p/T/x/Tx/xT)
C-----------------------------------------
       INDIC_ECRIT_NPT = 0
        IF (
     $      INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $ .OR. INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $ ) THEN
         INDIC_ECRIT_NPT = 1
        END IF
C-----En NPT(Tx/xT), seuls sont écrits les résultats
C-----à la température finale (la plus basse)
        IF (
     $       INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT'
     $ .AND. I_PAS_BOUCLE_INTERNE .EQ. N_PAS_BOUCLE_INTERNE ) THEN
         INDIC_ECRIT_NPT = 1
        END IF
        IF (
     $       INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $ .AND. I_PAS_BOUCLE_EXTERNE .EQ. N_PAS_BOUCLE_EXTERNE ) THEN
         INDIC_ECRIT_NPT = 1
        END IF
C---------------------------------------------------------------------
C-----Test pour écriture conditionnelle, suivant le cas NPT(T/x/Tx/xT)
C---------------------------------------------------------------------
       IF ( INDIC_ECRIT_NPT .EQ. 1 ) THEN
C==================================================================
C=====Ecriture des quantités de défauts ponctuels
C=====par sous-réseau r et par défaut ponctuel :
C=====(i)  si r <= N_R_1,
C=====			 DP = L(r), 2(r), 3(r)  ou i(r) (i > 3) ;
C=====(ii) si N_R_1 < r <= N_R_1 + N_R_2,
C=====			 DP = L(r), 1(r), 3(r), i(r) (i > 3) ;
C=====(iii) si N_R_1 + N_R_2 < r <= N_R_1 + N_R_2 + N_R_3,
C=====			 DP = L(r), 1(r), 2(r), i(r) (i > 3) ;
C=====(iv) si N_R_1 + N_R_2 + N_R_3 < r,
C=====			 DP = i(r) (i > 3)
C==================================================================
C=====(i.e. les DP tels que IND_D_R_TYP # 0)
C==================================================================
C--------------------------------
C-----Boucle sur les sous-réseaux
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------
C-----Boucle sur les types
C-------------------------
        DO I_TYP = 0 , N_TYP
C----------------------
C-----Numéro du fichier
C----------------------
          I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'écriture de l'interstitiel
C-----pour le sous-réseau et l'espèce donnés
C--------------------------------------------
       INDIC_ECRIT_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR ) THEN
          INDIC_ECRIT_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o') )
     $  THEN
          INDIC_ECRIT_INTER = 1
        END IF
       END IF
C--------------------------------------------
C-----Toutes les combinaisons ( I_R , I_TYP )
C-----ne correspondent pas à un DP
C--------------------------------------------
          IF ( ( I_R .LE. N_R_1
     $     .AND. I_TYP .NE. 1 )
     $    .OR. ( I_R .GT. N_R_1
     $     .AND. I_R .LE. N_R_1 + N_R_2
     $     .AND. I_TYP .NE. 2 )
     $    .OR. ( I_R .GT. N_R_1 + N_R_2
     $     .AND. I_R .LE. N_R_1 + N_R_2 + N_R_3
     $     .AND. I_TYP .NE. 3 )
     $    .OR. ( I_R .GT. N_R_1 + N_R_2 + N_R_3
     $     .AND. INDIC_ECRIT_INTER .EQ. 1 ) ) THEN
C----------------------------------
C-----Ecriture de la quantité de DP
C----------------------------------
             WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_R ( I_TYP , I_R ) , H_FORM_D_R ( I_TYP , I_R )
C--------------------------------------
C-----Fin du test d'existence de défaut
C--------------------------------------
        END IF
C----------------------------------------------
C-----Fin des boucles sur types et sous-réseaux
C----------------------------------------------
       END DO
      END DO
C===========================================================
C=====Ecriture de l'énergie, du volume et de l'énergie libre
C===========================================================
C-------------------------------------------------
C-----Cas d'écriture de l'énergie libre par maille
C-------------------------------------------------
            IF 
     $    ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm' )
     $      THEN
             WRITE ( 110 , CAR_COL_VAR_E_L )
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE , 
     $       E_MAILLE , V_MAILLE , S_CONF_MAILLE , G_MAILLE
C------------------------------------------------
C-----Cas d'écriture de l'énergie libre par atome
C------------------------------------------------
            ELSE
C- - - - - - - - - - - - - - - - - - - - - - -
C- - -Dans le cas G/atome, deux possibilités :
C- - -G totale ou G de formation
C- - - - - - - - - - - - - - - - - - - - - - -
             IF ( INDIC_G .EQ. 'T' .OR. INDIC_G .EQ. 't' )
     $       THEN
              WRITE ( 110 , CAR_COL_VAR_E_L )
     $      ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $        TEMPERATURE ,
     $        E_AT , V_AT , S_CONF_AT , G_AT
             ELSE
              WRITE ( 110 , CAR_COL_VAR_E_L )
     $      ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $        TEMPERATURE ,
     $        E_AT , V_AT , S_CONF_AT , G_AT_FORM
C- - - - - - - - - - - - - - - - - -
C- - -Fin du test "totale/formation"
C- - - - - - - - - - - - - - - - - -
             END IF
C-------------------------------
C-----Fin du test "atome/maille"
C-------------------------------
            END IF
C======================================
C=====Ecriture des potentiels chimiques
C======================================
          WRITE ( 120 , CAR_COL_VAR_POT_CHIM )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE
C---------------------------------------------------------------------
C-----Fin du test pour écriture conditionnelle, suivant NPT(T/x/Tx/xT)
C---------------------------------------------------------------------
       END IF
C========================================
C========================================
C=====Fin des écritures dans les fichiers
C========================================
C========================================
C-----------------------------------------------------------------------
C-----Pourcentage du calcul effectué
C-----(avec indicateur d'écriture pour éviter les répétitions à l'écran)
C-----------------------------------------------------------------------
          I_POUR_CENT_NOUV = 100 * I_COMPTEUR_PAS
     $                     / ( N_PAS_BOUCLE_EXTERNE
     $                       * N_PAS_BOUCLE_INTERNE )
          I_ECRIT_POUR_CENT = 0
          IF ( I_POUR_CENT_NOUV .NE. I_POUR_CENT ) THEN
            I_POUR_CENT = I_POUR_CENT_NOUV
            I_ECRIT_POUR_CENT = 1
          END IF
C--------------------------------------
C-----Ecriture du pourcentage à l'écran
C--------------------------------------
          IF ( I_ECRIT_POUR_CENT .EQ. 1 ) THEN
           WRITE ( * , * )
     $      I_POUR_CENT , ' % du calcul NPT effectués'
          END IF
C============================================================
C=====Mise à jour des valeurs initiales entre deux itérations
C===== --> utile seulement en NPT(Tx)
C============================================================
C---------------------------------------------------------------------
C-----Cas NPT(Tx) : tableau pour la conservation des résultats 
C-----(aux diverses compo. = boucle interne) à la température courante
C-----pour l'initialisation à chaque compo. à la température suivante
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----L'indice I_PAS_BOUCLE_INTERNE dans X_NPT_COMPO_SVG
C-----repère les différentes compo. à la température courante T(N),
C----- --> val. init. f(x) utilisées à la température suivante T(N+1).
C---------------------------------------------------------------------
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $     ) THEN
         DO I = 1 , N_NPT
          X_NPT_COMPO_SVG ( I_PAS_BOUCLE_INTERNE , I ) = X_NPT ( I ) 
         END DO
        END IF
C#######################################################
C#####Fin de la double boucle sur les pas de température
C#####et de fraction atomique pour l'élément spécifié
C#######################################################
       END DO
      END DO
C##########################
C##########################
C#####Fin du traitement NPT
C##########################
C##########################
      END IF
      WRITE ( * , * ) 
     $ '----------------------------------------------------------'
      WRITE ( * , * )
     $ "REMARQUE concernant l'utilisation des fichiers de DP :"
      WRITE ( * , * ) 
     $ 'pour éviter des déconvenues lors des tracés de graphes'
      WRITE ( * , * ) 
     $ '(points joints de manière erronée), il peut être utile'
      WRITE ( * , * )
     $ 'de classer le contenu de chacun de ces fichiers'
      WRITE ( * , * ) 
     $ 'par valeurs croissantes de la composition mise en abscisse'
      WRITE ( * , * )
     $ '(la présente version du programme adpi'
      WRITE ( * , * )
     $ "n'effectue pas ce classement)"
      WRITE ( * , * ) 
     $ '----------------------------------------------------------'
      WRITE ( * , * )
      WRITE ( * , * ) '          ============================'
      WRITE ( * , * ) '          Calcul ADPI (DP non chargés)'
      WRITE ( * , * ) '          ----------------------------'
      WRITE ( * , * ) '                 FIN DU PROGRAMME'
      WRITE ( * , * ) '          ============================'
C-----Instruction permettant de terminer le programme
C-----en évitant le cas "DP CHARGES" ci-après
      GOTO 8888
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&& Fin du cas "DP NON CHARGES"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&& 
C&&&&& Cas "DP CHARGES"
C&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&& 
C-----Instruction de branchement du cas "DP CHARGES"
 3333 CONTINUE
C---------------------------------------------------------------------
C-----Calcul des nombres d'atomes de chaque espèce par maille
C-----(utiles au calcul de la somme pondérée des potentiels chimiques)
C-----(déjà présent dans "adpi sans charge")
C---------------------------------------------------------------------
      N_1_MAILLE = 0
      N_2_MAILLE = 0
      N_3_MAILLE = 0
      DO I_R = 1 , N_R_1
       N_1_MAILLE = N_1_MAILLE + P_R ( I_R )
      END DO
      DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
       N_2_MAILLE = N_2_MAILLE + P_R ( I_R )
      END DO
      DO I_R = 1 + N_R_2 + N_R_1 , N_R_3 + N_R_2 + N_R_1
       N_3_MAILLE = N_3_MAILLE + P_R ( I_R )
      END DO
C------------------------------------------------------------------
C-----Calcul de l'énergie et du volume de référence par maille (eV)
C------------------------------------------------------------------
      E_REF_MAILLE = E_REF / DFLOAT ( N_MAILLE_REF )
      V_REF_MAILLE = V_REF / DFLOAT ( N_MAILLE_REF )
C---------------------------------------------------------------------
C-----Calcul des enthalpies totale et par maille de la cellule sans DP
C---------------------------------------------------------------------
      H_REF = E_REF - V_REF * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
      H_REF_MAILLE = H_REF / DFLOAT ( N_MAILLE_REF )
C	write ( * , * ) V_REF_MAILLE
C------------------------------------------------------------------
C-----Calcul approché de la somme pondérée des potentiels chimiques
C-----pour une maille à partir de la relation de Gibbs-Duhem
C------------------------------------------------------------------
      S_POT = V_REF_MAILLE * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
     $      + E_REF_MAILLE
C===========================================================
C=====Potentiel chimique de 1 (élément de référence)
C=====(expression valable pour 2 ou 3 éléments intrinsèques)
C===========================================================
      POT_1 = S_POT - DFLOAT ( N_2_MAILLE ) * POT_2
     $              - DFLOAT ( N_3_MAILLE ) * POT_3
      POT_1 = POT_1 / DFLOAT ( N_1_MAILLE )
      write ( * , * ) 'mu(1)=',POT_1
      write ( * , * ) 'mu(2)=',POT_2
      write ( * , * ) 'mu(3)=',POT_3
C================================================
C=====Calcul optionnel de grandeurs préliminaires
C=====relatives aux complexes chargés éventuels
C================================================
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C--------------------------------------------------------------------
C-----Calcul du nombre de sites du sous-réseau r dans chaque complexe
C-----et du nombre d'atomes de type i dans ce complexe
C--------------------------------------------------------------------
        ALLOCATE ( U_COMPLEXE_S_R_Q ( N_TYPES_COMPLEXES_Q , N_R ) )
        ALLOCATE ( V_COMPLEXE_TYPE_Q ( N_TYPES_COMPLEXES_Q , N_TYP ) )
        U_COMPLEXE_S_R_Q = O
        V_COMPLEXE_TYPE_Q = 0
C-------------------------------------
C-----Boucle sur les complexes chargés
C-------------------------------------
        DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
C-------------------------------------
C-----Boucle sur les sites du complexe
C-------------------------------------
          DO I_SITE = 1 , NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q )
           DO I_R = 1 , N_R
            IF ( I_S_R_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE ) .EQ. I_R )
     $      THEN
            U_COMPLEXE_S_R_Q ( I_COMPLEXE_Q , I_R )
     $    = U_COMPLEXE_S_R_Q ( I_COMPLEXE_Q , I_R ) + 1
            END IF
           END DO
           DO I_TYP = 1 , N_TYP
            IF ( I_TYPE_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE )
     $      .EQ. I_TYP )
     $      THEN
            V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , I_TYP )
     $    = V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , I_TYP ) + 1
            END IF
           END DO
C-----------------------------------------------
C-----Fin de la boucle sur les sites du complexe
C-----------------------------------------------
          END DO
C-----------------------------------------------
C-----Fin de la boucle sur les complexes chargés
C-----------------------------------------------
        END DO
C--------------------------------------------------------
C-----Paramètres indicateurs de types et de sous-réseaux
C-----utiles au calcul des termes de potentiels chimiques
C-----relatifs aux complexes éventuels
C--------------------------------------------------------
      ALLOCATE ( ALPHA_TYPE_COMPLEXE ( 0 : N_TYP ) )
      ALLOCATE ( BETA_S_R_COMPLEXE ( 1 : N_R ) ) 
      ALPHA_TYPE_COMPLEXE = 1
      ALPHA_TYPE_COMPLEXE ( 0 ) = 0
      BETA_S_R_COMPLEXE = 1
      DO I_R = N_R_1 + N_R_2 + N_R_3 + 1 , N_R
         BETA_S_R_COMPLEXE ( I_R ) = 0
      END DO
C--------------------------------------------------------------------
C-----Type chimique normal de chaque sous-réseau
C-----(0 pour sous-réseaux interstitiels)
C-----utile au calcul des termes de potentiels chimiques de complexes
C--------------------------------------------------------------------
      ALLOCATE ( I_TYPE_NORMAL_S_R ( N_R ) )
      I_TYPE_NORMAL_S_R = 0
      DO I_R = 1 , N_R_1
        I_TYPE_NORMAL_S_R ( I_R ) = 1
      END DO
      DO I_R = N_R_1 + 1 , N_R_1 + N_R_2
        I_TYPE_NORMAL_S_R ( I_R ) = 2
      END DO
      DO I_R = N_R_1 + N_R_2 + 1 , N_R_1 + N_R_2 + N_R_3
        I_TYPE_NORMAL_S_R ( I_R ) = 3
      END DO
C---------------------------------------------------------
C-----Fin du test de prise en compte des complexes chargés
C---------------------------------------------------------
       END IF
C=======================================================
C=====Fin de calcul optionnel de grandeurs préliminaires
C=====relatives aux complexes chargés éventuels
C=======================================================
C###############################
C#####Calcul des Egc, Vgc et Hgc
C###############################
C--------------------------------------
C--------------------------------------
C-----Partie 1/3 : antisites et lacunes
C--------------------------------------
C--------------------------------------
C------------------------------------------------
C-----antisites des diverses espèces intrinsèques
C-----sur les divers sous-réseaux
C------------------------------------------------
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           DO I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP )
            DO I_Q = 1 , NQ_D_R_Q ( J_TYP , I_R ) 
            E_GC_D_R_Q ( I_Q , J_TYP , I_R )
     $    = E_B_D_R_Q ( I_Q , J_TYP , I_R )
     $    - E_REF
            V_GC_D_R_Q ( I_Q , J_TYP , I_R )
     $    = V_B_D_R_Q ( I_Q , J_TYP , I_R )
     $    - V_REF
            H_GC_D_R_Q ( I_Q , J_TYP , I_R )
     $    = E_GC_D_R_Q ( I_Q , J_TYP , I_R )
     $    + PRESSION * V_GC_D_R_Q ( I_Q , J_TYP , I_R )
           END DO
          END DO
         END IF
        END DO
       END DO
C---------------------------------------------
C-----Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C---------------------------------------------
      DO I_R = 1 , RHO_N ( N_TYP_INTR ) 
         DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R ) 
            E_GC_D_R_Q ( I_Q , 0 , I_R )
     $    = E_B_D_R_Q ( I_Q , 0 , I_R )
     $    - E_REF
            V_GC_D_R_Q ( I_Q , 0 , I_R )
     $    = V_B_D_R_Q ( I_Q , 0 , I_R )
     $    - V_REF
            H_GC_D_R_Q ( I_Q , 0 , I_R )
     $    = E_GC_D_R_Q ( I_Q , 0 , I_R )
     $    + PRESSION * V_GC_D_R_Q ( I_Q , 0 , I_R )
         END DO
      END DO
C----------------------------------------------------------
C----------------------------------------------------------
C-----Partie 2/3 (optionnelle) : interstitiels intrinsèques
C----------------------------------------------------------
C----------------------------------------------------------
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR 
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R 
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP_INTR , I_R ) 
            E_GC_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    = E_B_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    - E_REF
            V_GC_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    = V_B_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    - V_REF
            H_GC_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    = E_GC_D_R_Q ( I_Q , I_TYP_INTR , I_R )
     $    + PRESSION * V_GC_D_R_Q ( I_Q , I_TYP_INTR , I_R )
         END DO
        END DO
       END DO
      END IF
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----Partie 3/3 :
C-----lecture optionnelle des énergies d'éléments d'addition
C-----(substitutionnels et interstitiels)
C-----------------------------------------------------------
C-----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_R = 1 , RHO_N ( N_TYP_INTR )
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R ) 
            E_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = E_B_D_R_Q ( I_Q , I_TYP , I_R )
     $    - E_REF
            V_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = V_B_D_R_Q ( I_Q , I_TYP , I_R )
     $    - V_REF
            H_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = E_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    + PRESSION * V_GC_D_R_Q ( I_Q , I_TYP , I_R )
           END DO
          END DO
C-----------------------------------------------
C-----Cas de présence d'interstitiels d'addition
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R ) 
            E_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = E_B_D_R_Q ( I_Q , I_TYP , I_R )
     $    - E_REF
            V_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = V_B_D_R_Q ( I_Q , I_TYP , I_R )
     $    - V_REF
            H_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    = E_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $    + PRESSION * V_GC_D_R_Q ( I_Q , I_TYP , I_R )
           END DO
          END DO
        END IF
C---------------------------------------------
C---------------------------------------------
C-----Fin de boucle sur les types (partie 3/3)
C---------------------------------------------
C---------------------------------------------
      END DO
C-------------------------------------------
C-------------------------------------------
C-----Partie optionnelle : complexes chargés
C-------------------------------------------
C-------------------------------------------
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
       DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
        DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q ) 
            E_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    = E_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    - E_REF
            V_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    = V_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    - V_REF
            H_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    = E_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $    + PRESSION * V_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
         END DO
        END DO
      END IF
C	stop
C#################################################
C#################################################
C#####              DP CHARGES
C#####2 options : balayage en mu(électrons)
C#####ou recherche (par la méthode NR) de mu(élec) 
C#####tel que n+A=p+D (neutralité électrique)
C#################################################
C#################################################
C#################################################################
C#####Rappel : en présence de DP chargés :
C##### * mode muVT seul possible
C##### * pas de balayage en potentiels chimiques pour les éléments
C#####  --> potentiels chimiques fixés directement dans DATA.adpi
C#################################################################
C==============================
C=====Lecture du fichier de DdE
C===============================
       OPEN ( UNIT = 20 ,
     $ FILE = REP_DDE ( 1 : LONG_REP_DDE )
     $     // FICH_DDE ( 1 : LONG_FICH_DDE ) )
      DO I = 1 , 5
       READ ( 20 , * ) 
      END DO
      READ ( 20 , * ) E_MAX_DDE , E_MIN_DDE , N_PAS_DDE
      ALLOCATE ( TAB_DDE ( N_PAS_DDE , 2 ) )
      DO I_PAS = 1 , N_PAS_DDE
      READ ( 20 , * ) ( TAB_DDE ( I_PAS , I ) , I = 1 , 2 )
      END DO
C----------------------
C-----Ecriture de titre
C----------------------
      WRITE ( * , 2000 ) V_REF
      IF ( I_CALC_CHARGE .EQ. 1 ) THEN
       WRITE ( * , 2001 )
      END IF
      IF ( I_CALC_CHARGE .EQ. 2 ) THEN
       WRITE ( * , 2002 )
      END IF
C#######################################################################
C#######################################################################
C#####Boucle sur mu_e pour balayage ou algorithme NR (neutralité élec).
C#####Rappel : pas de balayage en potentiels chimiques pour les éléments
C#######################################################################
C#######################################################################
      F_COUR = 1.D100
      I_PAS_MU_ELEC = - 1
      I_ARRET = 0
      DO WHILE ( I_ARRET .EQ. 0 ) 
        I_PAS_MU_ELEC = I_PAS_MU_ELEC + 1
C=========================================
C=========================================
C=====Calcul de n et p, et de leur dérivée
C=========================================
C=========================================
         CALL
     $   CALC_N_P
     $ ( N_PAS_DDE , TAB_DDE ,
     $   E_MIN_DDE , E_MAX_DDE ,
     $   E_MAX_BV , E_MIN_BC ,
     $   POT_CHIM_ELEC , TEMPERATURE ,
     $   C_VOL_N , C_VOL_P ,
     $   D_C_VOL_N , D_C_VOL_P )
C===========================
C===========================
C=====Calcul des H_formation
C===========================
C===========================
C==================================
C=====Boucle sur les sous-réseaux 1
C==================================
      DO I_R = 1 , N_R_1
C------------
C-----Lacunes
C------------
       DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R )
           H_FORM_D_R_Q ( I_Q , 0 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 0 , I_R ) + POT_1
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C----------------
C-----Antisites 2
C----------------
       DO I_Q = 1 , NQ_D_R_Q ( 2 , I_R )
           H_FORM_D_R_Q ( I_Q , 2 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 2 , I_R ) + POT_1 - POT_2
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 2 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
       DO I_Q = 1 , NQ_D_R_Q ( 3 , I_R )
           H_FORM_D_R_Q ( I_Q , 3 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 3 , I_R ) + POT_1 - POT_3
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 3 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
     $   = H_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $   + POT_1 - POT_I_INIT ( I_TYP )
     $   + DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
         END DO
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 1
C=========================================
      END DO      
C============================================
C=====Boucle sur les sous-réseaux 2 éventuels
C============================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
       DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R )
           H_FORM_D_R_Q ( I_Q , 0 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 0 , I_R ) + POT_2
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C----------------
C-----Antisites 1
C----------------
       DO I_Q = 1 , NQ_D_R_Q ( 1 , I_R )
           H_FORM_D_R_Q ( I_Q , 1 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 1 , I_R ) - POT_1 + POT_2
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 1 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
       DO I_Q = 1 , NQ_D_R_Q ( 3 , I_R )
           H_FORM_D_R_Q ( I_Q , 3 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 3 , I_R ) + POT_2 - POT_3
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 3 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
     $   = H_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $   + POT_2 - POT_I_INIT ( I_TYP )
     $   + DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
         END DO
        END DO
C===================================================
C=====Fin de boucle sur les sous-réseaux 2 éventuels
C===================================================
      END DO      
C============================================
C=====Boucle sur les sous-réseaux 3 éventuels
C============================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
       DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R )
           H_FORM_D_R_Q ( I_Q , 0 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 0 , I_R ) + POT_3
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C----------------
C-----Antisites 1
C----------------
       DO I_Q = 1 , NQ_D_R_Q ( 1 , I_R )
           H_FORM_D_R_Q ( I_Q , 1 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 1 , I_R ) - POT_1 + POT_3
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 1 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C----------------
C-----Antisites 2
C----------------
       DO I_Q = 1 , NQ_D_R_Q ( 2 , I_R )
           H_FORM_D_R_Q ( I_Q , 2 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 2 , I_R ) + POT_3 - POT_2
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 2 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
       END DO
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
     $   = H_GC_D_R_Q ( I_Q , I_TYP , I_R )
     $   + POT_3 - POT_I_INIT ( I_TYP )
     $   + DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
         END DO
        END DO
C===================================================
C=====Fin de boucle sur les sous-réseaux 3 éventuels
C===================================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-réseaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C-----------------------------------------------------
C-----Eléments 1, 2 et 3 (ces deux derniers éventuels)
C-----------------------------------------------------
         DO I_Q = 1 , NQ_D_R_Q ( 1 , I_R )
           H_FORM_D_R_Q ( I_Q , 1 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 1 , I_R ) - POT_1
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 1 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
         END DO
         IF ( N_TYP_INTR .GT. 1 ) THEN
          DO I_Q = 1 , NQ_D_R_Q ( 2 , I_R )
           H_FORM_D_R_Q ( I_Q , 2 , I_R )
     $   = H_GC_D_R_Q ( I_Q , 2 , I_R ) - POT_2
     $   + DFLOAT ( Q_D_R_Q ( I_Q , 2 , I_R ) )
     $   * ( E_MAX_BV + POT_CHIM_ELEC )
          END DO
          IF ( N_TYP_INTR .GT. 2 ) THEN
           DO I_Q = 1 , NQ_D_R_Q ( 3 , I_R )
            H_FORM_D_R_Q ( I_Q , 3 , I_R )
     $    = H_GC_D_R_Q ( I_Q , 3 , I_R ) - POT_3
     $    + DFLOAT ( Q_D_R_Q ( I_Q , 3 , I_R ) )
     $    * ( E_MAX_BV + POT_CHIM_ELEC )
           END DO
          END IF
         END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
             H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
     $     = H_GC_D_R_Q ( I_Q , I_TYP , I_R ) - POT_I_INIT ( I_TYP )
     $     + DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $     * ( E_MAX_BV + POT_CHIM_ELEC )
          END DO
         END DO
C=================================================================
C=====Fin de boucle optionnelle sur les sous-réseaux interstitiels
C=================================================================
        END DO
       END IF
C        do i_typ = 0 , n_typ
C         do i_r = 1 , n_r
C        write ( * , * ) i_typ , i_r , nq_d_r_q ( i_typ , i_r )
C         end do
C        end do
C        do i_typ = 0 , n_typ
C         do i_r = 1 , n_r
C          do i_q = 1 , nq_d_r_q ( i_typ , i_r )
C        write ( * , * ) i_typ , i_r , i_q , 
C     $       h_form_d_r_q ( i_q , i_typ , i_r )
C          end do
C         end do
C        end do
C        write ( * , * ) pot_1 , pot_2 , pot_3
C		stop
C----------------------------------------------------
C----------------------------------------------------
C-----Partie optionnelle : Hform de complexes chargés
C----------------------------------------------------
C----------------------------------------------------
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C=====================================
C=====Boucle sur les complexes chargés
C=====================================
        DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
          DELTA_MU_COMPLEXE = 0.D0
C-------------------------------------
C-------------------------------------
C-----Boucle sur les sites du complexe
C-------------------------------------
C-------------------------------------
         DO I_SITE = 1 , NOMBRE_SITES_COMPLEXE_Q ( I_COMPLEXE_Q )
           I_S_R_COUR = I_S_R_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE )
           I_TYPE_COUR = I_TYPE_COMPLEXE_Q ( I_COMPLEXE_Q , I_SITE )
           I_TYPE_NORMAL_COUR = I_TYPE_NORMAL_S_R ( I_S_R_COUR )
C------------------------------------------------------------------
C-----Préliminaire :
C-----calcul du potentiel chimique de l'élément sur le site courant
C-----et du potentiel chimique "normal" sur le site courant
C------------------------------------------------------------------
C-----Cas 1 : un type intrinsèque
            IF ( N_TYP_INTR .EQ. 1 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
                POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
                POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
                POT_CHIM_NORMAL_SITE = 0.D0
              END IF
C-----Cas 2 : deux types intrinsèques
           ELSE IF ( N_TYP_INTR .EQ. 2 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .EQ. 2 ) THEN
                POT_CHIM_SITE = POT_2
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
               POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 2 ) THEN
               POT_CHIM_NORMAL_SITE = POT_2
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
               POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
               POT_CHIM_NORMAL_SITE = 0.D0
              END IF
C-----Cas 3 : trois types intrinsèques
           ELSE IF ( N_TYP_INTR .EQ. 3 ) THEN
              IF ( I_TYPE_COUR .EQ. 1 ) THEN
                POT_CHIM_SITE = POT_1
              ELSE IF ( I_TYPE_COUR .EQ. 2 ) THEN
                POT_CHIM_SITE = POT_2
              ELSE IF ( I_TYPE_COUR .EQ. 3 ) THEN
                POT_CHIM_SITE = POT_3
              ELSE IF ( I_TYPE_COUR .NE. 0 ) THEN
                POT_CHIM_SITE = POT_I ( I_TYPE_COUR )
              ELSE
                POT_CHIM_SITE = 0.D0
              END IF
              IF ( I_TYPE_NORMAL_COUR .EQ. 1 ) THEN
                POT_CHIM_NORMAL_SITE = POT_1
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 2 ) THEN
                POT_CHIM_NORMAL_SITE = POT_2
              ELSE IF ( I_TYPE_NORMAL_COUR .EQ. 3 ) THEN
                POT_CHIM_NORMAL_SITE = POT_3
              ELSE IF ( I_TYPE_NORMAL_COUR .NE. 0 ) THEN
                POT_CHIM_NORMAL_SITE = POT_I ( I_TYPE_NORMAL_COUR )
              ELSE
               POT_CHIM_NORMAL_SITE = 0.D0
              END IF
           END IF
C		stop
C-------------------------------------------------------------------
C-----Calcul du terme de potentiel chimique pour le complexe courant
C-----(le même pour tous les états de charge d'un complexe donné)
C-------------------------------------------------------------------
           DELTA_MU_COMPLEXE
     $   = DELTA_MU_COMPLEXE
     $   + DFLOAT ( BETA_S_R_COMPLEXE ( I_S_R_COUR ) )
     $   * POT_CHIM_NORMAL_SITE  
     $   - DFLOAT ( ALPHA_TYPE_COMPLEXE ( I_TYPE_COUR ) )
     $   * POT_CHIM_SITE
C-----------------------------------------------
C-----------------------------------------------
C-----Fin de la boucle sur les sites du complexe
C-----------------------------------------------
C-----------------------------------------------
          END DO
C	stop
C-------------------------------------------------------
C-----Boucle sur les états de charge du complexe courant
C-----et calcul des enthalpies de formation associées
C-------------------------------------------------------
       DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q ) 
         H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $ = H_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $ + DELTA_MU_COMPLEXE
     $ + DFLOAT ( Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) )
     $ * ( E_MAX_BV + POT_CHIM_ELEC )
      END DO
C===============================================
C=====Fin de la boucle sur les complexes chargés
C===============================================
      END DO
C------------------------------------------------
C------------------------------------------------
C-----Fin de Hform(complexes chargés) (optionnel)
C------------------------------------------------
C------------------------------------------------
       END IF
C	stop
C=========================================================
C=========================================================
C=====Concentrations de charges par donneurs et accepteurs
C=========================================================
C=========================================================
      C_VOL_ACC = 0.D0
      C_VOL_DON = 0.D0
      D_C_VOL_ACC = 0.D0
      D_C_VOL_DON = 0.D0
C======================================
C=====Partie 1/3 : antisites et lacunes
C======================================
C------------------------------------------------
C------------------------------------------------
C-----Antisites des diverses espèces intrinsèques
C-----sur les divers sous-réseaux
C------------------------------------------------
C------------------------------------------------
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           DO I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP )
            DO I_Q = 1 , NQ_D_R_Q ( J_TYP , I_R ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_D_R_Q ( I_Q , J_TYP , I_R ) .GT. 0 ) THEN
               XX
     $       = DFLOAT ( Q_D_R_Q ( I_Q , J_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , J_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , J_TYP , I_R ) ) / K_T
C              write (*,*)
C     $ 'type=',j_typ,'sr=',i_r ,'q=',Q_D_R_Q ( I_Q , J_TYP , I_R ),
C     $   '  Contribution aux donneurs = ' ,xx
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_D_R_Q ( I_Q , J_TYP , I_R ) .LT. 0 ) THEN
               XX
     $       = - DFLOAT ( Q_D_R_Q ( I_Q , J_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , J_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , J_TYP , I_R ) ) / K_T
C              write (*,*)
C     $   'type=',j_typ,'sr=',i_r ,'q=',Q_D_R_Q ( I_Q , J_TYP , I_R ),
C     $   '  Contribution aux accepteurs = ' ,xx
            END IF
           END DO
          END DO
         END IF
        END DO
       END DO
C---------------------------------------------
C---------------------------------------------
C-----Lacunes L(r) (pour r <= rho(N_TYP_INTR))
C---------------------------------------------
C---------------------------------------------
      DO I_R = 1 , RHO_N ( N_TYP_INTR ) 
         DO I_Q = 1 , NQ_D_R_Q ( 0 , I_R ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_D_R_Q ( I_Q , 0 , I_R ) .GT. 0 ) THEN
               XX
     $       = DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , 0 , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) ) / K_T
C              write (*,*)
C     $ 'type=0 ','sr=',i_r ,'q=',Q_D_R_Q ( I_Q , 0 , I_R ),
C     $   '  Contribution aux donneurs = ' ,xx
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_D_R_Q ( I_Q , 0 , I_R ) .LT. 0 ) THEN
               XX
     $       = - DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , 0 , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , 0 , I_R ) ) / K_T
C              write (*,*)
C     $ 'type=0 ','sr=',i_r ,'q=',Q_D_R_Q ( I_Q , 0 , I_R ),
C     $   '  Contribution aux accepteurs = ' ,xx
             END IF
         END DO
      END DO
C========================================================
C=====Partie 2/3 optionnelle : interstitiels intrinsèques
C========================================================
      IF ( ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' )
     $ .AND. ( INDIC_INTER_INTR .EQ. 'O'
     $    .OR. INDIC_INTER_INTR .EQ. 'o' ) )
     $ THEN
       DO I_TYP_INTR = 1 , N_TYP_INTR 
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R 
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP_INTR , I_R ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) .GT. 0 ) THEN
               XX
     $       = DFLOAT ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP_INTR , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) ) / K_T
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) .LT. 0 ) THEN
               XX
     $       = - DFLOAT ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP_INTR , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP_INTR , I_R ) ) / K_T
             END IF
         END DO
        END DO
       END DO
      END IF
C===================================================
C=====Partie 3/3 (optionnelle) : éléments d'addition
C=====(substitutionnels et interstitiels)
C===================================================
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_R = 1 , RHO_N ( N_TYP_INTR )
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP , I_R ) .GT. 0 ) THEN
               XX
     $       = DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) ) / K_T
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP , I_R ) .LT. 0 ) THEN
               XX
     $       = - DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) ) / K_T
             END IF
           END DO
          END DO
C-----------------------------------------------
C-----------------------------------------------
C-----Cas de présence d'interstitiels d'addition
C-----------------------------------------------
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
         DO I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R
          DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP , I_R ) .GT. 0 ) THEN
               XX
     $       = DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) ) / K_T
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_D_R_Q ( I_Q , I_TYP , I_R ) .LT. 0 ) THEN
               XX
     $       = - DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) )
     $       * DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * DFLOAT ( Q_D_R_Q ( I_Q , I_TYP , I_R ) ) / K_T
             END IF
           END DO
          END DO
C-------------------------------------------------------
C-----Fin du test de présence d'interstitiels d'addition
C-------------------------------------------------------
        END IF
C---------------------------------------------
C---------------------------------------------
C-----Fin de boucle sur les types (partie 3/3)
C---------------------------------------------
C---------------------------------------------
      END DO
C===========================================
C=====Partie optionnelle : complexes chargés
C===========================================
      DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
          I_R = I_S_R_MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q )
          DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q ) 
C---------------------------------
C-----Contribution des DP donneurs
C---------------------------------
             IF ( Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) .GT. 0 ) THEN
C			write ( * , * ) 'Qcmplx>0' , I_Q , I_COMPLEXE_Q
               XX
     $       = Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $       * DEXP ( - H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) / K_T )
     $       * DFLOAT ( MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
C			write ( * , * ) 'XX = ' , XX
               C_VOL_DON = C_VOL_DON + XX
               D_C_VOL_DON = D_C_VOL_DON
     $       - XX * Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) / K_T
             END IF
C-----------------------------------
C-----Contribution des DP accepteurs
C-----------------------------------
             IF ( Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) .LT. 0 ) THEN
C			write ( * , * ) 'Qcmplx<0' , I_Q , I_COMPLEXE_Q
               XX
     $       = - Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $       * DEXP ( - H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) / K_T )
     $       * DFLOAT ( MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) )
     $       * DFLOAT ( N_MAILLE_REF * P_R ( I_R ) )
C			write ( * , * ) 'XX = ' , XX
               C_VOL_ACC = C_VOL_ACC + XX
               D_C_VOL_ACC = D_C_VOL_ACC
     $       - XX * Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) / K_T
             END IF
           END DO
          END DO
C=======================================================
C=====Balayage en potentiel chimique (mode 1)
C===== -> mise à jour du potentiel chimique électronique
C=======================================================
       IF ( I_CALC_CHARGE .EQ. 1 ) THEN
        write ( * , 3001 ) I_PAS_MU_ELEC , POT_CHIM_ELEC ,
     $        C_VOL_N + C_VOL_ACC , C_VOL_P + C_VOL_DON 
        POT_CHIM_ELEC = POT_CHIM_ELEC + D_POT_CHIM_ELEC
       END IF 
C=======================================================
C=====Fonction-objectif n+A-p-D et sa dérivée (mode 2)
C===== -> mise à jour du potentiel chimique électronique
C=======================================================
       IF ( I_CALC_CHARGE .EQ. 2 ) THEN
        F_COUR = C_VOL_N - C_VOL_P + C_VOL_ACC - C_VOL_DON
C        write ( * , * ) 'F_COUR  = ' , F_COUR 
        D_F_COUR = D_C_VOL_N - D_C_VOL_P + D_C_VOL_ACC - D_C_VOL_DON
        write ( * , 3002 ) I_PAS_MU_ELEC , POT_CHIM_ELEC ,
     $                  C_VOL_N , C_VOL_ACC , C_VOL_P , C_VOL_DON ,
     $                  F_COUR , D_F_COUR 
        POT_CHIM_ELEC = POT_CHIM_ELEC - F_COUR / D_F_COUR
       END IF
C--------------------------------
C-----Test pour l'arrêt en mode 1
C--------------------------------
        IF ( I_CALC_CHARGE .EQ. 1 ) THEN
         IF ( POT_CHIM_ELEC .GE. POT_CHIM_ELEC_MAX )
     $    I_ARRET = 1 
        END IF
C--------------------------------
C-----Test pour l'arrêt en mode 2
C--------------------------------
        IF ( I_CALC_CHARGE .EQ. 2 ) THEN
         IF ( DABS ( F_COUR ) .LT. EPS
     $     .OR. I_PAS_MU_ELEC .GT. N_MAX_PAS_NR_ELEC )
     $    I_ARRET = 1 
        END IF
C#######################################################################
C#######################################################################
C#####Fin de la boucle de balayage ou d'algorithme NR
C#####pour le potentiel chimique électronique
C#####Rappel : pas de balayage en potentiels chimiques pour les éléments
C#######################################################################
C#######################################################################
      END DO
      IF ( I_CALC_CHARGE .EQ. 2 ) THEN
        IF ( I_PAS_MU_ELEC .GE. N_MAX_PAS_NR_ELEC ) THEN
          write ( * , * ) 'Nmax(pas NR) atteint'
         END IF
      END IF
C	stop
C===================================
C===================================
c=====En mode 2 (neutralité élec.),
C=====calcul des fractions atomiques
C===================================
C===================================
      IF ( I_CALC_CHARGE .EQ. 2 ) THEN
C-------------------------------------------------
C-----Termes au numérateur des fractions atomiques
C-------------------------------------------------
          SOMME_1_TYP_2 = 0.D0
          SOMME_2_TYP_2 = 0.D0
          SOMME_1_TYP_3 = 0.D0
          SOMME_2_TYP_3 = 0.D0
C-----SOMME_I déjà déclaré dans "adpi"
          ALLOCATE ( SOMME_I ( N_TYP_INTR + 1 : N_TYP ) )
          SOMME_I = 0.D0
C---------------------------------------------------
C-----Termes au dénominateur des fractions atomiques
C---------------------------------------------------
          SOMME_0 = 0.D0
          SOMME_MAILLE = 0.D0
C         write ( * , * ) 'P_R = ' , P_R
C==================================
C=====Boucle sur les sous-réseaux 1
C==================================
      DO I_R = 1 , N_R_1
         SOMME_0 = DFLOAT ( P_R ( I_R ) )
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_0
     $   = SOMME_0 - DFLOAT ( P_R ( I_R ) ) * X_DP
        END DO
C----------------
C-----Antisites 2
C----------------
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
        END DO
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        I_TYP = 3
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
        END DO
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
             SOMME_I ( I_TYP )
     $     = SOMME_I ( I_TYP )
     $     + DFLOAT ( P_R ( I_R ) ) * X_DP
         END DO
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 1
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-réseaux 2
C==================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
         SOMME_1_TYP_2_R = 1.D0
         SOMME_0 = SOMME_0 + DFLOAT ( P_R ( I_R ) )
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_0
     $   = SOMME_0 - DFLOAT ( P_R ( I_R ) ) * X_DP
          SOMME_1_TYP_2_R = SOMME_1_TYP_2_R - X_DP
        END DO
C----------------
C-----Antisites 1
C----------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_1_TYP_2_R = SOMME_1_TYP_2_R - X_DP
        END DO
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        I_TYP = 3
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
           SOMME_1_TYP_2_R = SOMME_1_TYP_2_R - X_DP
        END DO
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
            SOMME_1_TYP_2_R
     $    = SOMME_1_TYP_2_R - X_DP
            SOMME_I ( I_TYP )
     $    = SOMME_I ( I_TYP )
     $    + DFLOAT ( P_R ( I_R ) ) * X_DP
         END DO
        END DO
C-------------------------------------
C-----Terme de somme 1 pour l'espèce 2
C-------------------------------------
             SOMME_1_TYP_2 = SOMME_1_TYP_2
     $                     + SOMME_1_TYP_2_R * DFLOAT ( P_R ( I_R ) )
C=========================================
C=====Fin de boucle sur les sous-réseaux 2
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-réseaux 3
C==================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
         SOMME_1_TYP_3_R = 1.D0
         SOMME_0 = SOMME_0 + DFLOAT ( P_R ( I_R ) )
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_0
     $   = SOMME_0 - DFLOAT ( P_R ( I_R ) ) * X_DP
           SOMME_1_TYP_3_R = SOMME_1_TYP_3_R - X_DP
        END DO
C----------------
C-----Antisites 1
C----------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_1_TYP_3_R = SOMME_1_TYP_3_R - X_DP
        END DO
C----------------
C-----Antisites 2
C----------------
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
           SOMME_1_TYP_3_R
     $   = SOMME_1_TYP_3_R - X_DP
        END DO
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
            SOMME_1_TYP_3_R
     $    = SOMME_1_TYP_3_R - X_DP
            SOMME_I ( I_TYP )
     $    = SOMME_I ( I_TYP )
     $    + DFLOAT ( P_R ( I_R ) ) * X_DP
         END DO
        END DO
C-------------------------------------
C-----Terme de somme 1 pour l'espèce 3
C-------------------------------------
         SOMME_1_TYP_3 = SOMME_1_TYP_3
     $                 + SOMME_1_TYP_3_R * DFLOAT ( P_R ( I_R ) )
C=========================================
C=====Fin de boucle sur les sous-réseaux 3
C=========================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-réseaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C---------------------------------------------
C-----Eléments 1, 2 et 3 (ce dernier éventuel)
C---------------------------------------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_MAILLE_R = X_DP
        END DO
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_2
     $   = SOMME_2_TYP_2
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
           SOMME_MAILLE_R = SOMME_MAILLE_R + X_DP
        END DO
         IF ( N_TYP_INTR .GT. 2 ) THEN
          I_TYP = 3
           DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
            X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
           SOMME_2_TYP_3
     $   = SOMME_2_TYP_3
     $   + DFLOAT ( P_R ( I_R ) ) * X_DP
           SOMME_MAILLE_R = SOMME_MAILLE_R + X_DP
           END DO
         END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
           X_DP = DEXP ( - H_FORM_D_R_Q ( I_Q , I_TYP , I_R ) / K_T )
             SOMME_I ( I_TYP )
     $     = SOMME_I ( I_TYP )
     $     + DFLOAT ( P_R ( I_R ) ) * X_DP
             SOMME_MAILLE_R
     $     = SOMME_MAILLE_R + X_DP
         END DO
        END DO
          SOMME_MAILLE = SOMME_MAILLE
     $                 + SOMME_MAILLE_R * DFLOAT ( P_R ( I_R ) )
C=================================================================
C=====Fin de boucle optionnelle sur les sous-réseaux interstitiels
C=================================================================
        END DO
       END IF
C=====================================================
C=====Contributions optionnelles des complexes chargés
C=====================================================
C---------------------------------------------
C-----Initialisation des termes de sommes
C-----liés aux complexes chargés,
C-----(utiles au calcul des fractions atomiques,
C-----donc initialisation toujours faite,
C-----même sans complexes chargés)
C---------------------------------------------
        ALLOCATE ( SOMME_COMPLEXE_V_TYPE_I ( N_TYP_INTR + 1 : N_TYP ) )
        SOMME_COMPLEXE_U_TYPE_2 = 0.D0
        SOMME_COMPLEXE_V_TYPE_2 = 0.D0
        SOMME_COMPLEXE_U_TYPE_3 = 0.D0
        SOMME_COMPLEXE_V_TYPE_3 = 0.D0
        SOMME_COMPLEXE_V_TYPE_I = 0.D0
        SOMME_COMPLEXE_U_TOTALE = 0.D0
        SOMME_COMPLEXE_V_TOTALE = 0.D0
C	  stop
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' ) 
     $ THEN
C-------------------------------------------------
C-------------------------------------------------
C=====Boucle optionnelle sur les complexes chargés
C-------------------------------------------------
C-------------------------------------------------
       DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
C-------------------------------------------------------------------
C-----Fractions du complexe courant pour les divers états de charge,
C-----et cumul sur les charges pour calcul des fractions atomiques
C-------------------------------------------------------------------
       SOMME_Q_COUR = 0.D0
       DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
          X_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $  = DFLOAT ( MULTIPLICITE_COMPLEXE_Q ( I_COMPLEXE_Q ) )
     $  * DEXP ( - H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) / K_T )   
         SOMME_Q_COUR = SOMME_Q_COUR
     $ + X_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
       END DO
C	stop
C-------------------------------------------------------
C-----Termes de sommes relatifs à l'espèce 2 intrinsèque
C-----pour le calcul des fractions atomiques
C-------------------------------------------------------
       IF ( N_TYP_INTR .GT. 1 ) THEN
         SOMME_INTER_U = 0.D0
         DO I_R = N_R_1 + 1 , N_R_1 + N_R_2
          SOMME_INTER_U
     $  = SOMME_INTER_U 
     $  + DFLOAT ( U_COMPLEXE_S_R_Q ( I_COMPLEXE_Q , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $           ( I_COMPLEXE_Q ) ) )
         END DO
         SOMME_COMPLEXE_U_TYPE_2
     $ = SOMME_COMPLEXE_U_TYPE_2
     $ + SOMME_INTER_U
     $ * SOMME_Q_COUR
         SOMME_COMPLEXE_V_TYPE_2
     $ = SOMME_COMPLEXE_V_TYPE_2
     $  + DFLOAT ( V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , 2 ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $           ( I_COMPLEXE_Q ) ) )
     $ * SOMME_Q_COUR
       END IF
C	stop
C-------------------------------------------------------
C-----Termes de sommes relatifs à l'espèce 3 intrinsèque
C-----pour le calcul des fractions atomiques
C-------------------------------------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
         SOMME_INTER_U = 0.D0
         DO I_R = N_R_1 + N_R_2 + 1 , N_R_1 + N_R_2 + N_R_3
          SOMME_INTER_U
     $  = SOMME_INTER_U
     $  + DFLOAT ( U_COMPLEXE_S_R_Q ( I_COMPLEXE_Q , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $           ( I_COMPLEXE_Q ) ) )
         END DO
         SOMME_COMPLEXE_U_TYPE_3
     $ = SOMME_COMPLEXE_U_TYPE_3
     $ + SOMME_INTER_U
     $ * SOMME_Q_COUR
         SOMME_COMPLEXE_V_TYPE_3
     $ = SOMME_COMPLEXE_V_TYPE_3
     $  + DFLOAT ( V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , 3 ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $  ( I_COMPLEXE_Q ) ) )
     $ * SOMME_Q_COUR
       END IF
C------------------------------------------------------
C-----Termes de sommes relatifs aux éléments d'addition
C-----pour le calcul des fractions atomiques
C------------------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
         SOMME_COMPLEXE_V_TYPE_I ( I_TYP )
     $ = SOMME_COMPLEXE_V_TYPE_I ( I_TYP )
     $  + DFLOAT ( V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , I_TYP ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $  ( I_COMPLEXE_Q ) ) )
     $ * SOMME_Q_COUR
       END DO
C---------------------------------------------------------------
C-----Termes de sommes relatifs aux quantités de matière totales
C---------------------------------------------------------------
         SOMME_INTER_U = 0.D0
         DO I_R = 1 , N_R_1 + N_R_2 + N_R_3
          SOMME_INTER_U
     $  = SOMME_INTER_U
     $  + DFLOAT ( U_COMPLEXE_S_R_Q ( I_COMPLEXE_Q , I_R ) )
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $  ( I_COMPLEXE_Q ) ) )
         END DO
         SOMME_COMPLEXE_U_TOTALE
     $ = SOMME_COMPLEXE_U_TOTALE
     $ + SOMME_INTER_U
     $ * SOMME_Q_COUR
         SOMME_INTER_V = 0.D0
         DO I_TYP = 1 , N_TYP
          SOMME_INTER_V
     $  = SOMME_INTER_V
     $  + DFLOAT ( V_COMPLEXE_TYPE_Q ( I_COMPLEXE_Q , I_TYP ) )
         END DO
         SOMME_COMPLEXE_V_TOTALE
     $ = SOMME_COMPLEXE_V_TOTALE
     $  + SOMME_INTER_V
     $  * DFLOAT ( P_R ( I_S_R_MULTIPLICITE_COMPLEXE_Q
     $  ( I_COMPLEXE_Q ) ) )
     $ * SOMME_Q_COUR
C-----------------------------------------------
C-----------------------------------------------
C=====Fin de la boucle sur les complexes chargés
C-----------------------------------------------
C-----------------------------------------------
       END DO
C=================================================
C=====Fin du test de présence de complexes chargés
C=================================================
      END IF
C========================
C========================
C=====Fractions atomiques
C========================
C========================
      write (*,*) SOMME_COMPLEXE_U_TYPE_2
      write (*,*) SOMME_COMPLEXE_V_TYPE_2
      write (*,*) SOMME_COMPLEXE_U_TYPE_3
      write (*,*) SOMME_COMPLEXE_V_TYPE_3
      write (*,*) SOMME_COMPLEXE_U_TOTALE
      write (*,*) SOMME_COMPLEXE_V_TOTALE
      IF ( N_TYP_INTR .GT. 1 ) THEN
C---------------------------
C-----Fraction atomique de 2
C---------------------------
           X_AT ( 2 )
     $ = ( SOMME_2_TYP_2 + SOMME_1_TYP_2
     $   - SOMME_COMPLEXE_U_TYPE_2 + SOMME_COMPLEXE_V_TYPE_2 )
     $ / ( SOMME_0 + SOMME_MAILLE 
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
         IF ( X_AT ( 2 ) .LT. 0.D0 ) X_AT ( 2 ) = 0.D0
         IF ( N_TYP_INTR .GT. 2 ) THEN
C---------------------------
C-----Fraction atomique de 3
C---------------------------
           X_AT ( 3 )
     $ = ( SOMME_2_TYP_3 + SOMME_1_TYP_3  
     $   - SOMME_COMPLEXE_U_TYPE_3 + SOMME_COMPLEXE_V_TYPE_3 )
     $ / ( SOMME_0 + SOMME_MAILLE 
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
           IF ( X_AT ( 3 ) .LT. 0.D0 ) X_AT ( 3 ) = 0.D0
        END IF
      END IF
C------------------------------------------------
C-----Fractions atomiques des éléments d'addition  
C------------------------------------------------
       DO I_TYP = N_TYP_INTR + 1 , N_TYP
           X_AT ( I_TYP )
     $ = ( SOMME_I ( I_TYP ) + SOMME_COMPLEXE_V_TYPE_I ( I_TYP ) )
     $ / ( SOMME_0 + SOMME_MAILLE 
     $   - SOMME_COMPLEXE_U_TOTALE + SOMME_COMPLEXE_V_TOTALE )
        IF ( X_AT ( I_TYP ) .LT. 0.D0 ) X_AT ( I_TYP ) = 0.D0
       END DO
C---------------------------
C-----Fraction atomique de 1
C---------------------------
       X_AT ( 1 ) = 1.D0
       DO I_TYP = 2 , N_TYP
         X_AT ( 1 ) = X_AT ( 1 ) - X_AT ( I_TYP )
       END DO
       IF ( X_AT ( 1 ) .LT. 0.D0 ) X_AT ( 1 ) = 0.D0
C	write ( * , * ) ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
C----------------------------------------------------------------
C-----Fin du test "mode 2" pour le calcul des fractions atomiques
C----------------------------------------------------------------
      END IF
      WRITE ( * , * ) '==================='
      WRITE ( * , * ) 'Fractions atomiques :'
      WRITE ( * , * ) '==================='
      DO I_TYP = 1 , N_TYP
        WRITE ( * , * ) I_TYP , X_AT ( I_TYP )
      END DO
C      write (*,*) 'SOMME_1_TYP_2 = ' , SOMME_1_TYP_2
C      write (*,*) 'SOMME_2_TYP_2 = ' , SOMME_2_TYP_2
C      write (*,*) 'SOMME_1_TYP_3 = ' , SOMME_1_TYP_3
C      write (*,*) 'SOMME_2_TYP_3 = ' , SOMME_2_TYP_3
C      write (*,*) 'SOMME_0 = ' , SOMME_0
C      write (*,*) 'SOMME_MAILLE = ', SOMME_MAILLE
C======================================================
C======================================================
C=====Ecriture des H_form pour la neutralité électrique
C======================================================
C======================================================
      IF ( I_CALC_CHARGE .EQ. 2 ) THEN
       WRITE ( * , * ) '    ======================================'
       WRITE ( * , * ) '    Enthalpies de formation (eV) "n+A=p+D"'
       WRITE ( * , 3005 ) TEMPERATURE
       WRITE ( * , * ) '    ======================================'
       WRITE ( * , * ) '    type     s-r      q         Hf(eV)'
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C==================================
C=====Boucle sur les sous-réseaux 1
C==================================
      DO I_R = 1 , N_R_1
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C----------------
C-----Antisites 2
C----------------
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        I_TYP = 3
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 1
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-réseaux 2
C==================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C----------------
C-----Antisites 1
C----------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C--------------------------
C-----Antisites 3 éventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
        I_TYP = 3
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
      END IF
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 2
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-réseaux 3
C==================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
C------------
C-----Lacunes
C------------
        I_TYP = 0
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C----------------
C-----Antisites 1
C----------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C----------------
C-----Antisites 2
C----------------
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
C--------------------------
C-----Eléments > N_TYP_INTR
C--------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 3
C=========================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-réseaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C---------------------------------------------
C-----Eléments 1, 2 et 3 (ce dernier éventuel)
C---------------------------------------------
        I_TYP = 1
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
        I_TYP = 2
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
         IF ( N_TYP_INTR .GT. 2 ) THEN
        I_TYP = 3
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
         END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
        DO I_Q = 1 , NQ_D_R_Q ( I_TYP , I_R )
         write ( * , 3010 )
     $      I_TYP , I_R , Q_D_R_Q ( I_Q , I_TYP , I_R ) ,
     $      H_FORM_D_R_Q ( I_Q , I_TYP , I_R )
        END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
         END DO
C=================================================================
C=====Fin de boucle optionnelle sur les sous-réseaux interstitiels
C=================================================================
        END DO
       END IF
C====================================================
C=====Ecriture optionnelle pour les complexes chargés
C====================================================
       IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $ .OR. INDIC_COMPLEXES_Q .EQ. 'o' ) THEN
       WRITE ( * , * ) '      - -- - - - - - - - - - - -'
       WRITE ( * , * ) '    Complexes chargés :'
       WRITE ( * , * ) '    type        q         Hf(eV)'
       WRITE ( * , * ) '      - -- - - - - - - - - - - -'
        DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
         DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
         write ( * , 3020 )
     $      I_COMPLEXE_Q , Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) ,
     $      H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
          END DO
       WRITE ( * , * ) '      - - - - - - - - - - - - - - - - -'
        END DO
C==========================================================
C=====Fin d'écriture optionnelle pour les complexes chargés
C==========================================================
       END IF
C---------------------------------------
C-----Fin de l'écriture des Hf en mode 2
C---------------------------------------
      END IF
      WRITE ( * , * )
      WRITE ( * , * ) '          ========================'
      WRITE ( * , * ) '          Calcul ADPI (DP chargés)'
      WRITE ( * , * ) '          ------------------------'
      WRITE ( * , * ) '              FIN DU PROGRAMME'
      WRITE ( * , * ) '          ========================'
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&& Fin du cas "DP CHARGES"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C-----Instruction de branchement du cas "DP NON CHARGES"
 8888 CONTINUE
      END
C############################
C#####Inclusion de procédures
C############################
C----------------------------
C-----Références du programme
C----------------------------
C      INCLUDE
C     $ 'ref_prog.inc'
C     ==================================================================
C     I                                                                I
C     I Référence de programme (laboratoire, auteur)		       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 27/10/2005                              I
C     I Sainte Emeline						       I  
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   REF_PROG
C     ##################################################################
      WRITE ( * , * ) 
     $ '###############################################################'
      WRITE ( * , * ) 
     $ '#            Ce programme est la propriété du                 #'
      WRITE ( * , * ) 
     $ '# Laboratoire de Métallurgie Physique et Génie des Matériaux  #'
      WRITE ( * , * ) 
     $ "#     Université de Lille 1, Villeneuve d'Ascq - FRANCE       #"
      WRITE ( * , * ) 
     $ '#                ----------------------------                 #'
      WRITE ( * , * ) 
     $ '#    Auteur : Rémy Besson (Remy.Besson@univ-lille1.fr)        #'
      WRITE ( * , * ) 
     $ '###############################################################'
      WRITE ( * , * )
     $ "Veuillez presser la touche Entrée pour lancer l'exécution." 
      WRITE ( * , * )
     $ '----------------------------------------------------------'
      READ ( * , * )  
      RETURN
      END
C       INCLUDE
C     $ 'ref_prog_direct.inc'
C------------------------------
C-----Interruption du programme
C------------------------------
C      INCLUDE
C     $ 'interruption.inc'
C     ==================================================================
C     I          						       I
C     I Interruption du programme			   	       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 30/05/2003                              I
C     ==================================================================
C     ##################################################################
      SUBROUTINE INTERRUPTION
C     ##################################################################
      WRITE ( * , * ) '--------------------'
      WRITE ( * , * ) 'PROGRAMME INTERROMPU'
      WRITE ( * , * ) '--------------------'
      STOP
      RETURN
      END
C--------------------
C-----Format variable
C--------------------
C      INCLUDE
C     $ 'form_var_adpi.inc'
C     ==================================================================
C     I                                                                I
C     I Constitution d'un format variable d'écriture de réels          I
C     I et du titre associé pour le calcul ADPI			       I
C     I (potentiel chimique de chaque espèce,			       I
C     I fraction atomique de chaque espèce,			       I
C     I fraction de DP sur le sous-réseau considéré		       I
C     I et enthalpie de formation de ce DP)			       I
C     I et idem pour le fichier d'énergie libre		               I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 17/06/2005                              I
C     I Saint Hervé                                                    I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   FORMAT_VARIABLE_ADPI
     $ ( N_TYP_M , INDIC_AT_MAILLE_M , INDIC_G_M ,
     $   CAR_COL_VAR_X_DP_M , CAR_TITRE_VAR_X_DP_M ,
     $   CAR_COL_VAR_E_L_M , CAR_TITRE_VAR_E_L_M ,
     $   CAR_COL_VAR_POT_CHIM_M , CAR_TITRE_VAR_POT_CHIM_M ,
     $   CAR_COL_VAR_PB_CONV_BA_M , CAR_TITRE_VAR_PB_CONV_BA_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C-----------------------------------------------------
C-----Caractères utiles à la constitution d'un format   
C-----pour l'écriture d'un nombre variable de colonnes
C-----(potentiel chimique de chaque espèce,
C-----fraction atomique de chaque espèce,
C-----fraction de DP sur le sous-réseau considéré
C-----et enthalpie de formation de ce DP)
C-----------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_X_DP_M
C----------------------------
C-----Format du titre associé
C----------------------------
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP_M
C-----------------------------------------
C-----Idem pour le fichier d'énergie libre
C-----------------------------------------
      CHARACTER * 300 CAR_COL_VAR_E_L_M
      CHARACTER * 5000 CAR_TITRE_VAR_E_L_M
C-------------------------------------------------
C-----Idem pour le fichier de potentiels chimiques
C-------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_POT_CHIM_M
      CHARACTER * 5000 CAR_TITRE_VAR_POT_CHIM_M
C------------------------------------------
C-----Idem pour le fichier de "pb conv. BA"
C------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_PB_CONV_BA_M
      CHARACTER * 5000 CAR_TITRE_VAR_PB_CONV_BA_M
C---------------------------------------------------------
C-----Indicateur d'écriture des grandeurs thermodynamiques
C-----par atome (A/a) ou par maille (M/m)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_AT_MAILLE_M
C---------------------------------------------------------
C-----Indicateur d'écriture de l'énergie libre par atome :
C-----énergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_G_M
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C====================================
C=====Formats des colonnes de valeurs
C====================================
      CAR_COL_VAR_X_DP_M = '( 1X ,     ( G18.8E2 , 2X ) )'
      WRITE ( CAR_COL_VAR_X_DP_M ( 8 : 10 ) , '(I3)' ) 2 * N_TYP_M + 3
      CAR_COL_VAR_E_L_M = '( 1X ,     ( G14.7E3 , 3X ) )'
      WRITE ( CAR_COL_VAR_E_L_M ( 8 : 10 ) , '(I3)' )  N_TYP_M + 5
      CAR_COL_VAR_POT_CHIM_M = '( 1X ,     ( G18.8E2 , 2X ) )'
      WRITE ( CAR_COL_VAR_POT_CHIM_M ( 8 : 10 ) , '(I3)' )
     $                                             2 * N_TYP_M + 1
      CAR_COL_VAR_PB_CONV_BA_M
     $ = '( 1X ,     ( G18.8E2 , 2X ) , 12X , I1 )'
      WRITE ( CAR_COL_VAR_PB_CONV_BA_M ( 8 : 10 ) , '(I3)' )
     $                                                 N_TYP_M
C=======================
C=====Formats des titres
C=======================
      CAR_TITRE_VAR_X_DP_M = ' '
      CAR_TITRE_VAR_E_L_M = ' '
      CAR_TITRE_VAR_POT_CHIM_M = ' '
      CAR_TITRE_VAR_PB_CONV_BA_M = ' '
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fichiers de défauts ponctuels et de potentiels chimiques
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C------------------------------------
C------------------------------------
C-----Titres = "potentiels chimiques"
C------------------------------------
C------------------------------------
      DO I_TYP = 1 , N_TYP_M
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
          L_2 = LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA_M )
          IF ( I_TYP .EQ. 1 ) THEN
             WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 9 : L + 10 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 3 : L + 4 ) , '(A)' )
     $   'mu'
             WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 9 : L_1 + 10 ) ,
C           WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 3 : L_1 + 4 ) ,
     $ '(A)' ) 'mu'
            WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 4 : L_2 + 5 ) ,
     $ '(A)' ) 'mu'
          ELSE IF ( I_TYP .LE. 9 ) THEN
          WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 13 : L + 14 ) , '(A)' )
C        WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 3 : L + 4 ) , '(A)' )
     $  'mu'
C	    IF ( I_TYP .EQ. 2 ) THEN
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 13 : L_1 + 14 ) ,
     $ '(A)' ) 'mu'
            WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 13 : L_2 + 14 ) ,
C           WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 3 : L_1 + 4 ) ,
     $ '(A)' ) 'mu'
C	    ELSE
C            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 3 : L_1 + 4 ) ,
C     $ '(A)' ) 'mu'
C	    END IF
          ELSE IF ( I_TYP .LE. 99 ) THEN
             WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 7 : L + 8 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 3 : L + 4 ) , '(A)' )
     $    'mu'
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 3 : L_1 + 4 ) ,
     $ '(A)' ) 'mu'
            WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 3 : L_2 + 4 ) ,
     $ '(A)' ) 'mu'
          ELSE
           WRITE ( * , * ) '---------------------------'
           WRITE ( * , * ) 'Nombre de types limité à 99'
           WRITE ( * , * ) '---------------------------'
           CALL INTERRUPTION
          END IF
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
          L_2 = LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA_M )
          IF ( I_TYP .LE. 9 )  THEN
            WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 1 : L + 1 ) , '(I1)' )
     $              I_TYP
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 1 : L_1 + 1 ) ,
     $ '(I1)' ) I_TYP
            WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 1 : L_2 + 1 ) ,
     $ '(I1)' ) I_TYP
          ELSE 
            WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 1 : L + 2 ) , '(I2)' )
     $              I_TYP
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 1 : L_1 + 2 ) ,
     $ '(I2)' ) I_TYP
           WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 1 : L_2 + 2 ) ,
     $ '(I2)' ) I_TYP
          END IF
C-----Ecriture ci-dessous seulement pour le fichier "pb conv. BA"
          L_2 = LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA_M )
          IF ( I_TYP .EQ. 1 ) THEN
            WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 1 : L_2 + 6 ) ,
     $ '(A)' ) '(init)'
          END IF
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
          L_2 = LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA_M )
          WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 1 : L + 4 ) , '(A)' )
     $    '(eV)'
          WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 1 : L_1 + 4 ) ,
     $    '(A)' ) '(eV)'
         WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 1 : L_2 + 4 ) ,
     $    '(A)' ) '(eV)'
      END DO
C-----------------------------------
C-----------------------------------
C-----Titres = "fractions atomiques"
C-----------------------------------
C-----------------------------------
      DO I_TYP = 1 , N_TYP_M
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
C----------------------------------------------------------------------
C-----Le cas I_TYP = 1 doit être distingué car le texte qui le précède
C-----mu(N_TYP) est différent de celui des autres types x_at(I_TYP - 1)
C----------------------------------------------------------------------
           IF ( I_TYP .EQ. 1 ) THEN
             WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 15 : L + 18 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 4 : L + 7 ) , '(A)' )
     $             'x_at'
             WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 15 : L_1 + 18 ) ,
     $ '(A)' ) 'x_at'
C           WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 4 : L_1 + 7 ) ,
C    $ '(A)' ) 'x_at'
          ELSE IF ( I_TYP .LE. 9 ) THEN
             WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 16 : L + 19 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 4 : L + 7 ) , '(A)' )
     $             'x_at'
             WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 16 : L_1 + 19 ) ,
     $ '(A)' ) 'x_at'
C          WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 4 : L_1 + 7 ) ,
C    $ '(A)' ) 'x_at'
           ELSE
             WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 16 : L + 19 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 4 : L + 7 ) , '(A)' )
     $             'x_at'
             WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 16 : L_1 + 19 ) ,
     $ '(A)' ) 'x_at'
C           WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 4 : L_1 + 7 ) ,
C    $ '(A)' ) 'x_at'
           END IF
C         END IF
C         write ( * , * )
C    $    CAR_TITRE_VAR_X_DP_M ( 1 : LEN_TRIM ( CAR_TITRE_VAR_X_DP_M ) )
C---------------------------------------------
C-----Indices de types des fractions atomiques
C---------------------------------------------
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
C         write ( * , * ) L
          IF ( I_TYP .LE. 9 )  THEN
            WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 1 : L + 1 ) , '(I1)' )
     $              I_TYP
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 1 : L_1 + 1 ) ,
     $ '(I1)' ) I_TYP
          ELSE
            WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 1 : L + 2 ) , '(I2)' )
     $              I_TYP
            WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 1 : L_1 + 2 ) ,
     $ '(I2)' ) I_TYP
          END IF
          L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
          L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
C         write ( * , * ) L
      END DO
C     write ( * , * )
C     CAR_TITRE_VAR_X_DP_M ( 1 : LEN_TRIM ( CAR_TITRE_VAR_X_DP_M ) )
C--------------------------
C--------------------------
C-----Titre = "température"
C--------------------------
C--------------------------
      L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
      WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 12 : L + 25 ) , '(A)' )
C     WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 2 : L + 15 ) , '(A)' )
     $       'Température(K)'
      L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
      WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 12 : L_1 + 25 ) ,
C     WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 2 : L_1 + 15 ) ,
     $  '(A)' ) 'Température(K)'
C----------------------------
C----------------------------
C-----Titre = "indic. pb. BA"
C----------------------------
C----------------------------
      L_2 = LEN_TRIM ( CAR_TITRE_VAR_PB_CONV_BA_M )
      WRITE ( CAR_TITRE_VAR_PB_CONV_BA_M ( L_2 + 15 : L_2 + 32 ) ,
     $  '(A)' ) 'Indic. pb conv. BA'
C-----------------------------------------------------------
C-----------------------------------------------------------
C-----Titres = "fractions et enthalpies de formation des DP"
C-----------------------------------------------------------
C-----------------------------------------------------------
      L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
       WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 10 : L + 13 ) , '(A)' ) 'x_DP'
C     WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 2 : L + 5 ) , '(A)' ) 'x_DP'
      L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
       WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 12 : L + 22 ) , '(A)' )
C     WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 7 : L + 16 ) , '(A)' )
     $       'H_f(eV)_DP'
C      write ( * , * )
C     $CAR_TITRE_VAR_X_DP_M ( 1 : LEN_TRIM ( CAR_TITRE_VAR_X_DP_M ) )
C= = = = = = = = = = = = = = 
C= = =Fichier d'énergie libre
C= = = = = = = = = = = = = = 
C------------------------
C-----Fractions atomiques
C------------------------
      DO I_TYP = 1 , N_TYP_M
          L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
           IF ( I_TYP .EQ. 1 ) THEN
            WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 5 ) , '(A)' )
     $             'x_at'
           ELSE IF ( I_TYP .LE. 9 ) THEN
             WRITE ( CAR_TITRE_VAR_E_L_M ( L + 11 : L + 14 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_E_L_M ( L + 3 : L + 6 ) , '(A)' )
     $             'x_at'
           ELSE
            WRITE ( CAR_TITRE_VAR_E_L_M ( L + 10 : L + 13 ) , '(A)' )
C           WRITE ( CAR_TITRE_VAR_E_L_M ( L + 3 : L + 6 ) , '(A)' )
     $             'x_at'
           END IF
          L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
          IF ( I_TYP .LE. 9 )  THEN
            WRITE ( CAR_TITRE_VAR_E_L_M ( L + 1 : L + 1 ) , '(I1)' )
     $              I_TYP
          ELSE
            WRITE ( CAR_TITRE_VAR_E_L_M ( L + 1 : L + 2 ) , '(I2)' )
     $              I_TYP
          END IF
      END DO
C----------------
C-----Température
C----------------
      L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
      WRITE ( CAR_TITRE_VAR_E_L_M ( L + 12 : L + 25 ) , '(A)' )
C     WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 15 ) , '(A)' )
     $       'Température(K)'
C----------------------------------------------------------------
C-----Energie, volume, entropie de configuration et énergie libre
C----------------------------------------------------------------
      L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
      IF ( INDIC_AT_MAILLE_M .EQ. 'M' .OR. INDIC_AT_MAILLE_M .EQ. 'm' )
     $ THEN
C       WRITE ( CAR_TITRE_VAR_E_L_M ( L + 10 : L + 21 ) , '(A)' )
       WRITE ( CAR_TITRE_VAR_E_L_M ( L + 5 : L + 16 ) , '(A)' )
     $        'E(eV/maille)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
       WRITE ( CAR_TITRE_VAR_E_L_M ( L + 7 : L + 20 ) , '(A)' )
     $        'V(A**3/maille)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
        WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 20 ) , '(A)' ) 
     $        'S_conf(eV/K/maille)'
C      WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 20 ) , '(A)' ) 
C    $        'S_conf(eV/K/maille)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
        WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 13 ) , '(A)' ) 
     $        'G(eV/maille)'
C      WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 13 ) , '(A)' ) 
C    $        'G(eV/m)'
      ELSE
        WRITE ( CAR_TITRE_VAR_E_L_M ( L + 8 : L + 15 ) , '(A)' )
C      WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 9 ) , '(A)' )
     $        'E(eV/at)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
       WRITE ( CAR_TITRE_VAR_E_L_M ( L + 8 : L + 17 ) , '(A)' )
     $        'V(A**3/at)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
       WRITE ( CAR_TITRE_VAR_E_L_M ( L + 8 : L + 22 ) , '(A)' )
     $        'S_conf(eV/K/at)'
       L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
       IF ( INDIC_G_M .EQ. 'T' .OR. INDIC_G_M .EQ. 't' ) THEN
C         WRITE ( CAR_TITRE_VAR_E_L_M ( L + 10 : L + 21 ) , '(A)' )
         WRITE ( CAR_TITRE_VAR_E_L_M ( L + 3 : L + 14 ) , '(A)' )
     $          'G_tot(eV/at)'
        ELSE
C         WRITE ( CAR_TITRE_VAR_E_L_M ( L + 5 : L + 17 ) , '(A)' )
          WRITE ( CAR_TITRE_VAR_E_L_M ( L + 3 : L + 15 ) , '(A)' )
     $          'G_form(eV/at)'
        END IF
      END IF
      RETURN
      END
C------------------------------------------------
C-----Factorisation LU d'une matrice inversible
C-----(déclaré également dans sys_non_lin_NR.inc)
C------------------------------------------------
C     INCLUDE
C    $ 'fact_lu.inc'
C------------------------------------------------
C-----Inversion d'une matrice par la méthode LU
C-----(déclaré également dans sys_non_lin_NR.inc)
C------------------------------------------------
C     INCLUDE
C    $ 'inverse_lu.inc'
C-----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de sous sous-réseau et de son type
C-----------------------------------------------------------------------
C      INCLUDE
C     $ 'ind_D_R_typ.inc'
C     ==================================================================
C     I                                                                I
C     I Attribution d'un indice à chaque DP d'un alliage ordonné       I
C     I unaire, binaire ou ternaire avec éléments d'addition   	       I	
C     I en fonction de son type et de son sous-réseau                  I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 07/06/2005                              I
C     I Saint Gilbert						       I
C     ==================================================================
C     ##################################################################
         SUBROUTINE 
     $   INDICE_D_R_TYP
     $ ( N_TYP_M , N_TYP_INTR_M , 
     $   N_R_M , N_R_1_M , N_R_2_M , N_R_3_M ,
     $   I_R_INTER_M , I_R_INTER_INTR_M ,
     $   E_B_D_R_M ,
     $   N_TYP_DP_M , IND_DP_TYP_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C-------------------------------------------------------------------
C-----Indice du défaut en fonction de son sous-réseau et de son type
C-------------------------------------------------------------------
      INTEGER * 4 IND_DP_TYP_M ( 0 : N_TYP_M , N_R_M )
C---------------------------------------
C-----Energie de SC de chaque type de DP
C---------------------------------------
      REAL * 8 E_B_D_R_M ( 0 : N_TYP_M , N_R_M )
C---------------------------------------------------------
C-----Indicateur de présence de sous-réseaux interstitiels
C---------------------------------------------------------
      CHARACTER * 1 I_R_INTER_M
C-----------------------------------------------------------------
C-----Indicateur de prise en compte des interstitiels intrinsèques
C-----------------------------------------------------------------
      CHARACTER * 1 I_R_INTER_INTR_M
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C-------------------
C-----Initialisation
C-------------------
      IND_DP_TYP_M  = 0
C-----------------------------------------
C-----Nombre de sous-réseaux interstitiels
C-----------------------------------------
      N_R_INTER = N_R_M - N_R_1_M - N_R_2_M - N_R_3_M
C--------------------------------
C-----Indice de chaque type de DP
C--------------------------------
      N_TYP_DP_M = 0
      DO I_R = 1 , N_R_M
        DO I_TYP = 0 , N_TYP_M
C----------------------------------------------------
C-----Indicateur de prise en compte de l'interstitiel
C-----pour le sous-réseau et le type courant
C----------------------------------------------------
       INDIC_INTER = 0
       IF ( I_TYP .NE. 0
     $ .AND. ( I_R_INTER_M .EQ. 'O' .OR. I_R_INTER_M .EQ. 'o ') )
     $  THEN
        IF ( I_TYP .GT. N_TYP_INTR_M ) THEN
          INDIC_INTER = 1
        ELSE IF ( ( I_TYP .LE. N_TYP_INTR_M )
     $ .AND. ( I_R_INTER_INTR_M .EQ. 'O'
     $    .OR. I_R_INTER_INTR_M .EQ. 'o') )
     $  THEN
          INDIC_INTER = 1
          END IF
         END IF
        IF ( ( I_R .LE. N_R_1_M
     $   .AND. I_TYP .NE. 1 )
     $  .OR. ( I_R .GT. N_R_1_M
     $   .AND. I_R .LE. N_R_1_M + N_R_2_M
     $   .AND. I_TYP .NE. 2 )
     $  .OR. ( I_R .GT. N_R_1_M + N_R_2_M
     $   .AND. I_R .LE. N_R_1_M + N_R_2_M + N_R_3_M
     $   .AND. I_TYP .NE. 3 )
     $  .OR. ( I_R .GT. N_R_1_M + N_R_2_M + N_R_3_M
     $   .AND. INDIC_INTER .EQ. 1 ) ) THEN
C-----Test ci-dessous pour écarter les DP "non utiles"
C-----(e.g. certains interstitiels coûteux)
           IF ( E_B_D_R_M ( I_TYP , I_R ) .LT. 0.D0 ) THEN
                N_TYP_DP_M = N_TYP_DP_M + 1
                IND_DP_TYP_M ( I_TYP , I_R ) = N_TYP_DP_M
           END IF
         END IF
        END DO
      END DO
      RETURN
      END
C------------------------------------
C-----Procédures du simplexe linéaire
C------------------------------------
C      INCLUDE
C     $ 'simplexe_lin.inc'
C     ==================================================================
C     I                                                                I
C     I Procédures du simplexe linéaire				       I
C     I (tirées des "recettes numériques" 			       I
C     I version 2005)						       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 02/12/2005                              I
C     I Sainte Viviane 						       I
C     ==================================================================
C     ##################################################################
      SUBROUTINE SIMPLX ( A , M , N ,
     $                    MP , NP ,
     $                    M1 , M2 , M3 ,
     $                    ICASE , IZROV , IPOSV )
C     ##################################################################
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 IPOSV ( M ) , IZROV ( N )
      REAL * 8 EPS   
      INTEGER * 4 ICASE , M , M1 , M2 , M3 , MP , N , NP
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      INTEGER * 4 MMAX , NMAX
      PARAMETER ( MMAX = 100 , NMAX = 100 , EPS = 1.D-6 )
      REAL * 8 BMAX , Q1
      INTEGER * 4 L1 ( NMAX ) , L3 ( MMAX ) 
      INTEGER * 4 I , IP , IS , K , KH , KP , NL1
      IF ( M .NE. M1 + M2 + M3 ) THEN
       write ( * , * ) 'Erreur nombres de contraintes'
        CALL INTERRUPTION
      END IF
      NL1 = N
      DO K = 1 , N
        L1 ( K ) = K
        IZROV ( K ) = K
      END DO
      DO I = 1 , M
        IF ( A ( I + 1 , 1 ) .LT. 0.D0 ) THEN
         write ( * , * )  'Erreur tableau initial'
         CALL INTERRUPTION
        END IF
        IPOSV ( I ) = N + I
      END DO
      IF ( M2 + M3 .EQ. 0 ) GO TO 30
      DO I = 1 , M2
        L3 ( I ) = 1
      END DO
      DO K = 1 , N + 1
        Q1 = 0.D0
        DO I = M1 + 1 , M
          Q1 = Q1 + A ( I + 1 , K )
        END DO
        A ( M + 2 , K ) = - Q1
      END DO
10    CALL SIMP1 ( A , MP , NP , M + 1 ,
     $             L1 , NL1 , 0 , KP , BMAX )
      IF ( BMAX .LE. EPS .AND. A ( M + 2 , 1 ) .LT. - EPS ) THEN
        ICASE = - 1
        RETURN
      ELSE IF ( BMAX .LE. EPS .AND. A ( M + 2 , 1 ) .LE. EPS ) THEN
          DO IP = M1 + M2 + 1 , M
            IF ( IPOSV ( IP ) .EQ. IP + N ) THEN
              CALL SIMP1 ( A , MP , NP , IP ,
     $                     L1 , NL1 , 1 , KP , BMAX )
              IF ( BMAX .GT. EPS ) GO TO 1
            ENDIF
          END DO
        DO I = M1 + 1 , M1 + M2
          IF ( L3 ( I - M1 ) .EQ. 1 ) THEN
            DO K = 1 , N + 1
              A ( I + 1 , K ) = - A ( I + 1 , K )
           END DO 
          ENDIF
        END DO
        GO TO 30
      ENDIF
      CALL SIMP2 ( A , M , N , MP , NP , IP , KP )
      IF ( IP .EQ. 0 ) THEN
        ICASE = - 1
        RETURN
      ENDIF
1     CALL SIMP3 ( A , MP , NP , M + 1 , N , IP , KP )
      IF ( IPOSV ( IP ) .GE. N + M1 + M2 + 1 ) THEN
        DO K = 1 , NL1
          IF ( L1 ( K ) .EQ. KP ) GO TO 2
        END DO
2       NL1 = NL1 - 1
        DO IS = K , NL1
          L1 ( IS ) = L1 ( IS + 1 )
        END DO
      ELSE
        KH = IPOSV ( IP ) - M1 - N
        IF ( KH .GE. 1 ) THEN
          IF ( L3 ( KH ) .NE. 0 ) THEN
            L3 ( KH ) = 0
            A ( M + 2 , KP + 1 ) = A ( M + 2 , KP + 1 ) + 1.D0
            DO I = 1 , M + 2
             A ( I , KP + 1 ) = - A ( I , KP + 1 )
            END DO
          END IF  
        END IF
      ENDIF
      IS = IZROV ( KP )
      IZROV ( KP ) = IPOSV ( IP )
      IPOSV ( IP ) = IS
      GO TO 10
30    CALL SIMP1 ( A , MP , NP , 0 ,
     $             L1 , NL1 , 0 , KP , BMAX )
      IF ( BMAX .LE. EPS ) THEN
        ICASE = 0
        RETURN
      ENDIF
      CALL SIMP2 ( A , M , N , MP , NP , IP , KP )
      IF ( IP .EQ. 0 ) THEN
        ICASE = 1
        RETURN
      ENDIF
      CALL SIMP3 ( A , MP , NP , M , N , IP , KP )
      IS = IZROV ( KP )
      IZROV ( KP ) = IPOSV ( IP )
      IPOSV ( IP ) = IS
      GO TO 30
      END
C     ##################################################################
      SUBROUTINE SIMP1 ( A , MP , NP , MM ,
     $                   LL , NLL , IABF , KP , BMAX )
C     ##################################################################
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 LL ( NP )
      REAL * 8 BMAX
      INTEGER * 4 IABF , KP , MM , MP , NLL , NP
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      REAL * 8 TEST
      INTEGER * 4 K
      IF ( NLL .LE. 0 ) THEN
        BMAX = 0.D0
      ELSE 
        KP = LL ( 1 )
        BMAX = A ( MM + 1 , KP + 1 )
        DO K = 2 , NLL
          IF ( IABF .EQ. 0 ) THEN
            TEST = A ( MM + 1 , LL ( K ) + 1 ) - BMAX
          ELSE
            TEST = DABS ( A ( MM + 1 , LL ( K ) + 1 ) ) - DABS ( BMAX )
          ENDIF
          IF ( TEST .GT. 0.D0 ) THEN
            BMAX = A ( MM + 1 , LL ( K ) + 1 )
            KP = LL ( K )
          ENDIF
        END DO
      END IF
      RETURN
      END
C     ##################################################################
      SUBROUTINE SIMP2 ( A , M , N , MP , NP , IP , KP )
C     ##################################################################
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 IP , KP , M , MP , N , NP
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      REAL * 8 EPS , Q , Q0 , Q1 , QP
      INTEGER * 4 I , K
      PARAMETER ( EPS = 1.D-6 )
      IP = 0
      DO I = 1 , M
        IF ( A ( I + 1 , KP + 1 ) .LT. - EPS ) GOTO 1
      END DO
      RETURN
1     Q1 = - A ( I + 1 , 1 ) / A ( I + 1 , KP + 1 )
      IP = I
      DO I = IP + 1 , M
        IF ( A ( I + 1 , KP + 1 ) .LT. - EPS ) THEN
         Q = - A ( I + 1 , 1 ) / A ( I + 1 , KP + 1 )
         IF ( Q .LT. Q1 ) THEN
          IP = I
          Q1 = Q
         ELSE IF ( Q .EQ. Q1 ) THEN
          DO K = 1 , N
           QP = - A ( IP + 1 , K + 1 ) / A ( IP + 1 , KP + 1 )
           Q0 = - A ( I + 1 , K + 1 ) / A ( I + 1 , KP + 1 )
           IF ( Q0 .NE. QP ) GOTO 2
          END DO
2         IF ( Q0 .LT. QP ) IP = I
         END IF
        END IF 
      END DO
      RETURN
      END
C     ##################################################################
      SUBROUTINE SIMP3 ( A , MP , NP , I1 , K1 , IP , KP )
C     ##################################################################
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 I1 , IP , K1 , KP , MP , NP
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      REAL * 8 PIV
      INTEGER * 4 II , KK
      PIV = 1.D0 / A ( IP + 1 , KP + 1 )
      DO II = 1 , I1 + 1
          IF ( II - 1 .NE. IP ) THEN
            A ( II , KP + 1 ) = A ( II , KP + 1 ) * PIV
            DO KK = 1 , K1 + 1
              IF ( KK - 1 .NE. KP ) THEN
                A ( II , KK ) = A ( II , KK )
     $                        - A ( IP + 1 , KK ) * A ( II , KP + 1 )
              ENDIF
            END DO
          ENDIF
      END DO
      DO KK = 1 , K1 + 1
        IF ( KK - 1 .NE. KP )
     $  A ( IP + 1 , KK ) = - A ( IP + 1 , KK ) * PIV
      END DO
      A ( IP + 1 , KP + 1 ) = PIV
      RETURN
      END
C----------------------------------------------------------
C-----Analyse d'une ligne de caractères
C-----pour distinction des chaînes séparées par des espaces
C----------------------------------------------------------
C      INCLUDE
C     $ 'analyse_ligne.inc'
C     ==================================================================
C     I                                                                I
C     I Analyse d'une ligne de caractères			       I
C     I pour distinction des chaînes séparées par des espaces          I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 21/12/2004                              I
C     I Saint Thomas  	    	                                       I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   ANALYSE_LIGNE
     $ ( LIGNE_M ,
     $   N_CHAINES_M ,
     $   I_DEBUT_M , I_FIN_M ,
     $   LONG_CHAINE_M , CHAINE_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      CHARACTER * 1000 LIGNE_M
      INTEGER * 4 I_DEBUT_M ( 1000 )
      INTEGER * 4 I_FIN_M ( 1000 )
      INTEGER * 4 LONG_CHAINE_M ( 1000 )
      CHARACTER * 1000 CHAINE_M ( 1000 )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      CHARACTER I_L_PREC
      CHARACTER I_L
C----------------------------------------------------------
C-----Recherche du nombre de chaînes distinctes de la ligne
C-----et du caractère de début de chaque chaîne
C----------------------------------------------------------
      LONG_LIGNE = LEN_TRIM ( LIGNE_M )
      L = 1
      N_CHAINES_M = 0
      I_L_PREC = ' '
      DO WHILE ( L .LE. LONG_LIGNE ) 
        I_L = LIGNE_M ( L : L )
	IF ( I_L .NE. ' ' ) THEN
	 IF ( I_L_PREC .EQ. ' ' ) THEN
	   N_CHAINES_M = N_CHAINES_M + 1
	   I_DEBUT_M ( N_CHAINES_M ) = L
	 END IF 
        ELSE
         IF ( I_L_PREC .NE. ' ' ) THEN
	   I_FIN_M ( N_CHAINES_M ) = L - 1
	 END IF
        END IF
        I_L_PREC = I_L
	L = L + 1
      END DO
      I_FIN_M ( N_CHAINES_M ) = LONG_LIGNE
C-----------------------------------------------
C-----Longueur et contenu des chaînes distinctes
C-----------------------------------------------
      DO I = 1 , N_CHAINES_M
  	 LONG_CHAINE_M ( I )
     $ = I_FIN_M ( I ) - I_DEBUT_M ( I ) + 1
	 DO K = 0 , LONG_CHAINE_M ( I ) - 1
	   CHAINE_M ( I ) ( K + 1 : K + 1 )
     $   = LIGNE_M ( I_DEBUT_M ( I ) + K : I_DEBUT_M ( I ) + K )
	 END DO 
      END DO
      RETURN
      END
C---------------------------------------------------------
C-----Classement d'une liste d'entiers par ordre croissant
C---------------------------------------------------------
C      INCLUDE
C     $'tri_cr_entier.inc'
C     ==================================================================
C     I								       I
C     I Classement d'une liste d'entiers par ordre croissant	       I
C     I (procédure par récurrence) 				       I
C     I          						       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 23/02/2005                              I
C     I Saint Pierre Damien					       I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   TRI_CROISSANT_ENTIER
     $ ( N_VAL_M , VAL_ENTIER_M , IND_INIT_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C-------------------
C-----Liste initiale
C-------------------
      INTEGER * 4 VAL_ENTIER_M ( N_VAL_M )
C---------------------------------------------
C-----Indice dans la liste initiale
C-----des éléments classés par ordre croissant
C---------------------------------------------
      DIMENSION IND_INIT_M ( N_VAL_M )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C-------------------     
C-----Initialisation
C-------------------     
      DO I = 1 , N_VAL_M
	IND_INIT_M ( I ) = I
      END DO
C-----------------------------------------------------
C-----Balayage de la liste et insertion par récurrence
C-----------------------------------------------------
 1000     FORMAT
     $    ( 2X , 1000 ( 'a(', I2 , ' )   <   ' ) )  
 1100 FORMAT ('Liste initiale : ' , 2X , 1000 ( F7.2 , 6X ) )  
      DO I = 2 , N_VAL_M
C 	write ( * , * ) '-------'
C	write ( * , * ) 'I = ' , I
C	write ( * , * ) '-------'
C----------------------------------------------------------------
C-----Recherche du rang croissant de la nouvelle valeur à insérer
C-----dans la liste partielle dont le classement est déjà connu
C----------------------------------------------------------------
	J_0 = 0
	DO J = 1 , I - 2
	  J_1 = IND_INIT_M ( J )
	  J_2 = IND_INIT_M ( J + 1 )
	  IF ( VAL_ENTIER_M ( I ) .GT. VAL_ENTIER_M ( J_1 ) 
     $   .AND. VAL_ENTIER_M ( I ) .LE. VAL_ENTIER_M ( J_2 ) )
     $         J_0 = J
	END DO
C---------------------------------------------------------------
C-----Mise à jour des liens entre indices initiaux et croissants
C---------------------------------------------------------------
C--------------------------------------------------------------
C-----Premier cas :
C-----il existe des éléments de la liste partielle déjà classée
C-----encadrant le nouvel élément
C--------------------------------------------------------------
	IF ( J_0 .NE. 0 ) THEN
C	  write ( * , * )
C    $   'Cas 1 : "il existe j tel que a(j) < a(i) < a(j+1)"' 
C	  write ( * , * ) ' J_0 = ' , J_0
	  DO K = I - 1 , J_0 + 1 , - 1
C----------------------------------------------------------
C-----Remarque : il faut balayer K en décroissant
C-----pour éviter de perdre la première valeur au recopiage
C----------------------------------------------------------
	    IND_INIT_M ( K + 1 ) = IND_INIT_M ( K )
	  END DO
	  IND_INIT_M ( J_0 + 1 ) = I
        ELSE
C---------------------------------------------------
C-----Deuxième cas :
C-----le nouvel élément est supérieur
C-----à tous ceux de la liste partielle déjà classée
C---------------------------------------------------
	  I_1 = IND_INIT_M ( I - 1 )
	  IF ( VAL_ENTIER_M ( I ) .GT. VAL_ENTIER_M ( I_1 ) ) THEN
C         write ( * , * ) 
C    $   'Cas 2 : "quel que soit j, a(i) > a(j)"'
	   IND_INIT_M ( I ) = I 
	  ELSE
C---------------------------------------------------
C-----Troisième cas :
C-----le nouvel élément est inférieur
C-----à tous ceux de la liste partielle déjà classée
C---------------------------------------------------
C         write ( * , * ) 
C    $   'Cas 3 : "quel que soit j, a(i) < a(j)"'
           DO K = I , 2 , - 1
C----------------------------------------------------------
C-----Remarque : il faut balayer K en décroissant
C-----pour éviter de perdre la première valeur au recopiage
C----------------------------------------------------------
            IND_INIT_M ( K ) = IND_INIT_M ( K - 1 )
           END DO
	   IND_INIT_M ( 1 ) = I
	  END IF
	END IF
C	  write ( * , 1100 ) ( val_m ( k ) , k = 1 , n_val_m )
C	  write ( * , 1000 )
C     $ ( IND_INIT_M ( k ) , k = 1 , n_val_m )
      END DO	
      RETURN
      END
C-----------------------------------------------
C-----Déterminant et inverse d'une matrice 3 x 3
C-----------------------------------------------
C      INCLUDE
C     $'det_inv_mat_3_3.inc'
C     ==================================================================
C     I                                                                I
C     I Inversion d'une matrice 3 x 3 à l'aide de la formule :	       I
C     I Inverse de M 						       I
C     I = transposée (matrice des cofacteurs de M) / déterminant (M)   I
C     I                                                                I
C     I Calcul du déterminant de trois vecteurs			       I
C     I par développement première colonne de la matrice correspondanteI
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 12/05/2003                              I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   INV_MAT_3_3
     $ ( D_MAT_M , D_MAT_INV_M ) 
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      DIMENSION D_MAT_M ( 3 , 3 )
      DIMENSION D_MAT_INV_M ( 3 , 3 )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      DIMENSION D_MAT_AUX ( 3 , 3 )
C-------------------------------------------------------
C-----Calcul de la comatrice = matrice des cofacteurs
C-----(cofacteur ( i , j ) = (-1)i+j x mineur ( i , j ))
C-------------------------------------------------------
      D_MAT_AUX ( 1 , 1 ) = D_MAT_M ( 2 , 2 ) * D_MAT_M ( 3 , 3 )
     $			  - D_MAT_M ( 2 , 3 ) * D_MAT_M ( 3 , 2 )
      D_MAT_AUX ( 2 , 1 ) = D_MAT_M ( 1 , 2 ) * D_MAT_M ( 3 , 3 )
     $                    - D_MAT_M ( 1 , 3 ) * D_MAT_M ( 3 , 2 )
      D_MAT_AUX ( 2 , 1 ) = - D_MAT_AUX ( 2 , 1 )
      D_MAT_AUX ( 3 , 1 ) = D_MAT_M ( 1 , 2 ) * D_MAT_M ( 2 , 3 )
     $                    - D_MAT_M ( 1 , 3 ) * D_MAT_M ( 2 , 2 )
      D_MAT_AUX ( 1 , 2 ) = D_MAT_M ( 2 , 1 ) * D_MAT_M ( 3 , 3 )
     $                    - D_MAT_M ( 3 , 1 ) * D_MAT_M ( 2 , 3 )
      D_MAT_AUX ( 1 , 2 ) = - D_MAT_AUX ( 1 , 2 )
      D_MAT_AUX ( 2 , 2 ) = D_MAT_M ( 1 , 1 ) * D_MAT_M ( 3 , 3 )
     $                    - D_MAT_M ( 1 , 3 ) * D_MAT_M ( 3 , 1 )
      D_MAT_AUX ( 3 , 2 ) = D_MAT_M ( 1 , 1 ) * D_MAT_M ( 2 , 3 )
     $                    - D_MAT_M ( 2 , 1 ) * D_MAT_M ( 1 , 3 )
      D_MAT_AUX ( 3 , 2 ) = - D_MAT_AUX ( 3 , 2 )
      D_MAT_AUX ( 1 , 3 ) = D_MAT_M ( 2 , 1 ) * D_MAT_M ( 3 , 2 )
     $                    - D_MAT_M ( 3 , 1 ) * D_MAT_M ( 2 , 2 )
      D_MAT_AUX ( 2 , 3 ) = D_MAT_M ( 1 , 1 ) * D_MAT_M ( 3 , 2 )
     $                    - D_MAT_M ( 3 , 1 ) * D_MAT_M ( 1 , 2 )
      D_MAT_AUX ( 2 , 3 ) = - D_MAT_AUX ( 2 , 3 )
      D_MAT_AUX ( 3 , 3 ) = D_MAT_M ( 1 , 1 ) * D_MAT_M ( 2 , 2 )
     $                    - D_MAT_M ( 1 , 2 ) * D_MAT_M ( 2 , 1 )
C---------------------------------------
C-----Calcul du déterminant de D_MAT_M
C-----par développement première colonne
C---------------------------------------
      DET_1 = D_MAT_M ( 2 , 2 ) * D_MAT_M ( 3 , 3 )
     $      - D_MAT_M ( 2 , 3 ) * D_MAT_M ( 3 , 2 )
      DET_1 = D_MAT_M ( 1 , 1 ) * DET_1
      DET_2 = D_MAT_M ( 1 , 2 ) * D_MAT_M ( 3 , 3 )
     $      - D_MAT_M ( 3 , 2 ) * D_MAT_M ( 1 , 3 )
      DET_2 = - DET_2
      DET_2 = D_MAT_M ( 2 , 1 ) * DET_2
      DET_3 = D_MAT_M ( 1 , 2 ) * D_MAT_M ( 2 , 3 )
     $      - D_MAT_M ( 2 , 2 ) * D_MAT_M ( 1 , 3 )
      DET_3 = D_MAT_M ( 3 , 1 ) * DET_3
      DET = DET_1 + DET_2 + DET_3
      IF ( DET .EQ. 0.D0 ) THEN
	WRITE ( * , * ) "Attention : la matrice n'est pas inversible"
        CALL INTERRUPTION
      END IF
C------------------------
C-----Calcul de l'inverse
C------------------------
      DO I = 1 , 3
	DO J = 1 , 3
	  D_MAT_INV_M ( I , J ) = D_MAT_AUX ( J , I ) / DET
	END DO
      END DO
      RETURN
      END
C     ##################################################################
      SUBROUTINE DETERMINANT ( B_1_M , B_2_M , B_3_M , DET_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      DIMENSION B_1_M ( 3 ) , B_2_M ( 3 ) , B_3_M ( 3 )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C--------------------------------------------------------------
C-----Développement première colonne de la matrice des vecteurs
C--------------------------------------------------------------
      DET_1 = B_2_M ( 2 ) * B_3_M ( 3 )
     $      - B_2_M ( 3 ) * B_3_M ( 2 )
      DET_1 = B_1_M ( 1 ) * DET_1
      DET_2 = B_2_M ( 1 ) * B_3_M ( 3 )
     $      - B_2_M ( 3 ) * B_3_M ( 1 )
      DET_2 = - DET_2
      DET_2 = B_1_M ( 2 ) * DET_2
      DET_3 = B_2_M ( 1 ) * B_3_M ( 2 )
     $      - B_2_M ( 2 ) * B_3_M ( 1 )
      DET_3 = B_1_M ( 3 ) * DET_3
      DET_M = DET_1 + DET_2 + DET_3
      RETURN
      END
C     ==================================================================
C     I                                                                I
C     I Factorisation d'une matrice carrée réelle		       I
C     I en un produit LU d'une matrice L triangulaire inférieure       I
C     I à diagonale unité et d'une matrice U triangulaire supérieure   I
C     I                                                                I
C     I La matrice à décomposer étant M, le programme donne L et U     I
C     I telles que P * M * Q = L * U				       I
C     I (P et Q = matrices de permutations sur lignes et colonnes)     I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 24/01/2005                              I
C     I Saint Timothée						       I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   FACT_LU
     $ ( N_M ,
     $   MAT_M ,
     $   P_MAT_M , Q_MAT_M ,
     $   L_MAT_M , U_MAT_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 MAT_M ( N_M , N_M ) ,
     $         L_MAT_M ( N_M , N_M ) ,
     $         U_MAT_M ( N_M , N_M ) 
      REAL * 8 P_MAT_M ( N_M , N_M ) ,
     $	       Q_MAT_M ( N_M , N_M )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      REAL * 8 VECT ( N_M ) , L_1_MAT ( N_M , N_M ) 
C======================================================
C=====Matrices de permutations lignes (P), colonnes (Q)
C=====et matrice L initialisées à l'identité
C======================================================
       DO I = 1 , N_M
         DO J = 1 , N_M
                P_MAT_M ( I , J ) = 0.D0
                Q_MAT_M ( I , J ) = 0.D0
C--------------------------------------
C-----Matrice triangulaire inférieure L
C--------------------------------------
                L_MAT_M ( I , J ) = 0.D0
                L_1_MAT ( I , J ) = 0.D0
                IF( I .EQ. J ) THEN
                        P_MAT_M ( I , J ) = 1.D0
                        Q_MAT_M ( I , J ) = 1.D0
                        L_1_MAT ( I , J ) = 1.D0
                END IF
C--------------------------------------
C-----Matrice triangulaire supérieure U
C--------------------------------------
                U_MAT_M ( I , J ) = MAT_M ( I , J )
         END DO
       END DO
C=================================================
C-----Factorisation LU de la matrice MAT_M 
C-----de dimension N_M x N_M en N_M - 1 itérations
C=================================================
      DO IT = 1 , N_M - 1
C	write ( * , * ) 'it = ' , it
C-------------------------------------------------------------------
C-----Recherche du pivot global = Max (|U ( I , J )|) pour I,J >= IT
C-----de la sous-matrice (N_M - IT + 1 , N_M - IT + 1) de U
C-------------------------------------------------------------------
          U_PIV = 0.D0
          DO I = IT , N_M
            DO J = IT , N_M
               IF ( DABS ( U_MAT_M ( I , J ) ) .GT. U_PIV ) THEN
                  U_PIV = DABS ( U_MAT_M ( I , J ) )
                  I_PIV = I
                  J_PIV = J
               END IF
C		write ( * , * ) u_piv
             END DO
           END DO           
C         write ( * , * ) ' I_PIV , J_PIV = ' , I_PIV , J_PIV
C---------------------------------------------------------------------
C-----Permutation des lignes ( IT , I_PIV ) et colonnes ( IT , J_PIV )
C-----pour obtenir P ( IT , I_PIV ) * A * P ( IT , J_PIV )
C-----et ( P_MAT_M , Q_MAT_M )
C---------------------------------------------------------------------
C-----------------------------------------------------------------
C-----(i) permutation des lignes ( IT , I_PIV )
C-----(revient à une multiplication à gauche par P ( IT , I_PIV ))
C-----------------------------------------------------------------
           DO J = 1 , N_M
                VECT ( J ) = P_MAT_M ( IT , J ) 
                P_MAT_M ( IT , J ) = P_MAT_M ( I_PIV , J )
                P_MAT_M ( I_PIV , J ) = VECT ( J )
                VECT ( J ) = U_MAT_M ( IT , J )
                U_MAT_M ( IT , J ) = U_MAT_M ( I_PIV , J )
                U_MAT_M ( I_PIV , J ) = VECT ( J )
           END DO
C-----------------------------------------------------------------
C-----(ii) permutation des colonnes ( IT , J_PIV )
C-----(revient à une multiplication à droite par P ( IT , J_PIV ))
C-----------------------------------------------------------------
        DO I = 1 , N_M
                VECT ( I ) = Q_MAT_M ( I , IT )
                Q_MAT_M ( I , IT ) = Q_MAT_M ( I , J_PIV )
                Q_MAT_M ( I , J_PIV ) = VECT ( I )
                VECT ( I ) = U_MAT_M ( I , IT )
                U_MAT_M ( I , IT ) = U_MAT_M ( I , J_PIV )
                U_MAT_M ( I , J_PIV ) = VECT ( I )
        END DO
C	write ( * , * ) 'U provisoire après permutation'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( U_MAT_M ( I , J ) , J = 1 , N_M )
        END DO 
C------------------------------------------------------------------
C-----Permutation des COLONNES (pas des lignes) ( IT , I_PIV ) de L
C------------------------------------------------------------------
         DO I = 1 , N_M
                VECT ( I ) = L_1_MAT ( I , IT )
                L_1_MAT ( I , IT ) = L_1_MAT ( I , I_PIV )
                L_1_MAT ( I , I_PIV ) = VECT ( I )
        END DO
C       write ( * , * ) 'L provisoire après permutation'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( L_1_MAT ( I , J ) , J = 1 , N_M )
        END DO 
C--------------------------------------------------
C-----Opérations sur L
C-----(multiplication à droite par M-1)
C-----correspondant à la triangularisation de U
C-----(celle-ci est effectuée seulement ensuite,
C-----car elle change les éléments de matrice  de U
C-----utiles pour les opérations sur L)
C--------------------------------------------------
        DO I = 1 , N_M
          DO J = IT + 1 , N_M
            L_1_MAT ( I , IT ) = L_1_MAT ( I , IT )
     $                         + L_1_MAT ( I , J )
     $                         / U_MAT_M ( IT , IT )
     $                         * U_MAT_M ( J , IT )
          END DO
        END DO
C       write ( * , * ) 'L provisoire après opération'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( L_1_MAT ( I , J ) , J = 1 , N_M )
        END DO   
C-------------------------------------------------
C-----Triangularisation de U (multiplication par M
C-----de P ( IT , I_PIV ) * A * P ( IT , J_PIV ) )
C-----Le pivot est en position ( IT , IT )
C-------------------------------------------------
         DO I = IT + 1 , N_M
C-------------------------------------------------------
C-----Cet élément de matrice de la ligne I change
C-----lors de la combinaison linéaire des lignes I et IT
C-----=> on le mémorise
C-------------------------------------------------------
            UU = U_MAT_M ( I , IT )
C		write ( * , * ) ' UU = ' , UU
            DO J = 1 , N_M
              U_MAT_M ( I , J )
     $      = U_MAT_M ( I , J )
     $      - UU / U_MAT_M ( IT , IT ) * U_MAT_M ( IT , J )
            END DO
        END DO
C        write ( * , * ) 'U provisoire après triangulation'
        DO I = 1 , N_M
C         WRITE ( * , * ) ( U_MAT_M ( I , J ) , J = 1 , N_M )
        END DO
C-----------------------------------------
C---Fin de la boucle de N_M - 1 itérations
C-----------------------------------------
        END DO
C       write ( * , * ) 'P_MAT_M pour multiplication par L'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( P_MAT_M ( I , J ) , J = 1 , N_M )
        END DO     
C------------------------------------------------------------
C-----On termine le calcul de L en multipliant à gauche par P
C------------------------------------------------------------
        DO I = 1 , N_M
          DO J = 1 , N_M
            DO K = 1 , N_M
                L_MAT_M ( I , J )
     $        = L_MAT_M ( I , J )
     $        + P_MAT_M ( I , K ) * L_1_MAT ( K , J )
            END DO
          END DO
        END DO
        RETURN
        END
C      INCLUDE
C     $ 'inverse_lu.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul de l'inverse d'une matrice préalablement décomposée     I
C     I en produit L*U avec matrices de permutations P et Q	       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 24/01/2005                              I
C     I Saint Timothée						       I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   INVERSE_LU
     $ ( N_M ,
     $   L_MAT_M , U_MAT_M ,
     $   P_MAT_M , Q_MAT_M ,
     $   MAT_INV_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 L_MAT_M ( N_M , N_M ) , U_MAT_M ( N_M , N_M ) ,
     $         P_MAT_M ( N_M , N_M ) , Q_MAT_M ( N_M , N_M ) ,
     $         MAT_INV_M ( N_M , N_M )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      REAL * 8 C_MAT ( N_M , N_M )
      REAL * 8 MAT_INV_LU ( N_M , N_M )
      REAL * 8 B_MAT ( N_M , N_M )
C-----------------------------------------
C-----Matrice inverse de L * U = P * M * Q
C-----------------------------------------
       DO I = 1 , N_M
         DO J = 1 , N_M
		C_MAT ( I , J ) = 0.D0
         END DO
       END DO 
C=======================================================
C-----Calcul de MAT_INV_LU inverse de L*U en deux étapes
C-----(L * U * MAT_INV_LU = I)
C=======================================================
C---------------------------------------
C-----(i) Résolution de L * C_MAT = I
C-----(où C_MAT = U * MAT_INV_LU)
C-----avec C_MAT triangulaire inférieure
C-----à diagonale unité
C---------------------------------------
	DO J = 1 , N_M
	  C_MAT ( J , J ) = 1.D0
	  DO I = J + 1 , N_M
	   DO K = 1 , I - 1
		C_MAT ( I , J )
     $        = C_MAT ( I , J )
     $        - L_MAT_M ( I , K ) * C_MAT ( K , J )
	   END DO
          END DO
	END DO	
C       write ( * , * ) 'N = ', N_M
C       write ( * , * ) 'Matrice C_MAT'
        DO I = 1 , N_M
C               WRITE ( * , * )
C    $        ( C_MAT ( I , J ) , J = 1 , N_M )
        END DO
C       write ( * , * ) 'Matrice U'
        DO I = 1 , N_M
C               WRITE ( * , * )
C    $        ( U_MAT_M ( I , J ) , J = 1 , N_M )
        END DO          	
C----------------------------------------------
C-----(ii) Résolution de U * MAT_INV_LU = C_MAT
C----- avec ikl = 0 si k > l
C----------------------------------------------
        DO J = 1 , N_M
         MAT_INV_LU ( N_M , J )
     $ = C_MAT ( N_M , J ) / U_MAT_M ( N_M , N_M )
	  DO I = N_M - 1 , 1 , - 1
C----------------------------------------------------------
C-----ATTENTION : LA SYNTAXE "STEP -1" PROVOQUE DES ERREURS
C----------------------------------------------------------
	   MAT_INV_LU ( I , J ) = C_MAT ( I , J )
           DO K = N_M , I + 1 , - 1
                MAT_INV_LU ( I , J )
     $        = MAT_INV_LU ( I , J )
     $        - U_MAT_M ( I , K ) * MAT_INV_LU ( K , J )
           END DO
C         write ( * , * ) 'MAT_INV_LU ( I , J )' ,
C    $                     I , J , MAT_INV_LU ( I , J )
	  MAT_INV_LU ( I , J ) = MAT_INV_LU ( I , J )
     $                           / U_MAT_M ( I , I )
         END DO
       END DO       
C      write ( * , * )
C    $ 'matrice MAT_INV_LU (telle que U * MAT_INV_LU = C_MAT)'
       DO I = 1 , N_M
C           WRITE ( * , * )
C    $    ( MAT_INV_LU ( I , J ) , J = 1 , N_M )
        END DO            
C=================================================
C=====Calcul de l'inverse de la matrice initiale
C=====à l'aide des matrices de permutations P et Q
C=====MAT_INV = Q_MAT*MAT_INV_LU*P_MAT
C=================================================
       DO I = 1 , N_M
          DO J = 1 , N_M
            B_MAT ( I , J ) = 0.D0
            DO K = 1 , N_M
                 B_MAT ( I , J )
     $         = B_MAT ( I , J )
     $         + MAT_INV_LU ( I , K ) * P_MAT_M ( K , J )
            END DO
          END DO
        END DO
        DO I = 1 , N_M
          DO J = 1 , N_M
            MAT_INV_M ( I , J ) = 0.D0
            DO K = 1 , N_M
               MAT_INV_M ( I , J )
     $       = MAT_INV_M ( I , J )
     $       + Q_MAT_M ( I , K ) * B_MAT ( K , J )
            END DO
          END DO
        END DO
	RETURN
	END
C      INCLUDE
C     $ 'prod_mat_vect.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul du produit M.V d'un vecteur par une matrice 	       I
C     I								       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 06/07/2005                              I
C     I Sainte Maria Goretti                                           I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   PROD_MAT_VECT
     $ ( N_M , MAT_M , VECT_M , VECT_PROD_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 MAT_M ( N_M , N_M ) ,
     $         VECT_M ( N_M ) ,
     $         VECT_PROD_M ( N_M )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      VECT_PROD_M = 0.D0
      DO I_C = 1 , N_M
	DO J_C = 1 , N_M
	  VECT_PROD_M ( I_C ) 
     $  = VECT_PROD_M ( I_C )
     $  + MAT_M ( I_C , J_C ) * VECT_M ( J_C )
        END DO
      END DO
      RETURN
      END
C     ==================================================================
C     I                                                                I
C     I Recherche du facteur lambda (entre 0 et 1) correctif          I
C     I pour le pas de Newton-Raphson			 	       I
C     I (méthode de convergence globale) tel que		       I
C     I x_1(nouveau) = x_0(ancien) + lambda * ( x_1_NR - x_0 )         I
C     I avec x_1_NR - x_0 = pas_NR (pas NR complet)		       I
C     I								       I
C     I Procédure :						       I
C     I à partir du pas NR complet (lambda = 1, x_1 = x_0 + pas_NR)    I
C     I lambda est réduit progressivement			       I
C     I (x_1 = x_0 + lambda * pas_NR)				       I
C     I jusqu'à ce que soit vérifié le critère			       I
C     I phi(x_1) <= phi(x_0) + alpha grad(phi).( x_1 - x_0 )	       I
C     I	où phi(x) fonction scalaire de x est définie par :	       I
C     I	phi(x) = (1/2) F(x).F(x) = (1/2) somme_i (F_i)^2               I
C     I (avec F(x) la fonction vectorielle dont on cherche une racine) I
C     I	=> grad(phi)_j = somme_i F_i(dF_i/dx_j) = somme_i F_i.J_ij     I
C     I (J = matrice des dérivées partielles)				       I
C     I 							       I
C     I Le pas initial (réduction éventuelle de lambda = 1 )           I
C     I est traité différemment des autres			       I
C     I (approximation quadratique au lieu de cubique)		       I
C     I 							       I
C     I lambda est forcé à être entre deux valeurs minimale et maximaleI
C     I 							       I
C     I La procédure renvoie le nouveau vecteur x_1 accepté            I
C     I (et la valeur de phi en ce point)			       I
C     I 							       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 26/07/2006                              I
C     I Sainte Anne                                                    I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   CALC_PAS_NRCG 
C----------------------------------------------------------
C----Paramètres spécifiques ADPI de la fonction vectorielle
C----------------------------------------------------------
     $ ( N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C----------------------------------------------------
C-----Paramètres déterminés en dehors de la procédure
C----------------------------------------------------
     $   N_M ,
     $   X_0_M , PHI_0_M , GRAD_PHI_0_M ,
     $   ALPHA_NRCG_NPT_M ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT_M ,
     $   INDIC_TYPE_REDUC_NRCG_NPT_M , COEF_REDUC_NRCG_NPT_M ,
C--------------------------------------------------------------------
C-----Paramètre de contrôle de l'écriture à l'écran par CALC_PAS_NRCG
C--------------------------------------------------------------------
     $   I_ECRIT_CALC_PAS_NRCG_M ,
C-----------------------------------------------
C-----Paramètres déterminés dans cette procédure
C-----------------------------------------------
     $   PAS_NR_M , I_FONCTION_NON_DEFINIE_M ,
     $   X_1_M , PHI_1_M ) 
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C      PARAMETER ( TOLERANCE_X = 1.D-7 )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      REAL * 8 X_0_M ( N_M ) 
      REAL * 8 X_1_M ( N_M ) 
      REAL * 8 GRAD_PHI_0_M ( N_M ) 
      REAL * 8 PAS_NR_M ( N_M ) 
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-réseau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-réseau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-réseau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantité de matière totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de référence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C---------------------------------------
C-----Facteurs de réduction du pas de NR
C---------------------------------------
      REAL * 8 LAMBDA
      REAL * 8 LAMBDA_MIN 
      REAL * 8 LAMBDA_2
      REAL * 8 LAMBDA_TEMP
C--------------------------------------
C-----Valeur de la fonction vectorielle
C-----au point x_1 estimé
C--------------------------------------
      REAL * 8 F_1 ( N_M )
C-----------------------------------------------------
C-----Indicateur de signe négatif pour chaque variable
C-----------------------------------------------------
      INTEGER * 4 INDIC_COMP_NEG ( N_M )
C-----Matrice des dérivées partielles :
C-----inutile ici, mais présente en argument de F_NR_AVEC_DERIV
      REAL * 8 MAT_DERIV ( N_M )
C======================================
C=====Calcul de grad(phi).( x_1 - x_0 )
C======================================
C$$$$$ ATTENTION : vérifier si x_1 ici = x_1_NR complet
C$$$$$ ou si x_1 doit être réduit au fur et à mesure des pas NRCG
C$$$$$ (c'est le x_1 apparaissant dans le critère
C$$$$$ phi(x_1) <= phi(x_0) + alpha grad(phi).( x_1 - x_0 ))
      GRAD_PHI_PAS_NR = 0.D0
      DO I = 1 , N_M
	  GRAD_PHI_PAS_NR = GRAD_PHI_PAS_NR 
     $  	  	      + GRAD_PHI_0_M ( I ) * PAS_NR_M ( I )
      END DO	
C================================================================
C=====Calcul de lambda minimum autorisé.
C=====Il s'agit d'une borne indiquant une convergence "factice" :
C=====dans ce cas, la procédure s'arrête et transmet x_1 à sys_NR
C================================================================
C      X_TEST = 0.D0
C      DO I = 1 , N_M
C        YY = X_0_M ( I )
C        IF ( I .LT. N_M ) THEN
C         YY = DEXP ( DLOG ( 10.D0 ) * X_0_M ( I ) )
C        END IF
C        X_TEMP = DABS ( PAS_NR_M ( I ) )
C     $         / MAX ( DABS ( YY ) , 1.D0 )
C        IF ( X_TEMP .GT. X_TEST ) X_TEST = X_TEMP
C      END DO
C      LAMBDA_MIN = TOLERANCE_X / X_TEST
C-----Après divers essais, on utilise finalement pour lambda_min
C-----une valeur fixe (faible) spécifiée dans DATA.adpi.
      LAMBDA_MIN = VALEUR_LAMBDA_MIN_NRCG_NPT_M 
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
      write ( * , * ) 'Dans NRCG : LAMBDA_MIN = ' , LAMBDA_MIN	
      write ( * , * ) 'Dans NRCG : ALPHA = ' , ALPHA_NRCG_NPT_M 
      end if
C----------------------------------------------
C-----Initialisation de lambda (pas NR complet)
C----------------------------------------------
      LAMBDA = 1.D0
C##############################################
C#####Début de la boucle de réduction de lambda
C##############################################
      I_PAS_CG = 0
      DO
	I_PAS_CG = I_PAS_CG + 1
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
C       write ( * , 1 )
 	 write ( * , * ) '   =========='
	 write ( * , * ) '   I_PAS_CG = ' , I_PAS_CG
	 write ( * , * ) '   =========='
         write ( * , * ) 'PHI_2 = ' , PHI_2
C       write ( * , 1 )
        end if
C============================================================
C=====Calcul du x_1 estimé pour le lambda courant
C=====(l'initialisation de lambda à 1 fait que
C=====le pas NR complet est toujours essayé en premier lieu).
C============================================================
C-----Remarque : en cas de X_1(I)<0, la démarche suit deux étapes :
C-----1) la variable PAS_NR, initialement le pas NR complet,
C-----est modifiée, composante par composante, pour remédier à X_1(I)<0.
C-----2) une fois tous les X_1(I)>0, PAS_NR n'évolue plus,
C-----et c'est la réduction de lambda qui est effectuée.
C-----Rappel : en fonctionnement "normal",
C-----la variable I_FONCTION_NON_DEFINIE vaut 0.
C-----Elle peut passer à 1 dans deux cas :
C-----1) ci-dessous, lorsque l'une ou l'autre des composantes
C-----du x_1 estimé devient <0 ;
C-----2) en début de F_NR_AVEC_DERIV, si une quantité z est trouvée <0.
C-----Dans les deux cas, la valeur I_FONCTION_NON_DEFINIE = 1 est requise
C-----pour que f_NR effectue l'estimation de la fonction.
C-----Il y a deux appels à f_NR :
C-----* le premier dans sys_NR, pour évaluer F(x_0) en x_0 point courant ;
C-----* le second dans calc_pas_NRCG, pour évaluer F(x_1) (x_1 estimé). 
      I_FONCTION_NON_DEFINIE_M = 0
      INDIC_COMP_NEG = 0
      INDIC_EXISTENCE_COMP_NEG = 0 
	DO I = 1 , N_M - 1
C-----Rappel : excepté la dernière (nombre de mailles), les variables 
C-----x_NPT manipulées par les procédures sont les log_10 des x_DP.
C-----Le passage à xDP = 10^x_NPT est fait en début de F_NR.
	  X_1_M ( I ) = DEXP ( DLOG ( 10.D0 ) * X_0_M ( I ) )
     $              + LAMBDA * PAS_NR_M ( I )
        IF ( X_1_M ( I ) .LE. 0.D0 ) THEN
         I_FONCTION_NON_DEFINIE_M = 1
         INDIC_EXISTENCE_COMP_NEG = 1
         INDIC_COMP_NEG ( I ) = 1
         if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
          write(*,*)
     $   'Pas NR -> valeur négative pour variable' , I 
          end if
        END IF
        X_1_M ( I ) = DLOG ( X_1_M ( I ) ) / DLOG ( 10.D0 )
C-----Fin de la boucle sur les composantes "log_10"
	END DO
C-----La dernière composante est spécifiée directement (et non son log).
      X_1_M ( N_M ) = X_0_M ( N_M ) + LAMBDA * PAS_NR_M ( N_M )
      IF ( X_1_M ( N_M ) .LE. 0.D0 ) THEN
         I_FONCTION_NON_DEFINIE_M = 1
         INDIC_EXISTENCE_COMP_NEG = 1
         INDIC_COMP_NEG ( N_M ) = 1
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
          write(*,*)
     $   'Pas NR -> valeur négative pour variable' , N_M
       end if
      END IF
C----------------------------------------------------------------
C-----Réduction préliminaire du pas pour "élimination de X_1<0" :
C-----* option 1 : réduction suivant les seules comp. X_1(I)<0 
C-----* option 2 : réduction de toutes les composantes du pas
C----------------------------------------------------------------
      IF ( INDIC_EXISTENCE_COMP_NEG .EQ. 1 ) THEN
       IF ( INDIC_TYPE_REDUC_NRCG_NPT_M .EQ. 1 ) THEN
        DO I = 1 , N_M
         IF ( INDIC_COMP_NEG ( I ) .EQ. 1 ) THEN
          PAS_NR_M ( I ) = COEF_REDUC_NRCG_NPT_M * PAS_NR_M ( I ) 
          if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
          write(*,*)
     $   'composante ' , I , ' : Pas_NR -> Pas_NR * ' ,
     $    COEF_REDUC_NRCG_NPT_M 
          end if
         END IF
        END DO
       ELSE IF ( INDIC_TYPE_REDUC_NRCG_NPT_M .EQ. 2 ) THEN
          PAS_NR_M = COEF_REDUC_NRCG_NPT_M * PAS_NR_M
         if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then 
          write(*,*) ' : Pas_NR -> Pas_NR * ' ,
     $    COEF_REDUC_NRCG_NPT_M , ' pour toutes les composantes' 
         end if
       ELSE
        if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
         write(*,*) 'INDIC_TYPE_REDUC_NRCG_NPT_M mal choisi'
        end if
        stop
       END IF
      END IF
C       write ( * , * ) '------------------'
C       write ( * , * ) 'Dans calc_pas_NRCG :'
C       write ( * , * ) '------------------'
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
       write ( * , * ) 'LAMBDA = ' , LAMBDA
       write ( * , * ) 'PAS_NR_M = ' , PAS_NR_M
       write ( * , * ) 'X_0 (log_10 sauf dernière comp.) = ' , X_0_M
       write ( * , * ) 'X_1 (log_10 sauf dernière comp.) = ' , X_1_M
       end if
C-----Remarque : la composante ci-dessus vaut "NaN" si x_1<0,
C-----mais ce n'est pas gênant, car ce "NaN" ne sera pas utilisé
C-----grâce à l'information I_FONCTION_NON_DEFINIE = 1.
C=====================================================================
C=====Calcul de phi(x_1) = (1/2) F(x_1).F(x_1) = (1/2) somme_i (F_i)^2
C=====avec x_1 déduit de la valeur courante de lambda
C=====================================================================
C-----Ici, le calcul des dérivées partielles est inutile.
       INDIC_CALC_DERIV = 0
         CALL
     $   F_NR_AVEC_DERIV
     $ (
C---------------------------------------
C-----Paramètres dont dépend la fonction
C---------------------------------------
     $   N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C--------------------------------
C-----Paramètres généraux de f_NR
C--------------------------------
     $   I_FONCTION_NON_DEFINIE_M ,
     $   N_M , X_1_M , F_1 ,
     $   INDIC_CALC_DERIV , MAT_DERIV )
        PHI_1_M = 0.D0 
        DO I = 1 , N_M
          PHI_1_M 
     $  = PHI_1_M + F_1 ( I ) * F_1 ( I )
        END DO
	  PHI_1_M = 0.5D0 * PHI_1_M
C==============================
C==============================
C=====Test sur la valeur du pas
C==============================
C==============================
C=====================================================================
C=====Remarque : par rapport à la procédure "classique",
C=====un cas supplémentaire (traité en réduisant lambda) a été ajouté,
C=====qui correspond à la non-définition de la fonction
C=====(ce qui peut se produire car celle-ci contient des logarithmes)
C=====================================================================
C===============
C=====Cas "zéro"
C===============
       IF ( LAMBDA .LT. LAMBDA_MIN ) THEN
C $$$$$ Instruction ci-dessous à mettre en option ?
C 	  X_1_M = X_0_M
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
	  write ( * , * )
     $ '-----------------------------------------------------------'
	  write ( * , * )
     $ 'Convergence sur lambda "factice" = seuil lambda_min atteint'
C	  write ( * , * )
C     $ '    -> NRCG s'arrête et transmet la valeur courante de x_1.'
	  write ( * , * )
     $ '-----------------------------------------------------------'
        end if
 	RETURN
C================
C=====Premier cas : la fonction est définie et le pas est convenable
C=====(i.e. le lambda courant vérifie le critère de décroissance
C=====de phi(x) = (1/2) F(x).F(x) = (1/2) somme_i (F_i)^2
C===== => la procédure calc_pas_NRCG s'arrête et transmet à sys_NR
C=====le point X_1 qui devient le nouveau point courant.
C================
       ELSE IF ( I_FONCTION_NON_DEFINIE_M .EQ. 0 
     $  .AND. PHI_1_M .LE.
     $ ( PHI_0_M + ALPHA_NRCG_NPT_M * LAMBDA * GRAD_PHI_PAS_NR ) )
     $ THEN 
         if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
         write ( * , * ) '-------------------------------'
         write ( * , * ) 'Convergence sur lambda atteinte'
         write ( * , * ) '-------------------------------'
         end if
	  RETURN
C=================
C=====Deuxième cas : la fonction est définie et le pas doit être réduit
C=================	
       ELSE IF ( I_FONCTION_NON_DEFINIE_M .EQ. 0 ) THEN
C-----------------------------------------------------
C-----------------------------------------------------
C-----Première réduction (depuis lambda = 1 forcément)
C-----------------------------------------------------
C-----------------------------------------------------
	  IF ( LAMBDA .EQ. 1.D0 ) THEN
	    LAMBDA_TEMP
     $  = - GRAD_PHI_PAS_NR
     $	/ ( 2.D0 * ( PHI_1_M - PHI_0_M - GRAD_PHI_PAS_NR ) )
C------------------------------------
C------------------------------------
C-----Deuxième réduction et suivantes
C------------------------------------
C------------------------------------
	  ELSE
	   TERME_1 = PHI_1_M - PHI_0_M - LAMBDA * GRAD_PHI_PAS_NR
         TERME_2 = PHI_2 - PHI_0_M - LAMBDA_2 * GRAD_PHI_PAS_NR 
	   A = ( TERME_1 / LAMBDA * * 2 - TERME_2 / LAMBDA_2 * * 2 )
     $       / ( LAMBDA - LAMBDA_2 )
           B = ( - LAMBDA_2 * TERME_1 / LAMBDA * * 2
     $           + LAMBDA * TERME_2 / LAMBDA_2 * * 2 )
     $       / ( LAMBDA - LAMBDA_2 )
C---------------------------------------
C-----Calcul de la valeur de lambda
C-----(le cas A = 0 doit être distingué)
C---------------------------------------
          IF ( A .EQ. 0.D0 ) THEN
  	      LAMBDA_TEMP = - GRAD_PHI_PAS_NR / ( 2.D0 * B )
	      I_CAS = 1
	    ELSE
            DISCRIMINANT = B * B - 3.D0 * A * GRAD_PHI_PAS_NR	    
	      IF ( DISCRIMINANT .LT. 0.D0 ) THEN
	         LAMBDA_TEMP = 0.5D0 * LAMBDA	
	         I_CAS = 2
	      ELSE IF ( B .LE. 0.D0 ) THEN
	       	LAMBDA_TEMP = ( - B + DSQRT ( DISCRIMINANT ) )
     $		            / ( 3.D0 * A )	
	            I_CAS = 3
	      ELSE
		      LAMBDA_TEMP = - GRAD_PHI_PAS_NR
     $                      / ( B + DSQRT ( DISCRIMINANT ) )
	            I_CAS = 4
	      END IF
	  END IF
C------------------------------
C-----Valeur maximale de lambda
C------------------------------
          IF ( LAMBDA_TEMP .GT. 0.5D0 * LAMBDA )
     $         LAMBDA_TEMP = 0.5D0 * LAMBDA
C-------------------------------------------------
C-------------------------------------------------
C-----Fin du test sur réductions 1 ou (2 et suiv.)
C-------------------------------------------------
C-------------------------------------------------
	  END IF
C-------------------------------------------
C-----Mise à jour des variables incrémentées
C-------------------------------------------
       LAMBDA_2 = LAMBDA
       PHI_2 = PHI_1_M
C------------------------------
C-----Valeur minimale de lambda
C------------------------------
       LAMBDA = MAX ( LAMBDA_TEMP , 0.1D0 * LAMBDA )
C------------------------------
C-----Ecritures de vérification
C------------------------------
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
C       write ( * , * ) '******************'
C       write ( * , * ) 'Dans calc_pas_NRCG :'
C       write ( * , * ) '******************'
       write ( * , * ) 'I_CAS = ' , I_CAS
       write ( * , * ) 'GRAD_PHI_PAS_NR = ' , GRAD_PHI_PAS_NR
       write ( * , * ) 'PHI_0 = ' , PHI_0_M 
       write ( * , * ) 'PHI_1 = ' , PHI_1_M 
       write ( * , * ) 'TERME_1 = ' , TERME_1
       write ( * , * ) 'TERME_2 = ' , TERME_2
       write ( * , * ) 'A = ' , A
       write ( * , * ) 'B = ' , B
C       write ( * , * ) 'LAMBDA_TEMP = ' , LAMBDA_TEMP
C       write ( * , * ) '0.1D0 * LAMBDA =' , 0.1D0 * LAMBDA
C       write ( * , * ) 'Max des 2 = ' ,
C    $           MAX ( LAMBDA_TEMP , 0.1D0 * LAMBDA )
       write ( * , * ) 'LAMBDA_2 = ' , LAMBDA_2
       write ( * , * ) 'LAMBDA = ' , LAMBDA
       write ( * , * ) 'DISCRIMINANT = ' , DISCRIMINANT
       end if
C==================
C=====Troisième cas : la fonction n'est pas définie
C==================
      ELSE
C----------------------------------------
C-----Le pas est réduit dans ce cas aussi
C----------------------------------------
C        PAS_NR_M = PAS_NR_M * 0.5D0
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
C       write ( * , * ) '****************************'
	write ( * , * ) 'Dans calc_pas_NRCG :'
	write ( * , * ) ' Fonction non définie en X_1' 
C       write ( * , * ) '****************************'
C	write ( * , * ) '=> nouveau pas (la moitié du précédent) :'
C       write ( * , * ) PAS_NR_M
      write ( * , * )
     $ "=> aucune action pour ce pas NRCG (sauf 0,5*PAS_NR(I))"
C       write ( * , * ) '****************************'
       end if
C=====================================
C=====================================
C=====Fin du test sur la valeur du pas
C=====================================
C=====================================
      END IF
C############################################
C#####Fin de la boucle de réduction de lambda
C############################################
      END DO
      RETURN
      END
C     ##################################################################
C     ==================================================================
C     I								       I
C     I Déclaration de la procédure de Newton-Raphson      	       I
C     I (incluant les paramètres spécifiques à la fonction)	       I
C     I          						       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 26/07/2006                              I
C     I Sainte Anne                                                    I
C     ==================================================================
         SUBROUTINE
     $   SYS_NR
     $ (
C---------------------------------------
C-----Paramètres dont dépend la fonction
C---------------------------------------
     $   N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,  
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C----------------------------------
C-----Paramètres généraux de sys_NR
C----------------------------------
     $   N_NPT_M , X_NPT_M , F_NPT_M , J_NPT_M , P_J_NPT_M ,
     $   ALPHA_NRCG_NPT_M ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT_M ,
     $   INDIC_TYPE_REDUC_NRCG_NPT_M , COEF_REDUC_NRCG_NPT_M ,
     $   N_ITER_MAX_NPT_M , PRECISION_NPT_M )
C     ##################################################################	  
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-réseau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-réseau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-réseau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantité de matière totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de référence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C------------------------------------------------
C-----Point courant (vecteur), et en ce point :
C-----fonction et matrice des dérivées partielles
C------------------------------------------------
      REAL * 8 X_NPT_M ( N_NPT_M )
      REAL * 8 F_NPT_M ( N_NPT_M )
      REAL * 8 J_NPT_M ( N_NPT_M , N_NPT_M )
C---------------------------------------------------------
C------Fréquence de mise à jour de la matrice des dérivées
C---------------------------------------------------------
      INTEGER * 4 P_J_NPT_M
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C---------------------------------------------
C-----Matrices relatives à la décomposition LU
C---------------------------------------------
      REAL * 8 L_MAT ( N_NPT_M , N_NPT_M ) ,
     $         U_MAT ( N_NPT_M , N_NPT_M ) ,
     $         P_MAT ( N_NPT_M , N_NPT_M ) ,
     $         Q_MAT ( N_NPT_M , N_NPT_M )
C---------------------------------------
C-----Inverse de la matrice des dérivées
C---------------------------------------
      REAL * 8 J_INV ( N_NPT_M , N_NPT_M )
C----------------------
C-----Pas complet de NR
C----------------------
      REAL * 8 PAS_NR ( N_NPT_M )
C------------------
C-----Point courant
C------------------
      REAL * 8 X_0 ( N_NPT_M )
C--------------------------------------------------
C-----Gradient de la fonction auxiliaire phi en x_0
C--------------------------------------------------
      REAL * 8 GRAD_PHI_0 ( N_NPT_M )
C############
C#####Formats
C############
 1000 FORMAT ( 100 ( 2X , G18.13 ) )
 1100 FORMAT ( 100 ( 2X , G18.8 ) )
 1200 FORMAT ( 100 ( '-' ) )
C     write ( * , * ) '------------'
C     write ( * , * ) 'Début sys_NR'
C     write ( * , * ) '------------'
C-----Indicateur (interne à cette procédure) d'écriture de détails
      INDIC_ECRIT_DETAILS_SYS_NR = 0
C################################################################
C#####Algorithme de Newton-Raphson sur N_ITER_MAX_NPT itérations
C#####ou arrêt à convergence
C#####(critère sur l'écart entre deux points courants successifs)
C################################################################
      N_ITER = 0
      I_CONV = 0
      I_FONCTION_NON_DEFINIE = 0
      DO WHILE ( ( N_ITER .LE. N_ITER_MAX_NPT_M )
     $     .AND. ( I_CONV .EQ. 0 ) ) 
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
C      write ( * , 2 ) 
C      write ( * ,* ) '================'	
      write ( * ,* ) '                    *********************'	
      write ( * ,* ) '                    Itération (NR) numéro ' ,
     $                                    N_ITER
      write ( * ,* ) '                    *********************'	
C      write ( * ,* ) '================'	
C      write ( * , 2 ) 
C      write ( * , 1200 )
C      write ( * , * ) 'Point courant (log_10 sauf dernière comp.)'
C      write ( * , 1200 )
C      write ( * , 1000 ) X_NPT_M 
C      write ( * , 1200 )
        end if
       N_ITER = N_ITER + 1
C=============================================================
C=====Pour le point courant,
C=====calcul de la fonction (spécifiée analytiquement)
C=====ainsi que de dFi/dxj (analytiquement) et de son inverse,
C=====tous les P_J_NPT pas (et au premier pas)
C=============================================================
       INDIC_CALC_DERIV = 0
       IF ( MOD ( N_ITER , P_J_NPT_M ) .EQ. 0 .OR. N_ITER .EQ. 1 ) THEN
          INDIC_CALC_DERIV = 1
       END IF
C-----Remarque importante : f_NR renvoie I_FONCTION_NON_DEFINIE = 1 :
C-----1) si certaines variables x_1 sont négatives ;
C-----2) si certains termes complémentaires d'entropie z sont négatifs
C-----(se produit si le pas NR induit des variables x trop grandes).
C-----Dans ce 2ème cas, après avoir calculé les z et détecté z<0,
C-----f_NR s'arrête sans calculer la fonction F,
C-----et "passe la main" à calc_pas_NRCG pour modification du pas.
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
      write(*,*) "Avant F_NR, I_FONCTION_NON_DEFINIE = " ,
     $  I_FONCTION_NON_DEFINIE 
        end if
C	stop
         CALL
     $   F_NR_AVEC_DERIV
     $ (
C---------------------------------------
C-----Paramètres dont dépend la fonction
C---------------------------------------
     $   N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C--------------------------------
C-----Paramètres généraux de f_NR
C--------------------------------
     $   I_FONCTION_NON_DEFINIE ,
     $   N_NPT_M , X_NPT_M , F_NPT_M ,
     $   INDIC_CALC_DERIV , J_NPT_M )
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
        write ( * , * )
     $ ' X_NPT_M (log_10 sauf dernière comp.) = ' , X_NPT_M
        write ( * , * ) ' F_NPT_M = ' , F_NPT_M
        end if
C----------------------------------------------
C-----Tous les P_J_NPT pas (et au premier pas),
C-----factorisation LU et inversion
C-----de la matrice des dérivées partielles
C----------------------------------------------
      IF ( MOD ( N_ITER , P_J_NPT_M ) .EQ. 0 .OR. N_ITER .EQ. 1 ) THEN
C      write(*,*) "Début FACT_LU"
         CALL
     $   FACT_LU
     $ ( N_NPT_M ,
     $   J_NPT_M ,
     $   P_MAT , Q_MAT ,
     $   L_MAT , U_MAT )
C      write(*,*) "Fin FACT_LU"
         CALL
     $   INVERSE_LU
     $ ( N_NPT_M ,
     $   L_MAT , U_MAT ,
     $   P_MAT , Q_MAT ,
     $   J_INV )
      END IF
C-------------
C-----Ecriture
C-------------
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
C        write ( * , 1200 )
C        write ( * , * ) 'Fonction'
C        write ( * , 1200 )
C        write ( * , * ) F_NPT_M
C        write ( * , 1200 )
         write ( * , 1200 )
        write ( * , * ) 'Matrice des dérivées partielles'
        write ( * , 1200 )
        do i_c = 1 , N_NPT_M
          write ( * , 1100 )
     $  ( J_NPT_M ( i_c , j_c ) , j_c = 1 , N_NPT_M )
        end do
        write ( * , 1200 )
        write ( * , * ) 'Son inverse'
        write ( * , 1200 )
        do i_c = 1 , N_NPT_M
          write ( * , 1100 )
     $  ( J_INV ( i_c , j_c ) , j_c = 1 , N_NPT_M )
        end do
        end if
C=============================
C=====Calcul du pas NR complet
C=====PAS_NR = -(J-1).F
C=============================
         CALL
     $   PROD_MAT_VECT
     $ ( N_NPT_M , J_INV , F_NPT_M , PAS_NR )
         PAS_NR = - PAS_NR
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
       write ( * , * ) 'Dans sys_NR :'
       write ( * , * ) 'Calcul du pas NR complet :'
       write ( * , * ) 'PAS_NR = -(J-1).F'
       write ( * , * ) 'avec F_NPT = ' , F_NPT_M
       write ( * , * ) '-> PAS_NR complet (avant NRCG) = ' , PAS_NR
       end if
C====================================
C=====Réduction éventuelle de ce pas
C=====méthode de convergence globale)
C====================================
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Calcul de la fonction auxiliaire phi et de son gradient en x_0
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        PHI_0 = 0.D0
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
         write(*,*)'phi_0=',phi_0
        end if
        DO I = 1 , N_NPT_M
          PHI_0
     $  = PHI_0 + F_NPT_M ( I ) * F_NPT_M ( I )
        END DO
        PHI_0 = 0.5D0 * PHI_0
        GRAD_PHI_0 = 0.D0
	DO J = 1 , N_NPT_M
	  DO I = 1 , N_NPT_M
	      GRAD_PHI_0 ( J )
     $    = GRAD_PHI_0 ( J ) 
     $    + F_NPT_M ( I ) * J_NPT_M ( I , J )
	  END DO
	END DO
C- - - - - - - - - - - - - - - - - - - 
C-----Réduction éventuelle du pas de NR
C- - - - - - - - - - - - - - - - - - - 
      X_0 = X_NPT_M     ! Sauvegarde du point courant
                        ! avant appel à calc_pas_NRCG
                        ! (calc_pas_NRCG reçoit X_0 en entrée
                        ! et fournit X_NPT_M modifié en sortie)
C---------------------------------------------------------------------
C-----Paramètre de contrôle de l'écriture à l'écran par CALC_PAS_NRCG
C-----(pour limiter les écritures intermédiaires, lourdes si Tx ou xT)
C---------------------------------------------------------------------
         I_ECRIT_CALC_PAS_NRCG = 0
         CALL
     $   CALC_PAS_NRCG 
C----------------------------------------------------------
C----Paramètres spécifiques ADPI de la fonction vectorielle
C----------------------------------------------------------
     $ ( N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C----------------------------------------------------
C-----Paramètres déterminés en dehors de la procédure
C----------------------------------------------------
     $   N_NPT_M ,
     $   X_0 , PHI_0 , GRAD_PHI_0 ,
     $   ALPHA_NRCG_NPT_M ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT_M ,
     $   INDIC_TYPE_REDUC_NRCG_NPT_M , COEF_REDUC_NRCG_NPT_M ,
C--------------------------------------------------------------------
C-----Paramètre de contrôle de l'écriture à l'écran par CALC_PAS_NRCG
C--------------------------------------------------------------------
     $   I_ECRIT_CALC_PAS_NRCG ,
C-----------------------------------------------
C-----Paramètres déterminés dans cette procédure
C-----------------------------------------------
     $   PAS_NR , I_FONCTION_NON_DEFINIE ,
     $   X_NPT_M , PHI_1 ) 
C--------------------------------------------------------------
C-----Nouveau pas réduit après que X_NPT a été modifié par NRCG
C--------------------------------------------------------------
C Est-ce utile, puisque PAS_NR est également transmis par calc_pas_NRCG
C ci-dessus ? De plus, PAS_NR n'est plus utilisé ensuite dans sys_NR.
      PAS_NR = X_NPT_M - X_0
C-------------
C-----Ecriture
C-------------
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
      write ( * , 1200 )
      write ( * , * ) 'Nouveau point après NRCG'
      write ( * , 1200 )
      write ( * , * ) X_NPT_M
      write ( * , 1200 )
      write ( * , * ) 'Nouveau PAS_NR = ' , PAS_NR
       end if
C---------------------------------------------
C-----Test de convergence sur le vecteur-écart
C-----(maximum des composantes)
C---------------------------------------------
        VAL_TEST = 0.D0
        DO I_C = 1 , N_NPT_M
	  IF ( DABS ( PAS_NR ( I_C ) ) .GT. VAL_TEST )
     $    VAL_TEST = DABS ( PAS_NR ( I_C ) )
	END DO
        IF ( VAL_TEST .LE. PRECISION_NPT_M ) THEN
          I_CONV = 1
         if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then     
          write(*,*)
     $ '       *********************************************'
          write(*,*)
     $ '       ***** Critère de convergence NR atteint *****'
          write(*,*)
     $ '       *********************************************'
         end if
        END IF
C------------------------------------------------------
C-----Sortie de la boucle sans atteindre la convergence
C------------------------------------------------------
        IF ( ( N_ITER .EQ. N_ITER_MAX_NPT_M ) 
     $ .AND. ( I_CONV .EQ. 0 ) ) THEN
         if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
        WRITE ( * , * ) '-----------------------------------'
        WRITE ( * , * ) "Nombre maximal d'itérations atteint"
	  WRITE ( * , * ) 'sans convergence'
        WRITE ( * , * ) '-----------------------------------'
         end if
        END IF 
C##################################
C#####Fin de la boucle d'itérations
C##################################
      END DO
C     write ( * , * ) '----------'
C     write ( * , * ) 'Fin sys_NR'
C     write ( * , * ) '----------'
      RETURN
      END
C     ==================================================================
C     I                                                                I
C     I Interpolation linéaire a*x+b de 2 points	               I
C     I (en vue de l'estimation d'un troisième)			       I
C     I                   				               I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 18/10/2005                              I
C     I Saint Luc	                                               I
C     ==================================================================
C----------------------------------------
C-----Interruption du programme
C-----(commenté si déjà déclaré ailleurs)
C----------------------------------------
C     INCLUDE
C    $ '/home/besson/SOURCES/PROCEDURES/interruption.inc'
C     ##################################################################
         SUBROUTINE
     $   INTER_LIN
     $ ( X_0_M , X_1_M ,
     $   Y_0_M , Y_1_M ,
     $   A_M , B_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      DELTA_X = X_1_M - X_0_M 
      DELTA_Y = Y_1_M - Y_0_M
      DELTA = X_1_M * Y_0_M - X_0_M * Y_1_M
      IF ( DELTA_X .EQ. 0.D0 ) THEN
        WRITE ( * , * ) '----------------------'
	  WRITE ( * , * ) 'Interpolation linéaire :'
	  WRITE ( * , * ) 'dénominateur nul'
        WRITE ( * , * ) '----------------------'
	CALL INTERRUPTION
      END IF
      A_M = DELTA_Y / DELTA_X
      B_M = DELTA / DELTA_X
      RETURN
      END
C------------------------------------------
C-----Interpolation parabolique de 3 points
C-----(pour calcul d'un quatrième)
C------------------------------------------
C      INCLUDE
C     $'inter_para.inc'
C     ==================================================================
C     I                                                                I
C     I Interpolation parabolique a*x2+b*x+c de 3 points	       I
C     I (en vue de l'estimation d'un quatrième)			       I
C     I                   				               I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 18/10/2005                              I
C     I Saint Luc	                                               I
C     ==================================================================
C----------------------------------------
C-----Interruption du programme
C-----(commenté si déjà déclaré ailleurs)
C----------------------------------------
C     INCLUDE
C    $ '/home/besson/SOURCES/PROCEDURES/interruption.inc'
C     ##################################################################
         SUBROUTINE
     $   INTER_PARA
     $ ( X_0_M , X_1_M , X_2_M ,
     $   Y_0_M , Y_1_M , Y_2_M ,
     $   A_M , B_M , C_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
      DELTA_0 = ( X_2_M * * 2 - X_0_M * * 2 ) * ( X_1_M - X_0_M )
     $        - ( X_1_M * * 2 - X_0_M * * 2 ) * ( X_2_M - X_0_M )
      IF ( DELTA_0 .EQ. 0.D0 ) THEN
        WRITE ( * , * ) '-------------------------'
	WRITE ( * , * ) 'Interpolation parabolique :'
	WRITE ( * , * ) 'déterminant nul'
        WRITE ( * , * ) '-------------------------'
	CALL INTERRUPTION
      END IF
      DELTA_A = ( Y_2_M - Y_0_M ) * ( X_1_M - X_0_M )
     $        - ( Y_1_M - Y_0_M ) * ( X_2_M - X_0_M )
      DELTA_B = ( X_2_M * * 2 - X_0_M * * 2 ) * ( Y_1_M - X_0_M )
     $        - ( X_1_M * * 2 - X_0_M * * 2 ) * ( Y_2_M - X_0_M )
      A_M = DELTA_A / DELTA_0
      B_M = DELTA_B / DELTA_0
      C_M = Y_0_M - A_M * X_0_M * * 2 - B_M * X_0_M
      RETURN
      END
C------------------------------------------------------
C-----Calcul des quantités thermodynamiques dans l'ADPI
C------------------------------------------------------
C      INCLUDE
C     $'G_adpi.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul des quantités thermodynamiques dans l'ADPI              I
C     I (énergie et volume par maille, entropie de configuration,      I
C     I  énergie libre par maille, quantités par atome)		       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 28/03/2006                              I
C     I Saint Jean de Capistran                                        I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   G_ADPI
     $ ( N_TYP_M , N_TYP_INTR_M , 
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M , 
     $   P_R_M , 
     $   N_1_MAILLE_M , N_2_MAILLE_M , N_3_MAILLE_M , 
     $   E_REF_MAILLE_M , V_REF_MAILLE_M ,
     $   E_REF_TYP_M , 
     $   X_AT_M , 
     $   X_D_R_M , 
     $   E_GC_D_R_M , V_GC_D_R_M ,
     $   INDIC_COMPLEXES_M , N_TYPES_COMPLEXES_M ,
     $   MULTIPLICITE_COMPLEXE_M ,
     $   I_S_R_MULTIPLICITE_COMPLEXE_M ,
     $   E_GC_D_COMPLEXE_M , V_GC_D_COMPLEXE_M ,
     $   X_D_COMPLEXE_M ,
     $   TEMPERATURE_M , PRESSION_M ,
     $   N_AT_MAILLE_M ,
     $   E_MAILLE_M , V_MAILLE_M ,
     $   S_CONF_MAILLE_M , G_MAILLE_M , 
     $   E_AT_M , V_AT_M ,
     $   S_CONF_AT_M , G_AT_M ,
     $   G_AT_FORM_M )
C     ##################################################################
C########################################
C#####Utilisation du module de constantes
C########################################
      USE CONSTANTES
C###################################
C#####Types implicites des variables
C###################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
      INTEGER * 4 P_R_M ( 1 : N_R_M ) 
      REAL * 8 E_REF_TYP_M ( 1 : N_TYP_M )
      REAL * 8 X_AT_M ( 1 : N_TYP_M )
      REAL * 8 X_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
      REAL * 8 E_GC_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
      REAL * 8 V_GC_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
C--------------------------------------
C-----Nombre d'atomes par maille (réel)
C--------------------------------------
      REAL * 8 N_AT_MAILLE_M
C------------------------------------------------
C-----Indicateur de présence de défauts complexes
C------------------------------------------------
      CHARACTER * 1 INDIC_COMPLEXES_M
C--------------------------------------
C-----Quantités relatives aux complexes
C--------------------------------------
      INTEGER * 4 MULTIPLICITE_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      INTEGER * 4 I_S_R_MULTIPLICITE_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 E_GC_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 V_GC_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 X_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C-----------------------------------------------------
C-----Termes x et x.ln(x) pour les éléments d'addition
C-----------------------------------------------------
      REAL * 8 X_I ( N_TYP_INTR_M + 1 : N_TYP_M )
      REAL * 8 X_LOG_X_I ( N_TYP_INTR_M + 1 : N_TYP_M )
C        write ( * , * )
C    $   N_TYP_M , N_TYP_INTR_M ,
C    $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
C    $   P_R_M ,
C    $   N_1_MAILLE_M , N_2_MAILLE_M , N_3_MAILLE_M ,
C    $   E_REF_MAILLE_M , V_REF_MAILLE_M ,
C    $   E_REF_TYP_M ,
C    $   X_AT_M 
      write(*,*) "    Début de G_ADPI : taux de DP (1 ligne = 1 s-r) :"
      do I_R = 1 , N_R_M
       write (*,*) ( X_D_R_M ( I_TYP , I_R ) , I_TYP = 0 , N_TYP_M )
      end do
C===============================================================
C=====Energie et volume par maille, et entropie de configuration
C===============================================================
C----------------------
C-----Valeurs initiales
C----------------------
          E_MAILLE_M = E_REF_MAILLE_M
          V_MAILLE_M = V_REF_MAILLE_M
          S_CONF_MAILLE_M = 0.D0
C===================
C=====Sous-réseaux 1
C===================
C     write ( * , * ) 'Sous-réseaux 1'
      DO  I_R = 1 , N_R_1_M
C----------------------------
C-----Quantités préliminaires
C----------------------------
	D_P_R_COUR = DFLOAT ( P_R_M ( I_R ) )
	X_0 = X_D_R_M ( 0 , I_R )
	IF ( X_0 .LE. 0.D0 ) THEN
	  X_0 = 0.D0
	  X_LOG_X_0 = 0.D0
        ELSE
	  X_LOG_X_0 = X_0 * DLOG ( X_0 )
	END IF
	X_2 = X_D_R_M ( 2 , I_R )
        IF ( X_2 .LE. 0.D0 ) THEN
          X_2 = 0.D0
          X_LOG_X_2 = 0.D0
	ELSE
	  X_LOG_X_2 = X_2 * DLOG ( X_2 )
        END IF
        IF ( N_TYP_INTR_M .GT. 2 ) THEN
         X_3 = X_D_R_M ( 3 , I_R )
         IF ( X_3 .LE. 0.D0 ) THEN
           X_3 = 0.D0
           X_LOG_X_3 = 0.D0
         ELSE
	   X_LOG_X_3 = X_3 * DLOG ( X_3 )
         END IF
        END IF 
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          X_I ( I_TYP ) = X_D_R_M ( I_TYP , I_R )
          IF ( X_I ( I_TYP ) .LE. 0.D0 ) THEN
            X_LOG_X_I ( I_TYP ) = 0.D0
          ELSE
            X_LOG_X_I ( I_TYP ) = X_I ( I_TYP )
     $                          * DLOG ( X_I ( I_TYP ) )
          END IF
        END DO
C------------------------------------------------------------
C-----Energie et volume - contributions de 0, 2, 3 (éventuel)
C------------------------------------------------------------
 	E_MAILLE_M = E_MAILLE_M
     $             + D_P_R_COUR
     $             * ( E_GC_D_R_M ( 0 , I_R ) * X_0
     $               + E_GC_D_R_M ( 2 , I_R ) * X_2 )
        V_MAILLE_M = V_MAILLE_M
     $             + D_P_R_COUR
     $             * ( V_GC_D_R_M ( 0 , I_R ) * X_0
     $               + V_GC_D_R_M ( 2 , I_R ) * X_2 )
        IF ( N_TYP_INTR_M .GT. 2 ) THEN
        E_MAILLE_M = E_MAILLE_M
     $             + D_P_R_COUR
     $             * E_GC_D_R_M ( 3 , I_R ) * X_3 
        V_MAILLE_M = V_MAILLE_M
     $             + D_P_R_COUR
     $             * V_GC_D_R_M ( 3 , I_R ) * X_3
        END IF
C--------------------------------------------------------------------
C-----Entropie de configuration - contributions de 0, 2, 3 (éventuel)
C--------------------------------------------------------------------
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR
     $                  * ( X_LOG_X_0 + X_LOG_X_2 )
	TERME_S_CONF_MAILLE = 1.D0 - X_0 - X_2 
        IF ( N_TYP_INTR_M .GT. 2 ) THEN
         S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                   + D_P_R_COUR
     $                   * X_LOG_X_3
         TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_3
        END IF
C------------------------------------------
C-----Contributions des éléments d'addition
C------------------------------------------
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          E_MAILLE_M = E_MAILLE_M
     $               + D_P_R_COUR
     $               * E_GC_D_R_M ( I_TYP , I_R )
     $               * X_I ( I_TYP )
          S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                    + D_P_R_COUR
     $                    * X_LOG_X_I ( I_TYP )
C------------------------------------
C-----Terme complémentaire d'entropie
C------------------------------------
	  TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_I ( I_TYP )
	END DO
	IF ( TERME_S_CONF_MAILLE .LE. 0.D0 ) THEN
	  TERME_LOG_TERME = 0.D0
	ELSE
	  TERME_LOG_TERME = TERME_S_CONF_MAILLE
     $                    * DLOG ( TERME_S_CONF_MAILLE )
	END IF
	S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR * TERME_LOG_TERME
      END DO
C===================
C=====Sous-réseaux 2
C===================
C     write ( * , * ) 'Sous-réseaux 2'
      DO  I_R = 1 + N_R_1_M , N_R_2_M + N_R_1_M
C----------------------------
C-----Quantités préliminaires
C----------------------------
        D_P_R_COUR = DFLOAT ( P_R_M ( I_R ) )
        X_0 = X_D_R_M ( 0 , I_R )
        IF ( X_0 .LE. 0.D0 ) THEN
          X_0 = 0.D0
          X_LOG_X_0 = 0.D0
	ELSE
	  X_LOG_X_0 = X_0 * DLOG ( X_0 )
        END IF
        X_1 = X_D_R_M ( 1 , I_R )
        IF ( X_1 .LE. 0.D0 ) THEN
          X_1 = 0.D0
          X_LOG_X_1 = 0.D0
	ELSE
	  X_LOG_X_1 = X_1 * DLOG ( X_1 )
        END IF
        IF ( N_TYP_INTR_M .GT. 2 ) THEN
         X_3 = X_D_R_M ( 3 , I_R )
         IF ( X_3 .LE. 0.D0 ) THEN
          X_3 = 0.D0
          X_LOG_X_3 = 0.D0
         ELSE
          X_LOG_X_3 = X_3 * DLOG ( X_3 )
         END IF
        END IF
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          X_I ( I_TYP ) = X_D_R_M ( I_TYP , I_R )
          IF ( X_I ( I_TYP ) .LE. 0.D0 ) THEN
            X_LOG_X_I ( I_TYP ) = 0.D0
          ELSE
            X_LOG_X_I ( I_TYP ) = X_I ( I_TYP )
     $                          * DLOG ( X_I ( I_TYP ) )
          END IF
        END DO
C------------------------------------------------------------
C-----Energie et volume - contributions de 0, 1, 3 (éventuel)
C------------------------------------------------------------
        E_MAILLE_M = E_MAILLE_M
     $             + D_P_R_COUR
     $             * ( E_GC_D_R_M ( 0 , I_R ) * X_0
     $               + E_GC_D_R_M ( 1 , I_R ) * X_1 )
        V_MAILLE_M = V_MAILLE_M
     $             + D_P_R_COUR
     $             * ( V_GC_D_R_M ( 0 , I_R ) * X_0
     $               + V_GC_D_R_M ( 1 , I_R ) * X_1 )
       IF ( N_TYP_INTR_M .GT. 2 ) THEN
        E_MAILLE_M = E_MAILLE_M
     $             + D_P_R_COUR
     $              * E_GC_D_R_M ( 3 , I_R ) * X_3 
        V_MAILLE_M = V_MAILLE_M
     $             + D_P_R_COUR
     $             * V_GC_D_R_M ( 3 , I_R ) * X_3 
       END IF
C--------------------------------------------------------------------
C-----Entropie de configuration - contributions de 0, 1, 3 (éventuel)
C--------------------------------------------------------------------
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR
     $  	        * ( X_LOG_X_0 + X_LOG_X_1 )
        TERME_S_CONF_MAILLE = 1.D0 - X_0 - X_1
        IF ( N_TYP_INTR_M .GT. 2 ) THEN
         S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                   + D_P_R_COUR
     $                   * X_LOG_X_3 
         TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_3
        END IF
C------------------------------------------
C-----Contributions des éléments d'addition
C------------------------------------------
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          E_MAILLE_M = E_MAILLE_M
     $               + D_P_R_COUR 
     $               * E_GC_D_R_M ( I_TYP , I_R )
     $               * X_I ( I_TYP )
          S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                    + D_P_R_COUR
     $                    * X_LOG_X_I ( I_TYP )
C------------------------------------
C-----Terme complémentaire d'entropie
C------------------------------------
          TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_I ( I_TYP )
        END DO
        IF ( TERME_S_CONF_MAILLE .LE. 0.D0 ) THEN
          TERME_LOG_TERME = 0.D0
	ELSE
	  TERME_LOG_TERME = TERME_S_CONF_MAILLE
     $                    * DLOG ( TERME_S_CONF_MAILLE )
        END IF
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR * TERME_LOG_TERME
      END DO
C===================
C=====Sous-réseaux 3
C===================
C     write ( * , * ) 'Sous-réseaux 3'
      DO  I_R = 1 + N_R_1_M + N_R_2_M , N_R_3_M + N_R_2_M + N_R_1_M
C----------------------------
C-----Quantités préliminaires
C----------------------------
        D_P_R_COUR = DFLOAT ( P_R_M ( I_R ) )
        X_0 = X_D_R_M ( 0 , I_R )
        IF ( X_0 .LE. 0.D0 ) THEN
          X_0 = 0.D0
          X_LOG_X_0 = 0.D0
        ELSE
          X_LOG_X_0 = X_0 * DLOG ( X_0 )
        END IF
        X_1 = X_D_R_M ( 1 , I_R )
        IF ( X_1 .LE. 0.D0 ) THEN
          X_1 = 0.D0
          X_LOG_X_1 = 0.D0
        ELSE
          X_LOG_X_1 = X_1 * DLOG ( X_1 )
        END IF
        X_2 = X_D_R_M ( 2 , I_R )
        IF ( X_2 .LE. 0.D0 ) THEN
          X_2 = 0.D0
          X_LOG_X_2 = 0.D0
        ELSE
          X_LOG_X_2 = X_2 * DLOG ( X_2 )
        END IF
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          X_I ( I_TYP ) = X_D_R_M ( I_TYP , I_R )
          IF ( X_I ( I_TYP ) .LE. 0.D0 ) THEN
            X_LOG_X_I ( I_TYP ) = 0.D0
          ELSE
            X_LOG_X_I ( I_TYP ) = X_I ( I_TYP )
     $                          * DLOG ( X_I ( I_TYP ) )
          END IF
        END DO
C-------------------------------------------------
C-----Energie et volume - contributions de 0, 1, 2
C-------------------------------------------------
        E_MAILLE_M = E_MAILLE_M
     $             + D_P_R_COUR
     $             * ( E_GC_D_R_M ( 0 , I_R ) * X_0
     $               + E_GC_D_R_M ( 1 , I_R ) * X_1
     $               + E_GC_D_R_M ( 2 , I_R ) * X_2 )
        V_MAILLE_M = V_MAILLE_M
     $             + D_P_R_COUR
     $             * ( V_GC_D_R_M ( 0 , I_R ) * X_0
     $               + V_GC_D_R_M ( 1 , I_R ) * X_1
     $               + V_GC_D_R_M ( 2 , I_R ) * X_2 )
C---------------------------------------------------------
C-----Entropie de configuration - contributions de 0, 1, 2
C---------------------------------------------------------
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR
     $                  * ( X_LOG_X_0 + X_LOG_X_1 + X_LOG_X_2 )
        TERME_S_CONF_MAILLE = 1.D0 - X_0 - X_1 - X_2
C------------------------------------------
C-----Contributions des éléments d'addition
C------------------------------------------
        DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
          E_MAILLE_M = E_MAILLE_M
     $               + D_P_R_COUR
     $               * E_GC_D_R_M ( I_TYP , I_R )
     $               * X_I ( I_TYP )
          S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                    + D_P_R_COUR
     $                    * X_LOG_X_I ( I_TYP )
C------------------------------------
C-----Terme complémentaire d'entropie
C------------------------------------
          TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_I ( I_TYP ) 
        END DO
        IF ( TERME_S_CONF_MAILLE .LE. 0.D0 ) THEN
          TERME_LOG_TERME = 0.D0
        ELSE
	  TERME_LOG_TERME = TERME_S_CONF_MAILLE
     $                    * DLOG ( TERME_S_CONF_MAILLE )
        END IF
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR * TERME_LOG_TERME
      END DO
C===============================
C=====Sous-réseaux interstitiels
C===============================
C     write ( * , * ) 'Sous-réseaux interstitiels'
      DO  I_R = 1 + N_R_3_M + N_R_2_M + N_R_1_M , N_R_M
C----------------------------
C-----Quantités préliminaires
C----------------------------
        D_P_R_COUR = DFLOAT ( P_R_M ( I_R ) )
        TERME_S_CONF_MAILLE = 1.D0
        DO I_TYP = 1 , N_TYP_M
          X_COUR = X_D_R_M ( I_TYP , I_R )
          IF ( X_COUR .LE. 0.D0 ) THEN
            X_COUR = 0.D0
            X_LOG_X_COUR = 0.D0
	  ELSE
	    X_LOG_X_COUR = X_COUR * DLOG ( X_COUR )
          END IF
C----------------------
C-----Energie et volume
C----------------------
          E_MAILLE_M = E_MAILLE_M
     $               + D_P_R_COUR
     $               * E_GC_D_R_M ( I_TYP , I_R ) * X_COUR
          V_MAILLE_M = V_MAILLE_M
     $               + D_P_R_COUR
     $               * V_GC_D_R_M ( I_TYP , I_R ) * X_COUR
C------------------------------
C-----Entropie de configuration
C------------------------------
          S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                    + D_P_R_COUR * X_LOG_X_COUR
	  TERME_S_CONF_MAILLE = TERME_S_CONF_MAILLE - X_COUR
        END DO
        IF ( TERME_S_CONF_MAILLE .LE. 0.D0 ) THEN
          TERME_LOG_TERME = 0.D0
	ELSE
        TERME_LOG_TERME = TERME_S_CONF_MAILLE
     $                  * DLOG ( TERME_S_CONF_MAILLE )
        END IF
        S_CONF_MAILLE_M = S_CONF_MAILLE_M
     $                  + D_P_R_COUR * TERME_LOG_TERME
      END DO
C============================================================
C=====Energie, volume et entropie dus aux complexes éventuels
C============================================================
      E_MAILLE_COMPLEXES = 0.D0
      V_MAILLE_COMPLEXES = 0.D0
      S_CONF_MAILLE_COMPLEXES = 0.D0
         IF
     $ ( INDIC_COMPLEXES_M .EQ. 'O' .OR. INDIC_COMPLEXES_M .EQ. 'o' ) 
     $   THEN
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES_M
C------------
C-----Energie
C------------
         E_MAILLE_COMPLEXES
     $ = E_MAILLE_COMPLEXES
     $ + X_D_COMPLEXE_M ( I_COMPLEXE )
     $ * DFLOAT
     $ ( P_R_M ( I_S_R_MULTIPLICITE_COMPLEXE_M ( I_COMPLEXE ) ) )
     $ * E_GC_D_COMPLEXE_M ( I_COMPLEXE )	
C-----------
C-----Volume
C-----------
         V_MAILLE_COMPLEXES
     $ = V_MAILLE_COMPLEXES
     $ + X_D_COMPLEXE_M ( I_COMPLEXE )
     $ * DFLOAT
     $ ( P_R_M ( I_S_R_MULTIPLICITE_COMPLEXE_M ( I_COMPLEXE ) ) )
     $ * V_GC_D_COMPLEXE_M ( I_COMPLEXE )
C-------------
C-----Entropie	
C-------------
         X_COUR = X_D_COMPLEXE_M ( I_COMPLEXE )
         IF ( X_COUR .LE. 0.D0 ) THEN
            X_COUR = 0.D0
            X_LOG_X_COUR = 0.D0
         ELSE
            X_LOG_X_COUR = X_COUR * DLOG ( X_COUR )
          END IF
         S_CONF_MAILLE_COMPLEXES
     $ = S_CONF_MAILLE_COMPLEXES
     $ + X_COUR
     $ * DFLOAT
     $ ( P_R_M ( I_S_R_MULTIPLICITE_COMPLEXE_M ( I_COMPLEXE ) ) )
     $ * ( 1.D0
     $   + DLOG
     $   ( DFLOAT ( MULTIPLICITE_COMPLEXE_M ( I_COMPLEXE ) ) ) )
     $ - X_LOG_X_COUR
     $ * DFLOAT
     $ ( P_R_M ( I_S_R_MULTIPLICITE_COMPLEXE_M ( I_COMPLEXE ) ) )
        END DO
      END IF
C================================================================
C=====Energie, volume, entropie de configuration et énergie libre
C================================================================
         E_MAILLE_M
     $ = E_MAILLE_M + E_MAILLE_COMPLEXES
         V_MAILLE_M
     $ = V_MAILLE_M + V_MAILLE_COMPLEXES
           S_CONF_MAILLE_M
     $ = - K_B * S_CONF_MAILLE_M
     $   + K_B * S_CONF_MAILLE_COMPLEXES
         G_MAILLE_M
     $ = E_MAILLE_M
     $ - TEMPERATURE_M * S_CONF_MAILLE_M
     $ + PRESSION_M * V_MAILLE_M
C-----------------------------------------------------------------------
C-----Nombre d'atomes par maille et quantités thermodynamiques par atome
C-----------------------------------------------------------------------
         N_AT_MAILLE_M = DFLOAT ( N_1_MAILLE_M )
     $                 + DFLOAT ( N_2_MAILLE_M )
     $                 + DFLOAT ( N_3_MAILLE_M )
         DO I_R = 1 , N_R_1_M + N_R_2_M + N_R_3_M
	    N_AT_MAILLE_M = N_AT_MAILLE_M
     $                    - X_D_R_M ( 0 , I_R )
     $                    * DFLOAT ( P_R_M ( I_R ) )
         END DO 
         DO I_R = 1 + N_R_1_M + N_R_2_M + N_R_3_M , N_R_M
	  DO I_TYP = 1 , N_TYP_M
            N_AT_MAILLE_M = N_AT_MAILLE_M
     $                    + X_D_R_M ( I_TYP , I_R )
     $                    * DFLOAT ( P_R_M ( I_R ) )
	  END DO 
	 END DO
         E_AT_M = E_MAILLE_M / N_AT_MAILLE_M
         V_AT_M = V_MAILLE_M / N_AT_MAILLE_M
         S_CONF_AT_M = S_CONF_MAILLE_M / N_AT_MAILLE_M
         G_AT_M = G_MAILLE_M / N_AT_MAILLE_M
C-----------------------------------------
C-----Energie libre de formation par atome       
C-----------------------------------------
        G_AT_FORM_M = G_AT_M
        DO I_TYP = 1 , N_TYP_M
	 G_AT_FORM_M = G_AT_FORM_M
     $               - X_AT_M ( I_TYP ) * E_REF_TYP_M ( I_TYP )
        END DO
      RETURN
      END
C     ##################################################################
C----------------------------------------
C-----Fonction du système en mode NPT
C-----écrite pour N éléments intrinsèques
C----------------------------------------
C     INCLUDE
C    $'f_NR_N_intr.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul de la fonction vectorielle de D variables	       I
C     I spécifiée sous forme analytique 			       I
C     I intervenant dans l'ADPI NPT				       I
C     I pour l'algorithme de Newton-Raphson                            I
C     I								       I
C     I Pour un composé ave N éléments intrinsèques,                   I
C     I I éléments d'addition et R sous-réseaux		               I
C     I la dimension du système est :				       I
C     I D = 1 + R * ( N + I )			                       I
C     I								       I
C     I Les inconnues sont :					       I
C     I (i) le nombre de mailles				       I
C     I (ii) les logarithmes (log_10) des fractions de DP	       I
C     I								       I
C     I Remarque : la présente version est écrite pour N <= 3,         I
C     I mais le passage à N quelconque est possible		       I
C     I								       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 16/07/2020                              I
C     I Notre Dame du Mont Carmel                                      I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   F_NR_AVEC_DERIV
     $ (
C----------------------------------------------
C----Paramètres spécifiques ADPI de la fonction
C----------------------------------------------
     $   N_TYP_M , N_TYP_INTR_M ,
     $   N_R_1_M , N_R_2_M , N_R_3_M , N_R_M ,
     $   P_R_M ,
     $   IND_D_R_TYP_M ,
     $   H_GC_D_R_M ,
     $   K_T_M ,
     $   X_AT_M , N_AT_TOT_M ,
     $   H_REF_MAILLE_M ,
C--------------------------------------
C----Paramètres généraux de la fonction
C--------------------------------------
     $   I_FONCTION_NON_DEFINIE_M ,
     $   N_NPT_M , X_NPT_M , F_NPT_M ,
     $   INDIC_CALC_DERIV_M , MAT_DERIV_NPT_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C------------------------------------
C-----Vecteur et fonction vectorielle
C------------------------------------
      REAL * 8 X_NPT_M ( N_NPT_M ) ,
     $         F_NPT_M ( N_NPT_M )
C-----------------------------------------
C-----Matrice des dérivées partielles de F
C-----------------------------------------
      REAL * 8 MAT_DERIV_NPT_M ( N_NPT_M , N_NPT_M )
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-réseau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-réseau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-réseau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantité de matière totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de référence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C-----------------------------------------------------------------------
C-----Somme des nombres de sites par maille
C-----sur les sous-réseaux relatifs à chaque espèce chimique intrinsèque
C-----------------------------------------------------------------------
      INTEGER * 4 P_R_1 , P_R_2 , P_R_3
C----------------------------- 
C-----Nombre de mailles (réel)
C-----------------------------
      REAL * 8 M_MAILLES
C----------------------------------------------------------
C-----Fractions de DP en fonction des sous-réseaux et types
C----------------------------------------------------------
      REAL * 8 X_D_R ( 0 : N_TYP_M , N_R_M )
C-------------------------
C-----Quantités de matière
C-------------------------
      REAL * 8 N_AT ( N_TYP_M )
C--------------------------------------
C-----Termes complémentaires d'entropie
C--------------------------------------
      REAL * 8 Z_TYP_R ( N_R_M )
C-----Indicateur (interne à cette procédure) d'écriture de détails
      INDIC_ECRIT_DETAILS_F_NR = 0
C-----Notification du calcul analytique des dérivées
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
       if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
	  write ( * , * )
     $ " ==> F_NR en cours avec calcul analytique des dFi/dxj"  
       end if
      END IF
C====================
C=====Initialisations
C====================
C----------------------------------
C-----Initialisation de la fonction
C----------------------------------
      F_NPT_M = 0.D0
C----------------------------------------------------------------
C-----Initialisation de X_D_R (nécessaire pour le calcul des Z_R)
C----------------------------------------------------------------
      X_D_R = 0.D0
C-----------------------------------------------------------------------
C-----Somme des nombres de sites par maille
C-----sur les sous-réseaux relatifs à chaque espèce chimique intrinsèque
C-----------------------------------------------------------------------
      P_R_1 = 0
      P_R_2 = 0
      P_R_3 = 0
      DO I_R = 1 , N_R_1_M
	  P_R_1 = P_R_1 + P_R_M ( I_R )
      END DO
      DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
	  P_R_2 = P_R_2 + P_R_M ( I_R )
      END DO
      DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
	  P_R_3 = P_R_3 + P_R_M ( I_R )
      END DO
C-------------------------
C-----Quantités de matière
C-------------------------
      DO I_TYP = 1 , N_TYP_M
        N_AT ( I_TYP ) = N_AT_TOT_M * X_AT_M ( I_TYP )
      END DO
C------------------------------------
C-----Matrice des dérivées partielles
C-----mise à jour tous les P_J_M pas
C------------------------------------
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        MAT_DERIV_NPT_M = 0.D0
      END IF
C-------------
C-----Ecriture
C-------------
C	write ( * , * ) 'Dans f_NR, X_AT_M = ' , 	
C    $  X_AT_M
C	write ( * , * ) 'Dans f_NR, N_AT_TOT_M = ' ,
C    $  N_AT_TOT_M
C        write ( * , * ) 'Dans f_NR, N_AT = ' ,
C    $  N_AT
C=========================================
C=====Lien entre variables "muettes" X_NPT
C=====et variables du problème
C=========================================
C     write ( * , * ) 'N_R_1_M = ' , N_R_1_M
C     write ( * , * ) 'N_R_2_M = ' , N_R_2_M
      DO I_R = 1 , N_R_M
         DO I_TYP = 0 , N_TYP_M
C	  write ( * , * ) 'I_R = ' , I_R 
C	  write ( * , * ) 'I_TYP = ' , I_TYP
C	  write ( * , * ) 'IND_D_R_TYP = ' ,
C    $                     IND_D_R_TYP_M ( I_TYP , I_R )
	 END DO
      END DO	
C----------------------------------------------------------
C-----Variables "log_10 des fractions de défauts ponctuels"
C----------------------------------------------------------
      DO I_R = 1 , N_R_M
 	 DO I_TYP = 0 , N_TYP_M
        IF ( IND_D_R_TYP_M ( I_TYP , I_R ) .NE. 0 ) THEN
	   X_D_R ( I_TYP , I_R )
     $   = DEXP
     $   ( DLOG ( 10.D0 ) * X_NPT_M ( IND_D_R_TYP_M ( I_TYP , I_R ) ) )
        END IF
       END DO 
      END DO      
C	write ( * , * ) 'Dans f_NR :'
C	write ( * , * ) 'X_NPT_M : ' 
C	write ( * , * ) X_NPT_M
C       write ( * , * ) 'X_D_R (1 ligne = 1 s-r.) : '
C       do i_R = 1 , n_r_m 
C        write ( * , * ) ( X_D_R ( I_TYP , I_R ) , I_TYP = 0 , N_TYP_M )
C       end do
C       write ( * , * ) 'H_GC_D_R (1 ligne = 1 s-r.) : '
C       do i_R = 1 , n_r_m
C        write ( * , * )
C    $ ( H_GC_D_R_M ( I_TYP , I_R ) , I_TYP = 0 , N_TYP_M )
C       end do
C----------------------------------------
C-----Variable "nombre de mailles" (réel)
C----------------------------------------
      M_MAILLES = X_NPT_M ( N_NPT_M )
C=================================================
C=====Calcul des termes complémentaires d'entropie
C=================================================
      Z_TYP_R = 1.D0
C-------------------
C-----Sous-réseaux 1
C-------------------
      DO I_R = 1 , N_R_1_M
         DO I_TYP = 0 , N_TYP_M
	   IF ( I_TYP .NE. 1 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
	   END IF
         END DO
      END DO
C-----------------------------
C-----Sous-réseaux 2 éventuels
C-----------------------------
      DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
         DO I_TYP = 0 , N_TYP_M
           IF ( I_TYP .NE. 2 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C-----------------------------
C-----Sous-réseaux 3 éventuels
C-----------------------------
      DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
         DO I_TYP = 0 , N_TYP_M
           IF ( I_TYP .NE. 3 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C-----------------------------------------
C-----Sous-réseaux interstitiels éventuels
C-----------------------------------------
      DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M
         DO I_TYP = 0 , N_TYP_M
           IF ( I_TYP .NE. 0 ) THEN
                Z_TYP_R ( I_R )
     $        = Z_TYP_R ( I_R )
     $        - X_D_R ( I_TYP , I_R )
           END IF
         END DO
      END DO
C       write ( * , * ) '--------'
C	  write ( * , * ) 'Termes Z :'
C       write ( * , * ) '--------'
C	do I_R = 1 , N_R_M
C         write ( * , * ) 'I_R = ' , I_R
C         write ( * , * ) ' Z_TYP_R ( I_R ) = ' 
C	  write ( * , * )
C     $   Z_TYP_R ( I_R ) 
C	end do
C--------------------------------------------------------------------
C-----Test de positivité de ces termes Z
C-----(sinon la fonction n'est pas définie)
C-----Rque : I_FONCTION_NON_DEFINIE_M peut être modifié par F_NR,
C-----mais sa valeur en entrée de F_NR (déjà fixée dans calc_pas_NRCG
C-----par la présence de x_1(i)<0) a aussi de l'importance.
C--------------------------------------------------------------------
      DO I_R = 1 , N_R_M
	 IF ( Z_TYP_R ( I_R ) .LE. 0.D0 ) THEN
        WRITE ( * , * ) 
     $ '----------------------------------------------'
	  WRITE ( * , * )
     $ "Dans f_NR, le terme complémentaire d'entropie"
         WRITE ( * , * ) 
     $ 'Z_TYP_R ( I_R = ' , I_R , ' )'
         WRITE ( * , * ) 
     $ 'est négatif'
	 write ( * , * ) 'Z_TYP_R = ' , Z_TYP_R ( I_R )
	 I_FONCTION_NON_DEFINIE_M = 1
	 END IF
      END DO
C-----Si z<0 au point courant, la procédure s'arrête.
C-----C'est le cas aussi lorsque I_FONCTION_NON_DEFINIE_M 
C-----vaut déjà 1 à l'entrée de la procédure
C-----(valeur 1 affectée par calc_pas_NRCG). 
      IF ( I_FONCTION_NON_DEFINIE_M .EQ. 1 ) THEN
       if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
        WRITE ( * , * )
     $ 'Fonction non définie => le programme quitte f_NR'
        WRITE ( * , * )
     $ ' (pour réduction de PAS_NR dans calc_pas_NRCG)'
        end if
       RETURN
      END IF
C===================================
C=====Multiplicateurs de Lagrange
C=====pour les éléments intrinsèques
C===================================
           X_MULTI_1
     $ = - H_GC_D_R_M ( 0 , N_R_1_M )
     $   - K_T_M
     $   * DLOG ( X_D_R ( 0 , N_R_1_M )
     $          / Z_TYP_R ( N_R_1_M ) )
      IF ( N_TYP_INTR_M .GT. 1 ) THEN
             X_MULTI_2
     $   = - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $     - K_T_M
     $     * DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
     $            / Z_TYP_R ( N_R_1_M + N_R_2_M ) )
       IF ( N_TYP_INTR_M .GT. 2 ) THEN
            X_MULTI_3
     $  = - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $    - K_T_M
     $    * DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $           / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) )
       END IF
      END IF
C###########################################
C#####Composantes de la fonction vectorielle
C#####en fonction des variables du problème
C###########################################
C-------------------------
C-----Indice de composante
C-------------------------
      I_COMP = 0
C===========================================
C=====Composantes de la fonction vectorielle
C=====relatives aux sous-réseaux 1
C===========================================
      DO I_R = 1 , N_R_1_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-réseau 1, il y a une équation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. système d'équations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les éléments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 2 et 3, ainsi que pour les lacunes 
C-----(sauf pour le sous-réseau N_R_1)
C-----et pour les éléments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des éléments intrinsèques
C----------------------------------
C----------------------------------
C-----------------------  
C-----------------------  
C-----Elément 2 éventuel
C-----------------------
C-----------------------  
      IF ( N_TYP_INTR_M .GT. 1 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 2 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 2 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_1 - X_MULTI_2
C     $ * ( DLOG ( X_D_R ( 2 , I_R )
C     $          / Z_TYP_R ( I_R ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) ) )
C      write ( * , * ) 'Dans f_NR :'	
C      write ( * , * ) 'K_T_M = ' , K_T_M	
C      write ( * , * ) 'H_GC_D_R_M ( 2 , I_R ) = ' ,
C    $			H_GC_D_R_M ( 2 , I_R )
C      write ( * , * ) 'H_GC_D_R_M ( 0 , N_R_1_M ) = ' ,
C    $			H_GC_D_R_M ( 0 , N_R_1_M )
C      write ( * , * ) 'H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M ) = ' ,
C    $  		H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )	
C      write ( * , * ) 'X_D_R ( 2 , I_R ) = ' ,
C    $			X_D_R ( 2 , I_R )
C      write ( * , * ) 'X_D_R ( 0 , N_R_1_M ) = ' , 
C    $			X_D_R ( 0 , N_R_1_M )
C      write ( * , * ) 'X_D_R ( 0 , N_R_1_M + N_R_2_M ) = ' ,
C    $			X_D_R ( 0 , N_R_1_M + N_R_2_M )
C      write ( * , * ) ' Z_TYP_R ( I_R ) = ' , 
C    $			 Z_TYP_R ( I_R )
C      write ( * , * ) ' Z_TYP_R ( N_R_1_M ) = ' ,
C    $			 Z_TYP_R ( N_R_1_M )
C      write ( * , * ) ' Z_TYP_R ( N_R_1_M + N_R_2_M ) = ' ,
C    $			 Z_TYP_R ( N_R_1_M + N_R_2_M ) 	
C      write ( * , * )
C    $ 'log X_D_R(2 , I_R)/Z_TYP_R(I_R) = ' ,
C    $     DLOG ( X_D_R ( 2 , I_R )
C    $          / Z_TYP_R ( I_R ) )
C      write ( * , * )
C    $ 'logX_D_R(0 , N_R_1)/Z_TYP_R (N_R_1) = ' ,
C    $     DLOG ( X_D_R ( 0 , N_R_1_M )
C    $          / Z_TYP_R ( N_R_1_M ) )
C      write ( * , * )
C    $ 'logX_D_R(0 , N_R_1 + N_R_2 )/Z_TYP_R(N_R_1 + N_R_2) = ' ,
C    $     DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C    $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) )
C      write ( * , * )
C    $  '-----------------------'
C      write ( * , * )
C    $  'pour s-r 1 et élément 2 : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*****************************************************
C-----cf. système d'équations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les éléments d'addition
C-----(notes du 16/12/19).
C*****************************************************
C***********************************************************************
C-----L'indice J_COUR désigne une variable de la fonction vectorielle,
C-----correspondant à toutes les possibilités (I_R ; I_TYP) telles que
C-----IND_D_R_TYP_M ( I_TYP , I_R ) # 0.
C-----Le cas où IND_D_R_TYP ( I_TYP , I_R ) = 0, qui ne renvoie
C-----à aucune variable de cette fonction, recouvre deux situations :
C----- 1) les combinaisons (I_R ; I_TYP) qui ne renvoient à aucun DP,
C-----i.e. (i) un élément intrinsèque sur son sous-réseau "normal",
C-----    (ii) une lacune sur un sous-réseau interstitiel.
C----- 2) les combinaisons (I_R ; I_TYP) correspondant à un DP possible,
C-----    mais qui est écarté par une valeur nulle d'énergie de SC.
C-----Il faut alors noter que les termes complémentaires z^R
C-----présents dans les composantes Fi de la fonction
C-----ne dépendent que des variables autres que 1) et 2) ci-dessus.
C***********************************************************************
C***********************************************************************
C-----Dans chaque composante de F, le calcul des dérivées partielles
C-----de z^R nécessite d'éviter ces "fausses variables".
C-----Pour ce faire, on a d'abord eu recours à des tests adaptés
C-----aux divers cas (I_R ; I_TYP), e.g. "IF ( J_TYP .NE. 1 )" qui évite
C-----une "fausse variable" de type 1)(i).
C-----Cependant, un tel test n'évite pas le cas 2), e.g. un système
C-----avec certains antisites très coûteux (ou d'énergie incalculable),
C-----auxquels on associe donc une énergie de SC nulle.
C-----Dans ce cas, le programme affecte IND_D_R_TYP = 0 au DP,
C-----qui ne correspond alors plus à une variable du système NPT.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----[NB : la composante Fi (de la fonction vectorielle) associée
C-----à ce DP ne doit alors pas être créée par le programme
C-----(attention : ce n'est pas le cas pour l'instant, cf. supra :
C-----aucun test n'est présent pour écarter certains DP d'antisite
C-----I_TYP=2 sur s-r normalement "1", et Fi est toujours créée :
C-----le rejet de ces Fi inactifs - premier indice de MAT_DERIV_NPT_M -
C-----n'est donc pas implémenté,
C-----alors qu'il doit l'être au même titre que le rejet des dFi/dxj
C----- - deuxième indice de MAT_DERIV_NPT_M).]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Ainsi, pour ne pas inclure les fausses variables 1) ET 2), on peut
C-----remplacer le test précédent par un test "IF ( J_COUR .NE. 0 )",
C-----comme ci-dessous.
C***********************************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter le DP ( 2 , I_R ) (par E_SC = 0)
C----- = fausse variable 2).
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----(la composante Fi doit aussi être écartée,
C-----ce que ne fait pas la version actuelle :
C-----dans ce cas, le passage dans cette section
C-----"dérivées partielles" ne devrait même pas avoir lieu).
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 2 , I_R )
        END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
              if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
               write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
              end if
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
C----- --> si E_SC(lacune N_R_1) = 0, suggérer de permuter
C-----l'ordre des s-r de type 1 dans le fichier DATA.adpi
C-----(ce n'est pas possible si N_R_1 = 1 ! Que faire en NPT alors ?)
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
              if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
               write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
              end if
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M / Z_TYP_R ( N_R_1_M )
            END IF
           END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
C----- --> si E_SC(lacune N_R_1 + N_R_2) = 0, suggérer de permuter
C-----l'ordre des s-r de type 2 dans le fichier DATA.adpi
C-----(ce n'est pas possible si N_R_2 = 1 ! Que faire en NPT alors ?)
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : élément 2 éventuel
C------------------------------------
C------------------------------------
      END IF
C-----------------------
C-----------------------
C-----Elément 3 éventuel
C-----------------------
C-----------------------
      IF ( N_TYP_INTR_M .GT. 2 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 3 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 3 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_1 - X_MULTI_3
C     $ * ( DLOG ( X_D_R ( 3 , I_R )
C     $          / Z_TYP_R ( I_R ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C      write ( * , * )
C    $  '-----------------------'
C      write ( * , * )
C    $  'pour s-r 1 et élément 3 : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
        IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=',J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 3 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M )
            END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 + N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
          if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
           write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
          end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
            END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : élément 3 éventuel
C------------------------------------
C------------------------------------
      END IF
C--------------------------------------------------------
C--------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_1 > 1 
C-----(car la composante Fi pour les lacunes en N_R_1
C-----sert à éliminer le multiplicateur de Lagrange de 1)
C--------------------------------------------------------
C--------------------------------------------------------
        IF ( I_R .LT. N_R_1_M ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 0 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 0 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_1
C     $ * ( DLOG ( X_D_R ( 0 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) ) )
C      write ( * , * )
C    $  '-----------------------'
C      write ( * , * )
C    $  'pour s-r 1 et lacunes : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
             write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C-------------------------
C-------------------------
C-----Fin du cas : lacunes
C-------------------------
C-------------------------
       END IF
C------------------------------------------
C------------------------------------------
C-----Cas des éléments d'addition éventuels
C------------------------------------------
C------------------------------------------
C-------------------------------------
C-------------------------------------
C-----Boucle sur les types d'additions
C-------------------------------------
C-------------------------------------
       DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
C       WRITE(*,*) "s-r 1 : ADDITIONS ET DERIVEES"
C-----------------------------------------------------
C-----Indicateur d'existence d'une nouvelle composante
C-----de la fonction vectorielle étudiée
C-----------------------------------------------------
        I_NOUV_COMP = 0
          IF ( N_R_M .EQ. N_R_1_M ) THEN
           IF ( I_R .LT. N_R_1_M ) THEN 
              DELTA_1 = 0.D0
              I_NOUV_COMP = 1
           END IF
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $      .AND. ( N_R_2_M .NE. 0 ) ) THEN
           DELTA_1 = X_MULTI_1 - X_MULTI_2
           I_NOUV_COMP = 1
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
           DELTA_1 = X_MULTI_1 - X_MULTI_3
           I_NOUV_COMP = 1
        ELSE 
           DELTA_1 = X_MULTI_1
           I_NOUV_COMP = 1
        END IF
C------------------------------------------------------------------
C-----Si nouvelle composante, calcul de celle-ci et de ses dérivées
C------------------------------------------------------------------
        IF ( I_NOUV_COMP .EQ. 1 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( I_TYP , I_R )
     $ - H_GC_D_R_M ( I_TYP , N_R_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( I_TYP , I_R )
     $          / Z_TYP_R ( I_R ) )
     $   - DLOG ( X_D_R ( I_TYP , N_R_M )
     $          / Z_TYP_R ( N_R_M ) ) )
     $ + DELTA_1
C	write ( * , * ) 'Dans f_NR :'
C	write ( * , * ) 'X_D_R ( ' , I_TYP , ' , ' , I_R , ') = ' ,
C    $				X_D_R ( I_TYP , I_R )
C       write ( * , * ) 'Z_TYP_R ( I_R ) = ',
C    $				Z_TYP_R ( I_R )
C	write ( * , * ) 'X_D_R ( ' , I_TYP , ' ,  ' , N_R_M , ' = ' ,
C    $				X_D_R ( I_TYP , N_R_M )
C	write ( * , * ) 'Z_TYP_R ( N_R_M ) = ' ,
C    $				Z_TYP_R ( N_R_M )
C      write ( * , * )
C    $  '-----------------------'
C      write ( * , * )
C    $  'pour s-r 1 et élément add. : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 1 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( I_TYP , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
          J_COUR = IND_D_R_TYP_M ( I_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
           IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( I_TYP , N_R_M )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
           IF ( J_COUR .NE. 0 ) THEN 
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
             write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_1"
C-----On ajoute ci-dessous les dérivées du terme DELTA_1,
C-----en distinguant les différents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas où DELTA_1 = 0 => dérivée nulle
           IF ( N_R_M .EQ. N_R_1_M ) THEN
           IF ( I_R .LT. N_R_1_M ) THEN 
           END IF
C-----Cas où DELTA_1 = X_MULTI_1 - X_MULTI_2 = mu_A - mu_B
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $      .AND. ( N_R_2_M .NE. 0 ) ) THEN
C- - - Dérivées de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1_M 
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M )
            END IF
          END DO
C- - - Dérivées de - mu_B (ici N_R_M = N_R_1_M + N_R_2_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1_M + N_R_2_M 
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas où DELTA_1 = X_MULTI_1 - X_MULTI_3 = mu_A - mu_C
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
C- - - Dérivées de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1_M
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C- - - Dérivées de - mu_C (ici N_R_M = N_R_1_M + N_R_2_M + N_R_3_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
         J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
          DO J_TYP = 0 , N_TYP_M
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1 + N_R_2 + N_R_3
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas où DELTA_1 = X_MULTI_1 = + mu_A
        ELSE 
C- - - Dérivées de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C-----Fin des divers cas
        END IF
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C---------------------------------------------------------------------
C-----Fin du test "existence d'une nouvelle composante de la fonction"
C---------------------------------------------------------------------
         END IF
C-----------------------------------------------
C-----------------------------------------------
C-----Fin de la boucle sur les types d'additions
C-----------------------------------------------
C-----------------------------------------------
       END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 1
C=========================================
      END DO
C===============================================
C=====Composantes de la fonction vectorielle
C=====(éventuelles) relatives aux sous-réseaux 2
C===============================================
      DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-réseau 2, il y a une équation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. système d'équations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les éléments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 1 et 3, ainsi que pour les lacunes 
C-----(sauf pour le sous-réseau N_R_1 + N_R_2)
C-----et pour les éléments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des éléments intrinsèques
C----------------------------------
C----------------------------------
C--------------
C--------------
C-----Elément 1
C--------------
C--------------
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 1 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 1 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_2 - X_MULTI_1
C     $ * ( DLOG ( X_D_R ( 1 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter le DP ( 1 , I_R ) (par E_SC = 0)
C----- = fausse variable 2).
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 1 , I_R )
         END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
              if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
               write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
              end if
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( I_R )
           END IF
        END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
           DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 1 ) THEN
             J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M )
            END IF
           END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 sont requises
             J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C-----------------------  
C-----------------------
C-----Elément 3 éventuel
C-----------------------
C-----------------------
      IF ( N_TYP_INTR_M .GT. 2 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 3 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 3 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_2 - X_MULTI_3
C     $ * ( DLOG ( X_D_R ( 3 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
        IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 3 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 2 ) THEN
             J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
            END IF
           END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
            END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 + N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR 
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : élément 3 éventuel
C------------------------------------
C------------------------------------
       END IF
C----------------------------------------------------------
C----------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_2 > 1
C-----(car la composante Fi pour les lacunes en N_R_1+N_R_2
C-----sert à éliminer le multiplicateur de Lagrange de 2)
C----------------------------------------------------------
C----------------------------------------------------------
       IF ( I_R .LT. N_R_1_M + N_R_2_M ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 0 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 0 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_2
C     $ * ( DLOG ( X_D_R ( 0 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (e.g. type=2 sur s-r "2").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C-------------------------
C-------------------------
C-----Fin du cas : lacunes
C-------------------------
C-------------------------
       END IF
C------------------------------------------
C------------------------------------------
C-----Cas des éléments d'addition éventuels
C------------------------------------------
C------------------------------------------
C-------------------------------------
C-------------------------------------
C-----Boucle sur les types d'additions
C-------------------------------------
C-------------------------------------
       DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
C       WRITE(*,*) "s-r 2 : ADDITIONS ET DERIVEES"
C-----------------------------------------------------
C-----Indicateur d'existence d'une nouvelle composante
C-----de la fonction vectorielle étudiée
C-----------------------------------------------------
        I_NOUV_COMP = 0
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $ .AND. ( N_R_2_M .NE. 0 ) ) THEN
           IF ( I_R .LT. N_R_1_M + N_R_2_M ) THEN 
              DELTA_2 = 0.D0
              I_NOUV_COMP = 1
           END IF
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
           DELTA_2 = X_MULTI_2 - X_MULTI_3
           I_NOUV_COMP = 1
        ELSE 
           DELTA_2 = X_MULTI_2
           I_NOUV_COMP = 1
        END IF
C------------------------------------------------------------------
C-----Si nouvelle composante, calcul de celle-ci et de ses dérivées
C------------------------------------------------------------------
        IF ( I_NOUV_COMP .EQ. 1 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( I_TYP , I_R )
     $ - H_GC_D_R_M ( I_TYP , N_R_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( I_TYP , I_R )
     $          / Z_TYP_R ( I_R ) )
     $   - DLOG ( X_D_R ( I_TYP , N_R_M )
     $          / Z_TYP_R ( N_R_M ) ) )
     $ + DELTA_2
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 2 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( I_TYP , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
          J_COUR = IND_D_R_TYP_M ( I_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( I_TYP , N_R_M )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_2"
C-----On ajoute ci-dessous les dérivées du terme DELTA_2,
C-----en distinguant les différents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas où DELTA_2 = 0 => dérivée nulle
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $ .AND. ( N_R_2_M .NE. 0 ) ) THEN
           IF ( I_R .LT. N_R_1_M + N_R_2_M ) THEN 
           END IF
C-----Cas où DELTA_2 = X_MULTI_2 - X_MULTI_3 = mu_B - mu_C
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
C- - - Dérivées de mu_B
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1+N_R_2
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C- - - Dérivées de - mu_C (ici N_R_M = N_R_1_M + N_R_2_M + N_R_3_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1 + N_R_2 + N_R_3
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas où DELTA_2 = X_MULTI_2 = + mu_B
        ELSE 
C- - - Dérivées de mu_B
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1+N_R_2
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C-----Fin des divers cas
        END IF
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C---------------------------------------------------------------------
C-----Fin du test "existence d'une nouvelle composante de la fonction"
C---------------------------------------------------------------------
        END IF
C-----------------------------------------------
C-----------------------------------------------
C-----Fin de la boucle sur les types d'additions
C-----------------------------------------------
C-----------------------------------------------
       END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 2
C=========================================
      END DO
C===============================================
C=====Composantes de la fonction vectorielle
C=====(éventuelles) relatives aux sous-réseaux 3
C===============================================
      DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-réseau 3, il y a une équation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. système d'équations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les éléments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 1 et 2, ainsi que pour les lacunes 
C-----(sauf pour le sous-réseau N_R_1 + N_R_2 + N_R_3)
C-----et pour les éléments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des éléments intrinsèques
C----------------------------------
C----------------------------------
C--------------  
C--------------  
C-----Elément 1
C--------------
C--------------  
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 1 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 1 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_3 - X_MULTI_1
C     $ * ( DLOG ( X_D_R ( 1 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter le DP ( 1 , I_R ) (par E_SC = 0)
C----- = fausse variable 2).
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 1 , I_R )
         END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 N_R_3 sont requises
             J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C--------------  
C--------------  
C-----Elément 2
C--------------
C--------------  
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 2 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 2 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_3 - X_MULTI_2
C     $ * ( DLOG ( X_D_R ( 2 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
        IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 2 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 + N_R_2 + N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
          J_COUR
     $  = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C---------------------------
C---------------------------
C-----Fin du cas : élément 2
C---------------------------
C---------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_3 > 1
C-----(car la composante Fi pour les lacunes en N_R_1+N_R_2+N_R_3
C-----sert à éliminer le multiplicateur de Lagrange de 3)
C----------------------------------------------------------------
C----------------------------------------------------------------
       IF ( I_R .LT. N_R_1_M + N_R_2_M + N_R_3_M ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 0 , I_R )
C     $ - H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * DLOG ( X_D_R ( 0 , I_R )
     $        / Z_TYP_R ( I_R ) )
     $ + X_MULTI_3
C     $ * ( DLOG ( X_D_R ( 0 , I_R )
C     $          / Z_TYP_R ( I_R ) ) )
C     $   - DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , I_R )
         END IF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé) induites par z^I_R :
C-----attention à bien voir de quelles variables dépend z^I_R
C-----(éviter les fausses variables 1) ET 2)) !
C-----Le test initialement utilisé "IF ( J_TYP .NE. 1 )"
C-----n'écarte que les fausses variables 1) (e.g. type=2 sur s-r "2").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2+N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
            J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C-------------------------
C-------------------------
C-----Fin du cas : lacunes
C-------------------------
C-------------------------
       END IF
C------------------------------------------
C------------------------------------------
C-----Cas des éléments d'addition éventuels
C------------------------------------------
C------------------------------------------
C-------------------------------------
C-------------------------------------
C-----Boucle sur les types d'additions
C-------------------------------------
C-------------------------------------
       DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
C       WRITE(*,*) "s-r 3 : ADDITIONS ET DERIVEES"
C-----------------------------------------------------
C-----Indicateur d'existence d'une nouvelle composante
C-----de la fonction vectorielle étudiée
C-----------------------------------------------------
        I_NOUV_COMP = 0
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $ .AND. ( N_R_3_M .NE. 0 ) ) THEN
          IF ( I_R .LT. N_R_1_M + N_R_2_M + N_R_3_M ) THEN 
               DELTA_3 = 0.D0
               I_NOUV_COMP = 1
           END IF
        ELSE 
           DELTA_3 = X_MULTI_3
           I_NOUV_COMP = 1
        END IF
C------------------------------------------------------------------
C-----Si nouvelle composante, calcul de celle-ci et de ses dérivées
C------------------------------------------------------------------
        IF ( I_NOUV_COMP .EQ. 1 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( I_TYP , I_R )
     $ - H_GC_D_R_M ( I_TYP , N_R_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( I_TYP , I_R )
     $          / Z_TYP_R ( I_R ) )
     $   - DLOG ( X_D_R ( I_TYP , N_R_M )
     $          / Z_TYP_R ( N_R_M ) ) )
     $ + DELTA_3
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 3 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( I_TYP , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
          J_COUR = IND_D_R_TYP_M ( I_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / X_D_R ( I_TYP , N_R_M )
           END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_3"
C-----On ajoute ci-dessous les dérivées du terme DELTA_3,
C-----en distinguant les différents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas où DELTA_3 = 0 => dérivée nulle
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $ .AND. ( N_R_3_M .NE. 0 ) ) THEN
          IF ( I_R .LT. N_R_1_M + N_R_2_M + N_R_3_M ) THEN 
           END IF
C-----Cas où DELTA_3 = X_MULTI_3 = + mu_C
        ELSE 
C- - - Dérivées de mu_C
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2+N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrinsèque du sous-réseau N_R_1+N_R_2+N_R_3
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M
     $             ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C-----Fin des divers cas
        END IF
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C---------------------------------------------------------------------
C-----Fin du test "existence d'une nouvelle composante de la fonction"
C---------------------------------------------------------------------
        END IF
C-----------------------------------------------
C-----------------------------------------------
C-----Fin de la boucle sur les types d'additions
C-----------------------------------------------
C-----------------------------------------------
        END DO
C=========================================
C=====Fin de boucle sur les sous-réseaux 3
C=========================================
      END DO
C===========================================================
C=====Composantes de la fonction vectorielle
C=====(éventuelles) relatives aux sous-réseaux interstitiels
C===========================================================
C      write(*,*) "N_R_1_M = " , N_R_1_M 
C      write(*,*) "N_R_2_M = " , N_R_2_M 
C      write(*,*) "N_R_3_M = " , N_R_3_M 
C      write(*,*) "N_R_M = " , N_R_M 
C	 stop
      DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M
C--------------------------------------------------------------
C--------------------------------------------------------------
C-----Pour chaque sous-réseau interstitiel, il y a une équation
C-----pour chaque espèce 1, 2, 3,
C-----ainsi que pour chaque élément d'addition
C-----(sauf pour le sous-réseau N_R)
C--------------------------------------------------------------
C--------------------------------------------------------------
C-------------------------------------------------------
C-----Remarque importante : pour les s-r interstitiels,
C-----on omet les DP "non utiles",
C-----i.e. ceux pour lesquels "E_DP >= 0" :
C-----* pas de composante de F_NPT pour ces DP
C-----* pas de variable NPT pour ces DP
C-------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des éléments intrinsèques
C----------------------------------
C----------------------------------
C--------------  
C--------------  
C-----Elément 1
C--------------
C--------------  
C------------------------------------------------
C=====Test "E_DP < 0" pour la composante de F_NPT
C------------------------------------------------
       IF ( IND_D_R_TYP_M ( 1 , I_R ) .NE. 0 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 1 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( 1 , I_R )
     $          / Z_TYP_R ( I_R ) ) )
     $ - X_MULTI_1
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M )
C     $          / Z_TYP_R ( N_R_1_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrinsèques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
         IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 1 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
            END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
       END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
      END IF
C-----------------------
C-----------------------
C-----Elément 2 éventuel
C-----------------------
C-----------------------
      IF ( N_TYP_INTR_M .GT. 1 ) THEN
C------------------------------------------------
C=====Test "E_DP < 0" pour la composante de F_NPT
C------------------------------------------------
       IF ( IND_D_R_TYP_M ( 2 , I_R ) .NE. 0 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 2 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( 2 , I_R )
     $          / Z_TYP_R ( I_R ) ) )
     $ - X_MULTI_2
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrinsèques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
         J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
         IF ( J_COUR .NE. 0 ) THEN
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 2 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
             IF ( J_COUR .NE. 0 ) THEN
              if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
               write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
              end if
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( I_R )
            END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
        END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : élément 2 éventuel
C------------------------------------
C------------------------------------
      END IF
C-----------------------
C-----------------------
C-----Elément 3 éventuel
C-----------------------
C-----------------------
      IF ( N_TYP_INTR_M .GT. 2 ) THEN
C------------------------------------------------
C=====Test "E_DP < 0" pour la composante de F_NPT
C------------------------------------------------
       IF ( IND_D_R_TYP_M ( 3 , I_R ) .NE. 0 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( 3 , I_R )
C     $ + H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( 3 , I_R )
     $          / Z_TYP_R ( I_R ) ) )
     $ - X_MULTI_3
C     $   + DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C     $          / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succès - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrinsèques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 3 , I_R )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2+N_R_3 sont requises
C-----pour l'élimination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j à i fixé)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
            IF ( J_COUR .NE. 0 ) THEN
             if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
             end if
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
        END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : élément 3 éventuel
C------------------------------------
C------------------------------------
      END IF
C------------------------------------------
C------------------------------------------
C-----Cas des éléments d'addition éventuels
C------------------------------------------
C------------------------------------------
       DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
         IF ( I_R .LT. N_R_M ) THEN
C------------------------------------------------
C=====Test "E_DP < 0" pour la composante de F_NPT
C------------------------------------------------
         IF ( IND_D_R_TYP_M ( I_TYP , I_R ) .NE. 0 ) THEN
         I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_GC_D_R_M ( I_TYP , I_R )
     $ - H_GC_D_R_M ( I_TYP , N_R_M )
     $ + K_T_M
     $ * ( DLOG ( X_D_R ( I_TYP , I_R )
     $          / Z_TYP_R ( I_R ) )
     $   - DLOG ( X_D_R ( I_TYP , N_R_M )
     $          / Z_TYP_R ( N_R_M ) ) )
C-----------------------------------------------------------
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour éviter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE A TESTER (par comparaison avec muVT)
C-----SUR FeAl-B ou Fe3Al-C ou TiAl-O
C-----(i.e. des cas avec additions interstitielles)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITIONS INTERSTITIELLES
C-----(FeAl-B ou Fe3Al-C ou TiAl-O)
C-----(comme déjà fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
           J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( I_TYP , I_R )
          END IF
          DO J_TYP = 1 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / Z_TYP_R ( I_R )
           END IF
          END DO
          J_COUR = IND_D_R_TYP_M ( I_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
          IF ( J_COUR .NE. 0 ) THEN
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( I_TYP , N_R_M )
          END IF
          DO J_TYP = 1 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajouté ici,
C-----afin de pouvoir écarter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul éventuel
C-----des dérivées partielles de cette composante
C------------------------------------------------
         END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
        END IF
C-----------------------------------------------
C-----------------------------------------------
C-----Fin du cas : éléments d'addition éventuels
C-----------------------------------------------
C-----------------------------------------------
        END IF
	 END DO
C=====================================================
C=====Fin de boucle sur les sous-réseaux interstitiels
C=====================================================
      END DO
C	stop
C###################################################
C#####N_TYP_M composantes de la fonction vectorielle
C#####relatives aux quantités de matière
C###################################################
C-------------------------
C-----Espèce intrinsèque 1
C-------------------------
       I_COMP = I_COMP + 1
       F_NPT_M ( I_COMP ) = 0.D0
       DO I_R = 1 , N_R_1_M  
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + Z_TYP_R ( I_R ) * DFLOAT ( P_R_M ( I_R ) )
C	write ( * , * ) 'i_r = ' , I_R
C	write ( * , * ) 'Z_TYP_R ( I_R ) = ' ,
C     $			Z_TYP_R ( I_R )
C	write ( * , * ) 'F_NPT_M ( I_COMP ) = ' , 
C     $  		F_NPT_M ( I_COMP )
       END DO
       DO I_R =  N_R_1_M + 1 , N_R_M
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( 1 , I_R ) * DFLOAT ( P_R_M ( I_R ) )
C       write ( * , * ) 'i_r = ' , I_R
C      write ( * , * ) 'X_D_R ( 1 , I_R ) = ',
C    $                  X_D_R ( 1 , I_R )
C       write ( * , * ) 'F_NPT_M ( I_COMP ) = ' , 
C    $                   F_NPT_M ( I_COMP )
       END DO
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP ) * M_MAILLES - N_AT ( 1 )
C	write ( * , * ) 'M_MAILLES = ' , M_MAILLES
C	write ( * , * ) 'N_AT ( 1 ) = ' , N_AT ( 1 )
C	write ( * , * ) 'F_NPT_M ( I_COMP ) = ' ,
C    $                 F_NPT_M ( I_COMP )
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
       DO I_R = 1 , N_R_1_M  
          DO J_TYP = 0 , N_TYP_M
           IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
           END IF
          END DO
       END DO
       DO I_R = N_R_1_M + 1 , N_R_M
            J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C-----Test "E_DP < 0" pour la variable "x_DP sur s-r interstitiel"
          IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
          END IF
       END DO
C-----Dérivée par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 1 ) ) / M_MAILLES
      END IF
C------------------------------------
C-----Espèce intrinsèque 2 éventuelle
C------------------------------------
      IF ( N_TYP_INTR_M .GT. 1 ) THEN
       I_COMP = I_COMP + 1
       F_NPT_M ( I_COMP ) = 0.D0
       DO I_R = 1 , N_R_1_M  
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( 2 , I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
       DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + Z_TYP_R ( I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
       DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_M  
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( 2 , I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP ) * M_MAILLES - N_AT ( 2 )
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
       DO I_R = 1 , N_R_1_M
            J_COUR = IND_D_R_TYP_M ( 2 , I_R )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
       END DO
       DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
          DO J_TYP = 0 , N_TYP_M
           IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
           END IF
          END DO
       END DO
       DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_M
            J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C-----Test "E_DP < 0" pour la variable "x_DP sur s-r interstitiel"
          IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
          END IF
       END DO
C-----Dérivée par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 2 ) ) / M_MAILLES
      END IF
C--------------------------------------
C-----Fin du cas : espèce intrinsèque 2
C--------------------------------------
      END IF
C------------------------------------
C-----Espèce intrinsèque 3 éventuelle
C------------------------------------
      IF ( N_TYP_INTR_M .GT. 2 ) THEN
       I_COMP = I_COMP + 1
       F_NPT_M ( I_COMP ) = 0.D0
       DO I_R = 1 , N_R_1_M + N_R_2_M 
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( 3 , I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
       DO I_R =  N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + Z_TYP_R ( I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
       DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M  
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( 3 , I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP ) * M_MAILLES - N_AT ( 3 )
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
       DO I_R = 1 , N_R_1_M + N_R_2_M
            J_COUR = IND_D_R_TYP_M ( 3 , I_R )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
       END DO
       DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
          DO J_TYP = 0 , N_TYP_M
           IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
           END IF
          END DO
       END DO
       DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M
            J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C-----Test "E_DP < 0" pour la variable "x_DP sur s-r interstitiel"
          IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
          END IF
       END DO
C-----Dérivée par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 3 ) ) / M_MAILLES
      END IF
C--------------------------------------
C-----Fin du cas : espèce intrinsèque 3
C--------------------------------------
      END IF
C----------------------------------
C-----Eléments d'addition éventuels
C----------------------------------
      DO I_TYP = N_TYP_INTR_M + 1 , N_TYP_M
       I_COMP = I_COMP + 1
       DO I_R = 1 , N_R_M 
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP )
     $ + X_D_R ( I_TYP , I_R ) * DFLOAT ( P_R_M ( I_R ) )
       END DO
         F_NPT_M ( I_COMP )
     $ = F_NPT_M ( I_COMP ) * M_MAILLES - N_AT ( I_TYP )
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
         DO I_R = 1 , N_R_M 
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C-----Test "E_DP < 0" pour la variable "x_DP sur s-r interstitiel"
          IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + DFLOAT ( P_R_M ( I_R ) ) * M_MAILLES
          END IF
         END DO
        END IF
C-----Dérivée par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( I_TYP ) ) / M_MAILLES
C-------------------------------------------------
C-----Fin de la boucle sur les éléments d'addition
C-------------------------------------------------
      END DO
C################################################### 
C#####Dernière composante de la fonction vectorielle
C#####(h0 + ..)
C################################################### 
C------------------------------------
C-----Terme auxiliaire "Somme_Log(z)"
C------------------------------------
      SOMME_LOGZ = 0.D0
      DO I_R = 1 , N_R_M
        SOMME_LOGZ = SOMME_LOGZ
     $             + DFLOAT ( P_R_M ( I_R ) )
     $             * DLOG ( Z_TYP_R ( I_R ) )
      END DO
C--------------------------------------------------
C-----Termes R_1, R_1+R_2 et R_1+R_2+R_3 de (2.104)
C--------------------------------------------------
       W1 = H_GC_D_R_M ( 0 , N_R_1_M )
     $    + K_T_M * DLOG ( X_D_R ( 0 , N_R_1_M )
     $                   / Z_TYP_R ( N_R_1_M ) )
	 W1 = W1 * DFLOAT ( P_R_1 )
	 W2 = 0.D0
       IF ( N_TYP_INTR_M .GT. 1 ) THEN
        W2 = H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M )
     $     + K_T_M * DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M )
     $                    / Z_TYP_R ( N_R_1_M + N_R_2_M ) )
	 END IF
	 W2 = W2 * DFLOAT ( P_R_2 )
	 W3 = 0.D0
       IF ( N_TYP_INTR_M .GT. 2 ) THEN
        W3 = H_GC_D_R_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $      + K_T_M
     $      * DLOG ( X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $             / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M ) )
	  W3 = W3 * DFLOAT ( P_R_3 )
	 END IF
C==========================================
C=====Composante de la fonction vectorielle
C=====Il s'agit de (2.104) du mémoire d'HDR.
C==========================================
      I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_REF_MAILLE_M + W1 + W2 + W3 + K_T_M * SOMME_LOGZ
C-----Dérivées partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C-----Terme W1
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M ) * DFLOAT ( P_R_1 )
           DO J_TYP = 0 , N_TYP_M
            IF ( J_TYP .NE. 1 ) THEN
             J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_1_M ) * DFLOAT ( P_R_1 )
            END IF
           END DO
C-----Terme W2
           IF ( N_TYP_INTR_M .GT. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M ) * DFLOAT ( P_R_2 )
            DO J_TYP = 0 , N_TYP_M
             IF ( J_TYP .NE. 2 ) THEN
              J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M ) * DFLOAT ( P_R_2 )
             END IF
            END DO
           END IF
C-----Terme W3
           IF ( N_TYP_INTR_M .GT. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
     $    * DFLOAT ( P_R_3 )
            DO J_TYP = 0 , N_TYP_M
             IF ( J_TYP .NE. 3 ) THEN
              J_COUR
     $     = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M + N_R_3_M )
     $      * DFLOAT ( P_R_3 )
             END IF
            END DO
           END IF
C-----Dérivées du terme "Somme_Log(z)" : sous-réseaux intrinsèques 1
           DO I_R = 1 , N_R_1_M
            DO J_TYP = 0 , N_TYP_M
             IF ( J_TYP .NE. 1 ) THEN
              J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M * DFLOAT ( P_R_M ( I_R ) ) / Z_TYP_R ( I_R )
             END IF
            END DO
           END DO
C-----Dérivées du terme "Somme_Log(z)" : sous-réseaux intrinsèques 2 éventuels :
C-----test ci-dessous sans doute requis, car N_R_2 et N_R_3
C-----ne sont pas initialisés à 0 par défaut dans le programme principal
C----- -> A VERIFIER
           IF ( N_TYP_INTR_M .GT. 1 ) THEN
            DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
             DO J_TYP = 0 , N_TYP_M
              IF ( J_TYP .NE. 2 ) THEN
              J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M * DFLOAT ( P_R_M ( I_R ) ) / Z_TYP_R ( I_R )
             END IF
            END DO
           END DO
          END IF
C-----Dérivées du terme "Somme_Log(z)" : sous-réseaux intrinsèques 3 éventuels :
C-----test ci-dessous sans doute requis, car N_R_2 et N_R_3
C-----ne sont pas initialisés à 0 par défaut dans le programme principal
C----- -> A VERIFIER
           IF ( N_TYP_INTR_M .GT. 2 ) THEN
            DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
             DO J_TYP = 0 , N_TYP_M
              IF ( J_TYP .NE. 3 ) THEN
              J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M * DFLOAT ( P_R_M ( I_R ) ) / Z_TYP_R ( I_R )
             END IF
            END DO
           END DO
          END IF
C-----Dérivées du terme "Somme_Log(z)" : sous-réseaux interstitiels éventuels
          DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M
           DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C-----Test "E_DP < 0" pour la variable "x_DP sur s-r interstitiel"
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      - K_T_M * DFLOAT ( P_R_M ( I_R ) ) / Z_TYP_R ( I_R )
            END IF
           END DO
          END DO
C-----Fin du test "calcul des dérivées partielles" pour cette composante
        END IF
C-----------------
C-----Vérification
C-----------------
      IF ( I_COMP .NE. N_NPT_M ) THEN
	write ( * , * ) '-----------------------------------------'
	write ( * , * ) 'En fin de f_NR :'
	write ( * , * ) "nombre d'équations = " , I_COMP
	write ( * , * ) 'différent du nombre attendu = ' , N_NPT_M
	write ( * , * ) '-----------------------------------------'
	CALL INTERRUPTION
      END IF
C-------------
C-----Ecriture
C-------------
      if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
 	write ( * , * ) 'En fin de procédure F_NR :' 
      do i_npt = 1 , n_npt_m
 	write ( * , * ) 'f_NR (' , i_npt , ') = ' ,
     $		F_NPT_M ( i_npt ) 
      end do
      write (* , * ) '-------------------'
       end if
C	stop
      RETURN
      END
C     ==================================================================
C     I                                                                I
C     I Calcul des concentrations en porteurs intrinsèques n et        I
C     I et de leur dérivée                                             I
C     I à partir de la densité d'états                                 I
C     I en fonction du potentiel chimique électronique                 I
C     I								       I
C     ==================================================================
C     ==================================================================
C     I Dernière mise à jour : 11/07/2013                              I
C     I Saint Benoît                                                   I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   CALC_N_P
     $ (
     $   N_PAS_M , TAB_DDE_M ,
     $   E_MIN_M , E_MAX_M ,
     $   E_MAX_BV_M , E_MIN_BC_M ,
     $   POT_CHIM_ELEC_M  , TEMP_M ,
     $   C_VOL_N_M , C_VOL_P_M ,
     $   D_C_VOL_N_M , D_C_VOL_P_M )
C     ##################################################################
	 USE CONSTANTES
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####Déclaration des tableaux arguments de la procédure
C#######################################################
C--------------------------------------------------
C-----Pour chaque point, valeur de l'énergie et DdE
C--------------------------------------------------
      REAL * 8 TAB_DDE_M ( N_PAS_M , 2 )
C#####################################################
C#####Déclaration des tableaux internes à la procédure
C#####################################################
C---------------------------------
C-----Facteur thermodynamique kB*T
C---------------------------------
      REAL * 8 K_T
C--------------------
C-----Initialisations
C--------------------
      K_T = K_B * TEMP_M
      C_VOL_N_M = 0.D0
      C_VOL_P_M = 0.D0
      D_C_VOL_N_M = 0.D0
      D_C_VOL_P_M = 0.D0
      DO I_PAS = 1 , N_PAS_M - 1
        E_COUR = TAB_DDE_M ( I_PAS , 1 )
        E_SUIV = TAB_DDE_M ( I_PAS + 1 , 1 )
        DDE_COUR = TAB_DDE_M ( I_PAS , 2 )
        FACT_COUR = DEXP ( ( E_COUR - POT_CHIM_ELEC_M ) / K_T )
C-----------------
C-----Terme pour n
C-----------------
        IF ( E_COUR .GE. E_MIN_BC_M ) THEN
          C_VOL_N_M = C_VOL_N_M 
     $   + DDE_COUR * ( E_SUIV - E_COUR )
     $   / ( 1.D0 + FACT_COUR )
        END IF
C-----------------
C-----Terme pour p
C-----------------
        IF ( E_COUR .LE. E_MAX_BV_M ) THEN
          C_VOL_P_M = C_VOL_P_M 
     $   + DDE_COUR * ( E_SUIV - E_COUR )
     $   * ( 1.D0 - 1.D0 / ( 1.D0 + FACT_COUR ) )
        END IF
C-------------------------------
C-----Terme pour la dérivée de n
C-------------------------------
        IF ( E_COUR .GE. E_MIN_BC_M ) THEN
          D_C_VOL_N_M = D_C_VOL_N_M 
     $   + DDE_COUR * ( E_SUIV - E_COUR )
     $   * FACT_COUR / ( 1.D0 + FACT_COUR ) * * 2
     $   / K_T
        END IF
C-------------------------------
C-----Terme pour la dérivée de p
C-------------------------------
        IF ( E_COUR .LE. E_MAX_BV_M ) THEN
          D_C_VOL_P_M = D_C_VOL_P_M 
     $   - DDE_COUR * ( E_SUIV - E_COUR )
     $   * FACT_COUR / ( 1.D0 + FACT_COUR ) * * 2
     $   / K_T
        END IF
      END DO
      RETURN
      END

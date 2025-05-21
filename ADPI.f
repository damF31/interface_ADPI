C     ==================================================================
C     I     Approximation des D�fauts Ponctuels Ind�pendants (ADPI)    I
C     I             pour la description thermodynamique                I
C     I                 d'un alliage/compos� ordonn�                   I
C     I             (compos� interm�tallique, min�ral,...)             I
C     I        avec un, deux ou trois �l�ments intrins�ques            I
C     I                  et des nombres quelconques                    I
C     I          de sous-r�seaux et d'�l�ments d'addition              I
C     ==================================================================
C     I                                                                I
C     I Grandeurs physiques calcul�es :                                I
C     I -----------------------------                                  I
C     I * quantit�s de d�fauts ponctuels                               I
C     I * quantit�s thermodynamiques                                   I
C     I (�nergie libre, entropie, potentiels chimiques)                I
C     I en fonction des variables thermodynamiques intensives :        I
C     I composition = �cart � la stoechiom�trie, temp�rature           I
C     I (et �ventuellement pression)                                   I
C     I                                                                I
C     I Les r1 premiers sous-r�seaux contiennent l'esp�ce 1            I
C     I les r2 suivants l'esp�ce 2                                     I
C     I les r3 suivants l'esp�ce 3                                     I
C     I dans l'�tat fondamental sans d�faut                            I
C     I                                                                I
C     I Aucun sous-r�seau ne contient les �l�ments d'addition          I
C     I (i>3) � l'�tat fondamental                                     I
C     I                                                                I
C     I Les sous-r�seaux d'indices > r1 + r2 + r3 sont interstitiels   I
C     I                                                                I
C     I On note p(r) le nombre de sites du sous-r�seau r par maille    I
C     I                                                                I
C     I                 ============================                   I
C     I                 Trois modes de calcul ADPI :                   I
C     I                 ============================                   I
C     I                                                                I
C     I --------------------------------------------------             I
C     I Mode muVT --> ensemble grand canonique (mu(i),V,T)             I
C     I --------------------------------------------------             I
C     I                                                                I
C     I Relation de Gibbs-Duhem approch�e pour la pression :           I
C     I � cause de cette relation P(mu(i),T) approch�e,                I
C     I le balayage en potentiels chimiques r�sultant induit           I
C     I  une (l�g�re) d�rive en pression par rapport � la consigne.    I
C     I -> ajout d'une possibilit� de boucle d'autocoh�rence sur mu(1) I
C     I pour �viter cette d�rive en P.                                 I
C     I                                                                I
C     I En mode muVT, le programme ADPI effectue :                     I
C     I A) � partir de la relation de Gibbs-Duhem, le calcul approch�  I
C     I    de la valeur de N(1)*mu(1) + N(2)*mu(2) + N(3)*mu(3)        I
C     I    correspondant � la pression prescrite                       I
C     I                                                                I
C     I    Attention : il s'agit de la pression externe                I
C     I                exerc�e sur le syst�me par l'ext�rieur          I
C     I                                                                I
C     I B) un balayage en �carts de potentiels chimiques               I
C     I d_mu_i = mu(i_r�f) - mu(i) (i = 1, 2 ou 3 <> i_r�f)	       I
C     I             et en mu(i>3)                                      I
C     I								       I
C     I C) pour chaque valeur de d_mu_i et mu(i>3)                     I
C     I       l'�criture dans les fichiers ".adpi" des grandeurs :     I
C     I     (i) composition                                            I
C     I    (ii) quantit�s de d�fauts ponctuels x_d�f, avec d�f =       I
C     I  1(r) pour r1 < r <= r1 + r2 + r3                              I
C     I         ou r > r1 + r2 + r3 (1 interstitiel),                  I
C     I  2(r) pour r <= r1 ou r1 + r2 < r <= r1 + r2 + r3              I
C     I         ou r > r1 + r2 + r3 (2 interstitiel),                  I
C     I  3(r) pour r <= r1 + r2                                        I
C     I         ou r > r1 + r2 + r3 (3 interstitiel),                  I
C     I  L(r) pour r <= r1 + r2 + r3                                   I
C     I  i(r) pour tout r                                              I
C     I   (iii) �nergie, volume, entropie de configuration             I
C     I         et enthalpie libre par atome ou par maille             I
C     I                                                                I
C     I Possibilit� de fen�tres en composition pour l'�criture         I
C     I                                                                I
C     I Il a �t� ajout� ult�rieurement                                 I
C     I la possibilit� d'une boucle d'autocoh�rence sur mu(1)          I
C     I pour compenser la d�rive en P                                  I
C     I induite par la relation approch�e P(mu(i),T).                  I
C     I                                                                I
C     I -------------------------------------                          I
C     I Mode NPT=0 --> ensemble (N(i),P,T=0K)                          I
C     I -------------------------------------                          I
C     I                                                                I
C     I Minimisation de H (m�thode du simplexe)                        I
C     I par rapport aux quantit�s de DP                                I
C     I sous les contraintes de quantit�s de mati�re constantes        I
C     I                                                                I
C     I Variantes disponibles pour le mode NPT=0 : P/p, x              I
C     I                                                                I
C     I --------------------------------                               I
C     I Mode NPT --> ensemble (N(i),P,T)                               I
C     I --------------------------------                               I
C     I                                                                I
C     I R�solution par la m�thode de Newton-Raphson (NR)               I
C     I du syst�me non lin�aire d'inconnues ( M , x_d )                I
C     I correspondant � la minimisation de l'enthalpie libre           I
C     I                                                                I
C     I Variantes disponibles pour le mode NPT : P/p, T, x             I
C     I                                                                I
C     I Remarques :                                                    I
C     I un calcul NPT peut diverger si les param�tres sont mal choisis I
C     I => penser � jouer sur les param�tres suivants :                I
C     I * ordre (croissant / d�croissant) du balayage en T ou x        I
C     I * finesse du balayage (nombre de points)                       I
C     I * valeurs initiales des inconnues                              I
C     I * fr�quence de calcul de la matrice jacobienne                 I
C     I et �galement (plus rare)                                       I
C     I * nombre maximal d'it�rations de NR                            I
C     I * pr�cision pour l'arr�t de l'algorithme de NR                 I
C     I                                                                I
C     I ====================================================	       I
C     I          Possibilit�s de prise en compte :                     I
C     I * de DP complexes (approximative, en muVT seulement)           I
C     I * de DP charg�s (en muVT"point" seulement)                     I
C     I ====================================================	       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 24/05/2022                              I
C     I Mardi des Rogations                                            I
C     ==================================================================
      PROGRAM ADPI
      USE CONSTANTES
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C############################
C############################
C#####D�clarations et formats
C############################
C############################
C###################################
C#####Tableaux de dimension variable
C###################################
C--------------------------------------------------------
C-----Commentaires �crits dans le fichier de compte-rendu
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
C-----aux compositions � �crire dans les fichiers de r�sultats
C-----et demi-largeurs de fen�tres correspondantes
C-----(ainsi qu'un tableau de valeurs initiales
C-----pour la commodit� de lecture en raison des indices
C----- - absence de I_TYP_0)
C-------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          COEF_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_X_AT_INIT
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_X_AT
C---------------------------------------------------
C-----Matrices de d�composition LU pour le calcul
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
C-----� �crire dans les fichiers de r�sultats
C-----------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_0_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_INF_CTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          X_AT_SUP_CTR
C---------------------------------------------------
C-----Potentiels chimiques des �l�ments intrins�ques
C---------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_INTR
C-------------------------------------------------
C-----Potentiels chimiques des �l�ments d'addition
C-------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I
C--------------------------------------------
C-----Termes relatifs aux �l�ments d'addition
C--------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          ALPHA_I_R
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          SOMME_I
C------------------------------------------------------------------
C-----Valeurs initiales, nombres de pas et incr�ments
C-----des �carts de potentiels chimiques
C-----pour les �l�ments intrins�ques autres que celui de r�f�rence
C-----------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          D_POT_REF_INTR_INIT
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_D_POT_REF_INTR
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          PAS_D_POT_REF_INTR
C----------------------------------------------------
C-----Valeurs initiales, nombres de pas et incr�ments
C-----de potentiels chimiques des �l�ments d'addition
C----------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          POT_I_INIT
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_POT_I
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          PAS_POT_I
C---------------------------------------------------------
C-----Potentiels chimiques moyens des �l�ments d'addition,
C-----leurs carr�s et les variances correspondantes
C-----pour les valeurs �crites
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
C-----(utile pour la simulation de plusieurs boucles imbriqu�es
C-----en potentiels chimiques avec un seul indice)
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: 
     $                                          PROD_PART_N_POT_I
C------------------------------------------------
C-----Fractions d'antisites et lacunes (indice 0)
C-----dans les divers sous-r�seaux
C------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_D_R
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-r�seau
C-------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          P_R
C-----------------------------------------------------------------
C-----Nombre de sous-r�seaux pour chaque esp�ce intrins�que
C-----(utile pour le cas d'une extension � N esp�ces intrins�ques)
C-----------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          N_R_N
C-----------------------------------
C-----Nombres cumul�s correspondants
C-----------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          RHO_N
C====================================================================
C====================================================================
C=====Param�tres GC pour des d�fauts simples neutres (INDIC_CHARGE=0)
C====================================================================
C====================================================================
C--------------------------------------------------------------
C-----Energies "brutes" des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          E_B_D_R
C--------------------------------------------------------
C-----Energies GC des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          E_GC_D_R
C----------------------------------
C-----Enthalpies GC correspondantes
C----------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          H_GC_D_R
C-------------------------------------
C-----M�mes quantit�s pour les volumes
C-------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          V_B_D_R
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          V_GC_D_R
C--------------------------------------------------------------------
C-----Enthalpies de formation des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          H_FORM_D_R
C--------------------------------------
C-----Termes compl�mentaires d'entropie
C--------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          Z_TYP_R
C================================================
C================================================
C=====Param�tres pour des d�fauts simples charg�s
C================================================
C================================================
C--------------------------------------------------
C-----Densit� d'�tats �lectroniques :
C-----pour chaque point, valeur de l'�nergie et DdE
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          TAB_DDE
C-------------------------------------------------------
C-----Tableau des nombres d'�tats de charges pour les DP
C-------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          NQ_D_R_Q
C---------------------------------------------------------------
C-----Charges (enti�res) des d�fauts sur les divers sous-r�seaux
C---------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          Q_D_R_Q
C--------------------------------------------------------------
C-----Energies "brutes" des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          E_B_D_R_Q
C--------------------------------------------------------
C-----Energies GC des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          E_GC_D_R_Q
C----------------------------------
C-----Enthalpies GC correspondantes
C----------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          H_GC_D_R_Q
C-------------------------------------
C-----M�mes quantit�s pour les volumes
C-------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          V_B_D_R_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          V_GC_D_R_Q
C--------------------------------------------------------------------
C-----Enthalpies de formation des d�fauts sur les divers sous-r�seaux
C--------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : , : ) ::
     $                                          H_FORM_D_R_Q
C===================================================================
C===================================================================
C=====Param�tres relatifs aux d�fauts complexes (pour ADPI sans chg)
C===================================================================
C===================================================================
C-----------------
C-----Multiplicit�
C-----------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                MULTIPLICITE_COMPLEXE
C-----------------------------------------------------------
C-----Sous-r�seau sur lequel est calcul�e cette multiplicit�
C-----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                I_S_R_MULTIPLICITE_COMPLEXE
C--------------------------------
C-----Nombre de sites du complexe
C--------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NOMBRE_SITES_COMPLEXE
C--------------------------------------------------------------
C-----Num�ros de sous-r�seaux des sites occup�s par le complexe
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_S_R_COMPLEXE
C-----------------------------------------------
C-----Types chimiques du complexe sur ces sites,
C-----dans le m�me ordre (0 = lacune)
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
C-----Param�tres indicateurs de types et de sous-r�seaux
C-----utiles au calcul des termes de potentiels chimiques  
C-----associ�s aux complexes
C--------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  ALPHA_TYPE_COMPLEXE
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  BETA_S_R_COMPLEXE
C--------------------------------------------------------------------
C-----Type chimique normal de chaque sous-r�seau
C-----(0 pour sous-r�seaux interstitiels)
C-----utile au calcul des termes de potentiels chimiques de complexes
C--------------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                  I_TYPE_NORMAL_S_R
C----------------------------------------------------------
C-----Nombre de sites du sous-r�seau r dans chaque complexe
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
C-----Enthalpies de formation des d�fauts complexes
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          H_FORM_D_COMPLEXE
C-------------------------------------------
C-----Termes de somme sur les complexes
C-----relatifs aux �l�ments d'addition
C-----dans le calcul des fractions atomiques
C-------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                          SOMME_COMPLEXE_V_TYPE_I
C===================================================================
C===================================================================
C=====Param�tres relatifs aux d�fauts complexes charg�s (suffixe _Q)
C===================================================================
C===================================================================
C---------------------------------------
C-----Indicateur de DP complexes charg�s
C---------------------------------------
       CHARACTER INDIC_COMPLEXES_Q
C--------------------------------------------------
C-----Nombre d'�tats de charge pour chaque complexe
C--------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NQ_COMPLEXE_Q
C--------------------------------------------------
C-----Nombre d'�tats de charge pour chaque complexe
C--------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                Q_COMPLEXE_Q
C-----------------
C-----Multiplicit�
C-----------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                MULTIPLICITE_COMPLEXE_Q
C-----------------------------------------------------------
C-----Sous-r�seau sur lequel est calcul�e cette multiplicit�
C-----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                           I_S_R_MULTIPLICITE_COMPLEXE_Q
C--------------------------------
C-----Nombre de sites du complexe
C--------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) ::
     $                                NOMBRE_SITES_COMPLEXE_Q
C--------------------------------------------------------------
C-----Num�ros de sous-r�seaux des sites occup�s par le complexe
C--------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_S_R_COMPLEXE_Q
C------------------------------------------------------
C-----Types chimiques sur chacun des sites du complexe,
C-----dans le m�me ordre que les s-r (0 = lacune)
C------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                I_TYPE_COMPLEXE_Q
C-------------------------------------------------
C-----Energie et volume "bruts" de chaque complexe
C-----pour chaque �tat de charge
C-------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  E_B_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                              V_B_COMPLEXE_Q
C-------------------------------------------------------
C-----Energie, volume et enthalpie GC de chaque complexe
C-----pour chaque �tat de charge
C-------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  E_GC_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  V_GC_COMPLEXE_Q
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  H_GC_COMPLEXE_Q
C----------------------------------------------------------
C-----Nombre de sites du sous-r�seau r dans chaque complexe
C-----et nombre d'atomes de type i dans ce complexe
C----------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  U_COMPLEXE_S_R_Q
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  V_COMPLEXE_TYPE_Q
C-------------------------------
C-----Fractions de complexes
C-----pour chaque �tat de charge
C-------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                  X_COMPLEXE_Q
C--------------------------------------------------
C-----Enthalpies de formation des d�fauts complexes
C-----pour chaque �tat de charge
C--------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                             H_FORM_COMPLEXE_Q
C============================================================
C============================================================
C=====Cha�nes de caract�res pour suffixes de noms de fichiers
C============================================================
C============================================================
C--------------------------------------------------------
C-----Cha�ne de caract�res pour les suffixes de DP
C-----(d�faut ponctuel et sous-r�seau) et leurs longueurs
C--------------------------------------------------------
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_R
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_TYP
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_R
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_TYP
C--------------------------------------------------------
C-----Cha�ne de caract�res pour les suffixes de complexes
C-----et leurs longueurs
C--------------------------------------------------------
      CHARACTER * 4 , ALLOCATABLE , DIMENSION ( : ) :: W_COMPLEXE
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : ) :: L_W_COMPLEXE
C-----------------------------------------------------------------
C-----Tableau des noms de DP en fonction des types et sous-r�seaux
C-----et longueurs correspondantes
C-----------------------------------------------------------------
      CHARACTER * 10 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          NOM_D_R_TYP
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          LONG_NOM_D_R_TYP
C---------------------------------------------------
C-----Noms des DP en fonction de leur indice unique,
C-----longueurs des noms correspondantes,
C-----sous-r�seaux et types,
C-----�nergies, volumes et enthalpies GC
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
C=====Tableaux relatifs � l'ADPI NPT � T = 0 K
C=============================================
C=============================================
C---------------------------------------------------------------
C-----Indice du DP en fonction de son sous-r�seau et de son type
C---------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          IND_D_R_TYP 
C---------------------------------------------------------------------
C-----Tableau du simplexe contenant :
C-----(i) la fonction H � minimiser en premi�re ligne
C-----(i) les contraintes dans les lignes 2 � N_TYP + 1
C-----(i) la fonction auxiliaire dans la ligne N_TYP + 2
C-----Ce tableau contient N_TYP_D_R + 2 colonnes :
C-----la premi�re colonne contient les valeurs des contraintes,
C-----les autres colonnes contiennent les coefficients des contraintes
C---------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) :: TAB_SMPLX_D_R
C--------------------------------------
C-----Tableaux de r�sultats du simplexe
C--------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_1_SMPLX
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_2_SMPLX
C-----------------------------------------------
C-----En mode NPT=0 "balayage",
C-----tableau de r�sultats du simplexe pr�c�dent
C-----pour d�tection de limite de zone
C-----------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: I_2_SMPLX_PREC
C-------------------------------------------------------------
C-----Tableau des indices de la liste initiale I_2_SMPLX de DP
C-----class�s par ordre croissant
C-------------------------------------------------------------
      INTEGER * 4 , ALLOCATABLE , DIMENSION  ( : ) :: IND_INIT
C--------------------------------------------------------------
C-----Tableau des sous-r�seaux et types des DP constitutionnels
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
C-----Energies de r�f�rence des diverses esp�ces chimiques (eV/at.)
C------------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           E_REF_TYP
C=====================================
C=====================================
C=====Tableaux relatifs � l'ADPI - NPT
C=====================================
C=====================================
C-----------------------------------------
C-----log_10 des fractions initiales de DP
C-----------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           LOG_X_D_R_INIT
C--------------------------------------------------------
C-----Vecteur, fonction vectorielle et matrice jacobienne
C-----du syst�me non lin�aire � r�soudre en NPT
C--------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           X_NPT
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           F_NPT
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : , : ) ::
     $                                           J_NPT
C----------------------------------------
C-----En NPT, vecteur initial � un indice
C----------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION  ( : ) ::
     $                                           X_NPT_INIT
C------------------------------------------------
C-----En NPT, ce m�me vecteur lu dans un fichier,
C-----pour chaque composition / temp�rature
C------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_NPT_INIT_FICH
C------------------------------------------------------------
C-----En NPT(Tx) : tableau pour la conservation des r�sultats 
C-----(aux diverses compo. = boucle interne)
C-----� la temp�rature courante en vue de l'initialisation
C-----� chaque compo. � la temp�rature suivante
C------------------------------------------------------------
      REAL * 8 , ALLOCATABLE , DIMENSION ( : , : ) ::
     $                                          X_NPT_COMPO_SVG
C###############################
C#####Tableaux de dimension fixe
C###############################
C-------------------------------------------------
C-----Cha�ne de lecture dans le fichier de donn�es
C-------------------------------------------------
      CHARACTER * 100 CHAINE_LECT
C---------------------------------------
C-----Base du nom des fichiers de sortie
C---------------------------------------
      CHARACTER * 150 FICH_FIN
C---------------------------------------------------------
C-----Indicateur de pr�sence de sous-r�seaux interstitiels
C---------------------------------------------------------
      CHARACTER * 1 INDIC_R_INTER
C-----------------------------------------------------------------
C-----Indicateur de prise en compte des interstitiels intrins�ques
C-----------------------------------------------------------------
      CHARACTER * 1 INDIC_INTER_INTR
C------------------------------------------------
C-----Indicateur de pr�sence de d�fauts complexes
C------------------------------------------------
      CHARACTER * 1 INDIC_COMPLEXES
C---------------------------------------------------------
C-----Indicateur d'�criture des grandeurs thermodynamiques
C-----par atome (A/a) ou par maille (M/m)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_AT_MAILLE
C---------------------------------------------------------
C-----Indicateur d'�criture de l'�nergie libre par atome :
C-----�nergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_G
C------------------------------------------------------
C-----Indicateur du mode de calcul (muVT, NPT ou NPT=0)
C------------------------------------------------------
      CHARACTER * 10 INDIC_TYP_CALC
C-------------------------------------------------------------------
C-----Indicateur de type de calcul NPT(=0)
C----- * P/p : point (T,x) fix�s (modes NPT et NPT=0)
C----- * T : balayage en temp�rature (mode NPT)
C----- * x : balayage en composition (modes NPT et NPT=0)
C----- * Tx et xT  : doubles balayages en temp�rature et composition
C-----  --> boucle externe sur T ou x pour Tx ou xT respectivement
C----- (mode NPT)
C-------------------------------------------------------------------
      CHARACTER * 2 INDIC_TYP_CALC_NPT
C---------------------------------------------------------------
C-----Produit courant des nombres de pas des �l�ments d'addition
C---------------------------------------------------------------
      INTEGER * 4 PROD_PART_COUR
C-----------------------------------------------------
C-----Caract�res utiles � la constitution d'un format 
C-----pour l'�criture d'un nombre variable de colonnes
C-----(potentiel chimique de chaque esp�ce,
C-----fraction atomique de chaque esp�ce,
C-----fraction de DP sur le sous-r�seau consid�r�
C-----et �nergie de formation de ce DP)
C-----------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_X_DP
C------------------------------------------
C-----Format g�n�rique du titre associ�
C-----et format particulier � chaque d�faut
C------------------------------------------
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP_TYPE
C-----------------------------------------
C-----Idem pour le fichier d'�nergie libre
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
C-----Indicateur d'�criture des seuls points contenus dans la fen�tre
C--------------------------------------------------------------------
      CHARACTER * 1 INDIC_ECRIT_FENETRE
C-------------------------------------------------
C-----Ligne des fractions atomiques dans DATA.adpi
C-------------------------------------------------
      CHARACTER * 1000 LIGNE_X_AT
C--------------------------------------------------------
C-----Tableaux d'analyse de cette ligne cha�ne par cha�ne
C--------------------------------------------------------
      INTEGER * 4 I_DEBUT ( 1000 )
      INTEGER * 4 I_FIN ( 1000 )
      INTEGER * 4 LONG_CHAINE_X_AT ( 1000 )
      CHARACTER * 1000 CHAINE_X_AT ( 1000 )
C----------------------------------------------------
C-----Cha�ne totale de suffixe de fractions atomiques
C----------------------------------------------------
      CHARACTER * 1000 CHAINE_TOTALE_X_AT
C--------------------------------
C-----Facteur thermodynamique k*T
C--------------------------------
      REAL * 8 K_T
C--------------------------------------
C-----Nombre d'atomes par maille (r�el)
C--------------------------------------
      REAL * 8 N_AT_MAILLE
C---------------------------------
C-----Nombre d'atomes total (r�el)
C-----(pour calcul NPT)
C---------------------------------
      REAL * 8 N_AT_TOT
C-------------------------------------
C-----Nombre de mailles initial (r�el)
C-----(pour calcul NPT)
C-------------------------------------
      REAL * 8 N_MAILLES_INIT
C-------------------------------------------------
C-----Fr�quence (nombre de pas) de calcul de J_NPT
C-------------------------------------------------
      INTEGER * 4 P_J_NPT
C----------------------------------------
C-----Fichier de valeurs initiales en NPT
C----------------------------------------
      CHARACTER * 200 FICH_VAL_INIT
C-----------------------------------------------------------
C-----Caract�re indicateur d'�criture d'un fichier
C-----contenant les compositions o� la minimisation a �chou�
C-----(cas NPT=0)
C-----------------------------------------------------------
      CHARACTER * 1 CAR_FICH_550
C----------------------------------------------
C-----DP charg�s : R�pertoire et fichier de DdE
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
C     I Derni�re mise � jour : 21/06/2005                              I
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
C-----Formats de r�els
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
     $         4X , 'Incr�ment (eV) = ' , G12.6 )
  600 FORMAT ( "Nombre de types chimiques : " , 30X , I4 )
  601 FORMAT ( " Nombre d'�l�ments d'addition : " , 30X , I4 )
  602 FORMAT ( "   Nombre de sous-r�seaux : " , 36X , I4 )
  603 FORMAT
     $ ( '      Nombre de sous-r�seaux interstitiels :' , 26X , I4 )
  604 FORMAT
     $ ( '        Nombre de types de complexes pris en compte :' ,
     $   3X , I4 )
  705 FORMAT ( 15X , 100 ( 2X , I4 ) )
C-----Formats relatifs aux DP charg�s
 2000 FORMAT
     $( 'Charges (unit�s |e-|) "n, A, p, D" pour le volume de SC =' ,
     $  F10.3 , ' A^3' )
 2001 FORMAT ( 5X , 'It�r.' , 4X , 'mu(e) (eV)' ,
     $ 6X , '-[Q(n)+Q(A)]' , 6X , '+Q(p)+Q(D)' )
 2002 FORMAT ( 5X , 'It�r.' , 4X , 'mu(e) (eV)' ,
     $ 12X , '-Q(n)', 12X , '-Q(A)' , 12X , '+Q(p)' , 12X , '+Q(D)' ,
     $ 12X , 'Q(totale)' , 6X , 'dQ(totale)/dmu_e' )
 3001 FORMAT ( 2X , I6 , 3 ( 2X , G16.7 ) )
 3002 FORMAT ( 2X , I6 , 7 ( 2X , G16.7 ) )
 3005 FORMAT ( '              T = ' , F10.2 , ' K' )
 3010 FORMAT ( 3 ( 5X , I3 ) , 4X , G16.7  )
 3020 FORMAT ( 2 ( 5X , I3 ) , 4X , G16.7  )
C####################################
C####################################
C#####Fin des d�clarations et formats
C####################################
C####################################
C-------------------------------------------------------------------
C-----Ecriture � l'�cran des r�f�rences du programme avant ex�cution
C-------------------------------------------------------------------
        CALL REF_PROG
C       CALL REF_PROG_DIRECT
C&&&&& NOTE : la partie "lecture du fichier de donn�es"
C&&&&& est commune aux cas sans et avec charges
C&&&&& => elle contient des tests sur la valeur de INDIC_CHARGE
C&&&&& A l'inverse, la suite du programme est scind�e en deux parties,
C&&&&& la premi�re pour DP non charg�s, la seconde pour DP charg�s,
C&&&&& et de tels tests y sont donc inutiles.
C##################################
C##################################
C#####Lecture du fichier de donn�es
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
C-----total et intrins�que (hors �l�ments d'addition)
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
C       WRITE ( * , * ) "�l�ments d'addition non disponibles"
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
C-----Nombre de sous-r�seaux
C-----(entendu au sens d'ensemble d'atomes de m�me environnement)
C----------------------------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) N_R
C--------------------------------------------------------
C-----Pr�sence de sous-r�seaux interstitiels (O/o ou N/n)
C--------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_R_INTER
C----------------------------------------------------------------
C-----Prise en compte des interstitiels intrins�ques (O/o ou N/n)
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
     $ 'Indicateur de pr�sence de sous-r�seaux interstitiels = ' ,
     $   INDIC_R_INTER
         WRITE ( * , * )
     $ '=> ne pas s�lectionner la prise en compte'
         WRITE ( * , * )
     $ 'des interstitiels intrins�ques'
         WRITE ( * , * )
     $ '----------------------------------------------------'
         CALL INTERRUPTION
        END IF
      END IF
C-----------------------------------------------------------------
C-----Nombres de sous-r�seaux occup�s par les esp�ces intrins�ques
C-----� l'�tat fondamental
C-----------------------------------------------------------------
      ALLOCATE ( N_R_N ( N_TYP_INTR ) )
        READ ( 10 , 4 )
        READ ( 10 , * ) 
     $ ( N_R_N ( I_TYP_INTR ) , I_TYP_INTR = 1 , N_TYP_INTR )
C-----------------------------------
C-----Calcul auxiliaire :
C-----nombres cumul�s correspondants
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
C-----Calcul du nombre de sous-r�seaux interstitiels
C---------------------------------------------------
      N_R_INTER = N_R
      DO I_TYP_INTR = 1 , N_TYP_INTR
        N_R_INTER = N_R_INTER - N_R_N ( I_TYP_INTR )
      END DO
      IF ( N_R_INTER .LT. 0 ) THEN
        WRITE ( * , * )
     $  '------------------------------------------------------------'
        WRITE ( * , * )
     $   'Le nombre total de sous-r�seaux doit �tre au moins �gal'
        WRITE ( * , * )
     $   '� la somme des nombres de sous-r�seaux intrins�ques'
        WRITE ( * , * )
     $   '(la diff�rence correspondant aux sous-r�seaux interstitiels)'
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
     $  'Incoh�rence entre le nombre de sous-r�seaux interstitiels'
        WRITE ( * , * )
     $  ' = ' , N_R_INTER
        WRITE ( * , * )
     $  "et l'indicateur de pr�sence de ces sous-r�seaux"
        WRITE ( * , * )
     $  '---------------------------------------------------------'
        CALL INTERRUPTION
      END IF
C-----------------------------------------------------------------
C-----Etape auxiliaire :
C-----recopiage des nombres de sous-r�seaux par esp�ce intrins�que
C-----(la suite du programme ne g�re pas encore N types intrins�ques)
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
C-----Lecture du nombre de sites par maille pour chaque sous-r�seau
C------------------------------------------------------------------
      ALLOCATE ( P_R ( N_R ) )
      READ ( 10 , 4 )
      READ ( 10 , * ) ( P_R ( I_R ) , I_R = 1 , N_R )
C-----------------------------------------------
C-----Pr�sence de d�fauts complexes (O/o ou N/n)
C-----------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , '(A)' ) INDIC_COMPLEXES
C-----------------------------------------------------------
C-----Mise � z�ro �ventuelle du nombre de types de complexes
C-----------------------------------------------------------
       IF ( INDIC_COMPLEXES .EQ. 'N' .OR. INDIC_COMPLEXES .EQ. 'n' ) 
     $ THEN
        N_TYPES_COMPLEXES = 0
      END IF
C-----------------------------------------------------------
C-----Temp�rature prescrite (K) (non utilis�e en mode NPT=0)
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
C-----Energies de r�f�rence des diverses esp�ces chimiques (eV/at.)
C------------------------------------------------------------------
      ALLOCATE ( E_REF_TYP ( N_TYP ) )
      READ ( 10 , 4 )
      READ ( 10 , * ) ( E_REF_TYP ( I_TYP ) , I_TYP = 1 , N_TYP )
C----------------------------------------------------
C-----Choix d'�criture des grandeurs thermodynamiques 
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
C-----Indicateur d'�criture de l'�nergie libre par atome :
C-----�nergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      READ ( 10 , 5 )
      READ ( 10 , * ) INDIC_G
      IF ( .NOT.
     $    ( INDIC_G .EQ. 'T' .OR. INDIC_G .EQ. 't'
     $ .OR. INDIC_G .EQ. 'F' .OR. INDIC_G .EQ. 'f' ) )
     $ THEN
        WRITE ( * , * ) '---------------------------------------'
        WRITE ( * , * ) "Ecriture de l'�nergie libre par atome  :"
        WRITE ( * , * ) 'choisir T/t (totale) ou F/f (formation)'
        WRITE ( * , * ) '---------------------------------------'
        CALL INTERRUPTION
      END IF
C----------------------------------------------------------------
C-----Le choix entre �nergie libre totale et de formation
C-----n'est permis que si l'�criture par atome a �t� s�lectionn�e
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
     $ "=> seule l'option d'�nergie libre totale est autoris�e"
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
C-----(pour produire une version "brid�e")
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
C-----Indicateur (0/1) de DP charg�s
C-----------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) INDIC_CHARGE
      IF ( .NOT. ( INDIC_CHARGE .EQ. 0 .OR. INDIC_CHARGE .EQ. 1 ) ) THEN
        write (* , * ) "V�rifier INDIC_CHARGE = 0/1"
        stop
      END IF
C=======================================================
C=====Sous-section "param�tres g�n�raux pour DP charg�s"
C=====(lecture optionnelle si INDIC_CHARGE = 1)
C=======================================================
       IF ( INDIC_CHARGE .EQ. 1 ) THEN
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 13 : 57 )
     $   .NE. 'Sous-section g�n�rale relative aux DP charg�s' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------------------
C-----Nombre maximum d'�tats de charg�s consid�r�s
C-----et type de calcul (1=balayage en mu_e, 2=neutralit� �lectrique)
C--------------------------------------------------------------------
      READ ( 10 , 6 )
      READ ( 10 , * ) N_MAX_CHARGES , I_CALC_CHARGE
C-------------------------------------------------------
C-----Pr�sence de d�fauts complexes charg�s (O/o ou N/n)
C-------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , '(A)' ) INDIC_COMPLEXES_Q
C-----INDIC_COMPLEXES ne contr�le que les complexes non charg�s
C-----et n'est pas utilis� en cas de DP charg�s.
C-----Cependant, pour la clart�, on impose la coh�rence
C-----entre INDIC_COMPLEXES et INDIC_COMPLEXES_Q
C-----(car INDIC_COMPLEXES est dans les "param�tres g�n�raux"
C-----et est donc requis, comme INDIC_COMPLEXES_Q, en cas de DP charg�s)
C----- INDIC_COMPLEXES_Q = O ==> INDIC_COMPLEXES = O
C----- INDIC_COMPLEXES_Q = N ==> INDIC_COMPLEXES = N
      IF ( INDIC_COMPLEXES_Q .EQ. 'o' .OR. INDIC_COMPLEXES_Q .EQ. 'O' )
     $ THEN
        IF ( INDIC_COMPLEXES .EQ. 'n' .OR. INDIC_COMPLEXES .EQ. 'N' )
     $  THEN
         write(*,*) "Pr�sence de complexes charg�s"
         write(*,*) "==> activer aussi INDIC_COMPLEXES (champ avant T)"
         WRITE(*,*)
         stop
        END IF
      END IF
      IF ( INDIC_COMPLEXES_Q .EQ. 'n' .OR. INDIC_COMPLEXES_Q .EQ. 'N' )
     $ THEN
        IF ( INDIC_COMPLEXES .EQ. 'o' .OR. INDIC_COMPLEXES .EQ. 'O' )
     $  THEN
         write(*,*) "Absence de complexes charg�s"
         write(*,*)
     $   "==> d�sactiver aussi INDIC_COMPLEXES (champ avant T)"
         WRITE(*,*)
         stop
        END IF
      END IF
C-------------------------------------------------------------------
C-----Mise � z�ro �ventuelle du nombre de types de complexes charg�s
C-------------------------------------------------------------------
       IF ( INDIC_COMPLEXES_Q .EQ. 'N' .OR. INDIC_COMPLEXES_Q .EQ. 'n' )
     $ THEN
        N_TYPES_COMPLEXES_Q = 0
      END IF
C-----------------------------------------------------
C-----Les DP charg�s ne sont pris en compte qu'en muVT
C-----------------------------------------------------
      IF ( INDIC_CHARGE .EQ. 1 ) THEN
        IF
     $ ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
     $ THEN
       WRITE ( * , * ) '-------------------------'
       WRITE ( * , * ) 'Traitement des DP charg�s :'
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
C-----Potentiel chimique �lectronique
C-----valeur initiale (mode 2) ou minimale (mode 1),
C-----puis incr�ment (mode 1) et valeur maximale (mode 1)
C--------------------------------------------------------
      READ ( 10 , 6 )
      READ ( 10 , * ) POT_CHIM_ELEC ,
     $                D_POT_CHIM_ELEC , POT_CHIM_ELEC_MAX
C---------------------------------
C-----R�pertoire et fichier de DdE
C---------------------------------
      READ ( 10 , 4 )
      READ ( 10 , ' ( A ) ' ) REP_DDE
      LONG_REP_DDE = INDEX ( REP_DDE , ' ' ) - 1
      READ ( 10 , 4 )
      READ ( 10 , ' ( A ) ' ) FICH_DDE
      LONG_FICH_DDE = INDEX ( FICH_DDE , ' ' ) - 1
C-----------------------------------------------
C-----Pr�cision d'arr�t et nombre maximum de pas
C-----pour l'algorithme NR "�lectronique"
C-----------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) EPS , N_MAX_PAS_NR_ELEC
C-----Fin de lecture optionnelle des param�tres g�n�raux pour DP charg�s
      END IF
C##########################################################
C#####Fin de la lecture de la section "PARAMETRES GENERAUX"
C##########################################################
C===================================================
C=====Ouverture pr�liminaire de tableaux n�cessaires 
C=====� la lecture des sections suivantes
C===================================================
C---------------------------------------------------------------------
C-----Fractions atomiques dans l'alliage
C-----telles que fix�es par les potentiels chimiques
C-----et fractions atomiques
C-----d�duites de la valeur pour l'�l�ment sp�cifi� et des contraintes
C---------------------------------------------------------------------
      ALLOCATE ( X_AT ( N_TYP ) )
      ALLOCATE ( X_AT_0_CTR ( N_TYP ) )
C----------------------------------------------------------------
C-----Param�tres pertinents seulement si N_TYP > 2,
C-----et m�me seulement en muVT pour X_AT_INF_CTR et X_AT_SUP_CTR
C-----(tableaux ouverts dans tous les cas n�anmoins)
C----------------------------------------------------------------
      ALLOCATE ( COEF_CTR ( N_TYP - 1 , N_TYP + 1 ) )
      ALLOCATE ( D_X_AT ( N_TYP ) )
      ALLOCATE ( D_X_AT_INIT ( N_TYP ) )
      ALLOCATE ( X_AT_INF_CTR ( N_TYP ) )
      ALLOCATE ( X_AT_SUP_CTR ( N_TYP ) )
C###########################################################
C#####PARAMETRES RELATIFS A UN CALCUL NPT / NPT=0
C#####Champs lus et utilis�s seulement si calcul NPT / NPT=0
C###########################################################
C-----Test sur le type de calcul : autre que muVT ?
      IF
     $ ( .NOT.
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' ) )
     $ THEN
C------------------------------------------------
C-----Recherche de l'en-t�te de la grande section
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
C-----* P/p : point (T,x) fix�s (modes NPT et NPT=0)
C-----* T : balayage en temp�rature (mode NPT)
C-----* x : balayage en composition (modes NPT et NPT=0)
C-----* Tx et xT  : doubles balayages en temp�rature et composition
C----- --> boucle externe sur T ou x pour Tx ou xT respectivement
C-----(mode NPT)
C------------------------------------------------------------------
C------------------------------------------------------------------
C----------------------------------------------
C-----Recherche de l'en-t�te de la sous-section
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
C- - -Test sur la validit� de cet indicateur lu
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
     $  'T : balayage en temp�rature'
         WRITE ( * , * ) 
     $  'x : balayage en composition'
         WRITE ( * , * )
     $  'Tx ou xT : doubles balayages'
         WRITE ( * , * ) 
     $ '------------------------------------'
         CALL INTERRUPTION
        END IF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Fin du test sur la validit� de cet indicateur lu
C- - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Avertissement concernant le mode NPT(Tx)
C     IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
C    $  .EQ. 'Tx' ) THEN
C       write(*,*) "************************************************"
C       write(*,*) "Vous avez s�lectionn� le mode NPT(Tx)."
C       write(*,*) "Cependant, � l'inverse de son homologue NPT(xT),"
C       write(*,*) "ce mode n'a pas �t� test� sur divers exemples."
C       write(*,*) "==> Il serait pr�f�rable de privil�gier NPT(xT)."
C       write(*,*) "Si vous souhaitez poursuivre tout de m�me"
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
C- - -Troisi�me test : dans le cas d'un calcul NPT=0 avec N_TYP # 3
C- - -(section PARAMETRES GENERAUX), l'option "x" n'est pas possible
C- - -(seule l'option "P/p" est possible).
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT )
     $  .EQ. 'x' ) THEN
       IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' 
     $ .AND. N_TYP .NE. 3 ) THEN 
        WRITE ( * , * )
        WRITE ( * , * )
     $ "Syst�me avec N_TYP # 3 types chimiques + mode NPT=0 :"
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
C-----Recherche de l'en-t�te de la sous-section
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
C- - -V�rification de somme(x_at)=1
C- - -dans les cas o� cette donn�e est utilis�e
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
C- - -Fin de v�rification de somme(x_at)=1
C- - -dans les cas o� cette donn�e est utilis�e
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
C-----autour de la stoechiom�trie :
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
C-----Recherche de l'en-t�te de la sous-section
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
C-----d'un balayage en T ou en fraction atomique pour l'esp�ce choisie
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
C-----Recherche de l'en-t�te de la sous-section
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
C-----Recherche de l'en-t�te de la sous-section
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
C-----Recherche de l'en-t�te de la sous-section
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
C-----Recherche de l'en-t�te de la sous-section
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
C-----d'un balayage en T ou en fraction atomique pour l'esp�ce choisie
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
C#####Champs lus et utilis�s seulement si calcul muVT ou NPT
C###########################################################
C-----Test du type de calcul : NPT ou muVT
      IF ( .NOT.
     $    ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' ) )
     $ THEN
C============================================
C============================================
C=====Section lue seulement si DP non charg�s
C============================================
C============================================
       IF ( INDIC_CHARGE .EQ. 0 ) THEN
C-------------------------------------------------------
C-----Recherche de l'en-t�te de section
C-----(devenue inutile car en-t�tes int�rieurs utilis�s) 
C-------------------------------------------------------
C       REWIND ( 10 )
C       CHAINE_LECT = ''
C       DO WHILE ( CHAINE_LECT ( 6 : 50 )
C    $   .NE. 'PARAMETRES RELATIFS A UN CALCUL muVT ou NPT' )
C         READ ( 10 , '(A)' ) CHAINE_LECT
C       END DO
C======================================================================
C=====1) Champs � inclure ssi N_TYP > 2, en NPT(x/Tx/xT) comme en muVT
C===== --> syst�me de N_TYP - 2 contraintes
C=====y compris en muVT lorsque le filtre "fen�tres"
C=====n'est pas s�lectionn� 
C=====(cf. section suivante "PARAMETRES SPECIFIQUES A UN CALCUL muVT"),
C=====bien que ces champs soient alors inutilis�s dans ce cas.
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
C-----Recherche de l'en-t�te
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
C-----La N_TYP - 1 �me contrainte est l'unit� de la somme des fractions
C----------------------------------------------------------------------
       DO I_TYP = 1 , N_TYP
        COEF_CTR ( N_TYP - 1 , I_TYP ) = 1.D0
       END DO
       COEF_CTR ( N_TYP - 1 , N_TYP + 1 ) = 0.D0
C===================================================
C=====Fin des donn�es 1) lues seulement si N_TYP > 2
C=====en NPT(x/Tx/xT) comme en muVT
C===================================================
      END IF
C======================================================================
C=====2) Champs � inclure (i) en NPT(x/Tx/xT),
C=====                    (ii) en muVT si N_TYP > 2,
C=====y compris en muVT lorsque le filtre "fen�tres"
C=====n'est pas s�lectionn�
C=====(cf. section suivante "PARAMETRES SPECIFIQUES A UN CALCUL muVT"),
C=====bien que ces champs soient alors inutilis�s dans ce cas.
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
C-----El�ment dont on �tudie l'effet de l'enrichissement
C-----(toutes les autres fractions atomiques
C-----sont fix�es dans des fen�tres � l'�criture)
C-----(toujours utilis� en NPT(x/Tx/xT),
C-----utilis� si N_TYP > 2 en muVT)
C-------------------------------------------------------
C-------------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $   .NE. "El�ment dont on �tudie l'effet de l'enrichissement" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 5 )
        READ ( 10 , * ) I_TYP_0
C--------------------------------------------------------------
C--------------------------------------------------------------
C-----Valeurs initiale et finale de fraction atomique �crite
C-----pour l'�l�ment sp�cifi� (dont on �tudie l'enrichissement)
C-----(toujours utilis� en NPT(x/Tx/xT), 
C-----utilis� si N_TYP > 2 en muVT)
C--------------------------------------------------------------
C--------------------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 50 )
     $   .NE. "Valeurs extr�males du domaine de fraction atomique" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 5 )
        READ ( 10 , * ) X_TYP_0_INIT , X_TYP_0_FIN
C-----Test pr�liminaire
        IF ( X_TYP_0_INIT .LT. 0.D0
     $  .OR. X_TYP_0_FIN .LT. 0.D0 ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Valeurs initiale et finale de fraction atomique �crite'
          WRITE ( * , * )
          WRITE ( * , * )
     $ "pour l'�l�ment sp�cifi� (dont on �tudie l'enrichissement) :"
          WRITE ( * , * )
     $ 'on doit avoir x_init > 0 et x_fin > 0'
          WRITE ( * , * )
     $ '---------------------------------------------------------'
        CALL INTERRUPTION
        END IF
C--------------------------------------------------------
C-----En plus de leur emploi en mode NPT(x/Tx/xT),
C-----X_TYP_0_INIT et X_TYP_0_FIN servent aussi � affecter
C-----X_AT_INF_CTR et X_AT_SUP_CTR, qui sont utilis�s
C-----seulement en mode muVT si N_TYP > 2.
C--------------------------------------------------------
         X_AT_INF_CTR ( I_TYP_0 ) = X_TYP_0_INIT
         X_AT_SUP_CTR ( I_TYP_0 ) = X_TYP_0_FIN
C=======================================================
C=====Fin de la lecture optionnelle
C=====des donn�es 2) toujours utilis�es en NPT(x/Tx/xT),
C=====ou utilis�es si N_TYP > 2 en muVT
C=======================================================
        END IF
C===================================
C===================================
C=====Fin du test "DP non charg�s ?"
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
C#####Champs lus et utilis�s seulement si calcul muVT
C####################################################
C-----Test sur le type de calcul : muVT ?
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ THEN  
C-------------------------------------------------------
C-----Recherche de l'en-t�te de section
C-----(devenue inutile car en-t�tes int�rieurs utilis�s) 
C-------------------------------------------------------
C       REWIND ( 10 )
C       CHAINE_LECT = ''
C       DO WHILE ( CHAINE_LECT ( 6 : 44 )
C    $   .NE. 'PARAMETRES SPECIFIQUES A UN CALCUL muVT' )
C         READ ( 10 , '(A)' ) CHAINE_LECT
C       END DO
C============================================
C============================================
C=====Premier cas de lecture : DP non charg�s
C===== --> balayage en mu_i
C============================================
C============================================
       IF ( INDIC_CHARGE .EQ. 0 ) THEN
C--------------------------------------------------------
C--------------------------------------------------------
C-----Type chimique intrins�que de r�f�rence
C-----pour les potentiels chimiques,
C-----pr�cision et nombre maximum d'it�rations
C-----pour l'arr�t de la boucle d'autocoh�rence sur POT_1
C--------------------------------------------------------
C--------------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. '(i) Type chimique intrins�que de r�f�rence' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
       READ ( 10 , 3 )
       READ ( 10 , * ) I_TYP_REF_MU , PRECISION_MU_1 , N_ITER_MAX_MU_1
C-----------------
C-----V�rification
C-----------------
      IF ( N_ITER_MAX_MU_1 .LT. 1 ) THEN
       write(*,*)
     $ " *** PROGRAMME INTERROMPU ***"
       write(*,*)
     $ "Autocoh�rence sur mu1 -> v�rifier que N_ITER_MAX_MU_1 > 0"
       write(*,*)
     $ "(N_ITER_MAX_MU_1 = 1 ==> pas d'autocoh�rence)"
        stop
      END IF
C----------------------------------------
C-----La suite du programme ne fonctionne
C-----que si le type de r�f�rence est 1
C----------------------------------------
      if ( I_TYP_REF_MU .ne. 1 ) then
        write ( * , * ) 'La pr�sente version ne fonctionne'
        write ( * , * ) 'que pour type_r�f.(mu) = 1'
        call interruption
      end if
C========================================
C=====Donn�es lues seulement si N_TYP > 2
C========================================
       IF ( N_TYP .GT. 2 ) THEN
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C-----Indicateur d'�criture des seuls points contenus dans la fen�tre
C-----(indicateur lu seulement pour plus de 2 types)
C--------------------------------------------------------------------
C--------------------------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 47 )
     $   .NE. "Indicateur d'�criture des seuls points contenus" )
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
     $          "Indicateur d'�criture des points  : O/o ou N/n" 
          WRITE ( * , * ) 
     $          '----------------------------------------------'
          CALL INTERRUPTION
        END IF
C------------------------------------------------
C------------------------------------------------
C-----Demi-largeurs des fen�tres en composition
C-----pour les �l�ments autres que celui sp�cifi�
C-----(seulement si N_TYP > 2)
C------------------------------------------------
C------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 42 )
     $   .NE. "Demi-largeurs des fen�tres en composition" )
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
C-----Dans le cas o� N_TYP = 2, rien n'est lu
C-----mais il faut tout de m�me initialiser I_ECRIT_FENETRE
C----------------------------------------------------------
      ELSE
        I_ECRIT_FENETRE = 0
C=============================================================
C=====Fin de la lecture des param�tres pertinents si N_TYP > 2
C=============================================================
      END IF
C-----------------------------------------------------------
C-----Recherche de l'en-t�te de la sous-section muVT
C-----sur les balayages en mu(intrins�ques) et mu(additions)
C-----(permet de sauter la section optionnelle pr�c�dente)
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
C-----Propri�t�s des s�ries d'�carts de potentiels chimiques
C-----pour les �l�ments intrins�ques autres que celui de r�f�rence
C-----------------------------------------------------------------
C-----------------------------------------------------------------
      ALLOCATE ( D_POT_REF_INTR_INIT ( N_TYP_INTR ) )
      ALLOCATE ( N_D_POT_REF_INTR ( N_TYP_INTR ) )
      ALLOCATE ( PAS_D_POT_REF_INTR ( N_TYP_INTR ) )
C--------------------------------------------------
C-----Initialisations (� 1 pour les nombres de pas)
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
C-----Connexion avec la suite du programme, dont la pr�sente version
C-----ne fonctionne pas pour plus de 3 types intrins�ques
C-----(m�me probl�me que plus haut pour N_R_1, N_R_2 et N_R_3)
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
C-----Propri�t�s des s�ries des mu_i (additions)
C-----La boucle de N_TYP_INTR + 1 � N_TYP ci-dessous
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
C=====Deuxi�me cas de lecture : DP charg�s
C===== --> pas de balayage en mu_i
C=========================================
C=========================================
      ELSE
C---------------------------------------------
C-----Recherche du d�but de cette sous-section
C---------------------------------------------
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 51 )
     $   .NE. 'Cas "muVT+charges" : sp�cification directe des mu_i' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
        DO I = 1 , 3
         READ ( 10 , * )
        END DO
C----------------------------------------------------------------
C-----mu(2) (eV) et �ventuellement mu(3) (si ternaire intrins�que)
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
C=====ou "muVT+point" (avec DP charg�s)
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
C#####Champs lus et utilis�s seulement si calcul NPT
C###################################################
C-----Test "mode NPT ?" pour la lecture de cette grande section
      IF
     $ ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
C---------------------------------------------
C-----Recherche de l'en-t�te de grande section
C---------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 43 )
     $   .NE. 'PARAMETRES SPECIFIQUES A UN CALCUL NPT' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-------------------------------------------------
C-------------------------------------------------
C-----Valeurs extr�males du domaine de temp�rature
C----- --> champs utilis�s en NPT(T/Tx/xT)
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
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 1 : 44 )
     $   .NE. "Valeurs extr�males du domaine de temp�rature" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C------------
C-----Lecture
C------------
        READ ( 10 , 1 )
        READ ( 10 , * ) T_INIT , T_FIN
C-----Test pr�liminaire 1
        IF ( T_INIT .LT. 0.D0
     $  .OR. T_FIN .LT. 0.D0 ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Domaine de temp�rature en NPT(T/Tx/xT) :'
          WRITE ( * , * )
     $ ' choisir T_min > 0 et T_max > 0'
          WRITE ( * , * )
     $ '---------------------------------------------------------'
        CALL INTERRUPTION
        END IF
C-----Test pr�liminaire 2
        IF ( T_INIT .LT. T_FIN ) THEN
          WRITE ( * , * )
     $ '---------------------------------------------------------'
          WRITE ( * , * )
     $ 'Domaine de temp�rature en NPT(T) :'
          WRITE ( * , * )
     $ ' veuillez choisir T_init > T_fin'
          WRITE ( * , * )
     $ ' (convergence tr�s sensible aux valeurs initiales sinon).'
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
C=====Mode NPT : param�tres de l'algorithme
C=====toujours lus en NPT (p, x ou T)
C==========================================
C==========================================
C-------------------------------------------
C-----Recherche de l'en-t�te de sous-section
C-------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 17 : 56 )
     $   .NE. "Sous-section NPT relative � l'algorithme" )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-----Lecture s�quentielle des champs (pas de recherche d'en-t�te)
C-----dans cette sous-section
C--------------------------------------------------------------------
C-----Nombre maximal d'it�rations pour l'algorithme de Newton-Raphson
C--------------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) N_ITER_MAX_NPT
C---------------------------------------------------------------------
C-----Pr�cision requise pour l'arr�t de l'algorithme de Newton-Raphson
C-----(maximum des valeurs absolues
C-----des composantes de l'�cart entre deux pas)
C---------------------------------------------------------------------
        READ ( 10 , 5 )
        READ ( 10 , * ) PRECISION_NPT
C-----------------------------------------------------------------
C-----Fr�quence (nombre de pas) de calcul de la matrice jacobienne
C-----------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) P_J_NPT
C----------------------------------------
C-----Valeur du param�tre alpha dans NRCG
C----------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) ALPHA_NRCG_NPT
C----------------------------------------------------------------------
C-----NRCG :
C----- 1) indicateur � l'�tape pr�liminaire pour "�limination de x_1<0"
C-----= 1 si r�duction du pas suivant les seules composantes x_1(i)<0,
C-----= 2 si r�duction de toutes les composantes du pas
C----- 2) valeur du coefficient de r�duction
C----------------------------------------------------------------------
        READ ( 10 , 7 )
        READ ( 10 , * ) INDIC_TYPE_REDUC_NRCG_NPT , COEF_REDUC_NRCG_NPT
C--------------------------------------------------------------------
C-----Valeur minimale "lambda_min" pour la r�duction du pas dans NRCG
C--------------------------------------------------------------------
        READ ( 10 , 4 )
        READ ( 10 , * ) VALEUR_LAMBDA_MIN_NRCG_NPT
C=================================
C=================================
C=====Mode NPT : valeurs initiales
C=================================
C=================================
C---------------------------------------------------------
C-----La variable "nombre de mailles" est initialis�e � 1,
C-----car sa valeur ne joue pas sur la r�solution
C-----(inconnues = variables intensives).
C---------------------------------------------------------
         N_MAILLES_INIT = 1.D0
C-----------------------------------------------------
C-----Ouverture du tableau des valeurs initiales de DP
C-----------------------------------------------------
      ALLOCATE ( LOG_X_D_R_INIT ( 0 : N_TYP , N_R ) )
C-------------------------------------------
C-----Recherche de l'en-t�te de sous-section
C-------------------------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 12 : 58 )
     $   .NE. 'Sous-section NPT relative aux valeurs initiales' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-----Lecture s�quentielle (pas de recherche d'en-t�te)
C-----pour le champ suivant, qui est lu dans tous les cas NPT,
C-----i.e. NPT(p), NPT(x/Tx/xT) et NPT(T)
        READ ( 10 , 4 )
        READ ( 10 , * ) INDIC_LECT_VAL_INIT_NPT
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Test de compatibilit� : la lecture des valeurs initiales
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
C-----indicateur de r�initialisation des variables NPT
C-----par les valeurs "DATA.adpi" � chaque d�but d'it�ration externe,
C-----i.e. entre les boucles externe (x) et interne (T)
C-----Valeur 1 : les valeurs "DATA.adpi" sont utilis�es
C-----� chaque incr�ment de x, avant le d�but de la boucle en T,
C-----� l'int�rieur de laquelle les val. init. sont affect�es
C-----proche en proche "PeP" de (x_n,T_p) � (x_n,T_p+1)
C-----Valeur 0 : les valeurs "DATA.adpi" sont utilis�es seulement
C-----avant le d�but des boucles, puis l'affectation se fait toujours
C-----en mode "PeP", m�me au changement de composition (x_n,T_P) ->
C----- (x_n+1,T_1) (N et P �tant les nombres de pas en x et T (resp.).
C-----Remarque : ce champ n'est lu (et ne doit donc �tre inclus en mode
C-----NPT(xT)), que si n'est pas s�lectionn�e ci-dessus l'option
C-----"lect. val. init. dans fich.", cette option prenant alors
C-----le pas sur les valeurs "DATA.adpi", et �tant utilis�e
C-----� chaque changement de x, avant le d�but de la boucle sur T
C-----(c'est-�-dire avec la m�me fr�quence que les valeurs "DATA.adpi"
C-----lorsque le champ pr�sent vaut 1).
C----------------------------------------------------------------------
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. "xT" ) THEN
C-----L'indicateur est pr�alablement initialis�, pour le cas o�
C-----"lect. val. init. fich." i.e. INDIC_LECT_VAL_INIT_NPT = 1
C-----(puisque dans ce cas, la lecture ci-dessous n'a pas lieu).
C-----(Cette initialisation n'est peut-�tre pas n�cessaire,
C-----mais elle est tout de m�me faite "par pr�caution").
         IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 0 ) THEN
C-----Recherche de l'en-t�te
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
C-----l'indicateur est tout de m�me initialis�,
C-----puisque dans ce cas, la lecture ci-dessus n'a pas lieu.
C-----(Cette initialisation n'est peut-�tre pas n�cessaire,
C-----mais elle est tout de m�me faite "par pr�caution").
         ELSE
           INDIC_LECT_VAL_INIT_XT = 0
         END IF
        END IF
C---------------------------------------------------------------
C---------------------------------------------------------------
C-----1) Cas de lecture dans un fichier de valeurs initiales :
C-----cette lecture est effectu�e seulement en mode NPT(x/Tx/xT),
C-----en raison du test pr�c�dent (sinon, le programme s'arr�te)
C---------------------------------------------------------------
C---------------------------------------------------------------
        IF ( INDIC_LECT_VAL_INIT_NPT .EQ. 1 ) THEN
C- - -D'apr�s la remarque ci-dessus, on se trouve forc�ment ici
C- - -dans l'un des cas NPT(x/Tx/xT).
C-----Recherche de l'en-t�te
         REWIND ( 10 )
         CHAINE_LECT = ''
         DO WHILE ( CHAINE_LECT ( 1 : 35 )
     $    .NE. 'Nom du fichier de valeurs initiales' )
           READ ( 10 , '(A)' ) CHAINE_LECT
         END DO
C-----Lecture s�quentielle (pas de recherche d'en-t�te interm�diaire)
C-----du nom de ce fichier et de son nombre de lignes
         READ ( 10 , '(A)' ) CHAINE_LECT
         READ ( 10 , '(A)' ) FICH_VAL_INIT
          LONG_FICH_VAL_INIT
     $  = INDEX ( FICH_VAL_INIT , ' ' ) - 1
         READ ( 10 , 4 )
         READ ( 10 , * ) N_LIGNES_FICH_VAL_INIT
C- - -Test sur la coh�rence entre la longueur du fichier
C- - -et le nombre de pas en composition :
C- - - le nombre de pas et le nombre de lignes
C- - -doivent �tre multiples ou sous-multiples l'un de l'autre.
C- - -D'apr�s la remarque ci-dessus, on se trouve forc�ment ici
C- - -dans l'un des cas NPT(x/Tx/xT)
C- - - ==> c'est N_PAS_X_TYP_0 lu qui est utilis�
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
C-----Recherche de l'en-t�te
         REWIND ( 10 )
         CHAINE_LECT = ''
         DO WHILE ( CHAINE_LECT ( 1 : 30 )
     $     .NE. 'log_10 des quantit�s initiales' )
          READ ( 10 , '(A)' ) CHAINE_LECT
         END DO
C--------------------------------------------------------------
C-----"log_10 des quantit�s initiales", partie 1/3 :
C-----lecture des quantit�s initiales d'antisites et de lacunes
C--------------------------------------------------------------
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Lecture "synth�tique" (adapt�e � N types intrins�ques)
C- - -des valeurs initiales pour les antisites
C- - -des diverses esp�ces intrins�ques sur les divers sous-r�seaux
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       DO I_TYP = 1 , N_TYP_INTR
         DO J_TYP = 1 , N_TYP_INTR
          IF ( J_TYP .NE. I_TYP ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( LOG_X_D_R_INIT ( J_TYP , I_R ) ,
     $     I_R = RHO_N ( I_TYP - 1 ) + 1 , RHO_N ( I_TYP ) )
C-----Test pour �viter des valeurs initiales aberrantes
C-----(induit erreur dans l'algorithme de r�solution)
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
C-----"log_10 des quantit�s initiales", partie 2/3 :
C-----lecture optionnelle des  quantit�s initiales
C-----relatives aux interstitiels intrins�ques
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
C-----"log_10 des quantit�s initiales", partie 3/3 :
C-----lecture optionnelle des quantit�s initiales d'�l�ments d'addition
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
C- - -Cas de pr�sence d'interstitiels d'addition
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
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 47 )
     $ .NE. 'PARAMETRES DE LA SUPERCELLULE DE REFERENCE' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------
C-----Energie de r�f�rence de la cellule sans d�faut (eV)
C--------------------------------------------------------
      READ ( 10 , 4 )
      READ ( 10 , * ) E_REF
C------------------------------------------------------------
C-----Volume de r�f�rence de la cellule sans d�faut (A * * 3)
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
C#####     Lecture des param�tres de SC des DP
C======================================================
C#####CAS DE DP NON CHARGES : DP simples (et complexes)
C======================================================
C######################################################
C-----Test "DP non charg�s ?"
      IF ( INDIC_CHARGE .EQ. 0 ) THEN
C###################################################
C#####PARAMETRES RELATIFS AUX DP SIMPLES NON CHARGES
C###################################################
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 59 )
     $ .NE. 'PARAMETRES DES SUPERCELLULES DE DP SIMPLES NON CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C-------------------------------------------------------
C-----Ouverture des tableaux
C-----d'�nergies "brutes" des cellules avec d�fauts (eV)
C-----et de m�me taille que la cellule de r�f�rence
C-----et �nergies GC correspondantes
C-------------------------------------------------------
      ALLOCATE ( E_B_D_R ( 0 : N_TYP , N_R ) )
      ALLOCATE ( E_GC_D_R ( 0 : N_TYP , N_R ) )
C-------------------------------------
C-----M�mes quantit�s pour les volumes
C-------------------------------------
      ALLOCATE ( V_B_D_R ( 0 : N_TYP , N_R ) )
      ALLOCATE ( V_GC_D_R ( 0 : N_TYP , N_R ) )
C----------------------------------------------------------
C-----Enthalpies GC des d�fauts sur les divers sous-r�seaux
C----------------------------------------------------------
      ALLOCATE ( H_GC_D_R ( 0 : N_TYP , N_R ) )
C---------------------------------------------------------------
C-----Enthalpies de formation des DP sur les divers sous-r�seaux
C---------------------------------------------------------------
      ALLOCATE ( H_FORM_D_R ( 0 : N_TYP , N_R ) )
C--------------------------------------
C-----Termes compl�mentaires d'entropie
C--------------------------------------
      ALLOCATE ( Z_TYP_R ( N_R ) )
C-----------------------------------------------------------
C-----Initialisation des �nergies et volumes de supercellule
C-----------------------------------------------------------
      E_B_D_R = 0.D0
      V_B_D_R = 0.D0
C===================================================
C=====Lecture des �nergies des cellules avec d�fauts
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
C-----lecture des �nergies "brutes" d'antisites et de lacunes
C------------------------------------------------------------
C------------------------------------------------------------
C---------------------------------------------------------------
C-----Lecture "synth�tique" (adapt�e � N types intrins�ques)
C-----des �nergies d'antisites des diverses esp�ces intrins�ques
C----- sur les divers sous-r�seaux
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
C-----lecture optionnelle des �nergies d'interstitiels intrins�ques
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
C-----lecture optionnelle des �nergies d'�l�ments d'addition
C-----(substitutionnels et interstitiels)
C-----------------------------------------------------------
C-----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         READ ( 10 , * )
     $   ( E_B_D_R ( I_TYP , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C-----------------------------------------------
C-----Cas de pr�sence d'interstitiels d'addition
C-----------------------------------------------
        IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
           READ ( 10 , 2 )
           READ ( 10 , * )
     $   ( E_B_D_R ( I_TYP , I_R ) ,
     $     I_R = 1 + RHO_N ( N_TYP_INTR ) , N_R )
        END IF
      END DO
C==================================================
C=====Lecture des volumes des cellules avec d�fauts
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
C-----Lecture "synth�tique" (adapt�e � N types intrins�ques)
C-----des volumes d'antisites des diverses esp�ces intrins�ques
C----- sur les divers sous-r�seaux
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
C-----lecture optionnelle des volumes d'interstitiels intrins�ques
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
C-----Lecture optionnelle des volumes d'�l�ments d'addition
C-----(substitutionnels et interstitiels)
C----------------------------------------------------------
C----------------------------------------------------------
      DO I_TYP = N_TYP_INTR + 1 , N_TYP
         READ ( 10 , 2 )
         READ ( 10 , * )
     $   ( V_B_D_R ( I_TYP , I_R ) , I_R = 1 , RHO_N ( N_TYP_INTR ) )
C-----------------------------------------------
C-----Cas de pr�sence d'interstitiels d'addition
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
C-----Recherche de l'en-t�te
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
C-----Nombre de types de d�fauts complexes pris en compte
C--------------------------------------------------------
       READ ( 10 , 4 )
       READ ( 10 , * ) N_TYPES_COMPLEXES
C--------------------------------------------------
C-----Ouverture des tableaux relatifs aux complexes
C-----dont les valeurs vont �tre lues
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
     $ .NE. 'Caract�ristiques des complexes' )
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
         WRITE ( * , * ) 'Num�roter les complexes de 1 � N'
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
C-----Fin de la lecture optionnelle des param�tres de complexes
C--------------------------------------------------------------
      ELSE
C----------------------------------------------------------------
C-----Sinon, allocation de tableaux pass�s en param�tres � G_ADPI
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
C#####     Lecture des param�tres de SC des DP
C==================================================
C#####CAS DE DP CHARGES : DP simples (et complexes)
C==================================================
C##################################################
      ELSE
C=========================================
C=========================================
C=====Premi�re partie : DP simples charg�s
C=========================================
C=========================================
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 47 )
     $   .NE. 'PARAMETRES RELATIFS AUX DP SIMPLES CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------------------
C-----Recherche de la section des param�tres (GC) de DP charg�s
C--------------------------------------------------------------
      CHAINE_LECT = ''
      DO WHILE ( CHAINE_LECT ( 1 : 21 )
     $  .NE. 'Pour chaque DP charg�' )
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO
      DO I = 1 , 3
        READ ( 10 , '(A)' ) CHAINE_LECT
      END DO    
C===========================================
C=====Ouverture des tableaux pour DP charg�s
C===========================================
C-------------------------------------------------------
C-----Tableau des nombres d'�tats de charges pour les DP
C-------------------------------------------------------
      ALLOCATE ( NQ_D_R_Q ( 0 : N_TYP , N_R ) )
C-------------------
C-----Charges des DP
C-------------------
      ALLOCATE ( Q_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C--------------------------------------------------
C-----Energies "brutes" des cellules DP charg�s
C-----et de m�me taille que la cellule de r�f�rence
C-----et �nergies GC correspondantes
C--------------------------------------------------
      ALLOCATE ( E_B_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
      ALLOCATE ( E_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C-------------------------------------
C-----M�mes quantit�s pour les volumes
C-------------------------------------
      ALLOCATE ( V_B_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
      ALLOCATE ( V_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C----------------------------------------------------------
C-----Enthalpies GC des d�fauts sur les divers sous-r�seaux
C----------------------------------------------------------
      ALLOCATE ( H_GC_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C---------------------------------------------------------------
C-----Enthalpies de formation des DP sur les divers sous-r�seaux
C---------------------------------------------------------------
      ALLOCATE ( H_FORM_D_R_Q ( N_MAX_CHARGES , 0 : N_TYP , N_R ) )
C-------------------------------------------------------
C-----Initialisation � 0 des param�tres de DP charg�s
C-----(inutile pour le calcul, mais permet d'identifier
C-----les triplets (i_q,i_typ,i_r) non affect�s � un DP)
C-------------------------------------------------------
      NQ_D_R_Q = 0
      E_B_D_R_Q = 0.D0
      V_B_D_R_Q = 0.D0
      Q_B_D_R_Q = 0.D0
C=====================================================
C=====Lecture des param�tres des cellules avec d�fauts
C=====================================================
C------------------------------------------------------------
C------------------------------------------------------------
C-----Partie 1/3 : 
C-----lecture des param�tres d'antisites et de lacunes
C------------------------------------------------------------
C------------------------------------------------------------
C---------------------------------------------------------------
C-----Lecture "synth�tique" (adapt�e � N types intrins�ques)
C-----des �nergies d'antisites des diverses esp�ces intrins�ques
C----- sur les divers sous-r�seaux
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
C-----lecture optionnelle des �nergies d'interstitiels intrins�ques
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
C-----lecture optionnelle des �nergies d'�l�ments d'addition
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
C-----Cas de pr�sence d'interstitiels d'addition
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
C=====Deuxi�me partie : DP complexes charg�s
C=====(optionnelle)
C===========================================
C===========================================
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C---------------------------
C-----Recherche de l'en-t�te
C---------------------------
        REWIND ( 10 )
        CHAINE_LECT = ''
        DO WHILE ( CHAINE_LECT ( 6 : 54 )
     $ .NE. 'PARAMETRES RELATIFS AUX DEFAUTS COMPLEXES CHARGES' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C--------------------------------------------------
C-----Nombre maximum de sites des complexes charg�s
C--------------------------------------------------
       READ ( 10 , 5 )
       READ ( 10 , * ) N_MAX_SITES_COMPLEXES_Q
C----------------------------------------------------------------
C-----Nombre de types de d�fauts complexes charg�s pris en compte
C----------------------------------------------------------------
       READ ( 10 , 4 )
       READ ( 10 , * ) N_TYPES_COMPLEXES_Q
C----------------------------------------------------------
C-----Ouverture des tableaux relatifs aux complexes charg�s
C-----dont les valeurs vont �tre lues
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
     $ .NE. 'Caract�ristiques des complexes charg�s' )
          READ ( 10 , '(A)' ) CHAINE_LECT
        END DO
C----------------------------------------------
C-----Boucle sur les types de complexes charg�s
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
          WRITE ( * , * ) 'Num�roter les complexes charg�s de 1 � N'
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
C-----Lecture des �nergies et volumes "bruts" pour chaque �tat de charge
        READ ( 10 , * ) NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
        DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q )
         READ ( 10 , * )
     $   Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) ,
     $   E_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) ,
     $   V_B_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
        END DO
C-----Fin de boucle sur complexes charg�s
       END DO   
C----------------------------------------------------------------------
C-----Fin de la lecture optionnelle des param�tres de complexes charg�s
C----------------------------------------------------------------------
       END IF
C-----Fin du test "DP non charg�s ?"
C-----(conditionne le choix des blocs de donn�es lus)
      END IF
C##################################################
C##### Fin de lecture des param�tres de SC des DP
C==================================================
C#####     CAS DE DP NON CHARGES OU CHARGES
C==================================================
C#####        DP simples (et complexes)
C##################################################
C############################################
C############################################
C#####Fin de la lecture du fichier de donn�es
C############################################
C############################################
C##############################################################
C&&&&& A partir d'ici, le programme est constitu�  de la mise 
C&&&&& bout � bout des deux cas (i) sans et (ii) avec charges,
C&&&&& avec simple "branchement" � l'un ou l'autre cas.
C&&&&& Bien que cela corresponde � des doublons d'instructions
C&&&&& (e.g pour ALPHA_TYPE_COMPLEXE ou I_TYPE_NORMAL_S_R),
C&&&&& cela permet aussi les doublons d'allocations
C&&&&& (ce qui serait source d'erreur si ces allocations
C&&&&& �taient r�ellement effectu�es deux fois).
C##############################################################
C&&&&& Une am�lioration ult�rieure pourrait consister
C&&&&& � "factoriser"' ces doublons.
C##############################################################
C&&&&& Pour l'instant, seule la lecture du fichier de donn�es
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
C-----Branchement � la section "DP charg�s"
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
C-----Ecritures � l'�cran
C------------------------
      WRITE ( * , 1100 )
      WRITE ( * , 600 ) N_TYP
      WRITE ( * , 601 ) N_TYP - N_TYP_INTR
      WRITE ( * , 602 ) N_R
      WRITE ( * , 603 ) N_R_INTER
      WRITE ( * , * )
     $ 'Nombre de sites par maille pour chaque sous-r�seau : '
      WRITE ( * , * )
     $ ( P_R ( I_R ) , I_R = 1 , N_R )
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' )
     $ THEN
       WRITE ( * , * )
     $ '--------------------------------------------------------------'
        WRITE ( * , * )
     $ 'Option "d�fauts complexes" s�lectionn�e'
       WRITE ( * , * )
     $ '- - - - - - - - - - - - - - - - - - - -'
       WRITE ( * , * )
     $ '   * La pr�sence de DP complexes peut induire'
       WRITE ( * , * )
     $ '     certaines fractions atomiques irr�alistes (> 1).'
       WRITE ( * , * )
     $ '     Ceci provient de taux de DP complexes > 1'
      WRITE ( * , * )
     $ "     pour certains potentiels chimiques, mais n'emp�che pas"
      WRITE ( * , * )
     $ "     l'utilisation pratique de l'ADPI+complexes"
      WRITE ( * , * )
     $ "     qui reste valable autour de la stoechiom�trie."
       WRITE ( * , * )
     $ '   * Un comportement irr�aliste (xDPcmplx, Gat(x))'
       WRITE ( * , * )
     $ "     peut aussi survenir pour des xi compatibles avec l'ADPI,"
       WRITE ( * , * )
     $ "     e.g. pour xB = 0 dans FeAl-B avec complexes incluant B." 
       WRITE ( * , * )
     $ "     Ces valeurs xB = 0 sont invalides,"
       WRITE ( * , * )
     $ "     car elles proviennent d'un for�age � 0 par le programme"
       WRITE ( * , * )
     $ "     de fractions atomiques trouv�es initialement < 0,"
       WRITE ( * , * )
     $ "     � cause du terme X (HDR, �q. 2.133) de d�nominateur de xi."
       WRITE ( * , * )
     $ "     L� encore, l'utilisation des r�sultats reste possible,"
       WRITE ( * , * )
     $ "     moyennant l'�limination de ces points xB = 0 artificiels"
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
     $ 'Pas de d�fauts complexes pris en compte'
        WRITE ( * , 1100 )
       END IF
C-------------------------------
C-----Format variable d'�criture
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
C-----Potentiels chimiques des �l�ments intrins�ques
C---------------------------------------------------
      ALLOCATE ( POT_INTR ( N_TYP_INTR ) )
C-------------------------------------------------
C-----Potentiels chimiques des �l�ments d'addition
C-------------------------------------------------
      ALLOCATE ( POT_I ( N_TYP_INTR + 1 : N_TYP ) )
C--------------------------------------------
C-----Termes relatifs aux �l�ments d'addition
C--------------------------------------------
      ALLOCATE ( ALPHA_I_R ( N_TYP_INTR + 1 : N_TYP ) )
      ALLOCATE ( SOMME_I ( N_TYP_INTR + 1 : N_TYP ) )
C---------------------------------------------------------------------
C-----Calcul des nombres d'atomes de chaque esp�ce par maille
C-----(utiles au calcul de la somme pond�r�e des potentiels chimiques)
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
C-----Calcul de l'�nergie et du volume de r�f�rence par maille (eV)
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
C-----Calcul approch� de la somme pond�r�e des potentiels chimiques
C-----pour une maille � partir de la relation de Gibbs-Duhem
C------------------------------------------------------------------
      S_POT = V_REF_MAILLE * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
     $      + E_REF_MAILLE
C==================================================
C==================================================
C=====Calcul optionnel de grandeurs pr�liminaires
C=====relatives aux complexes �ventuels NON CHARGES
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
C-----Calcul du nombre de sites du sous-r�seau r dans chaque complexe
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
     $'       Param�tres u_d(s-r) et v_d(type) des complexes :'
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
C=====Fin du calcul optionnel de grandeurs pr�liminaires
C=====relatives aux complexes �ventuels NON CHARGES
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
C-----Param�tres g�n�raux
C------------------------ 
      WRITE ( 100 , 1100 )
      DO I_LIGNE = 1 , N_LIGNES_COMMENTAIRE
        WRITE ( 100 , * ) COMMENTAIRE ( I_LIGNE )
      END DO
      WRITE ( 100 , 1100 )
      WRITE ( 100 , * ) 'Nombre de types : '
      WRITE ( 100 , 502 ) N_TYP
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) "Nombre d'�l�ments d'addition : "
      WRITE ( 100 , 502 ) N_TYP - N_TYP_INTR
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * ) 'Nombre de sous-r�seaux :'
      WRITE ( 100 , 502 ) N_R
      WRITE ( 100 , 1200 )
      IF ( N_R_3 .GT. 0 ) THEN
        WRITE ( 100 , * )
     $ 'Nombres de sous-r�seaux 1, 2, 3 et interstitiels :'
        WRITE ( 100 , 504 ) N_R_1 , N_R_2 , N_R_3 , N_R_INTER
        WRITE ( 100 , 1200 )
      ELSE
        WRITE ( 100 , * )
     $ 'Nombres de sous-r�seaux 1, 2 et interstitiels :'
        WRITE ( 100 , 503 ) N_R_1 , N_R_2 , N_R_INTER
        WRITE ( 100 , 1200 )
      END IF
      WRITE ( 100 , * )
     $   'Nombre de sites par maille pour chacun des ' , N_R ,
     $   ' sous-r�seaux : '
      WRITE ( 100 , * )
     $ ( P_R ( I_R ) , I_R = 1 , N_R )
      WRITE ( 100 , 1200 )
C-------------------------------
C-----Formule du compos� parfait
C-------------------------------
      IF ( N_TYP_INTR .GE. 2 ) THEN
 2212 FORMAT ( 1X , 'A(' , F4.2 , ')B(' , F4.2 , ')' )
 2213 FORMAT ( 1X , 'A(' , F4.2 , ')B(' , F4.2 , ')C(' , F4.2 , ')' )
        WRITE ( 100 , * )
     $     'Formule du compos� stoechiom�trique : '
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
     $   'D�fauts complexes :'
      WRITE ( 100 , * )
     $   '-----------------'
      WRITE ( 100 , * )
     $   '     Nombre de complexes pris en compte :'
      WRITE ( 100 , * ) N_TYPES_COMPLEXES
      WRITE ( 100 , * ) '         Num�ro :'
      WRITE ( 100 , 705 )
     $ ( I_COMPLEXE ,
     $   I_COMPLEXE = 1 , N_TYPES_COMPLEXES )
      WRITE ( 100 , * ) '         Multiplicit� :'
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
     $ '            Informations compl�mentaires :'
        WRITE ( 100 , * )
     $ '            ----------------------------'
       DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
        WRITE ( 100 , * )
     $ '            ------------------'
        WRITE ( 100 , * )
     $ '            Num�ro du complexe :'
        WRITE ( 100 , * )
     $ '            ------------------'
        WRITE ( 100 , * ) I_COMPLEXE
        WRITE ( 100 , * ) 
     $ '            Nombre de sites par sous-r�seau :'
        WRITE ( 100 , * )
     $ ( U_COMPLEXE_S_R ( I_COMPLEXE , I_R ) , I_R = 1 , N_R ) 
        WRITE ( 100 , * )
     $ "            Nombre d'atomes par type chimique :"
        WRITE ( 100 , * )
     $ ( V_COMPLEXE_TYPE ( I_COMPLEXE , I_TYP ) , I_TYP = 1 , N_TYP )
       END DO
       ELSE
        WRITE ( 100 , * )
     $ 'Pas de d�fauts complexes pris en compte'
       END IF
      WRITE ( 100 , 1200 )
C---------------------------------------------
C-----Fin d'�criture optionnelle des complexes
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
      WRITE ( 100 , * ) 'Temp�rature (K) : '
      WRITE ( 100 , 501 ) TEMPERATURE
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * )
     $ "Somme pond�r�e des potentiels chimiques (eV) (fonction de P)"
      WRITE ( 100 , * )
     $ '(sur les �l�ments intrins�ques) :'
      WRITE ( 100 , 501 ) S_POT
      WRITE ( 100 , 1200 )
      WRITE ( 100 , * )
     $  'Caract�ristiques des s�ries de potentiels chimiques :'
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
     $  "El�ment d'addition " , I_TYP , ':'
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
      WRITE ( 100 , * ) 'Temp�rature (K) : '
      WRITE ( 100 , 501 ) TEMPERATURE
      WRITE ( 100 , 1200 )
      END IF
C-----------------------------------------------------------
C-----Suffixes de sous-r�seau et de type pour les DP simples
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
C-----Suffixe de num�ro pour les DP complexes (optionnel)
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
C=====Matrices relatives � l'inversion LU du syt�me d'�quations
C=====(de dimension N_TYP - 1 ) donnant les fractions atomiques
C=====des �l�ments autres que I_TYP_0 en fonction des contraintes
C=====(en nombre N_TYP - 2) et de la fraction atomique de I_TYP_0
C=====pour le point courant en potentiels chimiques
C================================================================
C----------------------------------------
C-----Dimension de ce syst�me d'�quations
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
C=====Ouverture des fichiers de r�sultats par sous-r�seau et par DP
C=====ainsi que des fichiers d'�nergie libre et de potentiels chimiques
C=====(sauf pour le cas NPT=0)
C======================================================================
      IF ( .NOT.
     $   ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT=0' ) )
     $ THEN
C--------------------------------
C-----Boucle sur les sous-r�seaux
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------------------------
C-----Boucle sur les types (plus 0 = lacune)
C-------------------------------------------
      DO I_TYP = 0 , N_TYP
C----------------------
C-----Num�ro de fichier
C----------------------
        I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'�criture de l'interstitiel
C-----pour le sous-r�seau et l'esp�ce donn�s
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
C-----ne correspondent pas � un DP)
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
          WRITE ( I_FICH , * ) 'Sous-r�seau = ' , I_R 
          WRITE ( I_FICH , * ) '---------------'
          WRITE ( I_FICH , * ) 'D�faut ponctuel = ' ,
     $                       W_TYP ( I_TYP )
     $                     ( 1 : LEN_TRIM ( W_TYP ( I_TYP ) ) )
          WRITE ( I_FICH , * ) '---------------'
C-------------------------------------------
C-----Titre sp�cifique suivant le type de DP
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
C-----Fin de la boucle sur les sous-r�seaux
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
          WRITE ( I_FICH , * ) 'D�faut complexe = ' ,
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
C-----Ouverture du fichier d'�nergie libre
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
C=====Rep�rage des d�fauts et de leurs quantit�s GC
C=====par un indice unique (utile pour NPT=0 et NPT)
C===================================================
C-----------------------------------------------
C-----Indice de chaque DP
C-----en fonction de son sous-r�seau de son type
C-----------------------------------------------
C-------------------------------------------------
C-----Cet indice est mis � z�ro �galement
C-----pour les DP "r�els" mais "non utiles",
C-----i.e. tels que "E_DP >= 0".
C-----Le nombre de DP manipul�s N_TYP_D_R est aussi
C-----r�duit d'une unit� � chaque fois
C-----(ce nombre fixe la dimension du syst�me NPT :
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
C=====Noms de DP en fonction des types et sous-r�seaux
C=====et longueurs correspondantes
C=====================================================
      ALLOCATE ( NOM_D_R_TYP ( 0 : N_TYP , N_R ) )
      ALLOCATE ( LONG_NOM_D_R_TYP ( 0 : N_TYP , N_R ) )
C--------------------------------------
C-----Boucles de type et de sous-r�seau
C--------------------------------------
      DO I_TYP = 0 , N_TYP
       DO I_R = 1 , N_R
C--------------------------------------------
C-----Indicateur d'�criture de l'interstitiel
C-----pour le sous-r�seau et l'esp�ce donn�s
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
C-----Fin des boucles de type et de sous-r�seau
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
C#####Calcul des �nergies, volumes et enthalpies GC
C#####des DP simples (et �ventuellement complexes)
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
C	write ( * , * ) 'Fin du calcul des quantit�s GC'
C----------------------------------------------------------------
C-----Calcul optionnel des quantit�s GC des complexes non charg�s
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
C	write ( * , * ) 'Fin du calcul des quantit�s GC(cmplx)'
       END IF
C==================================================
C=====Noms de DP en fonction de leur indice unique,
C=====longueurs correspondantes,
C=====sous-r�seaux et types,
C=====et �nergies, volumes et enthalpies GC
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
C-----ne correspondent � aucun DP ( I_TYP , I_R )
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
C=====Fin du rep�rage des d�fauts et de leurs quantit�s GC
C=====par un indice unique (utile pour NPT=0 et NPT)
C=========================================================
C--------------------------------------------------------
C-----Param�tres indicateurs de types et de sous-r�seaux
C-----utiles au calcul des termes de potentiels chimiques
C-----relatifs aux complexes �ventuels
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
C-----Type chimique normal de chaque sous-r�seau
C-----(0 pour sous-r�seaux interstitiels)
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
C-----Ecriture de v�rification
C-----------------------------
C       write ( * , * ) 't_0(r) : ' ,
C    $ ( I_TYPE_NORMAL_S_R ( I_R ) , I_R = 1 , N_R )
C-----------------------------------------------
C-----Ouverture des tableaux de fractions de DP
C-----et initialisation � 0.
C-----En mode NPT, ces valeurs seront inchang�es
C-----pour les DP "non utiles", i.e. ceux pour
C-----lesquels est sp�cifi� "E_DP >= 0"
C-----(en particulier, certains interstitiels)
C-----------------------------------------------
      ALLOCATE ( X_D_R ( 0 : N_TYP , N_R ) )
      X_D_R = 0.D0
C--------------------------------------------------
C-----Ouverture des tableaux de complexes �ventuels
C--------------------------------------------------
C     write ( * , * ) 'N_TYPES_COMPLEXES = ' , N_TYPES_COMPLEXES
      ALLOCATE ( X_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
      ALLOCATE ( H_FORM_D_COMPLEXE ( N_TYPES_COMPLEXES ) )
      ALLOCATE ( SOMME_COMPLEXE_V_TYPE_I ( N_TYP_INTR + 1 : N_TYP ) )
      X_D_COMPLEXE = 0.D0
C#################################################
C#################################################
C#####Traitement du cas muVT (DP et �nergie libre)
C#################################################
C#################################################
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT' )
     $ THEN
C	write ( *  , * ) 'D�but traitement muVT'
C----------------------------------------------------
C-----Nombre total de pas sur les �l�ments d'addition
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
C-----Initialisation du nombre de points en composition �crits
C-----(contenus dans les fen�tres choisies),
C-----des valeurs moyennes de mu(1) - mu(2) et mu(1) - mu(3)
C-----ainsi que de leurs carr�s
C-----et idem pour les �l�ments d'addition
C-----(calcul de l'�cart-type) pour ces points
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
C	write ( * , * ) 'D�but des boucles muVT'
C-----------------------------------------------------------
C-----Ecriture concernant l'autocoh�rence �ventuelle sur mu1
C-----------------------------------------------------------
      IF ( N_ITER_MAX_MU_1 .EQ. 1 ) THEN
       write(*,*)
       write(*,*)
     $ "*** Vous avez choisi N_ITER_MAX_MU_1 = 1"
       write(*,*) " ==> pas d'autocoh�rence sur mu1 ***"
       write(*,*)
      ELSE
       write(*,*)
       write(*,*)
     $ "*** Vous avez choisi N_ITER_MAX_MU_1 > 1"
       write(*,*) " ==> boucle d'autocoh�rence (BA) sur mu1 ***"
       write(*,*) "  NB : la mise � jour de POT_1 est faite"
       write(*,*) "       en fin d'it�ration de BA,"
       write(*,*) "       et l'it�r. ITER_MU_1 = 1 correspond"
       write(*,*) "       aux propri�t�s avant BA"
       write(*,*) "       (e.g. pour x_at et G_AT �crits � l'�cran)."
       write(*,*) "  NB2 : les points du balayage en pot. chimiques"
       write(*,*) "        pour lesquels la BA a rencontr�"
       write(*,*) "        une difficult� de convergence sont �crits"
       write(*,*) "        eux aussi dans les fichiers de sortie."
       write(*,*) "        Ils sont �galement indiqu�s dans le fichier"
       write(*,*) '        "pb_conv_BA", en distinguant'
       write(*,*)
     $ " * nombre maxi. d'it�rations atteint (INDIC_PB_CONV_BA=1)"
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
C-----pour lesquels la BA ne s'est pas d�roul�e convenablement
C-------------------------------------------------------------
      NB_POINTS_PB_CONV_BA = 0
C###################################################
C#####D�but des boucles sur les potentiels chimiques
C###################################################
C---------------------------------------
C-----D�but du balayage en mu(1) - mu(2)
C---------------------------------------
        I_POUR_CENT = 0
        DO K_D_POT_2_1 = 1 , N_D_POT_2_1
C----------------------
C-----Balayage lin�aire
C----------------------
          D_POT_2_1 = D_POT_INIT_2_1
     $              + PAS_D_POT_2_1 * DFLOAT ( K_D_POT_2_1 )
C---------------------------
C-----Balayage logarithmique
C---------------------------
C         D_POT_2_1 = D_POT_INIT_2_1
C    $          + PAS_D_POT_2_1 * K_T * DLOG ( DFLOAT ( K_D_POT_2_1 ) )
C---------------------------------------
C-----D�but du balayage en mu(1) - mu(3)
C---------------------------------------
        DO K_D_POT_3_1 = 1 , N_D_POT_3_1
C----------------------
C-----Balayage lin�aire
C----------------------
          D_POT_3_1 = D_POT_INIT_3_1
     $              + PAS_D_POT_3_1 * DFLOAT ( K_D_POT_3_1 )
C-----------------------------------------------------------
C-----Calcul des potentiels chimiques des �l�ments 1, 2 et 3
C-----Pour 1, qui est l'�l�ment de r�f�rence, il s'agit de
C-----l'estimation initiale, avant boucle d'autocoh�rence.
C-----------------------------------------------------------
             POT_1_INIT = ( S_POT
     $                  + DFLOAT ( N_2_MAILLE ) * D_POT_2_1
     $                  + DFLOAT ( N_3_MAILLE ) * D_POT_3_1 )
     $             / DFLOAT ( N_1_MAILLE + N_2_MAILLE + N_3_MAILLE )
             POT_2 = POT_1_INIT - D_POT_2_1
             POT_3 = POT_1_INIT - D_POT_3_1
C---------------------------------------------------------------------
C---------------------------------------------------------------------
C-----D�but optionnel du balayage en mu ( i > N_TYP_INTR ) (additions)
C---------------------------------------------------------------------
C---------------------------------------------------------------------
        DO K_TOT_POT_I = 1 , N_TOT_POT_I
C----------------------------------------------------------------
C-----Indices des �l�ments d'addition � partir de l'indice unique
C-----et potentiels chimiques correspondants
C-----(simulation de plusieurs boucles de potentiels chimiques
C-----imbriqu�es � l'aide d'un seul indice K_TOT_POT_I)
C----------------------------------------------------------------
         J_COUR = K_TOT_POT_I - 1
C-------------------------------
C-----El�ments >= N_TYP_INTR + 2
C-------------------------------
         K_COUR_POT_I = MOD ( K_TOT_POT_I - 1 , N_TOT_POT_I )
         DO I_TYP = N_TYP , N_TYP_INTR + 2 , - 1
            PROD_PART_COUR = PROD_PART_N_POT_I ( I_TYP - 1 )
            K_COUR_POT_I
     $   =  J_COUR / PROD_PART_COUR + 1
            J_COUR = MOD ( J_COUR , PROD_PART_COUR ) 
C----------------------
C-----Balayage lin�aire
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
C-----El�ment N_TYP_INTR + 1
C---------------------------
       IF ( N_TYP .GT. N_TYP_INTR ) THEN
            K_COUR_POT_I = J_COUR + 1
C----------------------
C-----Balayage lin�aire
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
C-----Pr�cision et nombre maximum d'it�rations
C-----pour l'arr�t de la boucle d'autocoh�rence sur POT_1
C-----(lus dans le fichier DATA.adpi)
C--------------------------------------------------------
C      PRECISION_MU_1 = 1.D-8
C      N_ITER_MAX_MU_1 = 100
C-----Ligne ci-dessous uniquement pour test de la BA :
C-----lorsque celle-ci est divergente, cela est-il d�
C-----� la valeur initiale de POT_1, attribu�e syst�matiquement
C-----via l'estimation avec S_POT (cf. supra) ?
C-----d'apr�s les premiers tests sur Al_cfc(B,Ti),
C-----syst�me pour lequel S_POT => POT_1_INIT -3,7477 eV,
C-----la modification "manuelle" de POT_1_INIT comme ci-dessous
C-----(POT_1_INIT entre  -3,60 et -3,85 eV) ne change pas
C-----le comportement qui reste divergent.
C      POT_1_INIT = -3.60
C-------------------------------------------------------------
C-----Initialisation du potentiel chimique de r�f�rence POT_1,
C-----du nombre d'it�rations,
C-----de la valeur sauvegard�e (it�r. BA pr�c�dente) de POT_1,
C-----et de l'indicateur de "probl�me de convergence de la BA"
C-----avant la boucle d'autocoh�rence (BA)
C-------------------------------------------------------------
       ITER_MU_1 = 0
       POT_1 = POT_1_INIT
       POT_1_PREC = 1.D100
       INDIC_PB_CONV_BA = 0
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C&&&&&Boucle d'autocoh�rence sur POT_1
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
     $ "### Avant boucle d'autocoh�rence sur mu1 : " ,
     $         ' POT_1(initial) = ' , POT_1
         write(*,*) '----------'
        END IF
        write(*,*)
     $ '                ***** ITER_MU_1 = ' , ITER_MU_1 , ' *****'
       END IF
C--------------------------------------------------------
C--------------------------------------------------------
C-----Termes cumul�s sur les sous-r�seaux
C-----pour le calcul de la fraction atomique de l'alliage
C--------------------------------------------------------
C--------------------------------------------------------
C-------------------------------------------------
C-----Termes au num�rateur des fractions atomiques
C-------------------------------------------------
          SOMME_1_TYP_2 = 0.D0
          SOMME_2_TYP_2 = 0.D0
          SOMME_1_TYP_3 = 0.D0
          SOMME_2_TYP_3 = 0.D0
          SOMME_I = 0.D0
C---------------------------------------------------
C-----Termes au d�nominateur des fractions atomiques
C---------------------------------------------------
          SOMME_0 = 0.D0
          SOMME_MAILLE = 0.D0
C==================================
C==================================
C=====Boucle sur les sous-r�seaux 1
C==================================
C==================================
          DO I_R = 1 , N_R_1
C	     write ( * , * ) 'I_R = ' , I_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------
C-----Lacunes sur les sous-r�seaux 1
C-----------------------------------
           ALPHA_0_R
     $   = DEXP ( - ( H_GC_D_R ( 0 , I_R )
     $              + POT_1 ) / K_T )
C----------------------------------------------------------
C-----Antisites �ventuels 2(r), 3(r) sur les sous-r�seaux 1
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
C-----D�nominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_2_R + ALPHA_3_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
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
C=====Fractions de d�fauts ponctuels
C=====et �nergies de formation de ces DP
C=======================================
C------------
C-----Lacunes
C------------
           X_D_R ( 0 , I_R ) = ALPHA_0_R / DENOMINATEUR_X_D_R
           IF ( X_D_R ( 0 , I_R ) .LE. 0.D0 )
     $          X_D_R ( 0 , I_R ) = 1.D-100
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_1
C-----Rappel : les termes "SOMME" sont �galement calcul�s
C-----car utiles plus bas au calcul des fractions atomiques em muVT
           SOMME_0
     $   = SOMME_0
     $   + DFLOAT ( P_R ( I_R ) )
     $   * ( 1.D0 - X_D_R ( 0 , I_R ) )
C--------------------------
C-----Antisites 2 �ventuels
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Boucle sur les sous-r�seaux 2
C=====(lorsqu'ils existent)
C==================================
C==================================
        DO I_R = 1 + N_R_1 , N_R_2 + N_R_1
C	     write ( * , * ) 'I_R = ' , I_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------------------------
C-----Lacunes et antisites 1(r) sur les sous-r�seaux 2
C-----------------------------------------------------
           ALPHA_0_R
     $   = DEXP ( - ( H_GC_D_R ( 0 , I_R )
     $              + POT_2 ) / K_T )
           ALPHA_1_R
     $   = DEXP ( - ( H_GC_D_R ( 1 , I_R )
     $              - POT_1 + POT_2 ) / K_T )
C----------------------------------------------------
C-----Antisites �ventuels 3(r) sur les sous-r�seaux 2
C----------------------------------------------------
       ALPHA_3_R = 0.D0
       IF ( N_TYP_INTR .GT. 2 ) THEN
           ALPHA_3_R
     $   = DEXP ( - ( H_GC_D_R ( 3 , I_R )
     $              - POT_3 + POT_2 ) / K_T )
        END IF
C-----------------
C-----D�nominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_3_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
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
C=====Fractions de d�fauts ponctuels
C=====et �nergies de formation de ces DP
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
C-----Antisites 3 �ventuels
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
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
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
C=====Boucle sur les sous-r�seaux 3
C=====(lorsqu'ils existent)
C==================================
C==================================
        DO I_R = 1 + N_R_1 + N_R_2 , N_R_3 + N_R_2 + N_R_1
C=======================
C=====Coefficients alpha
C=======================
C-----------------------------------------------------------
C-----Antisites 1(r), 2(r) et lacunes sur les sous-r�seaux 3
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
C-----D�nominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_2_R + ALPHA_0_R
C-------------------------------------------------
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
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
C=====Fractions de d�fauts ponctuels
C=====et �nergies de formation de ces DP
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
C-----Traitement �ventuel des esp�ces > N_TYP_INTR 
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
C-----Terme de somme 1 pour l'esp�ce 3
C-------------------------------------
         SOMME_1_TYP_3 = SOMME_1_TYP_3
     $                 + SOMME_1_TYP_3_R * DFLOAT ( P_R ( I_R ) )
        END DO
C==========================================================
C==========================================================
C=====Boucle optionnelle sur les sous-r�seaux interstitiels
C==========================================================
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C=======================
C=====Coefficients alpha
C=======================
C-----------------
C-----D�fauts 1(r)
C-----------------
           ALPHA_1_R
     $   = DEXP ( - ( H_GC_D_R ( 1 , I_R ) - POT_1 ) / K_T )
C-----------------------------------
C-----D�fauts �ventuels 2(r) et 3(r)
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
C-----D�nominateur
C-----------------
           DENOMINATEUR_X_D_R
     $   = 1.D0 + ALPHA_1_R + ALPHA_2_R + ALPHA_3_R
C-------------------------------------------------
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
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
C=====Fractions de d�fauts ponctuels
C=====et �nergies de formation de ces DP
C=======================================
C---------------------------------------------
C-----El�ments 1, 2 et 3 (ce dernier �ventuel)
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
C-----El�ments d'addition �ventuels
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
C-----Fin de la boucle sur les sous-r�seaux interstitiels
C========================================================
C========================================================
        END DO
C----------------------------------------------------------
C-----Fin du test d'existence de sous-r�seaux interstitiels
C----------------------------------------------------------
       END IF
C---------------------------------------------
C-----Initialisation des termes de sommes
C-----utiles au calcul des fractions atomiques
C-----(toujours effectu�e, sinon erreur
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
C=====Calcul optionnel des quantit�s de d�fauts complexes
C=====(rapport�es au sous-r�seau sp�cifi� pour la multiplicit�)
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
C-----Calcul du potentiel chimique de l'�l�ment sur le site courant
C-----et du potentiel chimique "normal" sur le site courant
C------------------------------------------------------------------
C--------------------------------
C-----Cas 1 : un type intrins�que
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
C-----Cas 2 : deux types intrins�ques
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
C-----Cas 3 : trois types intrins�ques
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
C-----Ecriture de v�rification
C-----------------------------
C      if ( I_COMPLEXE .ge. 6 ) then
C         write ( * , * ) 'Complexe ' , I_COMPLEXE
C	  write ( * , * ) 'site : ' , I_SITE
C	  write ( * , * ) 'sous-r�seau : ' , I_S_R_COUR
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
C-----Ecriture de v�rification
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
C-----Facteur interm�diaire
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
C-----termes de sommes relatifs � l'esp�ce 2 intrins�que
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
C-----termes de sommes relatifs � l'esp�ce 3 intrins�que
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
C-----termes de sommes relatifs aux �l�ments d'addition
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
C-----termes de sommes relatifs aux quantit�s de mati�re totales
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
C-----Fractions atomiques des �l�ments d'addition  
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
C=====Calcul des quantit�s thermodynamiques dans l'ADPI           
C=====(�nergie et volume par maille, entropie de configuration, 
C=====�nergie libre par maille, quantit�s par atome)         
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
C-----Avant mise � jour de POT_1,
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
C-----(=> elle pourrait n'�tre faite que si N_ITER_MAX_MU_1 > 1)
C----------------------------------------------------------------
       POT_1_PREC = POT_1
C--------------------------------------------------------
C-----Si la BA est activ�e (i.e. si N_ITER_MAX_MU_1 > 1),
C-----mise � jour de POT_1, gr�ce � g_at=Somme(x_i*mu_i)
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
     $ " avant mise � jour de POT_1"
       write ( * , * ) 'ITER_MU_1 = ' , ITER_MU_1 , ' x_at : ' ,
     $ ( X_AT ( I_TYP ) , I_TYP = 1 , N_TYP )
       write(*,*) 'ITER_MU_1 = ' , ITER_MU_1 , ' G_AT = ' , G_AT
       write(*,*) 'ITER_MU_1 = ' , ITER_MU_1 ,
     $           ' POT_1 = ' , POT_1 , 'apr�s mise � jour (BA) via G_AT'
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
C&&&&&Fin de la boucle d'autocoh�rence sur POT_1
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      END DO
      IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
       write(*,*)
     $ "###### Fin de la boucle d'autocoh�rence sur mu_1 apr�s " ,
     $ ITER_MU_1 , " it�rations ######"
       write(*,*)
      END IF
C====================================================================
C====================================================================
C=====A partir des contraintes sur la composition
C=====et de la fraction atomique calcul�e pr�c�demment
C=====pour l'�l�ment sp�cifi� I_TYP_0,
C=====inversion LU du syst�me d'�quations de contraintes
C=====pour obtenir les valeurs des autres fractions atomiques
C=====autour desquelles doivent �tre centr�es les fen�tres d'�criture
C=====(toujours utilis� en NPT, seulement si N_TYP > 2 en muVT)
C====================================================================
C====================================================================
         IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT'
     $ .OR. ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'muVT'
     $  .AND. N_TYP .GT. 2 ) ) THEN
C----------------------------------------------------
C----------------------------------------------------
C-----Boucle sur les lignes du syst�me de contraintes
C----------------------------------------------------
C----------------------------------------------------
      DO I_TYP = 1 , N_SYS_CTR
C---------------------------------------------------------
C-----Boucle sur les colonnes du syst�me de contraintes
C-----en sautant la colonne I_TYP_0
C-----(�l�ment dont on �tudie l'effet de l'enrichissement)
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
C-----Fin de la boucle sur les lignes du syst�me de contraintes
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
C-----Inversion LU du syst�me de contraintes
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
C-----(pour l'�criture dans les fichiers de r�sultats)
C-----d�duites de la valeur courante pour l'�l�ment sp�cifi�
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
C-----Fin du test d'utilisation du syst�me de contraintes
C--------------------------------------------------------
       END IF
C===============================================================
C===============================================================
C=====Fin de l'utilisation optionnelle du syst�me de contraintes
C===============================================================
C===============================================================
C==================================================================
C==================================================================
C=====Ecriture des quantit�s de d�fauts ponctuels
C=====dans les fichiers de r�sultats
C=====par sous-r�seau r et par d�faut ponctuel :
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
C-----Seules sont �crites les valeurs de compositions
C-----comprises dans les fen�tres sp�cifi�es
C-----si N_TYP > 2
C----- -> indicateur de ce que le point de compo.  est
C----- dans la fen�tre (toujours vrai si 2 types) 
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
C-----Boucle sur les sous-r�seaux
C--------------------------------
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------
C-----Boucle sur les types
C-------------------------
        DO I_TYP = 0 , N_TYP
C----------------------
C-----Num�ro du fichier
C----------------------
          I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'�criture de l'interstitiel
C-----pour le sous-r�seau et l'esp�ce donn�s
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
C-----ne correspondent pas � un DP
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
C-----Potentiels chimiques des �l�ments intrins�ques
C---------------------------------------------------
        POT_INTR ( 1 ) = POT_1
         IF ( N _TYP_INTR .GT. 1 ) THEN
          POT_INTR ( 2 ) = POT_2
          IF ( N _TYP_INTR .GT. 2 ) THEN
           POT_INTR ( 3 ) = POT_3
          END IF
        END IF
C------------------------------------------------------------------
C-----Cas binaire (tous types confondus : intrins�ques + additions)
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
C-----(tous types confondus : intrins�ques + additions)
C----- -> fen�tres de composition �ventuelles
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
C-----Fin du test d'existence de d�faut
C--------------------------------------
        END IF
C----------------------------------------------
C----------------------------------------------
C-----Fin des boucles sur types et sous-r�seaux
C----------------------------------------------
C----------------------------------------------
       END DO
      END DO
C===================================================
C===================================================
C=====Ecriture �ventuelle des quantit�s de complexes
C===================================================
C===================================================
       IF ( INDIC_COMPLEXES .EQ. 'O' .OR. INDIC_COMPLEXES .EQ. 'o' ) 
     $ THEN
C-----------------------------
C-----Boucle sur les complexes
C-----------------------------
        DO I_COMPLEXE = 1 , N_TYPES_COMPLEXES
C----------------------
C-----Num�ro du fichier
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
C-----Cas ternaire et plus (fen�tres de composition)
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
C=====Ecriture de l'�nergie, du volume et de l'�nergie libre,
C=====en distinguant le cas binaire (fen�tres de composition sinon)
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
C-----Cas d'�criture de l'�nergie libre par maille
C-------------------------------------------------
            IF 
     $    ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm' )
     $      THEN
             WRITE ( 110 , CAR_COL_VAR_E_L )
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       E_MAILLE , V_MAILLE , S_CONF_MAILLE , G_MAILLE
C------------------------------------------------
C-----Cas d'�criture de l'�nergie libre par atome
C------------------------------------------------
            ELSE
C- - - - - - - - - - - - - - - - - - - - - - -
C- - -Dans le cas G/atome, deux possibilit�s :
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
C-----Fin du test sur "N_TYP et fen�tres"
C----------------------------------------
       END IF
C==================================================================
C==================================================================
C=====Ecriture des potentiels chimiques
C=====en distinguant le cas binaire (fen�tres de composition sinon)
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
C-----Cas ternaire et plus (fen�tres de composition)
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
C=====En cas de boucle d'autocoh�rence,
C=====�criture dans le fichier "pb_conv_BA" des points pour lesquels
C=====la BA ne s'est pas effectu�e convenablement
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
C-----Pourcentage du calcul effectu�
C-----et �criture correspondante
C-----si un seul type intrins�que
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
     $ ' % du calcul muVT effectu�s        # # # # #'
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
C-----Pourcentage du calcul effectu�
C-----et �criture correspondante
C-----si au moins deux types intrins�ques
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
     $ ' % du calcul muVT effectu�s        # # # # #'
          END IF
        END IF
C=======================================
C=====Fin de la boucle sur mu(1) - mu(2)
C=======================================
      END DO
C-------------------------------------------------------------------
C-----En cas de boucle d'autocoh�rence, �criture du nombre de points
C-----pours lesquels la BA ne s'est pas d�roul�e convenablement
C-------------------------------------------------------------------
       IF ( N_ITER_MAX_MU_1 .GT. 1 ) THEN
        write(*,*)
     $ "   ----------------------------------------------------------"
        write(*,*)
     $ "      Boucle d'autocoh�rence sur mu_1 : nombre de points"
        write(*,*)
     $ "     pour lesquels la BA a rencontr� un pb. = " ,
     $  NB_POINTS_PB_CONV_BA 
        write(*,*)
     $ '   Ces points sont r�pertori�s dans le fichier ""pb_conv_BA".'
        write(*,*)
     $ "   ----------------------------------------------------------"
       END IF
C-----------------------------------------------------------------------
C-----Valeurs moyennes de mu(1) - mu(2), mu(1) - mu(3)
C-----et de leurs carr�s (calcul de l'�cart-type) pour ces points
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
        WRITE ( 100 , * ) 'Nombre de points �crits :'
        WRITE ( 100 , * ) N_POINTS
C--------------------------------------
C-----Ecriture dans le fichier de liste
C--------------------------------------
        IF ( N_TYP_INTR .GE. 2 ) THEN
          WRITE ( 100 , 1200 )
          WRITE ( 100 , * )
     $   'Pour les ' , N_POINTS , 'points �crits :'
          WRITE ( 100 , * )
     $   'valeur moyenne de mu(1) - mu(2) = ' , D_POT_MOY_2_1 , 'eV'
          WRITE ( 100 , * )
     $   '�cart-type = ' , DSQRT ( VAR_D_POT_2_1 )
         IF ( N_TYP_INTR .GT. 2 ) THEN
          WRITE ( 100 , * )
     $   'valeur moyenne de mu(1) - mu(3) = ' , D_POT_MOY_3_1 , 'eV'
          WRITE ( 100 , * )
     $   '�cart-type = ' , DSQRT ( VAR_D_POT_3_1 )
         END IF
        END IF
C-----------------------
C-----Ecriture � l'�cran
C-----------------------
        IF ( N_TYP_INTR .GE. 2 ) THEN
          WRITE ( * , 1100 ) 
          WRITE ( * , * )
     $   'Pour les ' , N_POINTS , 'points �crits :'
          WRITE ( * , * )
     $   'valeur moyenne de mu(1) - mu(2) = ' , D_POT_MOY_2_1 , 'eV'
          WRITE ( * , * )
     $   '�cart-type = ' , DSQRT ( VAR_D_POT_2_1 )
         IF ( N_TYP_INTR .GT. 2 ) THEN
          WRITE ( * , * )
     $   'valeur moyenne de mu(1) - mu(3) = ' , D_POT_MOY_3_1 , 'eV'
          WRITE ( * , * )
     $   '�cart-type = ' , DSQRT ( VAR_D_POT_3_1 )
         END IF
        END IF
C-------------------------------------------
C-----Ecritures pour les �l�ments d'addition
C-------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
C	 IF ( I_TYP .NE. I_TYP_0 ) THEN
          WRITE ( * , * )
     $   'Pour les ' , N_POINTS , 'points �crits :'
          WRITE ( 100 , * )
     $   "valeur moyenne de mu pour l'�l�ment d'addition " , I_TYP ,
     $   '= ' , POT_I_MOY ( I_TYP )
        WRITE ( 100 , * )
     $   '�cart-type = ' , DSQRT ( VAR_POT_I ( I_TYP ) )
          WRITE ( * , * )
     $   "valeur moyenne de mu pour l'�l�ment d'addition " , I_TYP ,
     $   '= ' , POT_I_MOY ( I_TYP )
        WRITE ( * , * )
     $   '�cart-type = ' , DSQRT ( VAR_POT_I ( I_TYP ) )
C   	 END IF
        END DO
        WRITE ( 100 , * )
     $   'Utiliser ces valeurs pour affiner les s�ries �ventuelles'
        WRITE ( 100 , * )
     $   'mu(1) - mu(2), mu(1) - mu(3) et mu(additions)'
        WRITE ( 100 , 1100 )
        WRITE ( * , * )
     $   'Utiliser ces valeurs pour affiner les s�ries �ventuelles'
        WRITE ( * , * )
     $   'mu(1) - mu(2), mu(1) - mu(3) et mu(additions)'
        WRITE ( * , 1100 )
C--------------------------------------
C-----Cas o� aucun point n'a �t� trouv�
C--------------------------------------
      ELSE
        WRITE ( * , * )
     $   '------------------------------------------------------------'
        WRITE ( * , * )
     $   "Aucun point n'a �t� trouv�. Vous pouvez :"
        WRITE ( * , * )
     $   '* modifier les propri�t�s des s�ries de potentiels chimiques'
        WRITE ( * , * )
     $   '* choisir des fen�tres plus larges'
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
C-----Ouverture de tableaux de r�sultats du simplexe
C---------------------------------------------------
      ALLOCATE ( I_1_SMPLX ( N_TYP_D_R + 1 ) )
      ALLOCATE ( I_2_SMPLX ( N_TYP ) )
C--------------------------------------------------------------------
C-----Tableau des indices initiaux des DP class�s par ordre croissant
C--------------------------------------------------------------------
      ALLOCATE ( IND_INIT ( N_TYP ) )
C----------------------------------------
C-----Simplexe pr�c�dent (balayage NPT=0)
C----------------------------------------
      ALLOCATE ( I_2_SMPLX_PREC ( N_TYP ) )
C--------------------------------------------------
C-----Types et sous-r�seaux des DP constitutionnels
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
C-----Tableau des potentiels chimiques � T = 0 K
C-----------------------------------------------
      ALLOCATE ( POT_CHIM_0K ( N_TYP ) )
C===============================================================
C=====S�quence 1 d'op�rations relatives au cas de calcul "point"
C===============================================================
       IF ( INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'P'
     $ .OR. INDIC_TYP_CALC_NPT ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'p'
     $    ) THEN 
C-------------------------------------------------------
C-----Relecture dans DATA.adpi des fractions atomiques
C-----sous forme d'une cha�ne de caract�res
C-----(pour utilisation comme suffixe de nom de fichier)
C-------------------------------------------------------
C---------------------------
C-----Recherche de l'en-t�te
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
C-----sous forme de cha�nes de caract�res
C----------------------------------------
         CALL
     $   ANALYSE_LIGNE
     $ ( LIGNE_X_AT ,
     $   N_CHAINES_X_AT ,
     $   I_DEBUT , I_FIN ,
     $   LONG_CHAINE_X_AT , CHAINE_X_AT )
C        write ( * , * ) 'La ligne ' ,
C    $   ligne_x_at ( 1 : len_trim (ligne_x_at ) )
C        write ( * , * ) 'contient ' , N_CHAINES_x_at , 'cha�nes'
C        do i = 1 , N_CHAINES_X_AT
C         write ( * , * ) 'cha�ne, d�but, fin, longueur, contenu' ,
C    $                   i, I_DEBUT ( i ) , I_FIN ( i ) ,
C    $                   LONG_CHAINE_X_AT ( I ) ,
C    $            CHAINE_X_AT ( I ) ( 1 : LONG_CHAINE_X_AT ( I ) )
C        end do
C-----------------
C-----V�rification
C-----------------
      IF ( N_CHAINES_X_AT .NE. N_TYP ) THEN
        WRITE ( * , * ) '------------------------------------------'
        WRITE ( * , * ) 'Nombre de types chimiques = ' , N_TYP
        WRITE ( * , * ) '=> nombre de fractions atomiques incorrect'
        WRITE ( * , * ) "(v�rifier l'absence de tabulation)"
        WRITE ( * , * ) '------------------------------------------'
        CALL INTERRUPTION
      END IF
C-----------------------------------------------------
C-----Constitution d'une cha�ne unique de suffixe X_AT
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
C=====Fin de s�quence 1 d'op�rations relatives au cas de calcul "point"
C======================================================================
       ELSE
C====================================================================
C=====S�quence 1 d'op�rations relatives au cas de calcul "balayage x"
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
     $ 'Identification des DP (rotation dans le sens trigonom�trique)' )
  435  FORMAT ( 1X ,
     $ 24 ( '-' ) ,  'Limites de zones de DP' , 24 ( '-' ) ,
     $ '|' , 6 ( '-' ) ,
     $ 'DP et pot. chim. juste "apr�s" ces limites' , 6 ( '-' ) )
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
     $ "D�sirez-vous un fichier annexe contenant les compositions"
       WRITE ( * , * ) 
     $ "pour lesquelles la minimisation a �chou�"
       WRITE ( * , * ) 
     $ "(o pour oui, autre caract�re pour non) ?"
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
C-----pour lequel (i) soit la fonction est non born�e,
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
     $      I_POUR_CENT , ' % du calcul NPT=0 effectu�s'
          END IF
C-------------------------------------------------------
C-------------------------------------------------------
C-----Boucle sur theta (jusqu'� N_PAS_THETA + 1
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
C-----V�rification de la validit� de la composition
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
C#####(i) la fonction � minimiser en premi�re ligne
C#####(i) les contraintes dans les lignes 2 � N_TYP + 1
C#####(i) la fonction auxiliaire dans la ligne N_TYP + 2
C#####Ce tableau contient N_TYP_D_R + 2 colonnes :
C#####la premi�re colonne contient les valeurs des contraintes,
C#####les autres colonnes contiennent les coefficients des contraintes
C#####################################################################
C------------------------------------------
C-----Initialisation du tableau du simplexe
C------------------------------------------
      TAB_SMPLX_D_R = 0.D0
      TAB_SMPLX_D_R ( 1 , 2 ) = H_REF_MAILLE
C=================================
C=====Boucles sur les sous-r�seaux
C=================================
      DO I_R = 1 , N_R
C---------------------------------------------------------
C-----Partie du tableau correspondant au nombre de mailles
C-----(deuxi�me colonne)
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
C-----El�ment 2 �ventuel
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
C-----El�ment 3 �ventuel
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
C-----El�ments d'addition �ventuels
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
C-----pour le sous-r�seau et le type courant
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
C=====Remplissage de la premi�re ligne (fonction H)
C==================================================
         TAB_SMPLX_D_R ( 1 , IND_D_R_TYP ( I_TYP , I_R ) + 2 )
     $ = H_GC_D_R ( I_TYP , I_R )
C=========================================
C=====Remplissage des lignes 2 � N_TYP + 1
C=====(contraintes sous forme "f = 0")
C=========================================
C-----------------------------------------
C-----Ligne de la quantit� totale d'atomes
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
C-----Ligne �ventuelle de l'�l�ment 2
C-----(si au moins 2 types intrins�ques)
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
C-----Ligne de l'�l�ment 3 �ventuel
C-----(si 3 types intrins�ques)
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
C-----Lignes des �l�ments d'addition �ventuels
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
C-----Fin de la boucle sur les �l�ments d'addition
C-------------------------------------------------
         END DO
C---------------------------------- 
C-----Fin du test d'existence de DP
C---------------------------------- 
         END IF
C----------------------------------------------------------
C-----Fin des boucles de sous-r�seaux et de types chimiques
C----------------------------------------------------------
        END DO
      END DO
C==============================================
C=====Fin du remplissage du tableau du simplexe
C==============================================
C-----------------------------------------------------------
C-----La proc�dure du simplexe recherchant le maximum
C-----de la fonction lue en premi�re ligne de TAB_SMPLX_D_R,
C-----on transforme le probl�me du minimum d'enthalpie
C-----en maximum (- enthalpie)
C-----------------------------------------------------------
      DO I_VAR = 1 , N_TYP_D_R + 2
         TAB_SMPLX_D_R ( 1 , I_VAR ) =  - TAB_SMPLX_D_R ( 1 , I_VAR )
      END DO
C--------------------------------------------------------------
C-----Proc�dure du simplexe :
C-----il y a z�ro contraintes du type M1 ( < ) ou M2 ( > )
C-----et N_TYP contraintes du type M3 ( = )
C-----Les arguments de la proc�dures sont, dans l'ordre :
C----- * le tableau du simplexe ;
C----- * les nombres N_E de contraintes et N_VAR de variables ;
C----- * les nombres de lignes et de colonnes du tableau
C-----   (respectivement N_E + 2 , N_VAR + 1)
C----- * les nombres de contraintes du type (<),  (>) et (=)
C----- * l'indicateur d'absence de solution (IC)
C-----   et les deux tableaux de r�sultats
C--------------------------------------------------------------
         CALL
     $   SIMPLX
     $ ( TAB_SMPLX_D_R ,
     $   N_TYP , N_TYP_D_R + 1 ,
     $   N_TYP + 2 , N_TYP_D_R + 2 ,
     $   0 , 0 , N_TYP ,
     $   IC , I_1_SMPLX , I_2_SMPLX )
C-----------------------------------------------------------------------
C-----Une erreur de minimisation (fonction non born�e, pas de solution)
C-----est retenue pour �criture d'un avertissement en fin de calcul.
C-----Remarque : seuls les points de composition corrects sont concern�s
C-----------------------------------------------------------------------
         IF ( IC .NE. 0 .AND. I_COMPO .EQ. 1 ) THEN
           I_ERREUR_SIMPLEXE = 1
         END IF
C-----------------------------------------------------------------
C-----En mode "balayage x", comparaison du point avec le pr�c�dent
C-----(comparaison sur theta pour un m�me r)
C-----en excluant le premier point en theta,
C-----car il correspond � r diff�rent
C-----Mise � jour de I_2_SMPLX_PREC
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
C-----Classement par ordre croissant des indices des d�fauts trouv�s
C-------------------------------------------------------------------
         CALL
     $   TRI_CROISSANT_ENTIER
     $ ( N_TYP , I_2_SMPLX , IND_INIT )
C----------------------------------------------------------------
C-----Recherche des types et sous-r�seaux des DP constitutionnels
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
C=====Premi�re colonne (N_mailles)
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
C=====Deuxi�me colonne (n_d_1)
C=============================
C----------------------
C-----El�ment ( 1 , 2 )
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
C-----El�ment ( 2 , 2 ) pour au moins 2 types intrins�ques
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
C-----El�ment ( 3 , 2 ) pour 3 types intrins�ques
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
C-----El�ments ( 2 , 2 ) et/ou ( 3 , 2 )
C-----pour moins de 2 ou 3 types intrins�ques  
C--------------------------------------------
        DO I_TYP = N_TYP_INTR + 1 , N_TYP
          IF ( I_TYP_R_DP_CONST ( 1 , 1 ) .EQ. I_TYP ) THEN
              MAT_CALC_POT ( I_TYP , 2 ) = 1.D0
          END IF
        END DO
C==============================
C=====Troisi�me colonne (n_d_2)
C==============================
C----------------------
C-----El�ment ( 1 , 3 )
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
C-----El�ment ( 2 , 3 ) pour au moins 2 types intrins�ques
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
C-----El�ment ( 3 , 3 ) pour 3 types intrins�ques
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
C-----El�ments ( 2 , 3 ) et/ou ( 3 , 3 )
C-----pour moins de 2 ou 3 types intrins�ques
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
C-----Calcul des potentiels chimiques � T = 0 K
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
       WRITE ( 500 , * ) "La fonction n'est pas born�e"
       WRITE ( 500 , 1200 )
      ELSE IF ( IC .EQ. - 1 ) THEN
       WRITE ( 500 , * ) 'Aucune solution avec ces contraintes'
       WRITE ( 500 , 1200 )
      ELSE IF ( IC .EQ. 0 ) THEN
       WRITE ( 500 , * ) 'Une solution trouv�e :'
       WRITE ( 500 , 1200 ) 
       WRITE ( 500 , * ) 'Num�ros des variables non nulles :'
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
C-----Cas d'�criture du point de composition �  "erreur simplexe"
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
C-----�criture des r�sultats dans le fichier de liste
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
C-----Ecriture du compte-rendu d'erreur �ventuelle du simplexe
C-------------------------------------------------------------        
       IF ( I_ERREUR_SIMPLEXE .NE. 0 ) THEN
         WRITE ( 100 , 1200 )
         WRITE ( 100 , * )
     $ 'ATTENTION : lors de ce calcul, pour certaines compositions :'
         WRITE ( 100 , * )
     $ 'la minimisation de H �tait impossible :'
         WRITE ( 100 , * )
     $ "(i) soit la fonction n'est pas born�e"
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
C#####Traitement du cas NPT (DP et �nergie libre)
C################################################
C################################################
      IF ( INDIC_TYP_CALC ( 1 : LONG_INDIC_TYP_CALC ) .EQ. 'NPT' )
     $ THEN
C-----------------------------------
C-----Nombre d'inconnues du probl�me
C-----------------------------------
       N_NPT = N_TYP_D_R + 1
       write ( * , * )
     $ "Dimension du syst�me d'�quations NPT = " , N_NPT
       write ( * , * ) 
     $ '------------------------------------'
C--------------------------------------------
C-----Vecteur, fonction et matrice jacobienne
C-----du probl�me non lin�aire
C--------------------------------------------
       ALLOCATE ( X_NPT ( 1 : N_NPT ) )
       ALLOCATE ( F_NPT ( 1 : N_NPT ) ) 
       ALLOCATE ( J_NPT ( 1 : N_NPT , 1 : N_NPT ) ) 
C-------------------------------------------------
C-----Sauvegarde du vecteur initial � un indice
C-----(seulement pour comparaison au vecteur final
C----- --> crit�re de "convergence" NRCG)
C-------------------------------------------------
       ALLOCATE ( X_NPT_INIT ( 1 : N_NPT ) )
C- - - - - - - - - - - - - - - - - - - - - - - -
C- - -Le vecteur initial est conserv�
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
C- - -Variable "nombre de mailles" (format = r�el)
C- - -Rappel : N_MAILLES_INIT, qui ne joue pas sur les r�sultats,
C- - -(xDP = variables intensives) a �t� initialis� � 1.D0
C- - -au d�but du traitement NPT
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      X_NPT_INIT ( N_NPT ) = N_MAILLES_INIT
C---------------------------------------------------------------------
C-----Nombre initial d'atomes par maille (entier),
C-----dont est d�duit ci-dessous N_AT_TOT (fonction de N_MAILLES)
C-----N_AT_TOT = variable pass�e en argument de la proc�dure SYS_NR
C-----Rappel : N_AT_MAILLE_INIT est une constante,
C-----comme N_1_MAILLE, N_2_MAILLE, N_3_MAILLE calcul�s d�s la lecture
C-----du fichier de donn�es, et relatifs au syst�me sans DP.
C---------------------------------------------------------------------
       N_AT_MAILLE_INIT = N_1_MAILLE
     $                  + N_2_MAILLE
     $                  + N_3_MAILLE
C---------------------------------------------------------
C-----Nombre total d'atomes (format = r�el)
C-----= variable pass�e en argument de la proc�dure SYS_NR
C---------------------------------------------------------
      N_AT_TOT = N_MAILLES_INIT * DFLOAT ( N_AT_MAILLE_INIT )
C-----------------------------------------------------------------------
C-----Affectation des param�tres lus 
C-----(nombres de pas et valeurs extr�males pour les boucles sur T et x)
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
C-----Ouverture des tableaux de conservation des r�sultats
C-----d'une it�ration sur l'autre
C-----(pour affectation des valeurs initiales).
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Ce tableau n'est utile qu'en NPT(Tx).
C-----La boucle interne (premier indice du tableau) contient alors
C-----les val. init. aux diverses compo. � la temp�rature pr�c�dente.
C--------------------------------------------------------------------
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
         ALLOCATE ( X_NPT_COMPO_SVG ( N_PAS_BOUCLE_INTERNE , N_NPT ) )
       END IF
C=====================================================================
C=====Initialisation pr�alable du vecteur d'inconnues (x_d ; M)  
C---------------------------------------------------------------------
C=====En NPT(T/x/Tx/xT), l'initialisation est d'abord faite ici
C=====(i.e. avant la double boucle sur T et x)
C=====� l'aide des valeurs lues dans DATA.adpi.
C=====* En NPT(T) : cette initialisation concerne la T initiale.
C=====  Ensuite, initialisation "de proche en proche"
C===== (--> r�sultat du point T(N) = amorce du point T(N+1))
C=====* En NPT(x/Tx/xT) : cette initialisation de proche en proche
C=====  pourra �tre supplant�e par une autre initialisation
C=====  effectu�e  � l'int�rieur des boucles sur T et/ou x :
C=====    soit � l'aide des valeurs de DATA.adpi,
C=====    soit via un fichier de valeurs initiales (f(x)) (cf. infra).
C=====   [f(x) signifie que ces val. init. d�pendent de la compo.]
C=====================================================================
C----------------------------------------------------------------
C-----Variable "nombre de mailles" (format = r�el)
C-----Rappel : N_MAILLES_INIT, qui ne joue pas sur les r�sultats,
C-----(xDP = variables intensives) a �t� initialis� � 1.D0
C-----au d�but du traitement NPT
C----------------------------------------------------------------
      X_NPT ( N_NPT ) = N_MAILLES_INIT
C-----------------------------
C-----Variables "log_10(x_DP)"
C-----------------------------
      DO I_R = 1 , N_R
        DO I_TYP = 0 , N_TYP
C- - - - - - - - - - - - - - - - - - - - - - - - 
C- - -Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C- - -ne correspondent � aucun DP ( I_TYP , I_R )
C- - -ou � un DP non utile (E_DP > 0 ou = 0)
C- - -et ne sont donc pas initialis�es.
C- - - - - - - - - - - - - - - - - - - - - - - - 
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
                X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $        = LOG_X_D_R_INIT ( I_TYP , I_R )
          END IF
        END DO
      END DO
C===========================================================
C=====Fin de l'initialisation pr�alable du vecteur (x_d ; M)  
C===========================================================
C---------------------------------------------------------------------
C-----Fichier o� seront �crites pour chaque composition
C-----les valeurs converg�es des variables pour le calcul NPT en cours
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
C-----(dont le contenu est plac� dans X_NPT_INIT_FICH),
C-----puis ouverture du fichier o� sera �crit le "suivi"
C-----de l'utilisation de ce tableau (le num�ro de ligne pour chaque x).
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
C#####Double boucle sur la temp�rature
C#####et sur la fraction atomique pour l'�l�ment sp�cifi�
C########################################################
      I_POUR_CENT = 0
      I_COMPTEUR_PAS = 0
C# # # # # # # # # #
C# # #Boucle externe
C# # # # # # # # # #
      DO I_PAS_BOUCLE_EXTERNE = 1 , N_PAS_BOUCLE_EXTERNE
C       write(*,*) "I_PAS_BOUCLE_EXTERNE=",I_PAS_BOUCLE_EXTERNE
C--------------------------------------------------------------
C-----En NPT(xT) seul : � chaque d�but de composition,
C-----r�initialisation optionnelle avec les valeurs de DATA.adpi
C-----(optionnelle car seulement si INDIC_LECT_VAL_INIT_XT = 1)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Si INDIC_LECT_VAL_INIT_NPT = 1,
C-----cette r�initialisation sera supplant�e plus loin
C-----� l'aide des valeurs lues dans le fichier de val. init. :
C-----� l'int�rieur des boucles, juste avant la r�solution,
C-----et toujours en d�but de nouvelle compo. c'est-�-dire
C-----pour la premi�re temp�rature (I_PAS_BOUCLE_INTERNE = 1).
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
C=====Op�rations pr�liminaires � la r�solution
C=====et � faire pour chaque pas en composition et T
C===================================================
C===================================================
C=======================================================================
C=====1) Affectation de la composition courante :
C=====* cas NPT(x/Tx/xT) : calcul des fractions atomiques courantes X_AT
C=====� partir de la valeur courante de X_AT ( I_TYP_0 )
C=====et des contraintes en composition
C=====* cas NPT(p/P/T) : la composition X_AT n'est pas modifi�e
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
C-----Cas x/Tx : le balayage en x est trait� par la boucle interne
          IF ( INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'x'
     $    .OR. INDIC_TYP_CALC_NPT
     $       ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
          X_AT_0_COUR = X_TYP_0_INIT
     $                + DFLOAT ( I_PAS_BOUCLE_INTERNE )
     $                * PAS_BOUCLE_INTERNE
C-----Cas xT : le balayage en x est trait� par la boucle externe
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
C= = =� partir des contraintes sur la composition
C= = =et de la fraction atomique courante
C= = =pour l'�l�ment sp�cifi� I_TYP_0,
C= = =inversion LU du syst�me d'�quations de contraintes
C= = =pour obtenir les valeurs des autres fractions atomiques
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      IF ( N_TYP .GT. 2 ) THEN
C----------------------------------------------------
C-----Boucle sur les lignes du syst�me de contraintes
C----------------------------------------------------
       DO I_TYP = 1 , N_SYS_CTR
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -Boucle sur les colonnes du syst�me de contraintes
C- - -en sautant la colonne I_TYP_0
C- - -(�l�ment dont on �tudie l'effet de l'enrichissement)
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
C-----Fin de boucle sur les lignes du syst�me de contraintes
C-----------------------------------------------------------
      END DO
C---------------------------------------------------------
C-----Inversion LU de la matrice du syst�me de contraintes
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
C-----et � la valeur courante de X_AT ( I_TYP_0 )
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
C=====2) Affectation de la temp�rature courante :
C=====* cas NPT(T/Tx/xT) : calcul de la valeur courante
C=====* cas NPT(p/P/x) : on garde la valeur constante lue dans DATA.adpi
C=======================================================================
C= = = = = = = = = = = = 
C= = =* Cas NPT(T/Tx/xT)
C= = = = = = = = = = = = 
C-----Cas o� le balayage en T est trait� par la boucle interne
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'T'
     $  .OR. INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'xT' ) THEN
        TEMPERATURE = T_INIT
     $              + DFLOAT ( I_PAS_BOUCLE_INTERNE )
     $              * PAS_BOUCLE_INTERNE
         K_T = K_B * TEMPERATURE
C-----Cas o� le balayage en T est trait� par la boucle externe
        ELSE
     $  IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx' ) THEN
        TEMPERATURE = T_INIT
     $              + DFLOAT ( I_PAS_BOUCLE_EXTERNE )
     $              * PAS_BOUCLE_EXTERNE
         K_T = K_B * TEMPERATURE
C= = = = = = = = = = = = = = = = = = = = =
C= = =Fin de l'obtention de la temp�rature
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
C-----la temp�rature garde la valeur lue dans DATA.adpi
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
C= = =Fin du test "cas NPT(T/Tx/xT) ou NPT(p/P/x) ?"
C= = =pour l'affectation de la temp�rature courante
C= = = = = = = = = = = = = = = = = = = = = = = = = = 
      END IF
C=======================================================
C=====Fin de l'affectation 2) de la temp�rature courante
C=======================================================
C=====================================================================
C=====Si l'option INDIC_LECT_VAL_INIT_NPT = 1 est choisie,
C=====et donc dans les cas NPT(x/Tx/xT) seulement
C=====[puisque le programme s'arr�te si
C=====INDIC_LECT_VAL_INIT_NPT = 1 et NPT(p/P/T)],
C=====utilisation du fichier de valeurs initiales
C===== --> prend le pas sur l'initialisation par DATA.adpi
C=====r�alis�e plus haut en dehors des boucles sur T et/ou x.
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
C=====Cette op�ration a lieu (ssi INDIC_LECT_VAL_INIT_NPT = 1) :
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(x) � tous les pas en compo. (avec val. init. f(x))
C=====   [f(x) signifie que ces val. init. d�pendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(Tx) : � la premi�re temp�rature seulement
C=====               (--> premier pas de boucle externe)
C=====                pour toutes les compositions
C=====               (--> tous les pas de boucle interne)
C=====               (avec la bonne ligne du fichier).
C=====             Aux temp�ratures suivantes T(i),
C=====             les valeurs initiales sont celles (f(x)) � T(i-1)
C=====             conserv�es dans le tableau X_NPT_COMPO_SVG.
C=====   [f(x) signifie que ces val. init. d�pendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(xT) : � chaque composition
C=====               (--> chaque pas de boucle externe)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C=====Rappel : si INDIC_LECT_VAL_INIT_NPT = 0 :
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(T/x), l'initialisation se fait � l'aide de DATA.adpi
C=====� la premi�re T ou x, puis de proche en proche aux autres T ou x.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(Tx), l'initialisation se fait
C=====(i) � la premi�re T : par DATA.adpi (quelle que soit la compo.)
C=====(ii) aux autres T(i) : � l'aide des val. init. (=f(x))
C=====     conserv�es dans le tableau X_NPT_COMPO_SVG.
C=====(i.e. le sch�ma d'initialisation "Tx" est similaire dans les cas
C=====INDIC_LECT_VAL_INIT_NPT = 0 et INDIC_LECT_VAL_INIT_NPT = 1,
C=====en rempla�ant, � la premi�re, T "DATA.adpi" par "fich. val. init")
C=====   [f(x) signifie que ces val. init. d�pendent de la compo.]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C=====* en NPT(xT), l'initialisation est effectu�e
C=====(i) � l'aide des valeurs de DATA.adpi au d�but de chaque compo.
C=====(ii) de proche en proche � chaque temp�rature dans cette compo.
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
C- - -doivent �tre multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arr�te (test � la lecture 
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autoris�s sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou �gal �)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilis�e successivement
C- - -     pour N_FRAC points en composition cons�cutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou �gal �)
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
C- - -Test suppl�mentaire de pr�caution pour �viter
C- - -d'�ventuels d�passements dans le num�ro de ligne
C- - -(m�me si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui pr�c�dent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C----------------------------------------------------------------------
C-----Pour NPT(x) et si INDIC_LECT_VAL_INIT_NPT = 1, initialisation
C-----� l'aide des valeurs lues dans le fichier de valeurs initiales
C-----(prend alors le pas sur l'initialisation faite plus haut,
C-----en dehors de la boucle sur T ou x, � l'aide des valeurs initiales
C-----lues dans DATA.adpi)
C----------------------------------------------------------------------
          DO I_NPT = 1 , N_NPT
             X_NPT ( I_NPT ) = X_NPT_INIT_FICH ( I_COUR , I_NPT )
          END DO
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fin du double test pour v�rifier que l'on est en mode NPT(x)
C= = =et que INDIC_LECT_VAL_INIT_NPT = 1
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        END IF
      END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Deuxi�me cas : NPT(Tx) -> la boucle sur x est la boucle interne
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Le fichier de val. init. n'est utilis� qu'� la premi�re T.
C= = =Aux T suivantes, on utilise un tableau conservant les r�sultats
C= = =pour les diverses compositions � la T pr�c�dente.
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
       IF ( INDIC_TYP_CALC_NPT
     $    ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $    ) THEN
        IF ( I_PAS_BOUCLE_EXTERNE .EQ. 1 ) THEN
C----------------------------------------------------------------------
C-----Si NPT(Tx) et s'il s'agit du premier pas en T, r�initialisation :
C-----* si INDIC_LECT_VAL_INIT_NPT = 1 :
C-----par utilisation du fichier de valeurs initiales (= f(x))
C-----via l'indice de ligne I_COUR fonction de I_PAS_BOUCLE_INTERNE,
C-----* si INDIC_LECT_VAL_INIT_NPT = 0 :
C-----par utilisation de DATA.adpi (-> val. init. ind�pendantes de x).
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
C- - -doivent �tre multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arr�te (test � la lecture
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autoris�s sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou �gal �)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilis�e successivement
C- - -     pour N_FRAC points en composition cons�cutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou �gal �)
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
C- - -Test suppl�mentaire de pr�caution pour �viter
C- - -d'�ventuels d�passements dans le num�ro de ligne
C- - -(m�me si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui pr�c�dent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x/xT/Tx) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x/xT/Tx) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C- - -NPT(Tx) : au premier pas en T, si INDIC_LECT_VAL_INIT_NPT = 1,
C- - -initialisation (f(x)) � l'aide du fichier de valeurs initiales
C- - -via I_COUR d�duit de I_PAS_BOUCLE_INTERNE
C- - - --> initialisation (f(x)) qui a lieu � toutes les compo.
C- - -   [f(x) signifie que ces val. init. d�pendent de la compo.]
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
C- - -ne correspondent � aucun DP ( I_TYP , I_R )
C- - -ou � un DP non utile (E_DP > 0 ou = 0)
C- - -et ne sont donc pas initialis�es.
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
C-----utilisation des valeurs de T pr�c�dentes,
C-----conserv�es dans  X_NPT_COMPO_SVG
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
C= = =Fin du test pour v�rifier que l'on est en mode NPT(Tx)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
      END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Troisi�me cas : NPT(xT) -> la boucle sur x est la boucle externe
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Le fichier de val. init. est utilis� � chaque compo.
C= = =(i.e. � chaque pas de boucle externe),
C= = =avant le d�marrage de la boucle en T
C= = =(i.e. au premier pas de la boucle interne)
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = = Rque : cette r�initialisation "fich. val. init" propre � "xT"
C= = =n'�tant effectu�e qu'au premier pas de boucle interne,
C= = =elle pourrait donc �tre extraite de la boucle interne
C= = =et regroup�e avec l'initialisation "xT" par DATA.adpi
C= = =situ�e plus haut entre les deux boucles.
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
C- - -doivent �tre multiples ou sous-multiples l'un de l'autre,
C- - -sinon, le programme s'arr�te (test � la lecture
C- - -de INDIC_LECT_VAL_INIT_NPT si celui-ci vaut 1).
C- - -Les seuls cas autoris�s sont donc :
C- - -(i) N_PAS_X_TYP_0 multiple de (ou �gal �)
C- - -    N_LIGNES_FICH_VAL_INIT
C- - - (e.g. N_PAS_X_TYP_0 = 400
C- - - et N_LIGNES_FICH_VAL_INIT = 200 => N_FRAC = 2)
C- - - --> chaque ligne du fichier est utilis�e successivement
C- - -     pour N_FRAC points en composition cons�cutifs.
C- - -(i) N_PAS_X_TYP_0 sous-multiple de (ou �gal �)
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
C- - -Test suppl�mentaire de pr�caution pour �viter
C- - -d'�ventuels d�passements dans le num�ro de ligne
C- - -(m�me si ceux-ci sont a priori impossibles,
C- - -compte tenu des instructions qui pr�c�dent)
C- - - - - - - - - - - - - - - - - - - - - - - - - - -
          IF ( I_COUR .GT. N_LIGNES_FICH_VAL_INIT ) THEN
               I_COUR = N_LIGNES_FICH_VAL_INIT
          ELSE IF ( I_COUR .LT. 1 ) THEN
               I_COUR = 1
          END IF
           write(410,*)
     $    "Calcul NPT(x/xT/Tx) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $     " -> valeurs init. lues dans fichier, ligne " , I_COUR
           write(*,*)
     $    "Calcul NPT(x/xT/Tx) - it�r. " , I_PAS_BOUCLE_INTERNE ,
     $    " -> valeurs init. lues dans fichier, ligne " , I_COUR
C----------------------------------------------------------------------
C-----Pour NPT(xT) et si INDIC_LECT_VAL_INIT_NPT = 1, initialisation
C-----� l'aide des valeurs lues dans le fichier de valeurs initiales
C-----(prend alors le pas sur l'initialisation faite plus haut,
C-----en dehors de la boucle sur T ou x, � l'aide des valeurs initiales
C-----lues dans DATA.adpi)
C----------------------------------------------------------------------
          DO I_NPT = 1 , N_NPT
             X_NPT ( I_NPT ) = X_NPT_INIT_FICH ( I_COUR , I_NPT )
          END DO
C- - -Fin du test pour v�rifier que c'est le premier pas en T
         END IF
C= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
C= = =Fin du double test pour v�rifier que l'on est en mode NPT(xT)
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
C=====Fin des op�rations pr�liminaires � la r�solution
C=====================================================
C=====================================================
C======================================
C======================================
C=====R�solution du syst�me d'�quations
C=====� chaque pas de T et/ou de x
C======================================
C======================================
C	write ( * , * ) '###################'
C	write ( * , * ) 'D�but r�solution NR'	
C	write ( * , * ) '###################'
         CALL
     $   SYS_NR
     $ (
C---------------------------------------
C-----Param�tres dont d�pend la fonction
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
C-----Param�tres g�n�raux de sys_NR
C----------------------------------
     $   N_NPT , X_NPT , F_NPT , J_NPT , P_J_NPT ,
     $   ALPHA_NRCG_NPT ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT ,
     $   INDIC_TYPE_REDUC_NRCG_NPT , COEF_REDUC_NRCG_NPT ,
     $   N_ITER_MAX_NPT , PRECISION_NPT )
C	write ( * , * ) '#################'
C	write ( * , * ) 'Fin r�solution NR'	
C	write ( * , * ) '#################'
C================================================
C================================================
C=====Fin de la r�solution du syst�me d'�quations
C================================================
C================================================
C--------------------------------------------------------
C-----V�rification de la convergence de l'algorithme NRCG
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
     $ 'probl�me de convergence'
        WRITE ( * , * )
     $ '(le nouveau point est identique au point initial)'
        WRITE ( * , * )
     $ '==> modifier le point initial'
        WRITE ( * , * )
     $ '-------------------------------------------------'
        CALL INTERRUPTION
      END IF
C------------------------------------------------------------------
C-----Ecriture de la solution � l'�cran et dans le fichier "varNPT"
C----- --> une ligne = le vecteur converg� � chaque pas T ou x
C------------------------------------------------------------------
C       write ( * , 1200 )
C       write ( * , * ) 'Solution :'
C       write ( * , 1200 )
C       write ( * , * )
C    $ 'Quantit�s de DP (log_10) - Nombre de mailles'
C	write ( * , 1111 )
C    $ ( X_NPT ( I ) , I = 1 , N_NPT )
C       write ( * , 1200 )
C-------------------------------------------------------------------
C-----Ecriture dans le fichier "varNPT" :
C-----* en modes NPT(Tx/xT) -> pour les divers points de composition
C-----� la temp�rature vis�e (la derni�re du balayage)
C-----* en mode NPT(x) -> pour les divers points de composition
C-----� la temp�rature de consigne
C-------------------------------------------------------------------
C-----Rque : la variable I_LIGNE_FICH_VAR_NPT ne renvoie � rien !
C-----Ce n'est pas g�nant, car cet indice n'est pas utilis�
C-----� la lecture d'un fichier de val. init. du type "varNPT".
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
C=====Calcul des quantit�s thermodynamiques dans l'ADPI
C=====(�nergie et volume par maille, entropie de configuration,
C=====�nergie libre par maille, quantit�s par atome)
C==============================================================
C==============================================================
C-------------------------------------------------
C-----Passage des variables du syst�me d'�quations
C-----aux quantit�s de DP
C-------------------------------------------------
      DO I_R = 1 , N_R
        DO I_TYP = 0 , N_TYP
C------------------------------------------------
C-----Les valeurs IND_D_R_TYP ( I_TYP , I_R ) = 0
C-----ne correspondent � aucun DP ( I_TYP , I_R )
C-----ou � un DP non "utile" ("E_DP >= 0")
C------------------------------------------------
          IF ( IND_D_R_TYP ( I_TYP , I_R ) .NE. 0 ) THEN
            X_D_R ( I_TYP , I_R )
     $    = DEXP ( X_NPT ( IND_D_R_TYP ( I_TYP , I_R ) )
     $           * DLOG ( 10.D0 ) )
          END IF
        END DO
      END DO
C-----------------------------------------------------------------
C-----Proc�dure de calcul des quantit�s thermodynamiques de l'ADPI
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
C=====Fin du calcul des quantit�s thermodynamiques dans l'ADPI
C=============================================================
C=============================================================
C================================================================
C================================================================
C=====Potentiels chimiques (= multiplicateurs de Lagrange en NPT)
C================================================================
C================================================================
C=====Calcul pr�liminaire des termes compl�mentaires d'entropie
C=====(m�me calcul que dans F_NR_AVEC_DERIV)
      Z_TYP_R = 1.D0
C-------------------
C-----Sous-r�seaux 1
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
C-----Sous-r�seaux 2 �ventuels
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
C-----Sous-r�seaux 3 �ventuels
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
C-----Sous-r�seaux interstitiels �ventuels
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
C=====Potentiels chimiques des �l�ments intrins�ques
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
C=====Potentiels chimiques des �l�ments d'addition
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
C=====Boucle sur les sous-r�seaux 1
C==================================
          DO I_R = 1 , N_R_1
C------------
C-----Lacunes
C------------
           H_FORM_D_R ( 0 , I_R )
     $   = H_GC_D_R ( 0 , I_R ) + POT_INTR ( 1 )
C--------------------------
C-----Antisites 2 �ventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 1 ) THEN
           H_FORM_D_R ( 2 , I_R )
     $   = H_GC_D_R ( 2 , I_R ) + POT_INTR ( 1 ) - POT_INTR ( 2 )
       END IF
C--------------------------
C-----Antisites 3 �ventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) + POT_INTR ( 1 ) - POT_INTR ( 3 )
       END IF
C--------------------------
C-----El�ments > N_TYP_INTR
C--------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_INTR ( 1 ) - POT_I ( I_TYP )
           END DO
C=========================================
C=====Fin de boucle sur les sous-r�seaux 1
C=========================================
        END DO
C==================================
C=====Boucle sur les sous-r�seaux 2
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
C-----Antisites 3 �ventuels
C--------------------------
       IF ( N_TYP_INTR .GT. 2 ) THEN
           H_FORM_D_R ( 3 , I_R )
     $   = H_GC_D_R ( 3 , I_R ) - POT_INTR ( 3 ) + POT_INTR ( 2 )
      END IF
C-------------------------------------------------
C-----Traitement �ventuel des esp�ces > N_TYP_INTR
C-------------------------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) + POT_INTR ( 2 ) - POT_I ( I_TYP )
         END DO
C=========================================
C=====Fin de boucle sur les sous-r�seaux 2
C=========================================
        END DO
C==================================
C=====Boucle sur les sous-r�seaux 3
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
C-----Traitement �ventuel des esp�ces > N_TYP_INTR 
C-------------------------------------------------
         DO I_TYP = N_TYP_INTR + 1 , N_TYP
           H_FORM_D_R ( I_TYP , I_R )
     $   = H_GC_D_R ( I_TYP , I_R ) +  POT_INTR ( 3 ) - POT_I ( I_TYP )
         END DO
C=========================================
C=====Fin de boucle sur les sous-r�seaux 3
C=========================================
        END DO
C==========================================================
C=====Boucle optionnelle sur les sous-r�seaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C-----------------------------------------------------
C-----El�ments 1, 2 et 3 (ces deux derniers �ventuels)
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
C-----El�ments d'addition �ventuels
C----------------------------------
           DO I_TYP = N_TYP_INTR + 1 , N_TYP
             H_FORM_D_R ( I_TYP , I_R )
     $     = H_GC_D_R ( I_TYP , I_R ) - POT_I ( I_TYP )
          END DO
C========================================================
C=====Fin de la boucle sur les sous-r�seaux interstitiels
C========================================================
        END DO
C==========================================================
C=====Fin du test d'existence de sous-r�seaux interstitiels
C==========================================================
       END IF
C==============================================
C==============================================
C=====Fin du calcul des enthalpies de formation
C==============================================
C==============================================
C============================================
C============================================
C=====Ecriture dans les fichiers de r�sultats
C============================================
C============================================
C-----------------------------------------
C-----Indicateur d'�criture conditionnelle
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
C-----En NPT(Tx/xT), seuls sont �crits les r�sultats
C-----� la temp�rature finale (la plus basse)
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
C-----Test pour �criture conditionnelle, suivant le cas NPT(T/x/Tx/xT)
C---------------------------------------------------------------------
       IF ( INDIC_ECRIT_NPT .EQ. 1 ) THEN
C==================================================================
C=====Ecriture des quantit�s de d�fauts ponctuels
C=====par sous-r�seau r et par d�faut ponctuel :
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
C-----Boucle sur les sous-r�seaux
C--------------------------------
      DO I_R = 1 , N_R
C-------------------------
C-----Boucle sur les types
C-------------------------
        DO I_TYP = 0 , N_TYP
C----------------------
C-----Num�ro du fichier
C----------------------
          I_FICH = 1000 * I_R + I_TYP
C--------------------------------------------
C-----Indicateur d'�criture de l'interstitiel
C-----pour le sous-r�seau et l'esp�ce donn�s
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
C-----ne correspondent pas � un DP
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
C-----Ecriture de la quantit� de DP
C----------------------------------
             WRITE ( I_FICH , CAR_COL_VAR_X_DP )
     $     ( POT_INTR ( J_TYP ) , J_TYP = 1 , N_TYP_INTR ) ,
     $     ( POT_I ( J_TYP ) , J_TYP = N_TYP_INTR + 1 , N_TYP ) ,
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE ,
     $       X_D_R ( I_TYP , I_R ) , H_FORM_D_R ( I_TYP , I_R )
C--------------------------------------
C-----Fin du test d'existence de d�faut
C--------------------------------------
        END IF
C----------------------------------------------
C-----Fin des boucles sur types et sous-r�seaux
C----------------------------------------------
       END DO
      END DO
C===========================================================
C=====Ecriture de l'�nergie, du volume et de l'�nergie libre
C===========================================================
C-------------------------------------------------
C-----Cas d'�criture de l'�nergie libre par maille
C-------------------------------------------------
            IF 
     $    ( INDIC_AT_MAILLE .EQ. 'M' .OR. INDIC_AT_MAILLE .EQ. 'm' )
     $      THEN
             WRITE ( 110 , CAR_COL_VAR_E_L )
     $     ( X_AT ( J_TYP ) , J_TYP = 1 , N_TYP ) ,
     $       TEMPERATURE , 
     $       E_MAILLE , V_MAILLE , S_CONF_MAILLE , G_MAILLE
C------------------------------------------------
C-----Cas d'�criture de l'�nergie libre par atome
C------------------------------------------------
            ELSE
C- - - - - - - - - - - - - - - - - - - - - - -
C- - -Dans le cas G/atome, deux possibilit�s :
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
C-----Fin du test pour �criture conditionnelle, suivant NPT(T/x/Tx/xT)
C---------------------------------------------------------------------
       END IF
C========================================
C========================================
C=====Fin des �critures dans les fichiers
C========================================
C========================================
C-----------------------------------------------------------------------
C-----Pourcentage du calcul effectu�
C-----(avec indicateur d'�criture pour �viter les r�p�titions � l'�cran)
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
C-----Ecriture du pourcentage � l'�cran
C--------------------------------------
          IF ( I_ECRIT_POUR_CENT .EQ. 1 ) THEN
           WRITE ( * , * )
     $      I_POUR_CENT , ' % du calcul NPT effectu�s'
          END IF
C============================================================
C=====Mise � jour des valeurs initiales entre deux it�rations
C===== --> utile seulement en NPT(Tx)
C============================================================
C---------------------------------------------------------------------
C-----Cas NPT(Tx) : tableau pour la conservation des r�sultats 
C-----(aux diverses compo. = boucle interne) � la temp�rature courante
C-----pour l'initialisation � chaque compo. � la temp�rature suivante
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----L'indice I_PAS_BOUCLE_INTERNE dans X_NPT_COMPO_SVG
C-----rep�re les diff�rentes compo. � la temp�rature courante T(N),
C----- --> val. init. f(x) utilis�es � la temp�rature suivante T(N+1).
C---------------------------------------------------------------------
        IF ( INDIC_TYP_CALC_NPT
     $     ( 1 : LONG_INDIC_TYP_CALC_NPT ) .EQ. 'Tx'
     $     ) THEN
         DO I = 1 , N_NPT
          X_NPT_COMPO_SVG ( I_PAS_BOUCLE_INTERNE , I ) = X_NPT ( I ) 
         END DO
        END IF
C#######################################################
C#####Fin de la double boucle sur les pas de temp�rature
C#####et de fraction atomique pour l'�l�ment sp�cifi�
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
     $ 'pour �viter des d�convenues lors des trac�s de graphes'
      WRITE ( * , * ) 
     $ '(points joints de mani�re erron�e), il peut �tre utile'
      WRITE ( * , * )
     $ 'de classer le contenu de chacun de ces fichiers'
      WRITE ( * , * ) 
     $ 'par valeurs croissantes de la composition mise en abscisse'
      WRITE ( * , * )
     $ '(la pr�sente version du programme adpi'
      WRITE ( * , * )
     $ "n'effectue pas ce classement)"
      WRITE ( * , * ) 
     $ '----------------------------------------------------------'
      WRITE ( * , * )
      WRITE ( * , * ) '          ============================'
      WRITE ( * , * ) '          Calcul ADPI (DP non charg�s)'
      WRITE ( * , * ) '          ----------------------------'
      WRITE ( * , * ) '                 FIN DU PROGRAMME'
      WRITE ( * , * ) '          ============================'
C-----Instruction permettant de terminer le programme
C-----en �vitant le cas "DP CHARGES" ci-apr�s
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
C-----Calcul des nombres d'atomes de chaque esp�ce par maille
C-----(utiles au calcul de la somme pond�r�e des potentiels chimiques)
C-----(d�j� pr�sent dans "adpi sans charge")
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
C-----Calcul de l'�nergie et du volume de r�f�rence par maille (eV)
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
C-----Calcul approch� de la somme pond�r�e des potentiels chimiques
C-----pour une maille � partir de la relation de Gibbs-Duhem
C------------------------------------------------------------------
      S_POT = V_REF_MAILLE * PRESSION * FACT_CONV_KBAR_EV_SUR_A_3
     $      + E_REF_MAILLE
C===========================================================
C=====Potentiel chimique de 1 (�l�ment de r�f�rence)
C=====(expression valable pour 2 ou 3 �l�ments intrins�ques)
C===========================================================
      POT_1 = S_POT - DFLOAT ( N_2_MAILLE ) * POT_2
     $              - DFLOAT ( N_3_MAILLE ) * POT_3
      POT_1 = POT_1 / DFLOAT ( N_1_MAILLE )
      write ( * , * ) 'mu(1)=',POT_1
      write ( * , * ) 'mu(2)=',POT_2
      write ( * , * ) 'mu(3)=',POT_3
C================================================
C=====Calcul optionnel de grandeurs pr�liminaires
C=====relatives aux complexes charg�s �ventuels
C================================================
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C--------------------------------------------------------------------
C-----Calcul du nombre de sites du sous-r�seau r dans chaque complexe
C-----et du nombre d'atomes de type i dans ce complexe
C--------------------------------------------------------------------
        ALLOCATE ( U_COMPLEXE_S_R_Q ( N_TYPES_COMPLEXES_Q , N_R ) )
        ALLOCATE ( V_COMPLEXE_TYPE_Q ( N_TYPES_COMPLEXES_Q , N_TYP ) )
        U_COMPLEXE_S_R_Q = O
        V_COMPLEXE_TYPE_Q = 0
C-------------------------------------
C-----Boucle sur les complexes charg�s
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
C-----Fin de la boucle sur les complexes charg�s
C-----------------------------------------------
        END DO
C--------------------------------------------------------
C-----Param�tres indicateurs de types et de sous-r�seaux
C-----utiles au calcul des termes de potentiels chimiques
C-----relatifs aux complexes �ventuels
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
C-----Type chimique normal de chaque sous-r�seau
C-----(0 pour sous-r�seaux interstitiels)
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
C-----Fin du test de prise en compte des complexes charg�s
C---------------------------------------------------------
       END IF
C=======================================================
C=====Fin de calcul optionnel de grandeurs pr�liminaires
C=====relatives aux complexes charg�s �ventuels
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
C-----antisites des diverses esp�ces intrins�ques
C-----sur les divers sous-r�seaux
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
C-----Partie 2/3 (optionnelle) : interstitiels intrins�ques
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
C-----lecture optionnelle des �nergies d'�l�ments d'addition
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
C-----Cas de pr�sence d'interstitiels d'addition
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
C-----Partie optionnelle : complexes charg�s
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
C#####2 options : balayage en mu(�lectrons)
C#####ou recherche (par la m�thode NR) de mu(�lec) 
C#####tel que n+A=p+D (neutralit� �lectrique)
C#################################################
C#################################################
C#################################################################
C#####Rappel : en pr�sence de DP charg�s :
C##### * mode muVT seul possible
C##### * pas de balayage en potentiels chimiques pour les �l�ments
C#####  --> potentiels chimiques fix�s directement dans DATA.adpi
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
C#####Boucle sur mu_e pour balayage ou algorithme NR (neutralit� �lec).
C#####Rappel : pas de balayage en potentiels chimiques pour les �l�ments
C#######################################################################
C#######################################################################
      F_COUR = 1.D100
      I_PAS_MU_ELEC = - 1
      I_ARRET = 0
      DO WHILE ( I_ARRET .EQ. 0 ) 
        I_PAS_MU_ELEC = I_PAS_MU_ELEC + 1
C=========================================
C=========================================
C=====Calcul de n et p, et de leur d�riv�e
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
C=====Boucle sur les sous-r�seaux 1
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 1
C=========================================
      END DO      
C============================================
C=====Boucle sur les sous-r�seaux 2 �ventuels
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 2 �ventuels
C===================================================
      END DO      
C============================================
C=====Boucle sur les sous-r�seaux 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 3 �ventuels
C===================================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-r�seaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C-----------------------------------------------------
C-----El�ments 1, 2 et 3 (ces deux derniers �ventuels)
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
C-----El�ments d'addition �ventuels
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
C=====Fin de boucle optionnelle sur les sous-r�seaux interstitiels
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
C-----Partie optionnelle : Hform de complexes charg�s
C----------------------------------------------------
C----------------------------------------------------
        IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $  .OR. INDIC_COMPLEXES_Q .EQ. 'o' )
     $ THEN
C=====================================
C=====Boucle sur les complexes charg�s
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
C-----Pr�liminaire :
C-----calcul du potentiel chimique de l'�l�ment sur le site courant
C-----et du potentiel chimique "normal" sur le site courant
C------------------------------------------------------------------
C-----Cas 1 : un type intrins�que
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
C-----Cas 2 : deux types intrins�ques
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
C-----Cas 3 : trois types intrins�ques
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
C-----(le m�me pour tous les �tats de charge d'un complexe donn�)
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
C-----Boucle sur les �tats de charge du complexe courant
C-----et calcul des enthalpies de formation associ�es
C-------------------------------------------------------
       DO I_Q = 1 , NQ_COMPLEXE_Q ( I_COMPLEXE_Q ) 
         H_FORM_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $ = H_GC_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q )
     $ + DELTA_MU_COMPLEXE
     $ + DFLOAT ( Q_COMPLEXE_Q ( I_Q , I_COMPLEXE_Q ) )
     $ * ( E_MAX_BV + POT_CHIM_ELEC )
      END DO
C===============================================
C=====Fin de la boucle sur les complexes charg�s
C===============================================
      END DO
C------------------------------------------------
C------------------------------------------------
C-----Fin de Hform(complexes charg�s) (optionnel)
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
C-----Antisites des diverses esp�ces intrins�ques
C-----sur les divers sous-r�seaux
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
C=====Partie 2/3 optionnelle : interstitiels intrins�ques
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
C=====Partie 3/3 (optionnelle) : �l�ments d'addition
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
C-----Cas de pr�sence d'interstitiels d'addition
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
C-----Fin du test de pr�sence d'interstitiels d'addition
C-------------------------------------------------------
        END IF
C---------------------------------------------
C---------------------------------------------
C-----Fin de boucle sur les types (partie 3/3)
C---------------------------------------------
C---------------------------------------------
      END DO
C===========================================
C=====Partie optionnelle : complexes charg�s
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
C===== -> mise � jour du potentiel chimique �lectronique
C=======================================================
       IF ( I_CALC_CHARGE .EQ. 1 ) THEN
        write ( * , 3001 ) I_PAS_MU_ELEC , POT_CHIM_ELEC ,
     $        C_VOL_N + C_VOL_ACC , C_VOL_P + C_VOL_DON 
        POT_CHIM_ELEC = POT_CHIM_ELEC + D_POT_CHIM_ELEC
       END IF 
C=======================================================
C=====Fonction-objectif n+A-p-D et sa d�riv�e (mode 2)
C===== -> mise � jour du potentiel chimique �lectronique
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
C-----Test pour l'arr�t en mode 1
C--------------------------------
        IF ( I_CALC_CHARGE .EQ. 1 ) THEN
         IF ( POT_CHIM_ELEC .GE. POT_CHIM_ELEC_MAX )
     $    I_ARRET = 1 
        END IF
C--------------------------------
C-----Test pour l'arr�t en mode 2
C--------------------------------
        IF ( I_CALC_CHARGE .EQ. 2 ) THEN
         IF ( DABS ( F_COUR ) .LT. EPS
     $     .OR. I_PAS_MU_ELEC .GT. N_MAX_PAS_NR_ELEC )
     $    I_ARRET = 1 
        END IF
C#######################################################################
C#######################################################################
C#####Fin de la boucle de balayage ou d'algorithme NR
C#####pour le potentiel chimique �lectronique
C#####Rappel : pas de balayage en potentiels chimiques pour les �l�ments
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
c=====En mode 2 (neutralit� �lec.),
C=====calcul des fractions atomiques
C===================================
C===================================
      IF ( I_CALC_CHARGE .EQ. 2 ) THEN
C-------------------------------------------------
C-----Termes au num�rateur des fractions atomiques
C-------------------------------------------------
          SOMME_1_TYP_2 = 0.D0
          SOMME_2_TYP_2 = 0.D0
          SOMME_1_TYP_3 = 0.D0
          SOMME_2_TYP_3 = 0.D0
C-----SOMME_I d�j� d�clar� dans "adpi"
          ALLOCATE ( SOMME_I ( N_TYP_INTR + 1 : N_TYP ) )
          SOMME_I = 0.D0
C---------------------------------------------------
C-----Termes au d�nominateur des fractions atomiques
C---------------------------------------------------
          SOMME_0 = 0.D0
          SOMME_MAILLE = 0.D0
C         write ( * , * ) 'P_R = ' , P_R
C==================================
C=====Boucle sur les sous-r�seaux 1
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 1
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-r�seaux 2
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C-----Terme de somme 1 pour l'esp�ce 2
C-------------------------------------
             SOMME_1_TYP_2 = SOMME_1_TYP_2
     $                     + SOMME_1_TYP_2_R * DFLOAT ( P_R ( I_R ) )
C=========================================
C=====Fin de boucle sur les sous-r�seaux 2
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-r�seaux 3
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
C-----El�ments > N_TYP_INTR
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
C-----Terme de somme 1 pour l'esp�ce 3
C-------------------------------------
         SOMME_1_TYP_3 = SOMME_1_TYP_3
     $                 + SOMME_1_TYP_3_R * DFLOAT ( P_R ( I_R ) )
C=========================================
C=====Fin de boucle sur les sous-r�seaux 3
C=========================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-r�seaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C---------------------------------------------
C-----El�ments 1, 2 et 3 (ce dernier �ventuel)
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
C-----El�ments d'addition �ventuels
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
C=====Fin de boucle optionnelle sur les sous-r�seaux interstitiels
C=================================================================
        END DO
       END IF
C=====================================================
C=====Contributions optionnelles des complexes charg�s
C=====================================================
C---------------------------------------------
C-----Initialisation des termes de sommes
C-----li�s aux complexes charg�s,
C-----(utiles au calcul des fractions atomiques,
C-----donc initialisation toujours faite,
C-----m�me sans complexes charg�s)
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
C=====Boucle optionnelle sur les complexes charg�s
C-------------------------------------------------
C-------------------------------------------------
       DO I_COMPLEXE_Q = 1 , N_TYPES_COMPLEXES_Q
C-------------------------------------------------------------------
C-----Fractions du complexe courant pour les divers �tats de charge,
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
C-----Termes de sommes relatifs � l'esp�ce 2 intrins�que
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
C-----Termes de sommes relatifs � l'esp�ce 3 intrins�que
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
C-----Termes de sommes relatifs aux �l�ments d'addition
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
C-----Termes de sommes relatifs aux quantit�s de mati�re totales
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
C=====Fin de la boucle sur les complexes charg�s
C-----------------------------------------------
C-----------------------------------------------
       END DO
C=================================================
C=====Fin du test de pr�sence de complexes charg�s
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
C-----Fractions atomiques des �l�ments d'addition  
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
C=====Ecriture des H_form pour la neutralit� �lectrique
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
C=====Boucle sur les sous-r�seaux 1
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 1
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-r�seaux 2
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
C-----Antisites 3 �ventuels
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 2
C=========================================
      END DO      
C==================================
C=====Boucle sur les sous-r�seaux 3
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
C-----El�ments > N_TYP_INTR
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
C=====Fin de boucle sur les sous-r�seaux 3
C=========================================
      END DO      
C==========================================================
C=====Boucle optionnelle sur les sous-r�seaux interstitiels
C==========================================================
       IF ( INDIC_R_INTER .EQ. 'O' .OR. INDIC_R_INTER .EQ. 'o' ) THEN
        DO I_R = 1 + N_R_3 + N_R_2 + N_R_1 , N_R
C---------------------------------------------
C-----El�ments 1, 2 et 3 (ce dernier �ventuel)
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
C-----El�ments d'addition �ventuels
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
C=====Fin de boucle optionnelle sur les sous-r�seaux interstitiels
C=================================================================
        END DO
       END IF
C====================================================
C=====Ecriture optionnelle pour les complexes charg�s
C====================================================
       IF ( INDIC_COMPLEXES_Q .EQ. 'O'
     $ .OR. INDIC_COMPLEXES_Q .EQ. 'o' ) THEN
       WRITE ( * , * ) '      - -- - - - - - - - - - - -'
       WRITE ( * , * ) '    Complexes charg�s :'
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
C=====Fin d'�criture optionnelle pour les complexes charg�s
C==========================================================
       END IF
C---------------------------------------
C-----Fin de l'�criture des Hf en mode 2
C---------------------------------------
      END IF
      WRITE ( * , * )
      WRITE ( * , * ) '          ========================'
      WRITE ( * , * ) '          Calcul ADPI (DP charg�s)'
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
C#####Inclusion de proc�dures
C############################
C----------------------------
C-----R�f�rences du programme
C----------------------------
C      INCLUDE
C     $ 'ref_prog.inc'
C     ==================================================================
C     I                                                                I
C     I R�f�rence de programme (laboratoire, auteur)		       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 27/10/2005                              I
C     I Sainte Emeline						       I  
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   REF_PROG
C     ##################################################################
      WRITE ( * , * ) 
     $ '###############################################################'
      WRITE ( * , * ) 
     $ '#            Ce programme est la propri�t� du                 #'
      WRITE ( * , * ) 
     $ '# Laboratoire de M�tallurgie Physique et G�nie des Mat�riaux  #'
      WRITE ( * , * ) 
     $ "#     Universit� de Lille 1, Villeneuve d'Ascq - FRANCE       #"
      WRITE ( * , * ) 
     $ '#                ----------------------------                 #'
      WRITE ( * , * ) 
     $ '#    Auteur : R�my Besson (Remy.Besson@univ-lille1.fr)        #'
      WRITE ( * , * ) 
     $ '###############################################################'
      WRITE ( * , * )
     $ "Veuillez presser la touche Entr�e pour lancer l'ex�cution." 
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
C     I Derni�re mise � jour : 30/05/2003                              I
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
C     I Constitution d'un format variable d'�criture de r�els          I
C     I et du titre associ� pour le calcul ADPI			       I
C     I (potentiel chimique de chaque esp�ce,			       I
C     I fraction atomique de chaque esp�ce,			       I
C     I fraction de DP sur le sous-r�seau consid�r�		       I
C     I et enthalpie de formation de ce DP)			       I
C     I et idem pour le fichier d'�nergie libre		               I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 17/06/2005                              I
C     I Saint Herv�                                                    I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C-----------------------------------------------------
C-----Caract�res utiles � la constitution d'un format   
C-----pour l'�criture d'un nombre variable de colonnes
C-----(potentiel chimique de chaque esp�ce,
C-----fraction atomique de chaque esp�ce,
C-----fraction de DP sur le sous-r�seau consid�r�
C-----et enthalpie de formation de ce DP)
C-----------------------------------------------------
      CHARACTER * 300 CAR_COL_VAR_X_DP_M
C----------------------------
C-----Format du titre associ�
C----------------------------
      CHARACTER * 5000 CAR_TITRE_VAR_X_DP_M
C-----------------------------------------
C-----Idem pour le fichier d'�nergie libre
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
C-----Indicateur d'�criture des grandeurs thermodynamiques
C-----par atome (A/a) ou par maille (M/m)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_AT_MAILLE_M
C---------------------------------------------------------
C-----Indicateur d'�criture de l'�nergie libre par atome :
C-----�nergie libre totale (T/t) ou de formation (F/f)
C---------------------------------------------------------
      CHARACTER * 1 INDIC_G_M
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C= = =Fichiers de d�fauts ponctuels et de potentiels chimiques
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
           WRITE ( * , * ) 'Nombre de types limit� � 99'
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
C-----Le cas I_TYP = 1 doit �tre distingu� car le texte qui le pr�c�de
C-----mu(N_TYP) est diff�rent de celui des autres types x_at(I_TYP - 1)
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
C-----Titre = "temp�rature"
C--------------------------
C--------------------------
      L = LEN_TRIM ( CAR_TITRE_VAR_X_DP_M )
      WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 12 : L + 25 ) , '(A)' )
C     WRITE ( CAR_TITRE_VAR_X_DP_M ( L + 2 : L + 15 ) , '(A)' )
     $       'Temp�rature(K)'
      L_1 = LEN_TRIM ( CAR_TITRE_VAR_POT_CHIM_M )
      WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 12 : L_1 + 25 ) ,
C     WRITE ( CAR_TITRE_VAR_POT_CHIM_M ( L_1 + 2 : L_1 + 15 ) ,
     $  '(A)' ) 'Temp�rature(K)'
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
C= = =Fichier d'�nergie libre
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
C-----Temp�rature
C----------------
      L = LEN_TRIM ( CAR_TITRE_VAR_E_L_M )
      WRITE ( CAR_TITRE_VAR_E_L_M ( L + 12 : L + 25 ) , '(A)' )
C     WRITE ( CAR_TITRE_VAR_E_L_M ( L + 2 : L + 15 ) , '(A)' )
     $       'Temp�rature(K)'
C----------------------------------------------------------------
C-----Energie, volume, entropie de configuration et �nergie libre
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
C-----(d�clar� �galement dans sys_non_lin_NR.inc)
C------------------------------------------------
C     INCLUDE
C    $ 'fact_lu.inc'
C------------------------------------------------
C-----Inversion d'une matrice par la m�thode LU
C-----(d�clar� �galement dans sys_non_lin_NR.inc)
C------------------------------------------------
C     INCLUDE
C    $ 'inverse_lu.inc'
C-----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de sous sous-r�seau et de son type
C-----------------------------------------------------------------------
C      INCLUDE
C     $ 'ind_D_R_typ.inc'
C     ==================================================================
C     I                                                                I
C     I Attribution d'un indice � chaque DP d'un alliage ordonn�       I
C     I unaire, binaire ou ternaire avec �l�ments d'addition   	       I	
C     I en fonction de son type et de son sous-r�seau                  I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 07/06/2005                              I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C-------------------------------------------------------------------
C-----Indice du d�faut en fonction de son sous-r�seau et de son type
C-------------------------------------------------------------------
      INTEGER * 4 IND_DP_TYP_M ( 0 : N_TYP_M , N_R_M )
C---------------------------------------
C-----Energie de SC de chaque type de DP
C---------------------------------------
      REAL * 8 E_B_D_R_M ( 0 : N_TYP_M , N_R_M )
C---------------------------------------------------------
C-----Indicateur de pr�sence de sous-r�seaux interstitiels
C---------------------------------------------------------
      CHARACTER * 1 I_R_INTER_M
C-----------------------------------------------------------------
C-----Indicateur de prise en compte des interstitiels intrins�ques
C-----------------------------------------------------------------
      CHARACTER * 1 I_R_INTER_INTR_M
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C-------------------
C-----Initialisation
C-------------------
      IND_DP_TYP_M  = 0
C-----------------------------------------
C-----Nombre de sous-r�seaux interstitiels
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
C-----pour le sous-r�seau et le type courant
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
C-----Test ci-dessous pour �carter les DP "non utiles"
C-----(e.g. certains interstitiels co�teux)
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
C-----Proc�dures du simplexe lin�aire
C------------------------------------
C      INCLUDE
C     $ 'simplexe_lin.inc'
C     ==================================================================
C     I                                                                I
C     I Proc�dures du simplexe lin�aire				       I
C     I (tir�es des "recettes num�riques" 			       I
C     I version 2005)						       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 02/12/2005                              I
C     I Sainte Viviane 						       I
C     ==================================================================
C     ##################################################################
      SUBROUTINE SIMPLX ( A , M , N ,
     $                    MP , NP ,
     $                    M1 , M2 , M3 ,
     $                    ICASE , IZROV , IPOSV )
C     ##################################################################
C#######################################################
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 IPOSV ( M ) , IZROV ( N )
      REAL * 8 EPS   
      INTEGER * 4 ICASE , M , M1 , M2 , M3 , MP , N , NP
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 LL ( NP )
      REAL * 8 BMAX
      INTEGER * 4 IABF , KP , MM , MP , NLL , NP
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 IP , KP , M , MP , N , NP
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 A ( MP , NP )
      INTEGER * 4 I1 , IP , K1 , KP , MP , NP
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C-----Analyse d'une ligne de caract�res
C-----pour distinction des cha�nes s�par�es par des espaces
C----------------------------------------------------------
C      INCLUDE
C     $ 'analyse_ligne.inc'
C     ==================================================================
C     I                                                                I
C     I Analyse d'une ligne de caract�res			       I
C     I pour distinction des cha�nes s�par�es par des espaces          I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 21/12/2004                              I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      CHARACTER * 1000 LIGNE_M
      INTEGER * 4 I_DEBUT_M ( 1000 )
      INTEGER * 4 I_FIN_M ( 1000 )
      INTEGER * 4 LONG_CHAINE_M ( 1000 )
      CHARACTER * 1000 CHAINE_M ( 1000 )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
      CHARACTER I_L_PREC
      CHARACTER I_L
C----------------------------------------------------------
C-----Recherche du nombre de cha�nes distinctes de la ligne
C-----et du caract�re de d�but de chaque cha�ne
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
C-----Longueur et contenu des cha�nes distinctes
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
C     I (proc�dure par r�currence) 				       I
C     I          						       I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 23/02/2005                              I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C-------------------
C-----Liste initiale
C-------------------
      INTEGER * 4 VAL_ENTIER_M ( N_VAL_M )
C---------------------------------------------
C-----Indice dans la liste initiale
C-----des �l�ments class�s par ordre croissant
C---------------------------------------------
      DIMENSION IND_INIT_M ( N_VAL_M )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C-------------------     
C-----Initialisation
C-------------------     
      DO I = 1 , N_VAL_M
	IND_INIT_M ( I ) = I
      END DO
C-----------------------------------------------------
C-----Balayage de la liste et insertion par r�currence
C-----------------------------------------------------
 1000     FORMAT
     $    ( 2X , 1000 ( 'a(', I2 , ' )   <   ' ) )  
 1100 FORMAT ('Liste initiale : ' , 2X , 1000 ( F7.2 , 6X ) )  
      DO I = 2 , N_VAL_M
C 	write ( * , * ) '-------'
C	write ( * , * ) 'I = ' , I
C	write ( * , * ) '-------'
C----------------------------------------------------------------
C-----Recherche du rang croissant de la nouvelle valeur � ins�rer
C-----dans la liste partielle dont le classement est d�j� connu
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
C-----Mise � jour des liens entre indices initiaux et croissants
C---------------------------------------------------------------
C--------------------------------------------------------------
C-----Premier cas :
C-----il existe des �l�ments de la liste partielle d�j� class�e
C-----encadrant le nouvel �l�ment
C--------------------------------------------------------------
	IF ( J_0 .NE. 0 ) THEN
C	  write ( * , * )
C    $   'Cas 1 : "il existe j tel que a(j) < a(i) < a(j+1)"' 
C	  write ( * , * ) ' J_0 = ' , J_0
	  DO K = I - 1 , J_0 + 1 , - 1
C----------------------------------------------------------
C-----Remarque : il faut balayer K en d�croissant
C-----pour �viter de perdre la premi�re valeur au recopiage
C----------------------------------------------------------
	    IND_INIT_M ( K + 1 ) = IND_INIT_M ( K )
	  END DO
	  IND_INIT_M ( J_0 + 1 ) = I
        ELSE
C---------------------------------------------------
C-----Deuxi�me cas :
C-----le nouvel �l�ment est sup�rieur
C-----� tous ceux de la liste partielle d�j� class�e
C---------------------------------------------------
	  I_1 = IND_INIT_M ( I - 1 )
	  IF ( VAL_ENTIER_M ( I ) .GT. VAL_ENTIER_M ( I_1 ) ) THEN
C         write ( * , * ) 
C    $   'Cas 2 : "quel que soit j, a(i) > a(j)"'
	   IND_INIT_M ( I ) = I 
	  ELSE
C---------------------------------------------------
C-----Troisi�me cas :
C-----le nouvel �l�ment est inf�rieur
C-----� tous ceux de la liste partielle d�j� class�e
C---------------------------------------------------
C         write ( * , * ) 
C    $   'Cas 3 : "quel que soit j, a(i) < a(j)"'
           DO K = I , 2 , - 1
C----------------------------------------------------------
C-----Remarque : il faut balayer K en d�croissant
C-----pour �viter de perdre la premi�re valeur au recopiage
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
C-----D�terminant et inverse d'une matrice 3 x 3
C-----------------------------------------------
C      INCLUDE
C     $'det_inv_mat_3_3.inc'
C     ==================================================================
C     I                                                                I
C     I Inversion d'une matrice 3 x 3 � l'aide de la formule :	       I
C     I Inverse de M 						       I
C     I = transpos�e (matrice des cofacteurs de M) / d�terminant (M)   I
C     I                                                                I
C     I Calcul du d�terminant de trois vecteurs			       I
C     I par d�veloppement premi�re colonne de la matrice correspondanteI
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 12/05/2003                              I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   INV_MAT_3_3
     $ ( D_MAT_M , D_MAT_INV_M ) 
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      DIMENSION D_MAT_M ( 3 , 3 )
      DIMENSION D_MAT_INV_M ( 3 , 3 )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C-----Calcul du d�terminant de D_MAT_M
C-----par d�veloppement premi�re colonne
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      DIMENSION B_1_M ( 3 ) , B_2_M ( 3 ) , B_3_M ( 3 )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C--------------------------------------------------------------
C-----D�veloppement premi�re colonne de la matrice des vecteurs
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
C     I Factorisation d'une matrice carr�e r�elle		       I
C     I en un produit LU d'une matrice L triangulaire inf�rieure       I
C     I � diagonale unit� et d'une matrice U triangulaire sup�rieure   I
C     I                                                                I
C     I La matrice � d�composer �tant M, le programme donne L et U     I
C     I telles que P * M * Q = L * U				       I
C     I (P et Q = matrices de permutations sur lignes et colonnes)     I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 24/01/2005                              I
C     I Saint Timoth�e						       I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 MAT_M ( N_M , N_M ) ,
     $         L_MAT_M ( N_M , N_M ) ,
     $         U_MAT_M ( N_M , N_M ) 
      REAL * 8 P_MAT_M ( N_M , N_M ) ,
     $	       Q_MAT_M ( N_M , N_M )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
      REAL * 8 VECT ( N_M ) , L_1_MAT ( N_M , N_M ) 
C======================================================
C=====Matrices de permutations lignes (P), colonnes (Q)
C=====et matrice L initialis�es � l'identit�
C======================================================
       DO I = 1 , N_M
         DO J = 1 , N_M
                P_MAT_M ( I , J ) = 0.D0
                Q_MAT_M ( I , J ) = 0.D0
C--------------------------------------
C-----Matrice triangulaire inf�rieure L
C--------------------------------------
                L_MAT_M ( I , J ) = 0.D0
                L_1_MAT ( I , J ) = 0.D0
                IF( I .EQ. J ) THEN
                        P_MAT_M ( I , J ) = 1.D0
                        Q_MAT_M ( I , J ) = 1.D0
                        L_1_MAT ( I , J ) = 1.D0
                END IF
C--------------------------------------
C-----Matrice triangulaire sup�rieure U
C--------------------------------------
                U_MAT_M ( I , J ) = MAT_M ( I , J )
         END DO
       END DO
C=================================================
C-----Factorisation LU de la matrice MAT_M 
C-----de dimension N_M x N_M en N_M - 1 it�rations
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
C-----(revient � une multiplication � gauche par P ( IT , I_PIV ))
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
C-----(revient � une multiplication � droite par P ( IT , J_PIV ))
C-----------------------------------------------------------------
        DO I = 1 , N_M
                VECT ( I ) = Q_MAT_M ( I , IT )
                Q_MAT_M ( I , IT ) = Q_MAT_M ( I , J_PIV )
                Q_MAT_M ( I , J_PIV ) = VECT ( I )
                VECT ( I ) = U_MAT_M ( I , IT )
                U_MAT_M ( I , IT ) = U_MAT_M ( I , J_PIV )
                U_MAT_M ( I , J_PIV ) = VECT ( I )
        END DO
C	write ( * , * ) 'U provisoire apr�s permutation'
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
C       write ( * , * ) 'L provisoire apr�s permutation'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( L_1_MAT ( I , J ) , J = 1 , N_M )
        END DO 
C--------------------------------------------------
C-----Op�rations sur L
C-----(multiplication � droite par M-1)
C-----correspondant � la triangularisation de U
C-----(celle-ci est effectu�e seulement ensuite,
C-----car elle change les �l�ments de matrice  de U
C-----utiles pour les op�rations sur L)
C--------------------------------------------------
        DO I = 1 , N_M
          DO J = IT + 1 , N_M
            L_1_MAT ( I , IT ) = L_1_MAT ( I , IT )
     $                         + L_1_MAT ( I , J )
     $                         / U_MAT_M ( IT , IT )
     $                         * U_MAT_M ( J , IT )
          END DO
        END DO
C       write ( * , * ) 'L provisoire apr�s op�ration'
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
C-----Cet �l�ment de matrice de la ligne I change
C-----lors de la combinaison lin�aire des lignes I et IT
C-----=> on le m�morise
C-------------------------------------------------------
            UU = U_MAT_M ( I , IT )
C		write ( * , * ) ' UU = ' , UU
            DO J = 1 , N_M
              U_MAT_M ( I , J )
     $      = U_MAT_M ( I , J )
     $      - UU / U_MAT_M ( IT , IT ) * U_MAT_M ( IT , J )
            END DO
        END DO
C        write ( * , * ) 'U provisoire apr�s triangulation'
        DO I = 1 , N_M
C         WRITE ( * , * ) ( U_MAT_M ( I , J ) , J = 1 , N_M )
        END DO
C-----------------------------------------
C---Fin de la boucle de N_M - 1 it�rations
C-----------------------------------------
        END DO
C       write ( * , * ) 'P_MAT_M pour multiplication par L'
        DO I = 1 , N_M
C               WRITE ( * , * ) ( P_MAT_M ( I , J ) , J = 1 , N_M )
        END DO     
C------------------------------------------------------------
C-----On termine le calcul de L en multipliant � gauche par P
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
C     I Calcul de l'inverse d'une matrice pr�alablement d�compos�e     I
C     I en produit L*U avec matrices de permutations P et Q	       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 24/01/2005                              I
C     I Saint Timoth�e						       I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 L_MAT_M ( N_M , N_M ) , U_MAT_M ( N_M , N_M ) ,
     $         P_MAT_M ( N_M , N_M ) , Q_MAT_M ( N_M , N_M ) ,
     $         MAT_INV_M ( N_M , N_M )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C-----Calcul de MAT_INV_LU inverse de L*U en deux �tapes
C-----(L * U * MAT_INV_LU = I)
C=======================================================
C---------------------------------------
C-----(i) R�solution de L * C_MAT = I
C-----(o� C_MAT = U * MAT_INV_LU)
C-----avec C_MAT triangulaire inf�rieure
C-----� diagonale unit�
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
C-----(ii) R�solution de U * MAT_INV_LU = C_MAT
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
C=====� l'aide des matrices de permutations P et Q
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
C     I Derni�re mise � jour : 06/07/2005                              I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 MAT_M ( N_M , N_M ) ,
     $         VECT_M ( N_M ) ,
     $         VECT_PROD_M ( N_M )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C     I (m�thode de convergence globale) tel que		       I
C     I x_1(nouveau) = x_0(ancien) + lambda * ( x_1_NR - x_0 )         I
C     I avec x_1_NR - x_0 = pas_NR (pas NR complet)		       I
C     I								       I
C     I Proc�dure :						       I
C     I � partir du pas NR complet (lambda = 1, x_1 = x_0 + pas_NR)    I
C     I lambda est r�duit progressivement			       I
C     I (x_1 = x_0 + lambda * pas_NR)				       I
C     I jusqu'� ce que soit v�rifi� le crit�re			       I
C     I phi(x_1) <= phi(x_0) + alpha grad(phi).( x_1 - x_0 )	       I
C     I	o� phi(x) fonction scalaire de x est d�finie par :	       I
C     I	phi(x) = (1/2) F(x).F(x) = (1/2) somme_i (F_i)^2               I
C     I (avec F(x) la fonction vectorielle dont on cherche une racine) I
C     I	=> grad(phi)_j = somme_i F_i(dF_i/dx_j) = somme_i F_i.J_ij     I
C     I (J = matrice des d�riv�es partielles)				       I
C     I 							       I
C     I Le pas initial (r�duction �ventuelle de lambda = 1 )           I
C     I est trait� diff�remment des autres			       I
C     I (approximation quadratique au lieu de cubique)		       I
C     I 							       I
C     I lambda est forc� � �tre entre deux valeurs minimale et maximaleI
C     I 							       I
C     I La proc�dure renvoie le nouveau vecteur x_1 accept�            I
C     I (et la valeur de phi en ce point)			       I
C     I 							       I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 26/07/2006                              I
C     I Sainte Anne                                                    I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   CALC_PAS_NRCG 
C----------------------------------------------------------
C----Param�tres sp�cifiques ADPI de la fonction vectorielle
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
C-----Param�tres d�termin�s en dehors de la proc�dure
C----------------------------------------------------
     $   N_M ,
     $   X_0_M , PHI_0_M , GRAD_PHI_0_M ,
     $   ALPHA_NRCG_NPT_M ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT_M ,
     $   INDIC_TYPE_REDUC_NRCG_NPT_M , COEF_REDUC_NRCG_NPT_M ,
C--------------------------------------------------------------------
C-----Param�tre de contr�le de l'�criture � l'�cran par CALC_PAS_NRCG
C--------------------------------------------------------------------
     $   I_ECRIT_CALC_PAS_NRCG_M ,
C-----------------------------------------------
C-----Param�tres d�termin�s dans cette proc�dure
C-----------------------------------------------
     $   PAS_NR_M , I_FONCTION_NON_DEFINIE_M ,
     $   X_1_M , PHI_1_M ) 
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C      PARAMETER ( TOLERANCE_X = 1.D-7 )
C#######################################################
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      REAL * 8 X_0_M ( N_M ) 
      REAL * 8 X_1_M ( N_M ) 
      REAL * 8 GRAD_PHI_0_M ( N_M ) 
      REAL * 8 PAS_NR_M ( N_M ) 
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-r�seau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-r�seau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-r�seau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantit� de mati�re totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de r�f�rence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C---------------------------------------
C-----Facteurs de r�duction du pas de NR
C---------------------------------------
      REAL * 8 LAMBDA
      REAL * 8 LAMBDA_MIN 
      REAL * 8 LAMBDA_2
      REAL * 8 LAMBDA_TEMP
C--------------------------------------
C-----Valeur de la fonction vectorielle
C-----au point x_1 estim�
C--------------------------------------
      REAL * 8 F_1 ( N_M )
C-----------------------------------------------------
C-----Indicateur de signe n�gatif pour chaque variable
C-----------------------------------------------------
      INTEGER * 4 INDIC_COMP_NEG ( N_M )
C-----Matrice des d�riv�es partielles :
C-----inutile ici, mais pr�sente en argument de F_NR_AVEC_DERIV
      REAL * 8 MAT_DERIV ( N_M )
C======================================
C=====Calcul de grad(phi).( x_1 - x_0 )
C======================================
C$$$$$ ATTENTION : v�rifier si x_1 ici = x_1_NR complet
C$$$$$ ou si x_1 doit �tre r�duit au fur et � mesure des pas NRCG
C$$$$$ (c'est le x_1 apparaissant dans le crit�re
C$$$$$ phi(x_1) <= phi(x_0) + alpha grad(phi).( x_1 - x_0 ))
      GRAD_PHI_PAS_NR = 0.D0
      DO I = 1 , N_M
	  GRAD_PHI_PAS_NR = GRAD_PHI_PAS_NR 
     $  	  	      + GRAD_PHI_0_M ( I ) * PAS_NR_M ( I )
      END DO	
C================================================================
C=====Calcul de lambda minimum autoris�.
C=====Il s'agit d'une borne indiquant une convergence "factice" :
C=====dans ce cas, la proc�dure s'arr�te et transmet x_1 � sys_NR
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
C-----Apr�s divers essais, on utilise finalement pour lambda_min
C-----une valeur fixe (faible) sp�cifi�e dans DATA.adpi.
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
C#####D�but de la boucle de r�duction de lambda
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
C=====Calcul du x_1 estim� pour le lambda courant
C=====(l'initialisation de lambda � 1 fait que
C=====le pas NR complet est toujours essay� en premier lieu).
C============================================================
C-----Remarque : en cas de X_1(I)<0, la d�marche suit deux �tapes :
C-----1) la variable PAS_NR, initialement le pas NR complet,
C-----est modifi�e, composante par composante, pour rem�dier � X_1(I)<0.
C-----2) une fois tous les X_1(I)>0, PAS_NR n'�volue plus,
C-----et c'est la r�duction de lambda qui est effectu�e.
C-----Rappel : en fonctionnement "normal",
C-----la variable I_FONCTION_NON_DEFINIE vaut 0.
C-----Elle peut passer � 1 dans deux cas :
C-----1) ci-dessous, lorsque l'une ou l'autre des composantes
C-----du x_1 estim� devient <0 ;
C-----2) en d�but de F_NR_AVEC_DERIV, si une quantit� z est trouv�e <0.
C-----Dans les deux cas, la valeur I_FONCTION_NON_DEFINIE = 1 est requise
C-----pour que f_NR effectue l'estimation de la fonction.
C-----Il y a deux appels � f_NR :
C-----* le premier dans sys_NR, pour �valuer F(x_0) en x_0 point courant ;
C-----* le second dans calc_pas_NRCG, pour �valuer F(x_1) (x_1 estim�). 
      I_FONCTION_NON_DEFINIE_M = 0
      INDIC_COMP_NEG = 0
      INDIC_EXISTENCE_COMP_NEG = 0 
	DO I = 1 , N_M - 1
C-----Rappel : except� la derni�re (nombre de mailles), les variables 
C-----x_NPT manipul�es par les proc�dures sont les log_10 des x_DP.
C-----Le passage � xDP = 10^x_NPT est fait en d�but de F_NR.
	  X_1_M ( I ) = DEXP ( DLOG ( 10.D0 ) * X_0_M ( I ) )
     $              + LAMBDA * PAS_NR_M ( I )
        IF ( X_1_M ( I ) .LE. 0.D0 ) THEN
         I_FONCTION_NON_DEFINIE_M = 1
         INDIC_EXISTENCE_COMP_NEG = 1
         INDIC_COMP_NEG ( I ) = 1
         if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
          write(*,*)
     $   'Pas NR -> valeur n�gative pour variable' , I 
          end if
        END IF
        X_1_M ( I ) = DLOG ( X_1_M ( I ) ) / DLOG ( 10.D0 )
C-----Fin de la boucle sur les composantes "log_10"
	END DO
C-----La derni�re composante est sp�cifi�e directement (et non son log).
      X_1_M ( N_M ) = X_0_M ( N_M ) + LAMBDA * PAS_NR_M ( N_M )
      IF ( X_1_M ( N_M ) .LE. 0.D0 ) THEN
         I_FONCTION_NON_DEFINIE_M = 1
         INDIC_EXISTENCE_COMP_NEG = 1
         INDIC_COMP_NEG ( N_M ) = 1
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
          write(*,*)
     $   'Pas NR -> valeur n�gative pour variable' , N_M
       end if
      END IF
C----------------------------------------------------------------
C-----R�duction pr�liminaire du pas pour "�limination de X_1<0" :
C-----* option 1 : r�duction suivant les seules comp. X_1(I)<0 
C-----* option 2 : r�duction de toutes les composantes du pas
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
       write ( * , * ) 'X_0 (log_10 sauf derni�re comp.) = ' , X_0_M
       write ( * , * ) 'X_1 (log_10 sauf derni�re comp.) = ' , X_1_M
       end if
C-----Remarque : la composante ci-dessus vaut "NaN" si x_1<0,
C-----mais ce n'est pas g�nant, car ce "NaN" ne sera pas utilis�
C-----gr�ce � l'information I_FONCTION_NON_DEFINIE = 1.
C=====================================================================
C=====Calcul de phi(x_1) = (1/2) F(x_1).F(x_1) = (1/2) somme_i (F_i)^2
C=====avec x_1 d�duit de la valeur courante de lambda
C=====================================================================
C-----Ici, le calcul des d�riv�es partielles est inutile.
       INDIC_CALC_DERIV = 0
         CALL
     $   F_NR_AVEC_DERIV
     $ (
C---------------------------------------
C-----Param�tres dont d�pend la fonction
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
C-----Param�tres g�n�raux de f_NR
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
C=====Remarque : par rapport � la proc�dure "classique",
C=====un cas suppl�mentaire (trait� en r�duisant lambda) a �t� ajout�,
C=====qui correspond � la non-d�finition de la fonction
C=====(ce qui peut se produire car celle-ci contient des logarithmes)
C=====================================================================
C===============
C=====Cas "z�ro"
C===============
       IF ( LAMBDA .LT. LAMBDA_MIN ) THEN
C $$$$$ Instruction ci-dessous � mettre en option ?
C 	  X_1_M = X_0_M
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
	  write ( * , * )
     $ '-----------------------------------------------------------'
	  write ( * , * )
     $ 'Convergence sur lambda "factice" = seuil lambda_min atteint'
C	  write ( * , * )
C     $ '    -> NRCG s'arr�te et transmet la valeur courante de x_1.'
	  write ( * , * )
     $ '-----------------------------------------------------------'
        end if
 	RETURN
C================
C=====Premier cas : la fonction est d�finie et le pas est convenable
C=====(i.e. le lambda courant v�rifie le crit�re de d�croissance
C=====de phi(x) = (1/2) F(x).F(x) = (1/2) somme_i (F_i)^2
C===== => la proc�dure calc_pas_NRCG s'arr�te et transmet � sys_NR
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
C=====Deuxi�me cas : la fonction est d�finie et le pas doit �tre r�duit
C=================	
       ELSE IF ( I_FONCTION_NON_DEFINIE_M .EQ. 0 ) THEN
C-----------------------------------------------------
C-----------------------------------------------------
C-----Premi�re r�duction (depuis lambda = 1 forc�ment)
C-----------------------------------------------------
C-----------------------------------------------------
	  IF ( LAMBDA .EQ. 1.D0 ) THEN
	    LAMBDA_TEMP
     $  = - GRAD_PHI_PAS_NR
     $	/ ( 2.D0 * ( PHI_1_M - PHI_0_M - GRAD_PHI_PAS_NR ) )
C------------------------------------
C------------------------------------
C-----Deuxi�me r�duction et suivantes
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
C-----(le cas A = 0 doit �tre distingu�)
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
C-----Fin du test sur r�ductions 1 ou (2 et suiv.)
C-------------------------------------------------
C-------------------------------------------------
	  END IF
C-------------------------------------------
C-----Mise � jour des variables incr�ment�es
C-------------------------------------------
       LAMBDA_2 = LAMBDA
       PHI_2 = PHI_1_M
C------------------------------
C-----Valeur minimale de lambda
C------------------------------
       LAMBDA = MAX ( LAMBDA_TEMP , 0.1D0 * LAMBDA )
C------------------------------
C-----Ecritures de v�rification
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
C=====Troisi�me cas : la fonction n'est pas d�finie
C==================
      ELSE
C----------------------------------------
C-----Le pas est r�duit dans ce cas aussi
C----------------------------------------
C        PAS_NR_M = PAS_NR_M * 0.5D0
       if ( I_ECRIT_CALC_PAS_NRCG_M .EQ. 1 ) then
C       write ( * , * ) '****************************'
	write ( * , * ) 'Dans calc_pas_NRCG :'
	write ( * , * ) ' Fonction non d�finie en X_1' 
C       write ( * , * ) '****************************'
C	write ( * , * ) '=> nouveau pas (la moiti� du pr�c�dent) :'
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
C#####Fin de la boucle de r�duction de lambda
C############################################
      END DO
      RETURN
      END
C     ##################################################################
C     ==================================================================
C     I								       I
C     I D�claration de la proc�dure de Newton-Raphson      	       I
C     I (incluant les param�tres sp�cifiques � la fonction)	       I
C     I          						       I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 26/07/2006                              I
C     I Sainte Anne                                                    I
C     ==================================================================
         SUBROUTINE
     $   SYS_NR
     $ (
C---------------------------------------
C-----Param�tres dont d�pend la fonction
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
C-----Param�tres g�n�raux de sys_NR
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-r�seau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-r�seau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-r�seau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantit� de mati�re totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de r�f�rence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C------------------------------------------------
C-----Point courant (vecteur), et en ce point :
C-----fonction et matrice des d�riv�es partielles
C------------------------------------------------
      REAL * 8 X_NPT_M ( N_NPT_M )
      REAL * 8 F_NPT_M ( N_NPT_M )
      REAL * 8 J_NPT_M ( N_NPT_M , N_NPT_M )
C---------------------------------------------------------
C------Fr�quence de mise � jour de la matrice des d�riv�es
C---------------------------------------------------------
      INTEGER * 4 P_J_NPT_M
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C---------------------------------------------
C-----Matrices relatives � la d�composition LU
C---------------------------------------------
      REAL * 8 L_MAT ( N_NPT_M , N_NPT_M ) ,
     $         U_MAT ( N_NPT_M , N_NPT_M ) ,
     $         P_MAT ( N_NPT_M , N_NPT_M ) ,
     $         Q_MAT ( N_NPT_M , N_NPT_M )
C---------------------------------------
C-----Inverse de la matrice des d�riv�es
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
C     write ( * , * ) 'D�but sys_NR'
C     write ( * , * ) '------------'
C-----Indicateur (interne � cette proc�dure) d'�criture de d�tails
      INDIC_ECRIT_DETAILS_SYS_NR = 0
C################################################################
C#####Algorithme de Newton-Raphson sur N_ITER_MAX_NPT it�rations
C#####ou arr�t � convergence
C#####(crit�re sur l'�cart entre deux points courants successifs)
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
      write ( * ,* ) '                    It�ration (NR) num�ro ' ,
     $                                    N_ITER
      write ( * ,* ) '                    *********************'	
C      write ( * ,* ) '================'	
C      write ( * , 2 ) 
C      write ( * , 1200 )
C      write ( * , * ) 'Point courant (log_10 sauf derni�re comp.)'
C      write ( * , 1200 )
C      write ( * , 1000 ) X_NPT_M 
C      write ( * , 1200 )
        end if
       N_ITER = N_ITER + 1
C=============================================================
C=====Pour le point courant,
C=====calcul de la fonction (sp�cifi�e analytiquement)
C=====ainsi que de dFi/dxj (analytiquement) et de son inverse,
C=====tous les P_J_NPT pas (et au premier pas)
C=============================================================
       INDIC_CALC_DERIV = 0
       IF ( MOD ( N_ITER , P_J_NPT_M ) .EQ. 0 .OR. N_ITER .EQ. 1 ) THEN
          INDIC_CALC_DERIV = 1
       END IF
C-----Remarque importante : f_NR renvoie I_FONCTION_NON_DEFINIE = 1 :
C-----1) si certaines variables x_1 sont n�gatives ;
C-----2) si certains termes compl�mentaires d'entropie z sont n�gatifs
C-----(se produit si le pas NR induit des variables x trop grandes).
C-----Dans ce 2�me cas, apr�s avoir calcul� les z et d�tect� z<0,
C-----f_NR s'arr�te sans calculer la fonction F,
C-----et "passe la main" � calc_pas_NRCG pour modification du pas.
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
      write(*,*) "Avant F_NR, I_FONCTION_NON_DEFINIE = " ,
     $  I_FONCTION_NON_DEFINIE 
        end if
C	stop
         CALL
     $   F_NR_AVEC_DERIV
     $ (
C---------------------------------------
C-----Param�tres dont d�pend la fonction
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
C-----Param�tres g�n�raux de f_NR
C--------------------------------
     $   I_FONCTION_NON_DEFINIE ,
     $   N_NPT_M , X_NPT_M , F_NPT_M ,
     $   INDIC_CALC_DERIV , J_NPT_M )
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
        write ( * , * )
     $ ' X_NPT_M (log_10 sauf derni�re comp.) = ' , X_NPT_M
        write ( * , * ) ' F_NPT_M = ' , F_NPT_M
        end if
C----------------------------------------------
C-----Tous les P_J_NPT pas (et au premier pas),
C-----factorisation LU et inversion
C-----de la matrice des d�riv�es partielles
C----------------------------------------------
      IF ( MOD ( N_ITER , P_J_NPT_M ) .EQ. 0 .OR. N_ITER .EQ. 1 ) THEN
C      write(*,*) "D�but FACT_LU"
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
        write ( * , * ) 'Matrice des d�riv�es partielles'
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
C=====R�duction �ventuelle de ce pas
C=====m�thode de convergence globale)
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
C-----R�duction �ventuelle du pas de NR
C- - - - - - - - - - - - - - - - - - - 
      X_0 = X_NPT_M     ! Sauvegarde du point courant
                        ! avant appel � calc_pas_NRCG
                        ! (calc_pas_NRCG re�oit X_0 en entr�e
                        ! et fournit X_NPT_M modifi� en sortie)
C---------------------------------------------------------------------
C-----Param�tre de contr�le de l'�criture � l'�cran par CALC_PAS_NRCG
C-----(pour limiter les �critures interm�diaires, lourdes si Tx ou xT)
C---------------------------------------------------------------------
         I_ECRIT_CALC_PAS_NRCG = 0
         CALL
     $   CALC_PAS_NRCG 
C----------------------------------------------------------
C----Param�tres sp�cifiques ADPI de la fonction vectorielle
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
C-----Param�tres d�termin�s en dehors de la proc�dure
C----------------------------------------------------
     $   N_NPT_M ,
     $   X_0 , PHI_0 , GRAD_PHI_0 ,
     $   ALPHA_NRCG_NPT_M ,
     $   VALEUR_LAMBDA_MIN_NRCG_NPT_M ,
     $   INDIC_TYPE_REDUC_NRCG_NPT_M , COEF_REDUC_NRCG_NPT_M ,
C--------------------------------------------------------------------
C-----Param�tre de contr�le de l'�criture � l'�cran par CALC_PAS_NRCG
C--------------------------------------------------------------------
     $   I_ECRIT_CALC_PAS_NRCG ,
C-----------------------------------------------
C-----Param�tres d�termin�s dans cette proc�dure
C-----------------------------------------------
     $   PAS_NR , I_FONCTION_NON_DEFINIE ,
     $   X_NPT_M , PHI_1 ) 
C--------------------------------------------------------------
C-----Nouveau pas r�duit apr�s que X_NPT a �t� modifi� par NRCG
C--------------------------------------------------------------
C Est-ce utile, puisque PAS_NR est �galement transmis par calc_pas_NRCG
C ci-dessus ? De plus, PAS_NR n'est plus utilis� ensuite dans sys_NR.
      PAS_NR = X_NPT_M - X_0
C-------------
C-----Ecriture
C-------------
        if ( INDIC_ECRIT_DETAILS_SYS_NR .EQ. 1 ) then
      write ( * , 1200 )
      write ( * , * ) 'Nouveau point apr�s NRCG'
      write ( * , 1200 )
      write ( * , * ) X_NPT_M
      write ( * , 1200 )
      write ( * , * ) 'Nouveau PAS_NR = ' , PAS_NR
       end if
C---------------------------------------------
C-----Test de convergence sur le vecteur-�cart
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
     $ '       ***** Crit�re de convergence NR atteint *****'
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
        WRITE ( * , * ) "Nombre maximal d'it�rations atteint"
	  WRITE ( * , * ) 'sans convergence'
        WRITE ( * , * ) '-----------------------------------'
         end if
        END IF 
C##################################
C#####Fin de la boucle d'it�rations
C##################################
      END DO
C     write ( * , * ) '----------'
C     write ( * , * ) 'Fin sys_NR'
C     write ( * , * ) '----------'
      RETURN
      END
C     ==================================================================
C     I                                                                I
C     I Interpolation lin�aire a*x+b de 2 points	               I
C     I (en vue de l'estimation d'un troisi�me)			       I
C     I                   				               I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 18/10/2005                              I
C     I Saint Luc	                                               I
C     ==================================================================
C----------------------------------------
C-----Interruption du programme
C-----(comment� si d�j� d�clar� ailleurs)
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
      DELTA_X = X_1_M - X_0_M 
      DELTA_Y = Y_1_M - Y_0_M
      DELTA = X_1_M * Y_0_M - X_0_M * Y_1_M
      IF ( DELTA_X .EQ. 0.D0 ) THEN
        WRITE ( * , * ) '----------------------'
	  WRITE ( * , * ) 'Interpolation lin�aire :'
	  WRITE ( * , * ) 'd�nominateur nul'
        WRITE ( * , * ) '----------------------'
	CALL INTERRUPTION
      END IF
      A_M = DELTA_Y / DELTA_X
      B_M = DELTA / DELTA_X
      RETURN
      END
C------------------------------------------
C-----Interpolation parabolique de 3 points
C-----(pour calcul d'un quatri�me)
C------------------------------------------
C      INCLUDE
C     $'inter_para.inc'
C     ==================================================================
C     I                                                                I
C     I Interpolation parabolique a*x2+b*x+c de 3 points	       I
C     I (en vue de l'estimation d'un quatri�me)			       I
C     I                   				               I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 18/10/2005                              I
C     I Saint Luc	                                               I
C     ==================================================================
C----------------------------------------
C-----Interruption du programme
C-----(comment� si d�j� d�clar� ailleurs)
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
      DELTA_0 = ( X_2_M * * 2 - X_0_M * * 2 ) * ( X_1_M - X_0_M )
     $        - ( X_1_M * * 2 - X_0_M * * 2 ) * ( X_2_M - X_0_M )
      IF ( DELTA_0 .EQ. 0.D0 ) THEN
        WRITE ( * , * ) '-------------------------'
	WRITE ( * , * ) 'Interpolation parabolique :'
	WRITE ( * , * ) 'd�terminant nul'
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
C-----Calcul des quantit�s thermodynamiques dans l'ADPI
C------------------------------------------------------
C      INCLUDE
C     $'G_adpi.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul des quantit�s thermodynamiques dans l'ADPI              I
C     I (�nergie et volume par maille, entropie de configuration,      I
C     I  �nergie libre par maille, quantit�s par atome)		       I
C     I                                                                I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 28/03/2006                              I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
      INTEGER * 4 P_R_M ( 1 : N_R_M ) 
      REAL * 8 E_REF_TYP_M ( 1 : N_TYP_M )
      REAL * 8 X_AT_M ( 1 : N_TYP_M )
      REAL * 8 X_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
      REAL * 8 E_GC_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
      REAL * 8 V_GC_D_R_M ( 0 : N_TYP_M , 1 : N_R_M )
C--------------------------------------
C-----Nombre d'atomes par maille (r�el)
C--------------------------------------
      REAL * 8 N_AT_MAILLE_M
C------------------------------------------------
C-----Indicateur de pr�sence de d�fauts complexes
C------------------------------------------------
      CHARACTER * 1 INDIC_COMPLEXES_M
C--------------------------------------
C-----Quantit�s relatives aux complexes
C--------------------------------------
      INTEGER * 4 MULTIPLICITE_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      INTEGER * 4 I_S_R_MULTIPLICITE_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 E_GC_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 V_GC_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
      REAL * 8 X_D_COMPLEXE_M ( N_TYPES_COMPLEXES_M )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C-----------------------------------------------------
C-----Termes x et x.ln(x) pour les �l�ments d'addition
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
      write(*,*) "    D�but de G_ADPI : taux de DP (1 ligne = 1 s-r) :"
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
C=====Sous-r�seaux 1
C===================
C     write ( * , * ) 'Sous-r�seaux 1'
      DO  I_R = 1 , N_R_1_M
C----------------------------
C-----Quantit�s pr�liminaires
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
C-----Energie et volume - contributions de 0, 2, 3 (�ventuel)
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
C-----Entropie de configuration - contributions de 0, 2, 3 (�ventuel)
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
C-----Contributions des �l�ments d'addition
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
C-----Terme compl�mentaire d'entropie
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
C=====Sous-r�seaux 2
C===================
C     write ( * , * ) 'Sous-r�seaux 2'
      DO  I_R = 1 + N_R_1_M , N_R_2_M + N_R_1_M
C----------------------------
C-----Quantit�s pr�liminaires
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
C-----Energie et volume - contributions de 0, 1, 3 (�ventuel)
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
C-----Entropie de configuration - contributions de 0, 1, 3 (�ventuel)
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
C-----Contributions des �l�ments d'addition
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
C-----Terme compl�mentaire d'entropie
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
C=====Sous-r�seaux 3
C===================
C     write ( * , * ) 'Sous-r�seaux 3'
      DO  I_R = 1 + N_R_1_M + N_R_2_M , N_R_3_M + N_R_2_M + N_R_1_M
C----------------------------
C-----Quantit�s pr�liminaires
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
C-----Contributions des �l�ments d'addition
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
C-----Terme compl�mentaire d'entropie
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
C=====Sous-r�seaux interstitiels
C===============================
C     write ( * , * ) 'Sous-r�seaux interstitiels'
      DO  I_R = 1 + N_R_3_M + N_R_2_M + N_R_1_M , N_R_M
C----------------------------
C-----Quantit�s pr�liminaires
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
C=====Energie, volume et entropie dus aux complexes �ventuels
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
C=====Energie, volume, entropie de configuration et �nergie libre
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
C-----Nombre d'atomes par maille et quantit�s thermodynamiques par atome
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
C-----Fonction du syst�me en mode NPT
C-----�crite pour N �l�ments intrins�ques
C----------------------------------------
C     INCLUDE
C    $'f_NR_N_intr.inc'
C     ==================================================================
C     I                                                                I
C     I Calcul de la fonction vectorielle de D variables	       I
C     I sp�cifi�e sous forme analytique 			       I
C     I intervenant dans l'ADPI NPT				       I
C     I pour l'algorithme de Newton-Raphson                            I
C     I								       I
C     I Pour un compos� ave N �l�ments intrins�ques,                   I
C     I I �l�ments d'addition et R sous-r�seaux		               I
C     I la dimension du syst�me est :				       I
C     I D = 1 + R * ( N + I )			                       I
C     I								       I
C     I Les inconnues sont :					       I
C     I (i) le nombre de mailles				       I
C     I (ii) les logarithmes (log_10) des fractions de DP	       I
C     I								       I
C     I Remarque : la pr�sente version est �crite pour N <= 3,         I
C     I mais le passage � N quelconque est possible		       I
C     I								       I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 16/07/2020                              I
C     I Notre Dame du Mont Carmel                                      I
C     ==================================================================
C     ##################################################################
         SUBROUTINE
     $   F_NR_AVEC_DERIV
     $ (
C----------------------------------------------
C----Param�tres sp�cifiques ADPI de la fonction
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
C----Param�tres g�n�raux de la fonction
C--------------------------------------
     $   I_FONCTION_NON_DEFINIE_M ,
     $   N_NPT_M , X_NPT_M , F_NPT_M ,
     $   INDIC_CALC_DERIV_M , MAT_DERIV_NPT_M )
C     ##################################################################
      IMPLICIT REAL * 8 ( A - H , O - Z )
      IMPLICIT INTEGER * 4 ( I - N )
C#######################################################
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C------------------------------------
C-----Vecteur et fonction vectorielle
C------------------------------------
      REAL * 8 X_NPT_M ( N_NPT_M ) ,
     $         F_NPT_M ( N_NPT_M )
C-----------------------------------------
C-----Matrice des d�riv�es partielles de F
C-----------------------------------------
      REAL * 8 MAT_DERIV_NPT_M ( N_NPT_M , N_NPT_M )
C----------------------------------------------------------------------
C-----Indice de chaque DP en fonction de son type et de son sous-r�seau
C----------------------------------------------------------------------
      INTEGER * 4 IND_D_R_TYP_M ( 0 : N_TYP_M , N_R_M ) 
C--------------------------------------------------------------------
C-----Enthalpie GC de chaque DP en fonction du sous-r�seau et du type
C--------------------------------------------------------------------
      REAL * 8 H_GC_D_R_M ( 0 : N_TYP_M , N_R_M )
C----------------------------
C-----Facteur thermodynamique
C----------------------------
      REAL * 8 K_T_M
C-------------------------------------------------------
C-----Nombre de sites par maille pour chaque sous-r�seau
C-------------------------------------------------------
      INTEGER * 4 P_R_M ( N_R_M )
C------------------------
C-----Fractions atomiques
C------------------------
      REAL * 8 X_AT_M ( N_TYP_M )
C-------------------------------
C-----Quantit� de mati�re totale
C-------------------------------
      REAL * 8 N_AT_TOT_M
C--------------------------------------
C-----Enthalpie de r�f�rence par maille
C--------------------------------------
      REAL * 8 H_REF_MAILLE_M
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
C#####################################################
C-----------------------------------------------------------------------
C-----Somme des nombres de sites par maille
C-----sur les sous-r�seaux relatifs � chaque esp�ce chimique intrins�que
C-----------------------------------------------------------------------
      INTEGER * 4 P_R_1 , P_R_2 , P_R_3
C----------------------------- 
C-----Nombre de mailles (r�el)
C-----------------------------
      REAL * 8 M_MAILLES
C----------------------------------------------------------
C-----Fractions de DP en fonction des sous-r�seaux et types
C----------------------------------------------------------
      REAL * 8 X_D_R ( 0 : N_TYP_M , N_R_M )
C-------------------------
C-----Quantit�s de mati�re
C-------------------------
      REAL * 8 N_AT ( N_TYP_M )
C--------------------------------------
C-----Termes compl�mentaires d'entropie
C--------------------------------------
      REAL * 8 Z_TYP_R ( N_R_M )
C-----Indicateur (interne � cette proc�dure) d'�criture de d�tails
      INDIC_ECRIT_DETAILS_F_NR = 0
C-----Notification du calcul analytique des d�riv�es
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
C-----Initialisation de X_D_R (n�cessaire pour le calcul des Z_R)
C----------------------------------------------------------------
      X_D_R = 0.D0
C-----------------------------------------------------------------------
C-----Somme des nombres de sites par maille
C-----sur les sous-r�seaux relatifs � chaque esp�ce chimique intrins�que
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
C-----Quantit�s de mati�re
C-------------------------
      DO I_TYP = 1 , N_TYP_M
        N_AT ( I_TYP ) = N_AT_TOT_M * X_AT_M ( I_TYP )
      END DO
C------------------------------------
C-----Matrice des d�riv�es partielles
C-----mise � jour tous les P_J_M pas
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
C=====et variables du probl�me
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
C-----Variables "log_10 des fractions de d�fauts ponctuels"
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
C-----Variable "nombre de mailles" (r�el)
C----------------------------------------
      M_MAILLES = X_NPT_M ( N_NPT_M )
C=================================================
C=====Calcul des termes compl�mentaires d'entropie
C=================================================
      Z_TYP_R = 1.D0
C-------------------
C-----Sous-r�seaux 1
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
C-----Sous-r�seaux 2 �ventuels
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
C-----Sous-r�seaux 3 �ventuels
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
C-----Sous-r�seaux interstitiels �ventuels
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
C-----Test de positivit� de ces termes Z
C-----(sinon la fonction n'est pas d�finie)
C-----Rque : I_FONCTION_NON_DEFINIE_M peut �tre modifi� par F_NR,
C-----mais sa valeur en entr�e de F_NR (d�j� fix�e dans calc_pas_NRCG
C-----par la pr�sence de x_1(i)<0) a aussi de l'importance.
C--------------------------------------------------------------------
      DO I_R = 1 , N_R_M
	 IF ( Z_TYP_R ( I_R ) .LE. 0.D0 ) THEN
        WRITE ( * , * ) 
     $ '----------------------------------------------'
	  WRITE ( * , * )
     $ "Dans f_NR, le terme compl�mentaire d'entropie"
         WRITE ( * , * ) 
     $ 'Z_TYP_R ( I_R = ' , I_R , ' )'
         WRITE ( * , * ) 
     $ 'est n�gatif'
	 write ( * , * ) 'Z_TYP_R = ' , Z_TYP_R ( I_R )
	 I_FONCTION_NON_DEFINIE_M = 1
	 END IF
      END DO
C-----Si z<0 au point courant, la proc�dure s'arr�te.
C-----C'est le cas aussi lorsque I_FONCTION_NON_DEFINIE_M 
C-----vaut d�j� 1 � l'entr�e de la proc�dure
C-----(valeur 1 affect�e par calc_pas_NRCG). 
      IF ( I_FONCTION_NON_DEFINIE_M .EQ. 1 ) THEN
       if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
        WRITE ( * , * )
     $ 'Fonction non d�finie => le programme quitte f_NR'
        WRITE ( * , * )
     $ ' (pour r�duction de PAS_NR dans calc_pas_NRCG)'
        end if
       RETURN
      END IF
C===================================
C=====Multiplicateurs de Lagrange
C=====pour les �l�ments intrins�ques
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
C#####en fonction des variables du probl�me
C###########################################
C-------------------------
C-----Indice de composante
C-------------------------
      I_COMP = 0
C===========================================
C=====Composantes de la fonction vectorielle
C=====relatives aux sous-r�seaux 1
C===========================================
      DO I_R = 1 , N_R_1_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-r�seau 1, il y a une �quation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. syst�me d'�quations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les �l�ments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 2 et 3, ainsi que pour les lacunes 
C-----(sauf pour le sous-r�seau N_R_1)
C-----et pour les �l�ments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des �l�ments intrins�ques
C----------------------------------
C----------------------------------
C-----------------------  
C-----------------------  
C-----El�ment 2 �ventuel
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
C    $  'pour s-r 1 et �l�ment 2 : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*****************************************************
C-----cf. syst�me d'�quations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les �l�ments d'addition
C-----(notes du 16/12/19).
C*****************************************************
C***********************************************************************
C-----L'indice J_COUR d�signe une variable de la fonction vectorielle,
C-----correspondant � toutes les possibilit�s (I_R ; I_TYP) telles que
C-----IND_D_R_TYP_M ( I_TYP , I_R ) # 0.
C-----Le cas o� IND_D_R_TYP ( I_TYP , I_R ) = 0, qui ne renvoie
C-----� aucune variable de cette fonction, recouvre deux situations :
C----- 1) les combinaisons (I_R ; I_TYP) qui ne renvoient � aucun DP,
C-----i.e. (i) un �l�ment intrins�que sur son sous-r�seau "normal",
C-----    (ii) une lacune sur un sous-r�seau interstitiel.
C----- 2) les combinaisons (I_R ; I_TYP) correspondant � un DP possible,
C-----    mais qui est �cart� par une valeur nulle d'�nergie de SC.
C-----Il faut alors noter que les termes compl�mentaires z^R
C-----pr�sents dans les composantes Fi de la fonction
C-----ne d�pendent que des variables autres que 1) et 2) ci-dessus.
C***********************************************************************
C***********************************************************************
C-----Dans chaque composante de F, le calcul des d�riv�es partielles
C-----de z^R n�cessite d'�viter ces "fausses variables".
C-----Pour ce faire, on a d'abord eu recours � des tests adapt�s
C-----aux divers cas (I_R ; I_TYP), e.g. "IF ( J_TYP .NE. 1 )" qui �vite
C-----une "fausse variable" de type 1)(i).
C-----Cependant, un tel test n'�vite pas le cas 2), e.g. un syst�me
C-----avec certains antisites tr�s co�teux (ou d'�nergie incalculable),
C-----auxquels on associe donc une �nergie de SC nulle.
C-----Dans ce cas, le programme affecte IND_D_R_TYP = 0 au DP,
C-----qui ne correspond alors plus � une variable du syst�me NPT.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----[NB : la composante Fi (de la fonction vectorielle) associ�e
C-----� ce DP ne doit alors pas �tre cr��e par le programme
C-----(attention : ce n'est pas le cas pour l'instant, cf. supra :
C-----aucun test n'est pr�sent pour �carter certains DP d'antisite
C-----I_TYP=2 sur s-r normalement "1", et Fi est toujours cr��e :
C-----le rejet de ces Fi inactifs - premier indice de MAT_DERIV_NPT_M -
C-----n'est donc pas impl�ment�,
C-----alors qu'il doit l'�tre au m�me titre que le rejet des dFi/dxj
C----- - deuxi�me indice de MAT_DERIV_NPT_M).]
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Ainsi, pour ne pas inclure les fausses variables 1) ET 2), on peut
C-----remplacer le test pr�c�dent par un test "IF ( J_COUR .NE. 0 )",
C-----comme ci-dessous.
C***********************************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter le DP ( 2 , I_R ) (par E_SC = 0)
C----- = fausse variable 2).
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C-----(la composante Fi doit aussi �tre �cart�e,
C-----ce que ne fait pas la version actuelle :
C-----dans ce cas, le passage dans cette section
C-----"d�riv�es partielles" ne devrait m�me pas avoir lieu).
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
C----- --> si E_SC(lacune N_R_1) = 0, sugg�rer de permuter
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
C----- --> si E_SC(lacune N_R_1 + N_R_2) = 0, sugg�rer de permuter
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : �l�ment 2 �ventuel
C------------------------------------
C------------------------------------
      END IF
C-----------------------
C-----------------------
C-----El�ment 3 �ventuel
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
C    $  'pour s-r 1 et �l�ment 3 : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
          if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
           write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
          end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : �l�ment 3 �ventuel
C------------------------------------
C------------------------------------
      END IF
C--------------------------------------------------------
C--------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_1 > 1 
C-----(car la composante Fi pour les lacunes en N_R_1
C-----sert � �liminer le multiplicateur de Lagrange de 1)
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
            write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C-----Cas des �l�ments d'addition �ventuels
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
C-----de la fonction vectorielle �tudi�e
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
C-----Si nouvelle composante, calcul de celle-ci et de ses d�riv�es
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
C    $  'pour s-r 1 et �l�ment add. : F_NPT_M ( ' , I_COMP , ' ) = ' ,
C    $  F_NPT_M ( I_COMP )
C      write ( * , * )
C    $  '-----------------------'
C-----------------------------------------------------------
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 1 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - - 
           IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( I_TYP , N_R_M )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - 
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_1"
C-----On ajoute ci-dessous les d�riv�es du terme DELTA_1,
C-----en distinguant les diff�rents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas o� DELTA_1 = 0 => d�riv�e nulle
           IF ( N_R_M .EQ. N_R_1_M ) THEN
           IF ( I_R .LT. N_R_1_M ) THEN 
           END IF
C-----Cas o� DELTA_1 = X_MULTI_1 - X_MULTI_2 = mu_A - mu_B
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $      .AND. ( N_R_2_M .NE. 0 ) ) THEN
C- - - D�riv�es de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1 sont requises
C-----pour l'�limination des multiplicateurs de Lagrange
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1_M 
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M )
            END IF
          END DO
C- - - D�riv�es de - mu_B (ici N_R_M = N_R_1_M + N_R_2_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1_M + N_R_2_M 
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas o� DELTA_1 = X_MULTI_1 - X_MULTI_3 = mu_A - mu_C
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
C- - - D�riv�es de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1_M
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M )
           END IF
          END DO
C- - - D�riv�es de - mu_C (ici N_R_M = N_R_1_M + N_R_2_M + N_R_3_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
         J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
          DO J_TYP = 0 , N_TYP_M
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1 + N_R_2 + N_R_3
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas o� DELTA_1 = X_MULTI_1 = + mu_A
        ELSE 
C- - - D�riv�es de mu_A
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C=====Fin de boucle sur les sous-r�seaux 1
C=========================================
      END DO
C===============================================
C=====Composantes de la fonction vectorielle
C=====(�ventuelles) relatives aux sous-r�seaux 2
C===============================================
      DO I_R = N_R_1_M + 1 , N_R_1_M + N_R_2_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-r�seau 2, il y a une �quation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. syst�me d'�quations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les �l�ments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 1 et 3, ainsi que pour les lacunes 
C-----(sauf pour le sous-r�seau N_R_1 + N_R_2)
C-----et pour les �l�ments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des �l�ments intrins�ques
C----------------------------------
C----------------------------------
C--------------
C--------------
C-----El�ment 1
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter le DP ( 1 , I_R ) (par E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
           DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 1 ) THEN
             J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C-----------------------  
C-----------------------
C-----El�ment 3 �ventuel
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 2 ) THEN
             J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR 
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : �l�ment 3 �ventuel
C------------------------------------
C------------------------------------
       END IF
C----------------------------------------------------------
C----------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_2 > 1
C-----(car la composante Fi pour les lacunes en N_R_1+N_R_2
C-----sert � �liminer le multiplicateur de Lagrange de 2)
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (e.g. type=2 sur s-r "2").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C-----Cas des �l�ments d'addition �ventuels
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
C-----de la fonction vectorielle �tudi�e
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
C-----Si nouvelle composante, calcul de celle-ci et de ses d�riv�es
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 2 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / X_D_R ( I_TYP , N_R_M )
          END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_2"
C-----On ajoute ci-dessous les d�riv�es du terme DELTA_2,
C-----en distinguant les diff�rents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas o� DELTA_2 = 0 => d�riv�e nulle
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M )
     $ .AND. ( N_R_2_M .NE. 0 ) ) THEN
           IF ( I_R .LT. N_R_1_M + N_R_2_M ) THEN 
           END IF
C-----Cas o� DELTA_2 = X_MULTI_2 - X_MULTI_3 = mu_B - mu_C
        ELSE IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $      .AND. ( N_R_3_M .NE. 0 ) ) THEN
C- - - D�riv�es de mu_B
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2 sont requises
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1+N_R_2
C- - - et les DP tels que (E_SC = 0)
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
            IF ( J_COUR .NE. 0 ) THEN
             MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $     - K_T_M / Z_TYP_R ( N_R_1_M + N_R_2_M )
           END IF
          END DO
C- - - D�riv�es de - mu_C (ici N_R_M = N_R_1_M + N_R_2_M + N_R_3_M)
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1 + N_R_2 + N_R_3
C- - - et les DP tels que (E_SC = 0)
            IF ( J_COUR .NE. 0 ) THEN
              MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $      + K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C-----Cas o� DELTA_2 = X_MULTI_2 = + mu_B
        ELSE 
C- - - D�riv�es de mu_B
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici
        J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1+N_R_2
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C=====Fin de boucle sur les sous-r�seaux 2
C=========================================
      END DO
C===============================================
C=====Composantes de la fonction vectorielle
C=====(�ventuelles) relatives aux sous-r�seaux 3
C===============================================
      DO I_R = N_R_1_M + N_R_2_M + 1 , N_R_1_M + N_R_2_M + N_R_3_M
C----------------------------------------------------------
C----------------------------------------------------------
C-----Pour chaque sous-r�seau 3, il y a une �quation
C-----[= composante Fi de la fonction vectorielle,
C-----cf. syst�me d'�quations (2.100)-(2.103) de l'HDR
C-----et ses corrections pour les �l�ments d'addition
C-----(notes du 16/12/19)]
C-----pour les antisites 1 et 2, ainsi que pour les lacunes 
C-----(sauf pour le sous-r�seau N_R_1 + N_R_2 + N_R_3)
C-----et pour les �l�ments d'addition
C----------------------------------------------------------
C----------------------------------------------------------
C----------------------------------
C----------------------------------
C-----Cas des �l�ments intrins�ques
C----------------------------------
C----------------------------------
C--------------  
C--------------  
C-----El�ment 1
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter le DP ( 1 , I_R ) (par E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C           IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C--------------  
C--------------  
C-----El�ment 2
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (type=1 sur s-r "1").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
          J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
          J_COUR
     $  = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
          MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $  - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C---------------------------
C---------------------------
C-----Fin du cas : �l�ment 2
C---------------------------
C---------------------------
C----------------------------------------------------------------
C----------------------------------------------------------------
C-----Cas des lacunes, seulement si N_R_3 > 1
C-----(car la composante Fi pour les lacunes en N_R_1+N_R_2+N_R_3
C-----sert � �liminer le multiplicateur de Lagrange de 3)
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
          J_COUR = IND_D_R_TYP_M ( 0 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�) induites par z^I_R :
C-----attention � bien voir de quelles variables d�pend z^I_R
C-----(�viter les fausses variables 1) ET 2)) !
C-----Le test initialement utilis� "IF ( J_TYP .NE. 1 )"
C-----n'�carte que les fausses variables 1) (e.g. type=2 sur s-r "2").
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
            J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
            if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
            end if
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C-----Cas des �l�ments d'addition �ventuels
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
C-----de la fonction vectorielle �tudi�e
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
C-----Si nouvelle composante, calcul de celle-ci et de ses d�riv�es
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
C         WRITE(*,*) "s-r 3 : DERIVEES"
          J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / X_D_R ( I_TYP , N_R_M )
           END IF
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
           J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Nb3Sn-Cu-Ta
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITION EN SUBSTITUTION TEL QUE Nb3Sn-Cu
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C     WRITE(*,*) "DEBUT DES DERIVEES DE DELTA_3"
C-----On ajoute ci-dessous les d�riv�es du terme DELTA_3,
C-----en distinguant les diff�rents cas (cf. notes du 16/12/19)
C-----(on est ici dans la situation I_NOUV_COMP = 1)
C-----Cas o� DELTA_3 = 0 => d�riv�e nulle
        IF ( ( N_R_M .EQ. N_R_1_M + N_R_2_M + N_R_3_M )
     $ .AND. ( N_R_3_M .NE. 0 ) ) THEN
          IF ( I_R .LT. N_R_1_M + N_R_2_M + N_R_3_M ) THEN 
           END IF
C-----Cas o� DELTA_3 = X_MULTI_3 = + mu_C
        ELSE 
C- - - D�riv�es de mu_C
C-----Pas de test  "IF ( J_COUR .NE. 0 )" ici :
C-----les lacunes sur le s-r N_R_1+N_R_2+N_R_3 sont requises
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   - K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2+N_R_3)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C- - - Test  "IF ( J_COUR .NE. 0 )" ici :
C- - - on saute le type intrins�que du sous-r�seau N_R_1+N_R_2+N_R_3
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
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
C=====Fin de boucle sur les sous-r�seaux 3
C=========================================
      END DO
C===========================================================
C=====Composantes de la fonction vectorielle
C=====(�ventuelles) relatives aux sous-r�seaux interstitiels
C===========================================================
C      write(*,*) "N_R_1_M = " , N_R_1_M 
C      write(*,*) "N_R_2_M = " , N_R_2_M 
C      write(*,*) "N_R_3_M = " , N_R_3_M 
C      write(*,*) "N_R_M = " , N_R_M 
C	 stop
      DO I_R = N_R_1_M + N_R_2_M + N_R_3_M + 1 , N_R_M
C--------------------------------------------------------------
C--------------------------------------------------------------
C-----Pour chaque sous-r�seau interstitiel, il y a une �quation
C-----pour chaque esp�ce 1, 2, 3,
C-----ainsi que pour chaque �l�ment d'addition
C-----(sauf pour le sous-r�seau N_R)
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
C-----Cas des �l�ments intrins�ques
C----------------------------------
C----------------------------------
C--------------  
C--------------  
C-----El�ment 1
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrins�ques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 1 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 1 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
       END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
      END IF
C-----------------------
C-----------------------
C-----El�ment 2 �ventuel
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrins�ques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
         J_COUR = IND_D_R_TYP_M ( 2 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^(N_R_1+N_R_2)
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 2 ) THEN
            J_COUR = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
        END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : �l�ment 2 �ventuel
C------------------------------------
C------------------------------------
      END IF
C-----------------------
C-----------------------
C-----El�ment 3 �ventuel
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE TESTEE (avec succ�s - par comparaison avec muVT)
C-----SUR Cr23C6 et TiB2 (i.e. deux cas avec interstitiels intrins�ques)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC INTERSTITIELS INTRINSEQUES (Cr23C6 ou TiB2)
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
        J_COUR = IND_D_R_TYP_M ( 3 , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^I_R
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 1 , N_TYP_M
            J_COUR = IND_D_R_TYP_M ( J_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----pour l'�limination des multiplicateurs de Lagrange
           J_COUR = IND_D_R_TYP_M ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
           if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
              write(*,*)'Fi : i=', I_COMP,' -> xj : J_COUR=', J_COUR
           end if
           MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $   + K_T_M / X_D_R ( 0 , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - - - -
C-----Contributions aux dFi/dxj (divers j � i fix�)
C-----induites par z^N_R_1
C - - - - - - - - - - - - - - - - - - - - - - - - -
          DO J_TYP = 0 , N_TYP_M
C          IF ( J_TYP .NE. 3 ) THEN
            J_COUR
     $    = IND_D_R_TYP_M ( J_TYP , N_R_1_M + N_R_2_M + N_R_3_M )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
        END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
       END IF
C------------------------------------
C------------------------------------
C-----Fin du cas : �l�ment 3 �ventuel
C------------------------------------
C------------------------------------
      END IF
C------------------------------------------
C------------------------------------------
C-----Cas des �l�ments d'addition �ventuels
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
C-----D�riv�es partielles de cette composante de la fonction
C-----(tous les P_J_M pas)
C-----------------------------------------------------------
C*********************************************************
C*****Cf. plus haut les remarques sur le test "J_COUR = 0"
C*****pour �viter les "fausses variables" 1) et 2)
C*********************************************************
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C-----PARTIE A TESTER (par comparaison avec muVT)
C-----SUR FeAl-B ou Fe3Al-C ou TiAl-O
C-----(i.e. des cas avec additions interstitielles)
C-----PARTIE A VERIFIER EVENTUELLEMENT EN ECRIVANT EXPLICITEMENT
C-----LA MATRICE DES DERIVEES PARTIELLES
C-----POUR UN CAS AVEC ADDITIONS INTERSTITIELLES
C-----(FeAl-B ou Fe3Al-C ou TiAl-O)
C-----(comme d�j� fait pour Fe3AlC en l'absence d'additions)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF ( INDIC_CALC_DERIV_M .EQ. 1 ) THEN
           J_COUR = IND_D_R_TYP_M ( I_TYP , I_R )
C - - - - - - - - - - - - - - - - - - - - - - -
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
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
C-----Test "IF ( J_COUR .NE. 0 )" ajout� ici,
C-----afin de pouvoir �carter les DP (E_SC = 0)
C----- = fausses variables 2).
C - - - - - - - - - - - - - - - - - - - - - - -
           IF ( J_COUR .NE. 0 ) THEN
            MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    = MAT_DERIV_NPT_M ( I_COMP , J_COUR )
     $    - K_T_M / Z_TYP_R ( N_R_M )
           END IF
          END DO
C------------------------------------------------
C-----Fin du calcul �ventuel
C-----des d�riv�es partielles de cette composante
C------------------------------------------------
         END IF
C----------------------------------------------
C=====Fin du test "E_DP < 0" pour la composante
C----------------------------------------------
        END IF
C-----------------------------------------------
C-----------------------------------------------
C-----Fin du cas : �l�ments d'addition �ventuels
C-----------------------------------------------
C-----------------------------------------------
        END IF
	 END DO
C=====================================================
C=====Fin de boucle sur les sous-r�seaux interstitiels
C=====================================================
      END DO
C	stop
C###################################################
C#####N_TYP_M composantes de la fonction vectorielle
C#####relatives aux quantit�s de mati�re
C###################################################
C-------------------------
C-----Esp�ce intrins�que 1
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
C-----D�riv�es partielles de cette composante de la fonction
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
C-----D�riv�e par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 1 ) ) / M_MAILLES
      END IF
C------------------------------------
C-----Esp�ce intrins�que 2 �ventuelle
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
C-----D�riv�es partielles de cette composante de la fonction
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
C-----D�riv�e par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 2 ) ) / M_MAILLES
      END IF
C--------------------------------------
C-----Fin du cas : esp�ce intrins�que 2
C--------------------------------------
      END IF
C------------------------------------
C-----Esp�ce intrins�que 3 �ventuelle
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
C-----D�riv�es partielles de cette composante de la fonction
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
C-----D�riv�e par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( 3 ) ) / M_MAILLES
      END IF
C--------------------------------------
C-----Fin du cas : esp�ce intrins�que 3
C--------------------------------------
      END IF
C----------------------------------
C-----El�ments d'addition �ventuels
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
C-----D�riv�es partielles de cette composante de la fonction
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
C-----D�riv�e par rapport au nombre de mailles
           MAT_DERIV_NPT_M ( I_COMP , N_NPT_M )
     $ = ( F_NPT_M ( I_COMP ) + N_AT ( I_TYP ) ) / M_MAILLES
C-------------------------------------------------
C-----Fin de la boucle sur les �l�ments d'addition
C-------------------------------------------------
      END DO
C################################################### 
C#####Derni�re composante de la fonction vectorielle
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
C=====Il s'agit de (2.104) du m�moire d'HDR.
C==========================================
      I_COMP = I_COMP + 1
         F_NPT_M ( I_COMP )
     $ = H_REF_MAILLE_M + W1 + W2 + W3 + K_T_M * SOMME_LOGZ
C-----D�riv�es partielles de cette composante de la fonction
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
C-----D�riv�es du terme "Somme_Log(z)" : sous-r�seaux intrins�ques 1
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
C-----D�riv�es du terme "Somme_Log(z)" : sous-r�seaux intrins�ques 2 �ventuels :
C-----test ci-dessous sans doute requis, car N_R_2 et N_R_3
C-----ne sont pas initialis�s � 0 par d�faut dans le programme principal
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
C-----D�riv�es du terme "Somme_Log(z)" : sous-r�seaux intrins�ques 3 �ventuels :
C-----test ci-dessous sans doute requis, car N_R_2 et N_R_3
C-----ne sont pas initialis�s � 0 par d�faut dans le programme principal
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
C-----D�riv�es du terme "Somme_Log(z)" : sous-r�seaux interstitiels �ventuels
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
C-----Fin du test "calcul des d�riv�es partielles" pour cette composante
        END IF
C-----------------
C-----V�rification
C-----------------
      IF ( I_COMP .NE. N_NPT_M ) THEN
	write ( * , * ) '-----------------------------------------'
	write ( * , * ) 'En fin de f_NR :'
	write ( * , * ) "nombre d'�quations = " , I_COMP
	write ( * , * ) 'diff�rent du nombre attendu = ' , N_NPT_M
	write ( * , * ) '-----------------------------------------'
	CALL INTERRUPTION
      END IF
C-------------
C-----Ecriture
C-------------
      if ( INDIC_ECRIT_DETAILS_F_NR .EQ. 1 ) then
 	write ( * , * ) 'En fin de proc�dure F_NR :' 
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
C     I Calcul des concentrations en porteurs intrins�ques n et        I
C     I et de leur d�riv�e                                             I
C     I � partir de la densit� d'�tats                                 I
C     I en fonction du potentiel chimique �lectronique                 I
C     I								       I
C     ==================================================================
C     ==================================================================
C     I Derni�re mise � jour : 11/07/2013                              I
C     I Saint Beno�t                                                   I
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
C#####D�claration des tableaux arguments de la proc�dure
C#######################################################
C--------------------------------------------------
C-----Pour chaque point, valeur de l'�nergie et DdE
C--------------------------------------------------
      REAL * 8 TAB_DDE_M ( N_PAS_M , 2 )
C#####################################################
C#####D�claration des tableaux internes � la proc�dure
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
C-----Terme pour la d�riv�e de n
C-------------------------------
        IF ( E_COUR .GE. E_MIN_BC_M ) THEN
          D_C_VOL_N_M = D_C_VOL_N_M 
     $   + DDE_COUR * ( E_SUIV - E_COUR )
     $   * FACT_COUR / ( 1.D0 + FACT_COUR ) * * 2
     $   / K_T
        END IF
C-------------------------------
C-----Terme pour la d�riv�e de p
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

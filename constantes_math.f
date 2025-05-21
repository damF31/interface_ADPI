C     ##################################################################

      MODULE CONSTANTES_MATH

C     ##################################################################

C=========================
C=====Tenseur de Kronecker
C=========================

      REAL * 8 , DIMENSION ( 3 , 3 ) :: TENS_DELTA

      DATA TENS_DELTA

     $	   /   1.D0 ,   0.D0 ,   0.D0 ,

     $	       0.D0 ,   1.D0 ,   0.D0 ,

     $	       0.D0 ,   0.D0 ,   1.D0 /

C===========================
C=====Tenseur de permutation
C===========================

      REAL * 8 , DIMENSION ( 3 , 3 , 3 ) :: TENS_PERM

      DATA TENS_PERM

C--------------------------------
C-----TENS_PERM ( I , J , K = 1 )
C--------------------------------

         ! Col. I = 1		Col. I = 2	Col. I = 3

     $     /   0.D0 , 		0.D0 ,		0.D0 ,	! Ligne J = 1

     $         0.D0 ,   	0.D0 ,         -1.D0 ,	! Ligne J = 2

     $         0.D0 ,   	1.D0 ,		0.D0 ,	! Ligne J = 3

C--------------------------------
C-----TENS_PERM ( I , J , K = 2 )
C--------------------------------

         ! Col. I = 1           Col. I = 2      Col. I = 3

     $         0.D0 ,           0.D0 ,          1.D0 ,  ! Ligne J = 1

     $         0.D0 ,           0.D0 ,          0.D0 ,  ! Ligne J = 2

     $       - 1.D0 ,           0.D0 ,          0.D0 ,  ! Ligne J = 3

C--------------------------------
C-----TENS_PERM ( I , J , K = 3 )
C--------------------------------

         ! Col. I = 1           Col. I = 2      Col. I = 3

     $         0.D0 ,         - 1.D0 ,          0.D0 ,  ! Ligne J = 1

     $         1.D0 ,           0.D0 ,          0.D0 ,  ! Ligne J = 2

     $         0.D0 ,           0.D0 ,          0.D0 /  ! Ligne J = 3

      END MODULE CONSTANTES_MATH

      INTEGER MAX_N_SPIN_CORR_VECTORS, MAX_LEGS_WITH_SPIN_CORR
      PARAMETER (MAX_N_SPIN_CORR_VECTORS=20, MAX_LEGS_WITH_SPIN_CORR=2)
C     Lists the spin correlation vectors defined for each external leg
      COMPLEX*16 SPIN_CORR_VECTORS(NEXTERNAL,MAX_N_SPIN_CORR_VECTORS,4)
C     Indicates the number of spin correlations vectors defined for
C      each external leg
      INTEGER N_SPIN_CORR_VECTORS(NEXTERNAL)

      INTEGER MAX_SPIN_CORR_RUNS
C     Store the list of combination of spin_corr_vectors with which to
C      enhance the loop over helicity combinations
C     Also we use a temporary parameter const buffer to store the
C      value of max_spin_corr_runs because f2py
C     can't realize it is constant otherwise
C     INTEGER TMPCONST
C     PARAMETER(TMPCONST=MAX_N_SPIN_CORR_VECTORS**(MAX_LEGS_WITH_SPIN_C
C     ORR))
C     PARAMETER(MAX_SPIN_CORR_RUNS=TMPCONST)
      PARAMETER(MAX_SPIN_CORR_RUNS=MAX_N_SPIN_CORR_VECTORS**
     $ (MAX_LEGS_WITH_SPIN_CORR))
C     Store the number of spin-correlation runs defined by the user.
C     A run is just a pass through the helas calls for computing the
C      integrand for a specific helicity configuration
      INTEGER N_SPIN_CORR_RUNS
      INTEGER SPIN_CORR_RUNS(0:MAX_SPIN_CORR_RUNS,NEXTERNAL)
      COMMON/AA_WPWM_SPIN_CORRELATION_DATA/SPIN_CORR_VECTORS,
     $  N_SPIN_CORR_VECTORS, SPIN_CORR_RUNS, N_SPIN_CORR_RUNS

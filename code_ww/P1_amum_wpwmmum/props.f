      SUBROUTINE PROPS_BORN
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'genps.inc'
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'run.inc'
      INCLUDE 'coupl.inc'
      DOUBLE PRECISION PRMASS(-NEXTERNAL:0,LMAXCONFIGS)
      DOUBLE PRECISION PRWIDTH(-NEXTERNAL:0,LMAXCONFIGS)
      INTEGER POW(-NEXTERNAL:0,LMAXCONFIGS)
      DOUBLE PRECISION, PARAMETER :: ZERO = 0D0
      COMMON/PROPS/PRMASS,PRWIDTH,POW
      PRMASS(-1,1)  = ZERO
      PRWIDTH(-1,1) = ZERO
      POW(-1,1) = 2
      PRMASS(-2,1)  = ZERO
      PRWIDTH(-2,1) = ZERO
      POW(-2,1) = 1
      PRMASS(-1,2)  = ABS(MDL_MZ)
      PRWIDTH(-1,2) = ABS(MDL_WZ)
      POW(-1,2) = 2
      PRMASS(-2,2)  = ZERO
      PRWIDTH(-2,2) = ZERO
      POW(-2,2) = 1
      PRMASS(-1,3)  = ZERO
      PRWIDTH(-1,3) = ZERO
      POW(-1,3) = 1
      PRMASS(-2,3)  = ZERO
      PRWIDTH(-2,3) = ZERO
      POW(-2,3) = 1
      PRMASS(-1,4)  = ABS(MDL_MW)
      PRWIDTH(-1,4) = ABS(MDL_WW)
      POW(-1,4) = 2
      PRMASS(-2,4)  = ZERO
      PRWIDTH(-2,4) = ZERO
      POW(-2,4) = 1
      PRMASS(-1,5)  = ABS(MDL_MW)
      PRWIDTH(-1,5) = ABS(MDL_WW)
      POW(-1,5) = 2
      PRMASS(-2,5)  = ZERO
      PRWIDTH(-2,5) = ZERO
      POW(-2,5) = 2
      PRMASS(-1,6)  = ABS(MDL_MW)
      PRWIDTH(-1,6) = ABS(MDL_WW)
      POW(-1,6) = 2
      PRMASS(-2,6)  = ABS(MDL_MZ)
      PRWIDTH(-2,6) = ABS(MDL_WZ)
      POW(-2,6) = 2
      PRMASS(-1,7)  = ABS(MDL_MW)
      PRWIDTH(-1,7) = ABS(MDL_WW)
      POW(-1,7) = 2
      PRMASS(-2,7)  = ZERO
      PRWIDTH(-2,7) = ZERO
      POW(-2,7) = 2
      PRMASS(-1,8)  = ABS(MDL_MW)
      PRWIDTH(-1,8) = ABS(MDL_WW)
      POW(-1,8) = 2
      PRMASS(-2,8)  = ABS(MDL_MZ)
      PRWIDTH(-2,8) = ABS(MDL_WZ)
      POW(-2,8) = 2
      PRMASS(-1,9)  = ABS(MDL_MW)
      PRWIDTH(-1,9) = ABS(MDL_WW)
      POW(-1,9) = 2
      PRMASS(-2,9)  = ZERO
      PRWIDTH(-2,9) = ZERO
      POW(-2,9) = 2
      PRMASS(-1,10)  = ABS(MDL_MW)
      PRWIDTH(-1,10) = ABS(MDL_WW)
      POW(-1,10) = 2
      PRMASS(-2,10)  = ABS(MDL_MZ)
      PRWIDTH(-2,10) = ABS(MDL_WZ)
      POW(-2,10) = 2
      PRMASS(-1,11)  = ABS(MDL_MW)
      PRWIDTH(-1,11) = ABS(MDL_WW)
      POW(-1,11) = 2
      PRMASS(-2,11)  = ZERO
      PRWIDTH(-2,11) = ZERO
      POW(-2,11) = 2
      PRMASS(-1,12)  = ABS(MDL_MW)
      PRWIDTH(-1,12) = ABS(MDL_WW)
      POW(-1,12) = 2
      PRMASS(-2,12)  = ABS(MDL_MZ)
      PRWIDTH(-2,12) = ABS(MDL_WZ)
      POW(-2,12) = 2
      PRMASS(-1,13)  = ZERO
      PRWIDTH(-1,13) = ZERO
      POW(-1,13) = 1
      PRMASS(-2,13)  = ABS(MDL_MW)
      PRWIDTH(-2,13) = ABS(MDL_WW)
      POW(-2,13) = 2
      PRMASS(-1,14)  = ZERO
      PRWIDTH(-1,14) = ZERO
      POW(-1,14) = 1
      PRMASS(-2,14)  = ZERO
      PRWIDTH(-2,14) = ZERO
      POW(-2,14) = 1
      PRMASS(-1,15)  = ZERO
      PRWIDTH(-1,15) = ZERO
      POW(-1,15) = 2
      PRMASS(-2,15)  = ZERO
      PRWIDTH(-2,15) = ZERO
      POW(-2,15) = 1
      PRMASS(-1,16)  = ABS(MDL_MZ)
      PRWIDTH(-1,16) = ABS(MDL_WZ)
      POW(-1,16) = 2
      PRMASS(-2,16)  = ZERO
      PRWIDTH(-2,16) = ZERO
      POW(-2,16) = 1
      RETURN
      END

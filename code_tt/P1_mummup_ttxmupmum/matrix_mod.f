      SUBROUTINE ME_ACCESSOR_HOOK_4(P,HEL,USER_ALPHAS,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=7)
      REAL*8 PI
      PARAMETER (PI= 3.141592653589793D0)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQAMPSO)
      INTEGER HEL
      DOUBLE PRECISION USER_ALPHAS
CF2PY INTENT(IN)  :: P
CF2PY INTENT(IN)  :: HEL
CF2PY INTENT(IN)  :: USER_ALPHAS
CF2PY INTENT(OUT) :: ANS

      REAL*8 THIS_G

      INCLUDE 'coupl.inc'

C     ----------
C     BEGIN CODE
C     ----------

      IF (USER_ALPHAS.GT.0.0D0) THEN
        THIS_G = 2* DSQRT(USER_ALPHAS*PI)
        IF (THIS_G.NE.G) THEN
          G = THIS_G
          CALL UPDATE_AS_PARAM()
        ENDIF
      ENDIF

      CALL SMATRIXHEL_SPLITORDERS_4(P,HEL,ANS)

      END

      SUBROUTINE SMATRIX_4(P,ANS_SUMMED)
C     
C     Simple routine wrapper to provide the same interface for
C     backward compatibility for usage without split orders.
C     
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=7)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL), ANS_SUMMED
C     
C     VARIABLES
C     
      INTEGER I
      REAL*8 ANS(0:NSQAMPSO)
C     
C     BEGIN CODE
C     
      CALL SMATRIX_SPLITORDERS_4(P,ANS)
      ANS_SUMMED=ANS(0)

      END

      SUBROUTINE SMATRIXHEL_4(P,HEL,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER HEL
C     
C     GLOBAL VARIABLES
C     
      INTEGER USERHEL
      COMMON/HELUSERCHOICE_4/USERHEL
C     ----------
C     BEGIN CODE
C     ----------
      USERHEL=HEL
      CALL SMATRIX_4(P,ANS)
      USERHEL=-1

      END

C     Give access to the helicity definition to the f2py API.
      SUBROUTINE GET_HELICITY_DEFINITIONS_4(NHEL_OUT)
      IMPLICIT NONE

      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)
      INTEGER NHEL(NEXTERNAL,NCOMB)
      COMMON/BORN_HEL_CONFIGS_4/NHEL

      INTEGER NHEL_OUT(NCOMB,NEXTERNAL)
CF2PY INTENT(OUT) :: NHEL_OUT

      INTEGER I,J

      DO I=1,NEXTERNAL
        DO J=1,NCOMB
          NHEL_OUT(J,I) = NHEL(I,J)
        ENDDO
      ENDDO

      END

      SUBROUTINE SMATRIX_SPLITORDERS_4(P,ANS)
C     
C     Generated by MadGraph5_aMC@NLO v. 2.5.3, 2017-03-09
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     MadGraph StandAlone Version
C     
C     Returns amplitude squared summed/avg over colors
C     and helicities
C     for the point in phase space P(0:3,NEXTERNAL)
C     
C     Process: mu- mu+ > t t~ mu+ mu- QCD^2<=6 QED^2<=8 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NINITIAL
      PARAMETER (NINITIAL=2)
      INTEGER NPOLENTRIES
      PARAMETER (NPOLENTRIES=(NEXTERNAL+1)*6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=7)
      INTEGER HELAVGFACTOR
      PARAMETER (HELAVGFACTOR=4)
      LOGICAL CHOSEN_SO_CONFIGS(NSQAMPSO)
      DATA CHOSEN_SO_CONFIGS/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.
     $ ,.TRUE./
      COMMON/CHOSEN_BORN_SQSO_4/CHOSEN_SO_CONFIGS
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQAMPSO)
C     
C     LOCAL VARIABLES 
C     
      LOGICAL DO_INCLUDE_HEL_CONTRIB
      INTEGER NTRY
      REAL*8 T(NSQAMPSO), BUFF
      INTEGER IHEL,IDEN, I, J
C     For a 1>N process, them BEAMTWO_HELAVGFACTOR would be set to 1.
      INTEGER BEAMS_HELAVGFACTOR(2)
      DATA (BEAMS_HELAVGFACTOR(I),I=1,2)/2,2/
      INTEGER JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB)
      DATA NTRY/0/
      DATA GOODHEL/NCOMB*.FALSE./

      INTEGER HELICITIES(NEXTERNAL)

      DATA IDEN/ 4/
C     
C     GLOBAL VARIABLES
C     
      INTEGER NHEL(NEXTERNAL,NCOMB)
      DATA (NHEL(I,   1),I=1,6) / 1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,   2),I=1,6) / 1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,   3),I=1,6) / 1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,   4),I=1,6) / 1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,   5),I=1,6) / 1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,   6),I=1,6) / 1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,   7),I=1,6) / 1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,   8),I=1,6) / 1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,   9),I=1,6) / 1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  10),I=1,6) / 1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  11),I=1,6) / 1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  12),I=1,6) / 1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  13),I=1,6) / 1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  14),I=1,6) / 1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  15),I=1,6) / 1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  16),I=1,6) / 1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  17),I=1,6) / 1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  18),I=1,6) / 1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  19),I=1,6) / 1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  20),I=1,6) / 1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  21),I=1,6) / 1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  22),I=1,6) / 1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  23),I=1,6) / 1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  24),I=1,6) / 1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  25),I=1,6) / 1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  26),I=1,6) / 1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  27),I=1,6) / 1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  28),I=1,6) / 1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  29),I=1,6) / 1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  30),I=1,6) / 1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  31),I=1,6) / 1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  32),I=1,6) / 1, 1, 1,-1,-1, 1/
      DATA (NHEL(I,  33),I=1,6) /-1,-1,-1, 1, 1,-1/
      DATA (NHEL(I,  34),I=1,6) /-1,-1,-1, 1, 1, 1/
      DATA (NHEL(I,  35),I=1,6) /-1,-1,-1, 1,-1,-1/
      DATA (NHEL(I,  36),I=1,6) /-1,-1,-1, 1,-1, 1/
      DATA (NHEL(I,  37),I=1,6) /-1,-1,-1,-1, 1,-1/
      DATA (NHEL(I,  38),I=1,6) /-1,-1,-1,-1, 1, 1/
      DATA (NHEL(I,  39),I=1,6) /-1,-1,-1,-1,-1,-1/
      DATA (NHEL(I,  40),I=1,6) /-1,-1,-1,-1,-1, 1/
      DATA (NHEL(I,  41),I=1,6) /-1,-1, 1, 1, 1,-1/
      DATA (NHEL(I,  42),I=1,6) /-1,-1, 1, 1, 1, 1/
      DATA (NHEL(I,  43),I=1,6) /-1,-1, 1, 1,-1,-1/
      DATA (NHEL(I,  44),I=1,6) /-1,-1, 1, 1,-1, 1/
      DATA (NHEL(I,  45),I=1,6) /-1,-1, 1,-1, 1,-1/
      DATA (NHEL(I,  46),I=1,6) /-1,-1, 1,-1, 1, 1/
      DATA (NHEL(I,  47),I=1,6) /-1,-1, 1,-1,-1,-1/
      DATA (NHEL(I,  48),I=1,6) /-1,-1, 1,-1,-1, 1/
      DATA (NHEL(I,  49),I=1,6) /-1, 1,-1, 1, 1,-1/
      DATA (NHEL(I,  50),I=1,6) /-1, 1,-1, 1, 1, 1/
      DATA (NHEL(I,  51),I=1,6) /-1, 1,-1, 1,-1,-1/
      DATA (NHEL(I,  52),I=1,6) /-1, 1,-1, 1,-1, 1/
      DATA (NHEL(I,  53),I=1,6) /-1, 1,-1,-1, 1,-1/
      DATA (NHEL(I,  54),I=1,6) /-1, 1,-1,-1, 1, 1/
      DATA (NHEL(I,  55),I=1,6) /-1, 1,-1,-1,-1,-1/
      DATA (NHEL(I,  56),I=1,6) /-1, 1,-1,-1,-1, 1/
      DATA (NHEL(I,  57),I=1,6) /-1, 1, 1, 1, 1,-1/
      DATA (NHEL(I,  58),I=1,6) /-1, 1, 1, 1, 1, 1/
      DATA (NHEL(I,  59),I=1,6) /-1, 1, 1, 1,-1,-1/
      DATA (NHEL(I,  60),I=1,6) /-1, 1, 1, 1,-1, 1/
      DATA (NHEL(I,  61),I=1,6) /-1, 1, 1,-1, 1,-1/
      DATA (NHEL(I,  62),I=1,6) /-1, 1, 1,-1, 1, 1/
      DATA (NHEL(I,  63),I=1,6) /-1, 1, 1,-1,-1,-1/
      DATA (NHEL(I,  64),I=1,6) /-1, 1, 1,-1,-1, 1/
      COMMON/BORN_HEL_CONFIGS_4/NHEL

      INTEGER USERHEL
      DATA USERHEL/-1/
      COMMON/HELUSERCHOICE_4/USERHEL

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      DATA ((POLARIZATIONS(I,J),I=0,NEXTERNAL),J=0,5)/NPOLENTRIES*-1/
      COMMON/BORN_BEAM_POL_4/POLARIZATIONS



C     
C     FUNCTIONS
C     
      LOGICAL IS_BORN_HEL_SELECTED_4

C     ----------
C     BEGIN CODE
C     ----------


      NTRY=NTRY+1
      DO IHEL=1,NEXTERNAL
        JC(IHEL) = +1
      ENDDO
      DO I=1,NSQAMPSO
        ANS(I) = 0D0
      ENDDO
C     When spin-2 particles are involved, the Helicity filtering is
C      dangerous for the 2->1 topology.
C     This is because depending on the MC setup the initial PS points
C      have back-to-back initial states
C     for which some of the spin-2 helicity configurations are zero.
C      But they are no longer zero
C     if the point is boosted on the z-axis. Remember that HELAS
C      helicity amplitudes are no longer
C     lorentz invariant with expternal spin-2 particles (only the
C      helicity sum is).
C     For this reason, we simply remove the filterin when there is
C      only three external particles.
      IF (NEXTERNAL.LE.3) THEN
        DO IHEL=1,NCOMB
          GOODHEL(IHEL)=.TRUE.
        ENDDO
      ENDIF


      DO IHEL=1,NCOMB
        IF (USERHEL.EQ.-1.OR.USERHEL.EQ.IHEL) THEN
          IF (GOODHEL(IHEL) .OR. NTRY .LT. 2 .OR.USERHEL.NE.-1) THEN
            IF(NTRY.GE.2.AND.POLARIZATIONS(0,0).NE.-1.AND.
     $       (.NOT.IS_BORN_HEL_SELECTED_4(IHEL))) THEN
              CYCLE
            ENDIF
            DO_INCLUDE_HEL_CONTRIB = (POLARIZATIONS(0,0).EQ.
     $       -1.OR.IS_BORN_HEL_SELECTED_4(IHEL))
            DO I=1,NEXTERNAL
              HELICITIES(I) = NHEL(I,IHEL)
            ENDDO
            CALL MATRIX_4(P ,HELICITIES(1),JC(1), T)
            BUFF=0D0
            DO I=1,NSQAMPSO
              IF(DO_INCLUDE_HEL_CONTRIB) THEN
                ANS(I)=ANS(I)+T(I)
              ENDIF
              BUFF=BUFF+T(I)
            ENDDO
            IF (BUFF .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL)) THEN
              GOODHEL(IHEL)=.TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ANS(0)=0.0D0

      DO I=1,NSQAMPSO
        ANS(I)=ANS(I)/DBLE(IDEN)
        IF (CHOSEN_SO_CONFIGS(I)) THEN
          ANS(0)=ANS(0)+ANS(I)
        ENDIF
      ENDDO
      IF(USERHEL.NE.-1) THEN
        DO I=0,NSQAMPSO
          ANS(I)=ANS(I)*HELAVGFACTOR
        ENDDO
      ELSE
        DO J=1,NINITIAL
          IF (POLARIZATIONS(J,0).NE.-1) THEN
            DO I=0,NSQAMPSO
              ANS(I)=ANS(I)*BEAMS_HELAVGFACTOR(J)
              ANS(I)=ANS(I)/POLARIZATIONS(J,0)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      END

      SUBROUTINE SMATRIXHEL_SPLITORDERS_4(P,HEL,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER                 NCOMB
      PARAMETER (             NCOMB=64)
      INTEGER NSQAMPSO
      PARAMETER (NSQAMPSO=7)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS(0:NSQAMPSO)
      INTEGER HEL
CF2PY INTENT(IN)  :: P
CF2PY INTENT(IN)  :: HEL
CF2PY INTENT(OUT) :: ANS

C     
C     GLOBAL VARIABLES
C     
      INTEGER USERHEL
      COMMON/HELUSERCHOICE_4/USERHEL
C     ----------
C     BEGIN CODE
C     ----------
      USERHEL=HEL
      CALL SMATRIX_SPLITORDERS_4(P,ANS)
      USERHEL=-1

      END

      SUBROUTINE MATRIX_4(P,NHEL,IC,RES)

C     
C     Generated by MadGraph5_aMC@NLO v. 2.5.3, 2017-03-09
C     By the MadGraph5_aMC@NLO Development Team
C     Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
C     
C     Returns amplitude squared summed/avg over colors
C     for the point with external lines W(0:6,NEXTERNAL)
C     
C     Process: mu- mu+ > t t~ mu+ mu- QCD^2<=6 QED^2<=8 @1
C     
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NGRAPHS
      PARAMETER (NGRAPHS=50)
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER (NWAVEFUNCS=17, NCOLOR=1)
      INTEGER NAMPSO, NSQAMPSO
      PARAMETER (NAMPSO=4, NSQAMPSO=7)
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0D0,1D0))
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      REAL*8 RES(NSQAMPSO)
C     
C     LOCAL VARIABLES 
C     
      INTEGER I,J,M,N, SQSOIND
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS)
      COMPLEX*16 JAMP(NCOLOR,NAMPSO)
      COMPLEX*16 W(20,NWAVEFUNCS)
      COMPLEX*16 DUM0,DUM1
      DATA DUM0, DUM1/(0D0, 0D0), (1D0, 0D0)/
      logical sameflav_diags
      common /to_sameflav/sameflav_diags
C     
C     FUNCTION
C     
      INTEGER SQSOINDEX_4
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'coupl.inc'
C     
C     COLOR DATA
C     
      DATA DENOM(1)/1/
      DATA (CF(I,  1),I=  1,  1) /    3/
C     1 T(3,4)
C     ----------
C     BEGIN CODE
C     ----------


      CALL IXXXXX(P(0,1),ZERO,NHEL(1),+1*IC(1),W(1,1))
      CALL OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))
      CALL OXXXXX(P(0,3),MDL_MT,NHEL(3),+1*IC(3),W(1,3))
      CALL IXXXXX(P(0,4),MDL_MT,NHEL(4),-1*IC(4),W(1,4))
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))
      CALL OXXXXX(P(0,6),ZERO,NHEL(6),+1*IC(6),W(1,6))
      CALL FFV1P0_3(W(1,1),W(1,2),GC_3,ZERO,ZERO,W(1,7))
      CALL FFV1P0_3(W(1,5),W(1,6),GC_3,ZERO,ZERO,W(1,8))
      CALL FFV1_1(W(1,3),W(1,7),GC_2,MDL_MT,MDL_WT,W(1,9))
C     Amplitude(s) for diagram number 1
      CALL FFV1_0(W(1,4),W(1,9),W(1,8),GC_2,AMP(1))
      CALL FFV1_2(W(1,4),W(1,7),GC_2,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for diagram number 2
      CALL FFV1_0(W(1,10),W(1,3),W(1,8),GC_2,AMP(2))
      CALL FFV2_4P0_3(W(1,5),W(1,6),GC_68,GC_77,MDL_MZ,MDL_WZ,W(1,11))
C     Amplitude(s) for diagram number 3
      CALL FFV2_5_0(W(1,4),W(1,9),W(1,11),GC_69,GC_76,AMP(3))
C     Amplitude(s) for diagram number 4
      CALL FFV2_5_0(W(1,10),W(1,3),W(1,11),GC_69,GC_76,AMP(4))
      CALL FFV2_4P0_3(W(1,1),W(1,2),GC_68,GC_77,MDL_MZ,MDL_WZ,W(1,10))
      CALL FFV2_5_1(W(1,3),W(1,10),GC_69,GC_76,MDL_MT,MDL_WT,W(1,9))
C     Amplitude(s) for diagram number 5
      CALL FFV1_0(W(1,4),W(1,9),W(1,8),GC_2,AMP(5))
      CALL FFV2_5_2(W(1,4),W(1,10),GC_69,GC_76,MDL_MT,MDL_WT,W(1,12))
C     Amplitude(s) for diagram number 6
      CALL FFV1_0(W(1,12),W(1,3),W(1,8),GC_2,AMP(6))
C     Amplitude(s) for diagram number 7
      CALL FFV2_5_0(W(1,4),W(1,9),W(1,11),GC_69,GC_76,AMP(7))
C     Amplitude(s) for diagram number 8
      CALL FFV2_5_0(W(1,12),W(1,3),W(1,11),GC_69,GC_76,AMP(8))
      CALL FFS4_3(W(1,4),W(1,3),GC_116,MDL_MH,MDL_WH,W(1,12))
C     Amplitude(s) for diagram number 9
      CALL VVS1_0(W(1,10),W(1,11),W(1,12),GC_99,AMP(9))
      CALL FFV1P0_3(W(1,4),W(1,3),GC_2,ZERO,ZERO,W(1,9))
      CALL FFV1_2(W(1,5),W(1,7),GC_3,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 10
      CALL FFV1_0(W(1,13),W(1,6),W(1,9),GC_3,AMP(10))
      CALL FFV1_1(W(1,6),W(1,7),GC_3,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 11
      CALL FFV1_0(W(1,5),W(1,14),W(1,9),GC_3,AMP(11))
      CALL FFV2_5P0_3(W(1,4),W(1,3),GC_69,GC_76,MDL_MZ,MDL_WZ,W(1,7))
C     Amplitude(s) for diagram number 12
      CALL FFV2_4_0(W(1,13),W(1,6),W(1,7),GC_68,GC_77,AMP(12))
C     Amplitude(s) for diagram number 13
      CALL FFV2_4_0(W(1,5),W(1,14),W(1,7),GC_68,GC_77,AMP(13))
      CALL FFV2_4_2(W(1,5),W(1,10),GC_68,GC_77,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 14
      CALL FFV1_0(W(1,14),W(1,6),W(1,9),GC_3,AMP(14))
      CALL FFV2_4_1(W(1,6),W(1,10),GC_68,GC_77,ZERO,ZERO,W(1,13))
C     Amplitude(s) for diagram number 15
      CALL FFV1_0(W(1,5),W(1,13),W(1,9),GC_3,AMP(15))
C     Amplitude(s) for diagram number 16
      CALL FFV2_4_0(W(1,14),W(1,6),W(1,7),GC_68,GC_77,AMP(16))
C     Amplitude(s) for diagram number 17
      CALL FFV2_4_0(W(1,5),W(1,13),W(1,7),GC_68,GC_77,AMP(17))
      CALL FFV1P0_3(W(1,1),W(1,6),GC_3,ZERO,ZERO,W(1,13))
      CALL FFV1P0_3(W(1,5),W(1,2),GC_3,ZERO,ZERO,W(1,14))
      CALL FFV1_1(W(1,3),W(1,13),GC_2,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for diagram number 18
      CALL FFV1_0(W(1,4),W(1,10),W(1,14),GC_2,AMP(18))
      CALL FFV1_2(W(1,4),W(1,13),GC_2,MDL_MT,MDL_WT,W(1,15))
C     Amplitude(s) for diagram number 19
      CALL FFV1_0(W(1,15),W(1,3),W(1,14),GC_2,AMP(19))
      CALL FFV2_4P0_3(W(1,5),W(1,2),GC_68,GC_77,MDL_MZ,MDL_WZ,W(1,16))
C     Amplitude(s) for diagram number 20
      CALL FFV2_5_0(W(1,4),W(1,10),W(1,16),GC_69,GC_76,AMP(20))
C     Amplitude(s) for diagram number 21
      CALL FFV2_5_0(W(1,15),W(1,3),W(1,16),GC_69,GC_76,AMP(21))
      CALL FFV2_4P0_3(W(1,1),W(1,6),GC_68,GC_77,MDL_MZ,MDL_WZ,W(1,15))
      CALL FFV2_5_1(W(1,3),W(1,15),GC_69,GC_76,MDL_MT,MDL_WT,W(1,10))
C     Amplitude(s) for diagram number 22
      CALL FFV1_0(W(1,4),W(1,10),W(1,14),GC_2,AMP(22))
      CALL FFV2_5_2(W(1,4),W(1,15),GC_69,GC_76,MDL_MT,MDL_WT,W(1,17))
C     Amplitude(s) for diagram number 23
      CALL FFV1_0(W(1,17),W(1,3),W(1,14),GC_2,AMP(23))
C     Amplitude(s) for diagram number 24
      CALL FFV2_5_0(W(1,4),W(1,10),W(1,16),GC_69,GC_76,AMP(24))
C     Amplitude(s) for diagram number 25
      CALL FFV2_5_0(W(1,17),W(1,3),W(1,16),GC_69,GC_76,AMP(25))
C     Amplitude(s) for diagram number 26
      CALL VVS1_0(W(1,15),W(1,16),W(1,12),GC_99,AMP(26))
      CALL FFV1_2(W(1,5),W(1,13),GC_3,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 27
      CALL FFV1_0(W(1,12),W(1,2),W(1,9),GC_3,AMP(27))
      CALL FFV1_1(W(1,2),W(1,13),GC_3,ZERO,ZERO,W(1,17))
C     Amplitude(s) for diagram number 28
      CALL FFV1_0(W(1,5),W(1,17),W(1,9),GC_3,AMP(28))
C     Amplitude(s) for diagram number 29
      CALL FFV2_4_0(W(1,12),W(1,2),W(1,7),GC_68,GC_77,AMP(29))
C     Amplitude(s) for diagram number 30
      CALL FFV2_4_0(W(1,5),W(1,17),W(1,7),GC_68,GC_77,AMP(30))
      CALL FFV2_4_2(W(1,5),W(1,15),GC_68,GC_77,ZERO,ZERO,W(1,17))
C     Amplitude(s) for diagram number 31
      CALL FFV1_0(W(1,17),W(1,2),W(1,9),GC_3,AMP(31))
      CALL FFV2_4_1(W(1,2),W(1,15),GC_68,GC_77,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 32
      CALL FFV1_0(W(1,5),W(1,12),W(1,9),GC_3,AMP(32))
C     Amplitude(s) for diagram number 33
      CALL FFV2_4_0(W(1,17),W(1,2),W(1,7),GC_68,GC_77,AMP(33))
C     Amplitude(s) for diagram number 34
      CALL FFV2_4_0(W(1,5),W(1,12),W(1,7),GC_68,GC_77,AMP(34))
      CALL FFV1_2(W(1,1),W(1,14),GC_3,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 35
      CALL FFV1_0(W(1,12),W(1,6),W(1,9),GC_3,AMP(35))
      CALL FFV1_2(W(1,1),W(1,9),GC_3,ZERO,ZERO,W(1,5))
C     Amplitude(s) for diagram number 36
      CALL FFV1_0(W(1,5),W(1,6),W(1,14),GC_3,AMP(36))
C     Amplitude(s) for diagram number 37
      CALL FFV2_4_0(W(1,12),W(1,6),W(1,7),GC_68,GC_77,AMP(37))
      CALL FFV2_4_2(W(1,1),W(1,7),GC_68,GC_77,ZERO,ZERO,W(1,12))
C     Amplitude(s) for diagram number 38
      CALL FFV1_0(W(1,12),W(1,6),W(1,14),GC_3,AMP(38))
      CALL FFV2_4_2(W(1,1),W(1,16),GC_68,GC_77,ZERO,ZERO,W(1,14))
C     Amplitude(s) for diagram number 39
      CALL FFV1_0(W(1,14),W(1,6),W(1,9),GC_3,AMP(39))
C     Amplitude(s) for diagram number 40
      CALL FFV2_4_0(W(1,5),W(1,6),W(1,16),GC_68,GC_77,AMP(40))
C     Amplitude(s) for diagram number 41
      CALL FFV2_4_0(W(1,14),W(1,6),W(1,7),GC_68,GC_77,AMP(41))
C     Amplitude(s) for diagram number 42
      CALL FFV2_4_0(W(1,12),W(1,6),W(1,16),GC_68,GC_77,AMP(42))
      CALL FFV1_2(W(1,1),W(1,8),GC_3,ZERO,ZERO,W(1,16))
C     Amplitude(s) for diagram number 43
      CALL FFV1_0(W(1,16),W(1,2),W(1,9),GC_3,AMP(43))
C     Amplitude(s) for diagram number 44
      CALL FFV1_0(W(1,5),W(1,2),W(1,8),GC_3,AMP(44))
C     Amplitude(s) for diagram number 45
      CALL FFV2_4_0(W(1,16),W(1,2),W(1,7),GC_68,GC_77,AMP(45))
C     Amplitude(s) for diagram number 46
      CALL FFV1_0(W(1,12),W(1,2),W(1,8),GC_3,AMP(46))
      CALL FFV2_4_2(W(1,1),W(1,11),GC_68,GC_77,ZERO,ZERO,W(1,8))
C     Amplitude(s) for diagram number 47
      CALL FFV1_0(W(1,8),W(1,2),W(1,9),GC_3,AMP(47))
C     Amplitude(s) for diagram number 48
      CALL FFV2_4_0(W(1,5),W(1,2),W(1,11),GC_68,GC_77,AMP(48))
C     Amplitude(s) for diagram number 49
      CALL FFV2_4_0(W(1,8),W(1,2),W(1,7),GC_68,GC_77,AMP(49))
C     Amplitude(s) for diagram number 50
      CALL FFV2_4_0(W(1,12),W(1,2),W(1,11),GC_68,GC_77,AMP(50))
      if (.not.sameflav_diags) then
          amp(1:17)=dcmplx(0d0,0d0)
          amp(43:50)=dcmplx(0d0,0d0)
      endif
C     JAMPs contributing to orders QCD=3 QED=1
      JAMP(1,1)=-AMP(10)-AMP(11)+AMP(27)+AMP(28)+AMP(35)+AMP(36)-AMP
     $ (43)-AMP(44)
C     JAMPs contributing to orders QCD=2 QED=2
      JAMP(1,2)=-AMP(1)-AMP(2)-AMP(12)-AMP(13)+AMP(18)+AMP(19)+AMP(29)
     $ +AMP(30)+AMP(37)+AMP(38)-AMP(45)-AMP(46)
C     JAMPs contributing to orders QCD=1 QED=3
      JAMP(1,3)=-AMP(3)-AMP(4)-AMP(5)-AMP(6)-AMP(14)-AMP(15)+AMP(20)
     $ +AMP(21)+AMP(22)+AMP(23)+AMP(31)+AMP(32)+AMP(39)+AMP(40)-AMP(47)
     $ -AMP(48)
C     JAMPs contributing to orders QCD=0 QED=4
      JAMP(1,4)=-AMP(7)-AMP(8)-AMP(9)-AMP(16)-AMP(17)+AMP(24)+AMP(25)
     $ +AMP(26)+AMP(33)+AMP(34)+AMP(41)+AMP(42)-AMP(49)-AMP(50)

      RES = 0.D0
      DO M = 1, NAMPSO
        DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
            ZTEMP = ZTEMP + CF(J,I)*JAMP(J,M)
          ENDDO
          DO N = 1, NAMPSO
            SQSOIND = SQSOINDEX_4(M,N)
            RES(SQSOIND) = RES(SQSOIND) + ZTEMP*DCONJG(JAMP(I,N))
     $       /DENOM(I)
          ENDDO
        ENDDO
      ENDDO

      END

      SUBROUTINE GET_ME_4(P, ALPHAS, NHEL ,ANS)
      IMPLICIT NONE
C     
C     CONSTANT
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
C     
C     ARGUMENTS 
C     
      REAL*8 P(0:3,NEXTERNAL),ANS
      INTEGER NHEL
      DOUBLE PRECISION ALPHAS
      REAL*8 PI
CF2PY INTENT(OUT) :: ANS
CF2PY INTENT(IN) :: NHEL
CF2PY INTENT(IN) :: P(0:3,NEXTERNAL)
CF2PY INTENT(IN) :: ALPHAS
C     ROUTINE FOR F2PY to read the benchmark point.    
C     the include file with the values of the parameters and masses 
      INCLUDE 'coupl.inc'

      PI = 3.141592653589793D0
      G = 2* DSQRT(ALPHAS*PI)
      CALL UPDATE_AS_PARAM()
      IF (NHEL.NE.0) THEN
        CALL SMATRIXHEL_4(P, NHEL, ANS)
      ELSE
        CALL SMATRIX_4(P, ANS)
      ENDIF
      RETURN
      END

      SUBROUTINE INITIALISE_4(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.    
      IMPLICIT NONE
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH
C     USE SETPARA2 and not SETPARA so that ident_card.dat can be
C      automatically found
      CALL SETPARA2(PATH)  !first call to setup the paramaters    
      RETURN
      END

      LOGICAL FUNCTION IS_BORN_HEL_SELECTED_4(HELID)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=6)
      INTEGER    NCOMB
      PARAMETER (NCOMB=64)
C     
C     ARGUMENTS
C     
      INTEGER HELID
C     
C     LOCALS
C     
      INTEGER I,J
      LOGICAL FOUNDIT
C     
C     GLOBALS
C     
      INTEGER HELC(NEXTERNAL,NCOMB)
      COMMON/BORN_HEL_CONFIGS_4/HELC

      INTEGER POLARIZATIONS(0:NEXTERNAL,0:5)
      COMMON/BORN_BEAM_POL_4/POLARIZATIONS
C     ----------
C     BEGIN CODE
C     ----------

      IS_BORN_HEL_SELECTED_4 = .TRUE.
      IF (POLARIZATIONS(0,0).EQ.-1) THEN
        RETURN
      ENDIF

      DO I=1,NEXTERNAL
        IF (POLARIZATIONS(I,0).EQ.-1) THEN
          CYCLE
        ENDIF
        FOUNDIT = .FALSE.
        DO J=1,POLARIZATIONS(I,0)
          IF (HELC(I,HELID).EQ.POLARIZATIONS(I,J)) THEN
            FOUNDIT = .TRUE.
            EXIT
          ENDIF
        ENDDO
        IF(.NOT.FOUNDIT) THEN
          IS_BORN_HEL_SELECTED_4 = .FALSE.
          RETURN
        ENDIF
      ENDDO

      RETURN
      END





C     Set of functions to handle the array indices of the split orders


      INTEGER FUNCTION SQSOINDEX_4(ORDERINDEXA, ORDERINDEXB)
C     
C     This functions plays the role of the interference matrix. It can
C      be hardcoded or 
C     made more elegant using hashtables if its execution speed ever
C      becomes a relevant
C     factor. From two split order indices, it return the
C      corresponding index in the squared 
C     order canonical ordering.
C     
C     CONSTANTS
C     

      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=2, NSQUAREDSO=7, NAMPSO=4)
C     
C     ARGUMENTS
C     
      INTEGER ORDERINDEXA, ORDERINDEXB
C     
C     LOCAL VARIABLES
C     
      INTEGER I, SQORDERS(NSO)
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      DATA (AMPSPLITORDERS(  1,I),I=  1,  2) /    3,    1/
      DATA (AMPSPLITORDERS(  2,I),I=  1,  2) /    2,    2/
      DATA (AMPSPLITORDERS(  3,I),I=  1,  2) /    1,    3/
      DATA (AMPSPLITORDERS(  4,I),I=  1,  2) /    0,    4/
      COMMON/AMPSPLITORDERS_4/AMPSPLITORDERS
C     
C     FUNCTION
C     
      INTEGER SOINDEX_FOR_SQUARED_ORDERS_4
C     
C     BEGIN CODE
C     
      DO I=1,NSO
        SQORDERS(I)=AMPSPLITORDERS(ORDERINDEXA,I)+AMPSPLITORDERS
     $   (ORDERINDEXB,I)
      ENDDO
      SQSOINDEX_4=SOINDEX_FOR_SQUARED_ORDERS_4(SQORDERS)
      END

      INTEGER FUNCTION SOINDEX_FOR_SQUARED_ORDERS_4(ORDERS)
C     
C     This functions returns the integer index identifying the squared
C      split orders list passed in argument which corresponds to the
C      values of the following list of couplings (and in this order).
C     ['QCD', 'QED']
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO, NAMPSO
      PARAMETER (NSO=2, NSQSO=7, NAMPSO=4)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER SQSPLITORDERS(NSQSO,NSO)
      DATA (SQSPLITORDERS(  1,I),I=  1,  2) /    6,    2/
      DATA (SQSPLITORDERS(  2,I),I=  1,  2) /    5,    3/
      DATA (SQSPLITORDERS(  3,I),I=  1,  2) /    4,    4/
      DATA (SQSPLITORDERS(  4,I),I=  1,  2) /    3,    5/
      DATA (SQSPLITORDERS(  5,I),I=  1,  2) /    2,    6/
      DATA (SQSPLITORDERS(  6,I),I=  1,  2) /    1,    7/
      DATA (SQSPLITORDERS(  7,I),I=  1,  2) /    0,    8/
      COMMON/ALL_SQSPLITORDERS_4/SQSPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NSQSO
        DO J=1,NSO
          IF (ORDERS(J).NE.SQSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_SQUARED_ORDERS_4 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping in function'
      WRITE(*,*) 'SOINDEX_FOR_SQUARED_ORDERS_4'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

      SUBROUTINE GET_NSQSO_BORN_4(NSQSO)
C     
C     Simple subroutine returning the number of squared split order
C     contributions returned when calling smatrix_split_orders 
C     

      INTEGER    NSQUAREDSO
      PARAMETER  (NSQUAREDSO=7)

      INTEGER NSQSO
CF2PY INTENT(OUT) :: NSQSO

      NSQSO=NSQUAREDSO

      END

C     This is the inverse subroutine of SOINDEX_FOR_SQUARED_ORDERS_4.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_SQUARED_ORDERS_FOR_SOINDEX_4(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the squared
C      split order index in argument. Order values correspond to
C      following list of couplings (and in this order):
C     ['QCD', 'QED']
C     
C     CONSTANTS
C     
      INTEGER    NSO, NSQSO
      PARAMETER (NSO=2, NSQSO=7)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
CF2PY INTENT(IN) :: SOINDEX
CF2PY INTENT(OUT) :: ORDERS
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER SQSPLITORDERS(NSQSO,NSO)
      COMMON/ALL_SQSPLITORDERS_4/SQSPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.GT.0.AND.SOINDEX.LE.NSQSO) THEN
        DO I=1,NSO
          ORDERS(I) =  SQSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function GET_SQUARED_ORDERS_FOR_SOIN'
     $ //'DEX'
      WRITE(*,*) 'Could not find squared orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This is the inverse subroutine of getting amplitude SO orders.
C      Not directly useful, but provided nonetheless.
      SUBROUTINE GET_ORDERS_FOR_AMPSOINDEX_4(SOINDEX,ORDERS)
C     
C     This functions returns the orders identified by the split order
C      index in argument. Order values correspond to following list of
C      couplings (and in this order):
C     ['QCD', 'QED']
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=2, NAMPSO=4)
C     
C     ARGUMENTS
C     
      INTEGER SOINDEX, ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS_4/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      IF (SOINDEX.GT.0.AND.SOINDEX.LE.NAMPSO) THEN
        DO I=1,NSO
          ORDERS(I) =  AMPSPLITORDERS(SOINDEX,I)
        ENDDO
        RETURN
      ENDIF

      WRITE(*,*) 'ERROR:: Stopping function GET_ORDERS_FOR_AMPSOINDEX'
      WRITE(*,*) 'Could not find amplitude split orders index ',SOINDEX
      STOP

      END SUBROUTINE

C     This function is not directly useful, but included for
C      completeness
      INTEGER FUNCTION SOINDEX_FOR_AMPORDERS_4(ORDERS)
C     
C     This functions returns the integer index identifying the
C      amplitude split orders passed in argument which correspond to
C      the values of the following list of couplings (and in this
C      order):
C     ['QCD', 'QED']
C     
C     CONSTANTS
C     
      INTEGER    NSO, NAMPSO
      PARAMETER (NSO=2, NAMPSO=4)
C     
C     ARGUMENTS
C     
      INTEGER ORDERS(NSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J
      INTEGER AMPSPLITORDERS(NAMPSO,NSO)
      COMMON/AMPSPLITORDERS_4/AMPSPLITORDERS
C     
C     BEGIN CODE
C     
      DO I=1,NAMPSO
        DO J=1,NSO
          IF (ORDERS(J).NE.AMPSPLITORDERS(I,J)) GOTO 1009
        ENDDO
        SOINDEX_FOR_AMPORDERS_4 = I
        RETURN
 1009   CONTINUE
      ENDDO

      WRITE(*,*) 'ERROR:: Stopping function SOINDEX_FOR_AMPORDERS_4'
      WRITE(*,*) 'Could not find squared orders ',(ORDERS(I),I=1,NSO)
      STOP

      END

C     --------
C     Now defining F2PY hooks.
C     --------

      SUBROUTINE GET_SPLIT_ORDER_NAMES_4(SONAMES)
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NSO
      PARAMETER (NSO=2)
C     ARGUMENTS
C     
      CHARACTER*100 SONAMES(NSO)
      INTEGER IDX
CF2PY INTENT(OUT) :: SONAMES
C     
C     BEGIN CODE
C     
      SONAMES(1)='QCD'
      SONAMES(2)='QED'

      END


      SUBROUTINE CONFIGS_BORN
      IMPLICIT NONE
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'maxamps.inc'
      INCLUDE 'genps.inc'
      INTEGER IFOREST(2,-MAX_BRANCH:-1,LMAXCONFIGS)
      COMMON/TO_FOREST/ IFOREST
      INTEGER MAPCONFIG(0:LMAXCONFIGS), THIS_CONFIG
      COMMON/TO_MCONFIGS/MAPCONFIG, THIS_CONFIG
      INTEGER SPROP(MAXSPROC,-MAX_BRANCH:-1,LMAXCONFIGS)
      INTEGER TPRID(-MAX_BRANCH:-1,LMAXCONFIGS)
      COMMON/TO_SPROP/SPROP,TPRID
C     Diagram 2
      MAPCONFIG(1)=2
      IFOREST(1:2,-1,1)=[1,3]
      TPRID(-1,1)=24
      SPROP(1:1,-1,1)=0
      IFOREST(1:2,-2,1)=[-1,4]
C     Diagram 3
      MAPCONFIG(2)=3
      IFOREST(1:2,-1,2)=[1,3]
      TPRID(-1,2)=251
      SPROP(1:1,-1,2)=0
      IFOREST(1:2,-2,2)=[-1,4]
C     Diagram 4
      MAPCONFIG(3)=4
      IFOREST(1:2,-1,3)=[1,4]
      TPRID(-1,3)=24
      SPROP(1:1,-1,3)=0
      IFOREST(1:2,-2,3)=[-1,3]
C     Diagram 5
      MAPCONFIG(4)=5
      IFOREST(1:2,-1,4)=[1,4]
      TPRID(-1,4)=251
      SPROP(1:1,-1,4)=0
      IFOREST(1:2,-2,4)=[-1,3]
C     Number of configs
      MAPCONFIG(0)=4
      RETURN

      END


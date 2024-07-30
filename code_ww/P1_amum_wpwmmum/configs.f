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
C     Diagram 1
      MAPCONFIG(1)=1
      IFOREST(1:2,-1,1)=[4,3]
      SPROP(1:1,-1,1)=22
      TPRID(-1,1)=0
      IFOREST(1:2,-2,1)=[5,-1]
      SPROP(1:1,-2,1)=13
      TPRID(-2,1)=0
C     Diagram 2
      MAPCONFIG(2)=2
      IFOREST(1:2,-1,2)=[4,3]
      SPROP(1:1,-1,2)=23
      TPRID(-1,2)=0
      IFOREST(1:2,-2,2)=[5,-1]
      SPROP(1:1,-2,2)=13
      TPRID(-2,2)=0
C     Diagram 3
      MAPCONFIG(3)=3
      IFOREST(1:2,-1,3)=[5,3]
      SPROP(1:1,-1,3)=14
      TPRID(-1,3)=0
      IFOREST(1:2,-2,3)=[4,-1]
      SPROP(1:1,-2,3)=13
      TPRID(-2,3)=0
C     Diagram 4
      MAPCONFIG(4)=4
      IFOREST(1:2,-1,4)=[1,3]
      TPRID(-1,4)=24
      SPROP(1:1,-1,4)=0
      IFOREST(1:2,-2,4)=[-1,5]
      TPRID(-2,4)=14
      SPROP(1:1,-2,4)=0
      IFOREST(1:2,-3,4)=[-2,4]
C     Diagram 5
      MAPCONFIG(5)=5
      IFOREST(1:2,-1,5)=[1,3]
      TPRID(-1,5)=24
      SPROP(1:1,-1,5)=0
      IFOREST(1:2,-2,5)=[-1,4]
      TPRID(-2,5)=22
      SPROP(1:1,-2,5)=0
      IFOREST(1:2,-3,5)=[-2,5]
C     Diagram 6
      MAPCONFIG(6)=6
      IFOREST(1:2,-1,6)=[1,3]
      TPRID(-1,6)=24
      SPROP(1:1,-1,6)=0
      IFOREST(1:2,-2,6)=[-1,4]
      TPRID(-2,6)=23
      SPROP(1:1,-2,6)=0
      IFOREST(1:2,-3,6)=[-2,5]
C     Diagram 7
      MAPCONFIG(7)=7
      IFOREST(1:2,-1,7)=[1,3]
      TPRID(-1,7)=251
      SPROP(1:1,-1,7)=0
      IFOREST(1:2,-2,7)=[-1,4]
      TPRID(-2,7)=22
      SPROP(1:1,-2,7)=0
      IFOREST(1:2,-3,7)=[-2,5]
C     Diagram 8
      MAPCONFIG(8)=8
      IFOREST(1:2,-1,8)=[1,3]
      TPRID(-1,8)=251
      SPROP(1:1,-1,8)=0
      IFOREST(1:2,-2,8)=[-1,4]
      TPRID(-2,8)=23
      SPROP(1:1,-2,8)=0
      IFOREST(1:2,-3,8)=[-2,5]
C     Diagram 9
      MAPCONFIG(9)=9
      IFOREST(1:2,-1,9)=[1,4]
      TPRID(-1,9)=24
      SPROP(1:1,-1,9)=0
      IFOREST(1:2,-2,9)=[-1,3]
      TPRID(-2,9)=22
      SPROP(1:1,-2,9)=0
      IFOREST(1:2,-3,9)=[-2,5]
C     Diagram 10
      MAPCONFIG(10)=10
      IFOREST(1:2,-1,10)=[1,4]
      TPRID(-1,10)=24
      SPROP(1:1,-1,10)=0
      IFOREST(1:2,-2,10)=[-1,3]
      TPRID(-2,10)=23
      SPROP(1:1,-2,10)=0
      IFOREST(1:2,-3,10)=[-2,5]
C     Diagram 11
      MAPCONFIG(11)=11
      IFOREST(1:2,-1,11)=[1,4]
      TPRID(-1,11)=251
      SPROP(1:1,-1,11)=0
      IFOREST(1:2,-2,11)=[-1,3]
      TPRID(-2,11)=22
      SPROP(1:1,-2,11)=0
      IFOREST(1:2,-3,11)=[-2,5]
C     Diagram 12
      MAPCONFIG(12)=12
      IFOREST(1:2,-1,12)=[1,4]
      TPRID(-1,12)=251
      SPROP(1:1,-1,12)=0
      IFOREST(1:2,-2,12)=[-1,3]
      TPRID(-2,12)=23
      SPROP(1:1,-2,12)=0
      IFOREST(1:2,-3,12)=[-2,5]
C     Diagram 13
      MAPCONFIG(13)=13
      IFOREST(1:2,-1,13)=[5,3]
      SPROP(1:1,-1,13)=14
      TPRID(-1,13)=0
      IFOREST(1:2,-2,13)=[1,4]
      TPRID(-2,13)=24
      SPROP(1:1,-2,13)=0
      IFOREST(1:2,-3,13)=[-2,-1]
C     Diagram 14
      MAPCONFIG(14)=14
      IFOREST(1:2,-1,14)=[1,5]
      TPRID(-1,14)=13
      SPROP(1:1,-1,14)=0
      IFOREST(1:2,-2,14)=[-1,3]
      TPRID(-2,14)=14
      SPROP(1:1,-2,14)=0
      IFOREST(1:2,-3,14)=[-2,4]
C     Diagram 15
      MAPCONFIG(15)=15
      IFOREST(1:2,-1,15)=[4,3]
      SPROP(1:1,-1,15)=22
      TPRID(-1,15)=0
      IFOREST(1:2,-2,15)=[1,5]
      TPRID(-2,15)=13
      SPROP(1:1,-2,15)=0
      IFOREST(1:2,-3,15)=[-2,-1]
C     Diagram 16
      MAPCONFIG(16)=16
      IFOREST(1:2,-1,16)=[4,3]
      SPROP(1:1,-1,16)=23
      TPRID(-1,16)=0
      IFOREST(1:2,-2,16)=[1,5]
      TPRID(-2,16)=13
      SPROP(1:1,-2,16)=0
      IFOREST(1:2,-3,16)=[-2,-1]
C     Number of configs
      MAPCONFIG(0)=16
      RETURN

      END


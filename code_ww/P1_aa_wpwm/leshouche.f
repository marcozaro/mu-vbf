      SUBROUTINE GETLESHOUCHE_BORN
      IMPLICIT NONE
      INCLUDE 'nexternal.inc'
      INCLUDE 'maxamps.inc'
      INTEGER IDUP(NEXTERNAL,MAXPROC,MAXSPROC)
      INTEGER MOTHUP(2,NEXTERNAL)
      INTEGER ICOLUP(2,NEXTERNAL,MAXFLOW,MAXSPROC)
      COMMON/LESHOUCHE/IDUP,MOTHUP,ICOLUP
      IDUP(1:4,1,1)=[22,22,24,-24]
      MOTHUP(1,1: 4)=[  0,  0,  1,  1]
      MOTHUP(2,1: 4)=[  0,  0,  2,  2]
      RETURN
      END
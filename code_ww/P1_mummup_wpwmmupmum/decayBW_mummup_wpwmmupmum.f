      SUBROUTINE DECAYBW_MUMMUP_WPWMMUPMUM
      IMPLICIT NONE
      INCLUDE 'maxconfigs.inc'
      INCLUDE 'genps.inc'
      INTEGER GFORCEBW(-MAX_BRANCH:-1,LMAXCONFIGS)
      COMMON/INC_BW/GFORCEBW
      GFORCEBW(-1,5)=0
      GFORCEBW(-1,6)=0
      GFORCEBW(-1,7)=0
      GFORCEBW(-2,7)=0
      GFORCEBW(-1,8)=0
      GFORCEBW(-1,9)=0
      GFORCEBW(-2,9)=0
      GFORCEBW(-1,10)=0
      GFORCEBW(-1,12)=0
      GFORCEBW(-2,12)=0
      GFORCEBW(-1,13)=0
      GFORCEBW(-2,13)=0
      GFORCEBW(-3,13)=0
      GFORCEBW(-1,14)=0
      GFORCEBW(-2,14)=0
      GFORCEBW(-3,14)=0
      GFORCEBW(-1,15)=0
      GFORCEBW(-2,15)=0
      GFORCEBW(-3,15)=0
      GFORCEBW(-1,16)=0
      GFORCEBW(-2,16)=0
      GFORCEBW(-3,16)=0
      GFORCEBW(-1,17)=0
      GFORCEBW(-2,17)=0
      GFORCEBW(-3,17)=0
      GFORCEBW(-1,18)=0
      GFORCEBW(-2,18)=0
      GFORCEBW(-3,18)=0
      GFORCEBW(-1,19)=0
      GFORCEBW(-2,19)=0
      GFORCEBW(-3,19)=0
      GFORCEBW(-1,20)=0
      GFORCEBW(-2,20)=0
      GFORCEBW(-3,20)=0
      GFORCEBW(-1,21)=0
      GFORCEBW(-2,21)=0
      GFORCEBW(-3,21)=0
      GFORCEBW(-1,22)=0
      GFORCEBW(-2,22)=0
      GFORCEBW(-3,22)=0
      GFORCEBW(-1,23)=0
      GFORCEBW(-2,23)=0
      GFORCEBW(-3,23)=0
      GFORCEBW(-1,24)=0
      GFORCEBW(-2,24)=0
      GFORCEBW(-3,24)=0
      GFORCEBW(-1,25)=0
      GFORCEBW(-2,25)=0
      GFORCEBW(-3,25)=0
      GFORCEBW(-1,26)=0
      GFORCEBW(-2,26)=0
      GFORCEBW(-3,26)=0
      GFORCEBW(-1,27)=0
      GFORCEBW(-2,27)=0
      GFORCEBW(-3,27)=0
      GFORCEBW(-1,28)=0
      GFORCEBW(-2,28)=0
      GFORCEBW(-3,28)=0
      GFORCEBW(-1,29)=0
      GFORCEBW(-2,29)=0
      GFORCEBW(-3,29)=0
      GFORCEBW(-1,30)=0
      GFORCEBW(-2,30)=0
      GFORCEBW(-3,30)=0
      GFORCEBW(-1,31)=0
      GFORCEBW(-2,31)=0
      GFORCEBW(-3,31)=0
      GFORCEBW(-1,32)=0
      GFORCEBW(-2,32)=0
      GFORCEBW(-3,32)=0
      GFORCEBW(-1,33)=0
      GFORCEBW(-2,33)=0
      GFORCEBW(-3,33)=0
      GFORCEBW(-1,34)=0
      GFORCEBW(-2,34)=0
      GFORCEBW(-3,34)=0
      GFORCEBW(-1,35)=0
      GFORCEBW(-2,35)=0
      GFORCEBW(-3,35)=0
      GFORCEBW(-1,36)=0
      GFORCEBW(-2,36)=0
      GFORCEBW(-3,36)=0
      GFORCEBW(-1,37)=0
      GFORCEBW(-2,37)=0
      GFORCEBW(-3,37)=0
      GFORCEBW(-1,38)=0
      GFORCEBW(-2,38)=0
      GFORCEBW(-3,38)=0
      GFORCEBW(-1,39)=0
      GFORCEBW(-2,39)=0
      GFORCEBW(-3,39)=0
      GFORCEBW(-1,40)=0
      GFORCEBW(-2,40)=0
      GFORCEBW(-3,40)=0
      GFORCEBW(-1,41)=0
      GFORCEBW(-2,41)=0
      GFORCEBW(-3,41)=0
      GFORCEBW(-1,42)=0
      GFORCEBW(-2,42)=0
      GFORCEBW(-3,42)=0
      GFORCEBW(-1,43)=0
      GFORCEBW(-2,43)=0
      GFORCEBW(-3,43)=0
      GFORCEBW(-1,44)=0
      GFORCEBW(-2,44)=0
      GFORCEBW(-3,44)=0
      GFORCEBW(-1,45)=0
      GFORCEBW(-2,45)=0
      GFORCEBW(-3,45)=0
      GFORCEBW(-1,46)=0
      GFORCEBW(-2,46)=0
      GFORCEBW(-3,46)=0
      GFORCEBW(-1,47)=0
      GFORCEBW(-1,48)=0
      GFORCEBW(-2,48)=0
      GFORCEBW(-1,49)=0
      GFORCEBW(-1,50)=0
      GFORCEBW(-2,50)=0
      GFORCEBW(-1,51)=0
      GFORCEBW(-1,68)=0
      GFORCEBW(-1,69)=0
      GFORCEBW(-1,70)=0
      GFORCEBW(-2,70)=0
      GFORCEBW(-1,71)=0
      GFORCEBW(-1,72)=0
      GFORCEBW(-2,72)=0
      GFORCEBW(-1,73)=0
      GFORCEBW(-1,74)=0
      GFORCEBW(-2,74)=0
      GFORCEBW(-1,75)=0
      GFORCEBW(-1,76)=0
      GFORCEBW(-2,76)=0
      GFORCEBW(-1,81)=0
      GFORCEBW(-2,81)=0
      GFORCEBW(-1,82)=0
      GFORCEBW(-2,82)=0
      GFORCEBW(-1,83)=0
      GFORCEBW(-2,83)=0
      GFORCEBW(-1,84)=0
      GFORCEBW(-1,85)=0
      GFORCEBW(-2,85)=0
      GFORCEBW(-1,86)=0
      GFORCEBW(-1,87)=0
      GFORCEBW(-2,87)=0
      GFORCEBW(-1,88)=0
      GFORCEBW(-1,89)=0
      GFORCEBW(-2,89)=0
      GFORCEBW(-1,90)=0
      GFORCEBW(-1,91)=0
      GFORCEBW(-2,91)=0
      GFORCEBW(-1,92)=0
      GFORCEBW(-1,93)=0
      GFORCEBW(-2,93)=0
      GFORCEBW(-1,94)=0
      GFORCEBW(-1,95)=0
      GFORCEBW(-2,95)=0
      GFORCEBW(-1,96)=0
      GFORCEBW(-2,96)=0
      GFORCEBW(-1,97)=0
      GFORCEBW(-2,97)=0
      GFORCEBW(-1,98)=0
      GFORCEBW(-2,98)=0
      GFORCEBW(-1,99)=0
      GFORCEBW(-2,99)=0
      GFORCEBW(-1,100)=0
      GFORCEBW(-2,100)=0
      GFORCEBW(-1,101)=0
      GFORCEBW(-2,101)=0
      GFORCEBW(-1,102)=0
      GFORCEBW(-2,102)=0
      GFORCEBW(-1,103)=0
      GFORCEBW(-1,104)=0
      GFORCEBW(-2,104)=0
      GFORCEBW(-1,105)=0
      GFORCEBW(-1,106)=0
      GFORCEBW(-2,106)=0
      RETURN
      END


      ! matrix elements:
      ! 4: P1_mummup_ttxmupmum
      ! 3: P1_amum_ttxmum 
      ! 2: P1_amup_ttxmup 
      ! 1: P1_aa_ttx

      subroutine compute_me_doublereal(p,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
C returns the matrix element for the double real emission times
C (1-y1)*(1-y2), possibly approximated in the collinear limit(s)
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)
      ! 3-> single resolved collinear emission (y2=1)
      ! 4-> no resolved collinear emission (y1=y2=1)
      !!! note that xi are different in the various kinematics
      double precision y1,y2,omy1,omy2,xi1(4),xi2(4),ph1,ph2 
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      parameter (tiny=1d-4)
      double precision alp8pi

      double precision p_pass(0:3,6)
      double precision kp1(0:3), kp2(0:3), ksq1, ksq2
      double precision z1, z2

      double precision dot
      external dot 
      double precision pi
      parameter (pi=3.14159265359d0)

      double precision shat
      common /to_shat/shat

      double precision real2, real1c

      ! need to rearrange the momenta, as we assume mu+ mu-, while
      ! ME's have been generated accodring to what's written at the top
      ! of this file
      p_pass(:,:) = 0d0
      scvec(:,:) = 0d0

      alp8pi = dble(gal(1))**2*2

      ! k is the momentum entering the reduced matrix element
      ! kp is the direction in the orthogonal plane
      kp1(:) = (/0d0,cos(ph1),sin(ph1),0d0/)
      kp2(:) = (/0d0,cos(ph2),sin(ph2),0d0/)

      if (1d0-y1.gt.tiny.and.1d0-y2.gt.tiny) then
         !call write_momenta(p(0,1,1),6)
          p_pass(:,3:6) = p(:,3:6,1)
          p_pass(:,1) = p(:,2,1)
          p_pass(:,2) = p(:,1,1)
          !call write_momenta(p_pass,6)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_4(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy1*omy2
          !
      else if (1d0-y1.lt.tiny.and.1d0-y2.gt.tiny) then !collinear on leg 1 (mu+)
          p_pass(:,3:4) = p(:,3:4,2)
          p_pass(:,5) = p(:,6,2)
          p_pass(:,2) = p(:,2,2)
          p_pass(:,1) = p(:,1,2) * (1d0-xi1(2))
          p_pass(:,6) = 0d0
          !call write_momenta(p_pass,6)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_3(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)
          scvec(1,1) = kp1(0)
          scvec(1,2) = kp1(1)
          scvec(1,3) = kp1(2)
          scvec(1,4) = kp1(3)
          call set_spin_correlation_vectors_3(1,3,scvec)
          call SMATRIX_SPLITORDERS_3(p_pass,ANS_splitorders)
          ansk1 = ans_splitorders(0)
          call reset_spin_correlation_vectors_3()
          ksq1 = -xi1(2)*shat/2d0
          z1 = 1d0 - xi1(2)
          ans = alp8pi/-ksq1*(z1*ans+ansk1*4d0*(1d0-z1)/z1)/z1*omy2
          !
      else if (1d0-y1.gt.tiny.and.1d0-y2.lt.tiny) then !collinear on leg 2 (mu-)
          !call write_momenta(p(0,1,3),6)
          p_pass(:,3:4) = p(:,3:4,3)
          p_pass(:,5) = p(:,5,3)
          p_pass(:,2) = p(:,1,3)
          p_pass(:,1) = p(:,2,3) * (1d0-xi2(3))
          !call write_momenta(p_pass,6)
          call check_momenta(p_pass,5)
          call ME_ACCESSOR_HOOK_2(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)
          scvec(1,1) = kp2(0)
          scvec(1,2) = kp2(1)
          scvec(1,3) = kp2(2)
          scvec(1,4) = kp2(3)
          call set_spin_correlation_vectors_2(1,3,scvec)
          call SMATRIX_SPLITORDERS_2(p_pass,ANS_splitorders)
          ansk2 = ans_splitorders(0)
          call reset_spin_correlation_vectors_2()
          ksq2 = -xi2(3)*shat/2d0
          z2 = 1d0 - xi2(3)
          ans = alp8pi/-ksq2*(z2*ans+ansk2*4d0*(1d0-z2)/z2)/z2*omy1
          !
      else if (1d0-y1.lt.tiny.and.1d0-y2.lt.tiny) then !collinear on legs 1/2
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4) * (1d0-xi2(4))
          p_pass(:,1) = p(:,1,4) * (1d0-xi1(4))
          !call write_momenta(p(0,1,4),6)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_1(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)
          scvec(1,1) = kp1(0)
          scvec(1,2) = kp1(1)
          scvec(1,3) = kp1(2)
          scvec(1,4) = kp1(3)
          call set_spin_correlation_vectors_1(1,3,scvec)
          call SMATRIX_SPLITORDERS_1(p_pass,ANS_splitorders)
          ansk1 = ans_splitorders(0)
          call reset_spin_correlation_vectors_1()
          scvec(1,1) = kp2(0)
          scvec(1,2) = kp2(1)
          scvec(1,3) = kp2(2)
          scvec(1,4) = kp2(3)
          call set_spin_correlation_vectors_1(2,3,scvec)
          call SMATRIX_SPLITORDERS_1(p_pass,ANS_splitorders)
          ansk2 = ans_splitorders(0)
          call reset_spin_correlation_vectors_1()
          scvec(1,1) = kp1(0)
          scvec(1,2) = kp1(1)
          scvec(1,3) = kp1(2)
          scvec(1,4) = kp1(3)
          call set_spin_correlation_vectors_1(1,3,scvec)
          scvec(1,1) = kp2(0)
          scvec(1,2) = kp2(1)
          scvec(1,3) = kp2(2)
          scvec(1,4) = kp2(3)
          call set_spin_correlation_vectors_1(2,3,scvec)
          call SMATRIX_SPLITORDERS_1(p_pass,ANS_splitorders)
          ansk12 = ans_splitorders(0)
          call reset_spin_correlation_vectors_1()
          ksq1 = -xi1(4)*shat/2d0
          z1 = 1d0 - xi1(4)
          ksq2 = -xi2(4)*shat/2d0
          z2 = 1d0 - xi2(4)
          ans = alp8pi**2 / -ksq1 / -ksq2 * 
     #        (z1*z2*ans + z2*ansk1*4d0*(1d0-z1)/z1 + 
     #          z1*ansk2*4d0*(1d0-z2)/z2 +     
     #         ansk12*4d0*(1d0-z1)/z1*4d0*(1d0-z2)/z2) / z1 / z2

      endif

      return
      end


      subroutine compute_me_singlereal1a(p,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
C returns the matrix element for the mu+gamma single real emission times
C (1-y1), possibly approximated in the collinear limit(s)
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)
      ! 3-> single resolved collinear emission (y2=1)<-
      ! 4-> no resolved collinear emission (y1=y2=1)<-
      !!! note that xi are different in the various kinematics
      double precision y1,y2,omy1,omy2,xi1(4),xi2(4),ph1,ph2 
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      parameter (tiny=1d-4)
      double precision alp8pi

      double precision p_pass(0:3,6)
      double precision kp1(0:3), kp2(0:3), ksq1, ksq2
      double precision z1, z2

      double precision dot
      external dot 
      double precision pi
      parameter (pi=3.14159265359d0)

      double precision shat
      common /to_shat/shat

      double precision real2, real1c

      ! need to rearrange the momenta, as we assume mu+ mu-, while
      ! ME's have been generated accodring to what's written at the top
      ! of this file
      p_pass(:,:) = 0d0
      scvec(:,:) = 0d0

      alp8pi = dble(gal(1))**2*2

      ! k is the momentum entering the reduced matrix element
      ! kp is the direction in the orthogonal plane
      kp1(:) = (/0d0,cos(ph1),sin(ph1),0d0/)
      !kp2(:) = (/0d0,cos(ph2),sin(ph2),0d0/)

      if (1d0-y1.gt.tiny) then
          p_pass(:,3:4) = p(:,3:4,3)
          p_pass(:,5) = p(:,5,3)
          p_pass(:,2) = p(:,1,3)
          p_pass(:,1) = p(:,2,3)
          !call write_momenta(p_pass,6)
          call check_momenta(p_pass,5)
          call ME_ACCESSOR_HOOK_2(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy1
          !
      else if (1d0-y1.lt.tiny) then !collinear on leg 1 (mu+)
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4)
          p_pass(:,1) = p(:,1,4) * (1d0-xi1(4))
          !call write_momenta(p(0,1,4),6)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_1(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)
          scvec(1,1) = kp1(0)
          scvec(1,2) = kp1(1)
          scvec(1,3) = kp1(2)
          scvec(1,4) = kp1(3)
          call set_spin_correlation_vectors_1(1,3,scvec)
          call SMATRIX_SPLITORDERS_1(p_pass,ANS_splitorders)
          ansk1 = ans_splitorders(0)
          call reset_spin_correlation_vectors_1()
          ksq1 = -xi1(3)*shat/2d0
          z1 = 1d0 - xi1(3)
          ans = alp8pi/-ksq1*(z1*ans+ansk1*4d0*(1d0-z1)/z1)/z1
          !
      endif

      return
      end



      subroutine compute_me_singlereal1b(p,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
C returns the matrix element for the gamma-mu- single real emission times
C (1-y2), possibly approximated in the collinear limit(s)
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)
      ! 3-> single resolved collinear emission (y2=1)<-
      ! 4-> no resolved collinear emission (y1=y2=1)<-
      !!! note that xi are different in the various kinematics
      double precision y1,y2,omy1,omy2,xi1(4),xi2(4),ph1,ph2 
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      parameter (tiny=1d-4)
      double precision alp8pi

      double precision p_pass(0:3,6)
      double precision kp1(0:3), kp2(0:3), ksq1, ksq2
      double precision z1, z2

      double precision dot
      external dot 
      double precision pi
      parameter (pi=3.14159265359d0)

      double precision shat
      common /to_shat/shat

      double precision real2, real1c

      ! need to rearrange the momenta, as we assume mu+ mu-, while
      ! ME's have been generated accodring to what's written at the top
      ! of this file
      p_pass(:,:) = 0d0
      scvec(:,:) = 0d0

      alp8pi = dble(gal(1))**2*2

      ! k is the momentum entering the reduced matrix element
      ! kp is the direction in the orthogonal plane
      !kp1(:) = (/0d0,cos(ph1),sin(ph1),0d0/)
      kp2(:) = (/0d0,cos(ph2),sin(ph2),0d0/)

      if (1d0-y2.gt.tiny) then
          p_pass(:,1:5) = p(:,1:5,2)
          !call write_momenta(p_pass,6)
          call check_momenta(p_pass,5)
          call ME_ACCESSOR_HOOK_3(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy2
          !
      else if (1d0-y1.lt.tiny) then !collinear on leg 2 (mu-)
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4) * (1d0-xi2(4))
          p_pass(:,1) = p(:,1,4) 
          !call write_momenta(p(0,1,4),6)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_1(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)
          scvec(1,1) = kp2(0)
          scvec(1,2) = kp2(1)
          scvec(1,3) = kp2(2)
          scvec(1,4) = kp2(3)
          call set_spin_correlation_vectors_1(2,3,scvec)
          call SMATRIX_SPLITORDERS_1(p_pass,ANS_splitorders)
          ansk2 = ans_splitorders(0)
          call reset_spin_correlation_vectors_1()
          ksq2 = -xi2(2)*shat/2d0
          z2 = 1d0 - xi2(2)
          ans = alp8pi/-ksq2*(z2*ans+ansk2*4d0*(1d0-z2)/z2)/z2
          !
      endif

      return
      end





      subroutine compute_me_born_gaga(p,ans)
C returns the matrix element for the gamma-gamma born term
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)
      ! 3-> single resolved collinear emission (y2=1)
      ! 4-> no resolved collinear emission (y1=y2=1)
      ! in this case we use only 1
      double precision ans 
      double precision ans_splitorders(0:99)

      double precision p_pass(0:3,4)
      double precision shat
      common /to_shat/shat

      p_pass(:,:) = 0d0

      p_pass(:,1:4) = p(:,1:4,1)
      call check_momenta(p_pass,4)
      call ME_ACCESSOR_HOOK_1(p_pass,-1,0.118d0,ANS_splitorders)
      ans = ans_splitorders(0)

      return
      end

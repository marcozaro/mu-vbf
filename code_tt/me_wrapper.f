
      ! matrix elements:
      ! 4: P1_mummup_ttxmupmum
      ! 3: P1_amum_ttxmum 
      ! 2: P1_amup_ttxmup 
      ! 1: P1_aa_ttx

      subroutine compute_me_doublereal(p,icoll,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
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
      integer icoll
      double precision y1,y2,omy1,omy2,xi1(4),xi2(4),ph1,ph2 
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      common/to_coll_cutoff/tiny
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

      if (omy1.gt.tiny.and.omy2.gt.tiny) then
          p_pass(:,3:6) = p(:,3:6,1)
          p_pass(:,1) = p(:,2,1)
          p_pass(:,2) = p(:,1,1)
          call check_momenta(p_pass,6)
          call ME_ACCESSOR_HOOK_4(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy1*omy2
          !
      else if (omy1.lt.tiny.and.omy2.gt.tiny) then !collinear on leg 1 (mu+)
          p_pass(:,3:4) = p(:,3:4,2)
          p_pass(:,5) = p(:,6,2)
          p_pass(:,2) = p(:,2,2)
          p_pass(:,1) = p(:,1,2) * (1d0-xi1(2))
          p_pass(:,6) = 0d0
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
          z1 = 1d0 - xi1(icoll)
          ans = alp8pi/-ksq1*(z1*ans+ansk1*4d0*(1d0-z1)/z1)/z1*omy2
          !
      else if (omy1.gt.tiny.and.omy2.lt.tiny) then !collinear on leg 2 (mu-)
          p_pass(:,3:4) = p(:,3:4,3)
          p_pass(:,5) = p(:,5,3)
          p_pass(:,2) = p(:,1,3)
          p_pass(:,1) = p(:,2,3) * (1d0-xi2(3))
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
          z2 = 1d0 - xi2(icoll)
          ans = alp8pi/-ksq2*(z2*ans+ansk2*4d0*(1d0-z2)/z2)/z2*omy1
          !
      else if (omy1.lt.tiny.and.omy2.lt.tiny) then !collinear on legs 1/2
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4) * (1d0-xi2(4))
          p_pass(:,1) = p(:,1,4) * (1d0-xi1(4))
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
          z1 = 1d0 - xi1(icoll)
          ksq2 = -xi2(4)*shat/2d0
          z2 = 1d0 - xi2(icoll)
          ans = alp8pi**2 / -ksq1 / -ksq2 * 
     #        (z1*z2*ans + z2*ansk1*4d0*(1d0-z1)/z1 + 
     #          z1*ansk2*4d0*(1d0-z2)/z2 +     
     #         ansk12*4d0*(1d0-z1)/z1*4d0*(1d0-z2)/z2) / z1 / z2

      endif

      return
      end


      subroutine compute_me_singlereal1a(p,icoll,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
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
      integer icoll
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      common/to_coll_cutoff/tiny
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

      if (omy1.gt.tiny) then
          p_pass(:,3:4) = p(:,3:4,3)
          p_pass(:,5) = p(:,5,3)
          p_pass(:,2) = p(:,1,3)
          p_pass(:,1) = p(:,2,3)
          call check_momenta(p_pass,5)
          call ME_ACCESSOR_HOOK_2(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy1
          !
      else if (omy1.lt.tiny) then !collinear on leg 1 (mu+)
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4)
          p_pass(:,1) = p(:,1,4) * (1d0-xi1(4))
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
          ksq1 = -xi1(icoll)*shat/2d0
          z1 = 1d0 - xi1(icoll)
          ans = alp8pi/-ksq1*(z1*ans+ansk1*4d0*(1d0-z1)/z1)/z1
          !ans = alp8pi/-ksq1*((1+(1-z1)**2)/z1*ans)/z1
          !
      endif

      return
      end



      subroutine compute_me_singlereal1b(p,icoll,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,ans)
C returns the matrix element for the gamma-mu- single real emission times
C (1-y2), possibly approximated in the collinear limit(s)
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)<-
      ! 3-> single resolved collinear emission (y2=1)
      ! 4-> no resolved collinear emission (y1=y2=1)<-
      !!! note that xi are different in the various kinematics
      double precision y1,y2,omy1,omy2,xi1(4),xi2(4),ph1,ph2 
      integer icoll
      double precision ans, ansk1, ansk2, ansk12
      double precision ans_splitorders(0:99)
      integer max_sc_vectors
      parameter(max_sc_vectors=20)
      double precision scvec(max_sc_vectors,4)

      double precision tiny
      common/to_coll_cutoff/tiny
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

      if (omy2.gt.tiny) then
          p_pass(:,1:5) = p(:,1:5,2)
          call check_momenta(p_pass,5)
          call ME_ACCESSOR_HOOK_3(p_pass,-1,0.118d0,ANS_splitorders)
          ans = ans_splitorders(0)*omy2
          !
      else if (omy2.lt.tiny) then !collinear on leg 2 (mu-)
          p_pass(:,3:4) = p(:,3:4,4)
          p_pass(:,2) = p(:,2,4) * (1d0-xi2(4))
          p_pass(:,1) = p(:,1,4) 
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
          ksq2 = -xi2(icoll)*shat/2d0
          z2 = 1d0 - xi2(icoll)
          ans = alp8pi/-ksq2*(z2*ans+ansk2*4d0*(1d0-z2)/z2)/z2
          !
      endif

      return
      end





      subroutine compute_me_born_gaga(p,ans,icoll)
C returns the matrix element for the gamma-gamma born term
      implicit none
      include 'coupl.inc'
      double precision p(0:3,6,4)
      integer icoll
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

      p_pass(:,1:4) = p(:,1:4,icoll)
      call check_momenta(p_pass,4)
      call ME_ACCESSOR_HOOK_1(p_pass,-1,0.118d0,ANS_splitorders)
      ans = ans_splitorders(0)

      return
      end


! now we have functions for the subtracted contributions. 
! they also fill the analysis etc

      double precision function compute_subtracted_me_0(x,vegas_wgt,lum,tau,ycm,jac) 
      ! the photon-photon (born-like) contribution
      implicit none
      double precision x(8),vegas_wgt,lum,tau,ycm,jac

      double precision scoll
      common /to_scoll/scoll
      double precision shat
      common /to_shat/shat
      double precision mfin
      common /to_mfin/mfin
      double precision  mmin
      common /to_mmin/mmin
      double precision thresh

      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll
      logical passcuts
      external passcuts
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)
      ! stuff for the analysis
      integer pdgs(6), istatus(6)
      double precision p_an(0:3,6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos

      pdgs = (/22,22,6,-6,0,0/)
      istatus = (/-1,-1,1,1,1,1/)

      shat = tau * scoll
      thresh = mmin**2/shat

      icoll = 1

      jac0(icoll) = jac
      jac1a(icoll) = 0d0
      jac1b(icoll) = 0d0
      jac2(icoll) = 0d0
      call generate_kinematics(x, shat, thresh, icoll, 0, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
      call generate_momenta(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
      me(icoll) = 0d0
      ! boost the momenta to the lab frame. This is needed
      ! both for cuts and for the analysis
      call boost_to_lab_frame(p0(0,1,icoll),p_an,ycm)
      !! MZ need to understand how to deal with cuts
      if (passcuts(p_an,pdgs,istatus)) then 
        call compute_me_born_gaga(p0, me(icoll), icoll)

        if (fill_histos) then
            wgt_an(1) = jac0(icoll) * me(icoll) 
     &           * vegas_wgt * lum
            call analysis_fill(p_an,istatus,pdgs,wgt_an,icoll)
        endif
      endif
      compute_subtracted_me_0 = jac0(1) * me(1) * lum

      return
      end


      double precision function compute_subtracted_me_1a(x,vegas_wgt,lum,tau,ycm,jac,iqp) 
      ! the muon-gamma (single-real) contribution
      implicit none
      double precision x(8),vegas_wgt,lum,tau,ycm,jac
      integer iqp

      double precision scoll
      common /to_scoll/scoll
      double precision shat
      common /to_shat/shat
      double precision mfin
      common /to_mfin/mfin
      double precision  mmin
      common /to_mmin/mmin
      double precision thresh

      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision shat1a(4), shat1b(4), shat0(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll, isoft
      double precision mu2
      double precision getscale, qprime
      external getscale, qprime
      logical passcuts
      external passcuts
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)

      double precision exi1b,exi2b,exi1a,exi2a,exi12,exi22
      common/cexternal1/exi1b,exi2b,exi1a,exi2a,exi12,exi22

      ! stuff for the analysis
      integer pdgs(6), istatus(6)
      double precision p_an(0:3,6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos
      double precision delta_used
      common/to_delta_used/delta_used
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'input.inc'

      pdgs = (/-13,22,6,-6,-13,0/)
      istatus = (/-1,-1,1,1,1,1/)

      shat = tau * scoll
      thresh = mmin**2/shat

      if (iqp.ne.0) isoft = 0
      if (iqp.eq.0) isoft = 2
      ! generate the momenta for all kinematic configs
      do icoll = 3, 4
        jac1a(icoll) = jac
        jac2(icoll) = jac
        jac1b(icoll) = 0d0
        jac0(icoll) = 0d0
        if (iqp.eq.0) then

          call generate_kinematics(x, shat, thresh, icoll, isoft, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
          call generate_momenta(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
        else 
          call generate_kinematics3(x, shat, thresh, icoll, 0, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
           call generate_momenta3(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll),phi(icoll),icoll,
     &                      jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll),
     &                      shat1a(icoll), shat1b(icoll), shat0(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
          ! overwrite xi1/2
           xi1(icoll) = exi1a 
           xi2(icoll) = exi2a 

          jac1a(icoll) = jac1a(icoll)/2d0/shat1a(icoll)
        endif
      enddo

      do icoll = 3, 4
        mu2 = getscale(scoll, shat) 
        me(icoll) = 0d0
        ! check y against deltaI for the counterterm
        if (y1(3).lt.1d0-delta_used.and.icoll.eq.4)then
            cycle
        endif
        ! boost the momenta to the lab frame. This is needed
        ! both for cuts and for the analysis
        call boost_to_lab_frame(p1a(0,1,icoll),p_an,ycm)
        if (passcuts(p_an,pdgs,istatus)) then 
          call compute_me_singlereal1a(p1a,icoll,y1(icoll),y2(icoll),omy1(icoll),omy2(icoll),xi1,
     &                             xi2,ph1(icoll),ph2(icoll),me(icoll))

          if (fill_histos) then
            wgt_an(1) = jac1a(icoll) * me(icoll) / omy1(3) 
     &           * vegas_wgt * lum
            if (iqp.ne.0) wgt_an(1) = wgt_an(1) * qprime(1d0-xi2(icoll),shat,mu2,deltaI)
            if (icoll.eq.4) wgt_an(1) = - wgt_an(1) 
            call analysis_fill(p_an,istatus,pdgs,wgt_an,icoll)
          endif
        endif
      enddo

      if (iqp.eq.0) then
        compute_subtracted_me_1a =  
     &    (jac1a(3) * me(3) - jac1a(4) * me(4))
     &                       / omy1(3) * lum 
      else
        if (iqp.ne.2) write(*,*) 'ERROR, iqp-1a', iqp
        mu2 = getscale(scoll, shat)
        compute_subtracted_me_1a =  
     &    (jac1a(3) * me(3) * qprime(1-xi2(3),shat,mu2,deltaI) -
     &     jac1a(4) * me(4) * qprime(1-xi2(4),shat,mu2,deltaI)) 
     &                         / omy1(3) * lum 
      endif

      return
      end


      double precision function compute_subtracted_me_1b(x,vegas_wgt,lum,tau,ycm,jac, iqp) 
      ! the gamma-muon (single-real) contribution
      implicit none
      double precision x(8),vegas_wgt,lum,tau,ycm,jac
      integer iqp
      integer pdgs(6), istatus(6)

      double precision scoll
      common /to_scoll/scoll
      double precision shat
      common /to_shat/shat
      double precision mfin
      common /to_mfin/mfin
      double precision  mmin
      common /to_mmin/mmin
      double precision thresh

      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision shat1a(4), shat1b(4), shat0(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll, isoft
      double precision mu2
      double precision getscale, qprime
      external getscale, qprime
      logical passcuts
      external passcuts
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)

      double precision exi1b,exi2b,exi1a,exi2a,exi12,exi22
      common/cexternal1/exi1b,exi2b,exi1a,exi2a,exi12,exi22

      ! stuff for the analysis
      double precision p_an(0:3,6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos
      double precision delta_used
      common/to_delta_used/delta_used
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'input.inc'

      pdgs = (/22,13,6,-6,13,0/)
      istatus = (/-1,-1,1,1,1,1/)

      shat = tau * scoll
      thresh = mmin**2/shat

      if (iqp.ne.0) isoft = 0
      if (iqp.eq.0) isoft = 1
      ! generate the momenta for all kinematic configs
      do icoll = 2, 4, 2
        jac1b(icoll) = jac
        jac2(icoll) = jac
        jac1a(icoll) = 0d0
        jac0(icoll) = 0d0

        if (iqp.eq.0) then

          call generate_kinematics(x, shat, thresh, icoll, isoft, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
          call generate_momenta(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
        else 
          call generate_kinematics3(x, shat, thresh, icoll, 0, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
           call generate_momenta3(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll),phi(icoll),icoll,
     &                      jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll),
     &                      shat1a(icoll), shat1b(icoll), shat0(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
          ! overwrite xi1/2
           xi1(icoll) = exi1b 
           xi2(icoll) = exi2b 

          jac1b(icoll) = jac1b(icoll)/2d0/shat1b(icoll)
        endif
      enddo

      do icoll = 2, 4, 2
        mu2 = getscale(scoll, shat) 
        me(icoll) = 0d0
        ! check y against deltaI for the counterterm
        if (y2(2).lt.1d0-delta_used.and.icoll.eq.4) cycle
        ! boost the momenta to the lab frame. This is needed
        ! both for cuts and for the analysis
        call boost_to_lab_frame(p1b(0,1,icoll),p_an,ycm)
        if (passcuts(p_an,pdgs,istatus)) then 
          call compute_me_singlereal1b(p1b,icoll,y1(icoll),y2(icoll),omy1(icoll),omy2(icoll),xi1,
     &                             xi2,ph1(icoll),ph2(icoll),me(icoll))

          if (fill_histos) then
            wgt_an(1) = jac1b(icoll) * me(icoll) / omy2(2) 
     &           * vegas_wgt * lum
            if (iqp.ne.0) wgt_an(1) = wgt_an(1) * qprime(1d0-xi1(icoll),shat,mu2,deltaI)
            if (icoll.eq.4) wgt_an(1) = - wgt_an(1) 
            call analysis_fill(p_an,istatus,pdgs,wgt_an,icoll)
          endif
        endif
      enddo

      if (iqp.eq.0) then
        compute_subtracted_me_1b =  
     &    (jac1b(2) * me(2) - jac1b(4) * me(4))
     &                         / omy2(2) * lum 
      else
        if (iqp.ne.1) write(*,*) 'ERROR, iqp-1b', iqp
        mu2 = getscale(scoll, shat)
        compute_subtracted_me_1b =  
     &    (jac1b(2) * me(2) * qprime(1-xi1(2),shat,mu2,deltaI) - 
     &     jac1b(4) * me(4) * qprime(1-xi1(4),shat,mu2,deltaI))
     &                         / omy2(2) * lum

      endif

      return
      end


      double precision function compute_subtracted_me_2(x,vegas_wgt,lum,tau,ycm,jac) 
      ! the muon-muon (single-real) contribution
      implicit none
      double precision x(8),vegas_wgt,lum,tau,ycm,jac
      integer pdgs(6), istatus(6)

      double precision scoll
      common /to_scoll/scoll
      double precision shat
      common /to_shat/shat
      double precision mfin
      common /to_mfin/mfin
      double precision  mmin
      common /to_mmin/mmin
      double precision thresh

      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll
      logical passcuts
      external passcuts
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)
      ! stuff for the analysis
      double precision p_an(0:3,6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos
      double precision delta_used
      common/to_delta_used/delta_used

      pdgs = (/-13,13,6,-6,-13,13/)
      istatus = (/-1,-1,1,1,1,1/)

      shat = tau * scoll
      thresh = mmin**2/shat

      ! generate the momenta for all kinematic configs
      do icoll = 1, 4
        jac2(icoll) = jac
        jac1a(icoll) = 0d0
        jac1b(icoll) = 0d0
        jac0(icoll) = 0d0
        call generate_kinematics(x, shat, thresh, icoll, 0,  
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
        call generate_momenta(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
      enddo

C   icoll:
C   1-> doubly resolved collinear emissions
C   2-> single resolved collinear emission (y1=1)
C   3-> single resolved collinear emission (y2=1)
C   4-> no resolved collinear emission (y1=y2=1)
      do icoll = 1, 4
        me(icoll) = 0d0
        ! check y against deltaI for the counterterm
        if (y1(1).lt.1d0-delta_used.and.(icoll.eq.2.or.icoll.eq.4)) cycle
        if (y2(1).lt.1d0-delta_used.and.(icoll.eq.3.or.icoll.eq.4)) cycle
        !if ((y1(1).lt.1d0-deltaI.or.y2(1).lt.1d0-deltaI).and.icoll.eq.4) cycle
        ! boost the momenta to the lab frame. This is needed
        ! both for cuts and for the analysis
        call boost_to_lab_frame(p2(0,1,icoll),p_an,ycm)
        if (passcuts(p_an,pdgs,istatus)) then 
          call compute_me_doublereal(p2,icoll,y1(icoll),y2(icoll),omy1(icoll),omy2(icoll),xi1,
     &                             xi2,ph1(icoll),ph2(icoll),me(icoll))

          if (fill_histos) then
            wgt_an(1) = jac2(icoll) * me(icoll) / omy1(1) / omy2(1) 
     &           * vegas_wgt * lum
            if (icoll.eq.2.or.icoll.eq.3) wgt_an(1) = - wgt_an(1) 
            call analysis_fill(p_an,istatus,pdgs,wgt_an,icoll)
          endif
        endif
      enddo

      !write(*,*) 'YY', y1(1), y2(1), me, jac2
      compute_subtracted_me_2 = 
     &  (jac2(1) * me(1) - jac2(2) * me(2) - jac2(3) * me(3) + jac2(4) * me(4))
     &                       / omy1(1) / omy2(1) * lum

      return
      end




      double precision function compute_subtracted_me_0_qq(x,vegas_wgt,lum,tau,ycm,jac) 
      ! the photon-photon (born-like) contribution
      implicit none
      double precision x(8),vegas_wgt,lum,tau,ycm,jac

      double precision scoll
      common /to_scoll/scoll
      double precision shat
      common /to_shat/shat
      double precision mfin
      common /to_mfin/mfin
      double precision  mmin
      common /to_mmin/mmin
      double precision thresh

      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll
      logical passcuts
      external passcuts
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)
      double precision shat1a(4), shat1b(4), shat0(4)

      double precision exi1b,exi2b,exi1a,exi2a,exi12,exi22
      common/cexternal1/exi1b,exi2b,exi1a,exi2a,exi12,exi22

      ! stuff for the analysis
      integer pdgs(6), istatus(6)
      double precision p_an(0:3,6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos
      double precision getscale, qprime, mu2
      external getscale, qprime
      include 'input.inc'

      pdgs = (/22,22,6,-6,0,0/)
      istatus = (/-1,-1,1,1,1,1/)

      shat = tau * scoll
      thresh = mmin**2/shat

      icoll = 4

      jac0(icoll) = jac
      jac1a(icoll) = 0d0
      jac1b(icoll) = 0d0
      jac2(icoll) = 0d0
      call generate_kinematics3(x, shat, thresh, icoll, 0, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))

      call generate_momenta3(shat, mfin, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll),phi(icoll),icoll,
     &                      jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll),
     &                      shat1a(icoll), shat1b(icoll), shat0(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))

      ! overwrite xi1 and xi2
             xi1(icoll) = exi12 
             xi2(icoll) = exi22 
             
      me(icoll) = 0d0
      ! boost the momenta to the lab frame. This is needed
      ! both for cuts and for the analysis
      call boost_to_lab_frame(p0(0,1,icoll),p_an,ycm)
      !! MZ need to understand how to deal with cuts
      mu2 = getscale(scoll, shat) 
      if (passcuts(p_an,pdgs,istatus)) then 
        call compute_me_born_gaga(p0, me(icoll), icoll)

        if (fill_histos) then
            wgt_an(1) = jac0(icoll) * me(icoll) 
     &           * vegas_wgt * lum / (2*shat0(icoll))
     &         * qprime(1d0-xi1(icoll),shat,mu2,deltaI) * qprime(1d0-xi2(icoll),shat,mu2,deltaI)
            call analysis_fill(p_an,istatus,pdgs,wgt_an,icoll)
        endif
      endif
      compute_subtracted_me_0_qq = jac0(icoll) * me(icoll) * lum / (2*shat0(icoll))
     &         * qprime(1d0-xi1(icoll),shat,mu2,deltaI) * qprime(1d0-xi2(icoll),shat,mu2,deltaI)

      return
      end


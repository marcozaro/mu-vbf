      program driver
      implicit none

      double precision shat, mmin, thresh, jac, mtop
      common /to_shat/shat
      double precision x(10)
      integer i,j
      logical check
      parameter(check=.true.)

      double precision ans
      include 'coupl.inc'

      double precision integrand
      external integrand
      double precision integral,error,prob
      integer nprn
      logical fill_histos
      common /to_fill_histos/fill_histos

      shat = (1000d0)**2
      call setpara('Cards/param_card.dat')
      call printout()
      mtop = mdl_mt
      jac = 1d0
      mmin = 2d0*mtop
      thresh = mmin**2/shat
      !MZ leave this for now, to keep the RN sequence
      call fill_vegas_x(x)

      nprn = 0
      ! fill histogram only in the refine phase
      fill_histos = .false.
      call vegas01(10,integrand,0,10000,
     1        10,nprn,integral,error,prob)

      fill_histos = .true.
      call analysis_begin(1,"central value")
      call vegas01(10,integrand,1,50000,
     1        4,nprn,integral,error,prob)
      call analysis_end(1d0)

      return
      end


      double precision function integrand(x, vegas_wgt)
      implicit none
      double precision x(10), vegas_wgt
      include 'coupl.inc'
      double precision shat
      common /to_shat/shat
      double precision thresh
      double precision mtop, mmin, jac(4), me(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      integer icoll
      logical passcuts2
      external passcuts2
      double precision p2(0:3,6,4), p1a(0:3,5,4), p1b(0:3,5,4),p0(0:3,4,4)
      ! stuff for the analysis
      integer pdgs2(6), status2(6)
      double precision wgt_an(1)
      logical fill_histos
      common /to_fill_histos/fill_histos
      

      mtop = mdl_mt
      mmin = 2d0*mtop
      thresh = mmin**2/shat

      status2 = (/-1,-1,1,1,1,1/)
      pdgs2 = (/-13,13,6,-6,-13,13/)

      ! generate the momenta for all kinematic configs
      do icoll = 1, 4
        jac(icoll) = 1d0
        call generate_kinematics(x, shat, thresh, icoll, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll), jac(icoll))
        call generate_momenta(shat, mtop, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
      enddo

      do icoll = 1, 4
        me(icoll) = 0d0
        if (passcuts2(p2(0,1,icoll))) then 
          call compute_me_doublereal(p2,y1(icoll),y2(icoll),omy1(icoll),omy2(icoll),xi1,
     &                             xi2,ph1(icoll),ph2(icoll),me(icoll))

          if (fill_histos) then
            wgt_an(1) = jac(icoll) * me(icoll) / (1d0-y1(1)) / (1d0-y2(1)) * vegas_wgt
            if (icoll.eq.2.or.icoll.eq.3) wgt_an(1) = - wgt_an(1) 
            call analysis_fill(p2(0,1,icoll),status2,pdgs2,wgt_an,icoll)
          endif
        endif
      enddo
      integrand = jac(1) * me(1) - jac(2) * me(2) - jac(3) * me(3) + jac(4) * me(4)
      integrand = integrand / (1d0-y1(1)) / (1d0-y2(1))

      if (fill_histos) call HwU_add_points()

      return
      end


      subroutine fill_vegas_x(x)
C     fill the vegas x.
      implicit none
      double precision x(10)
      integer i
      double precision ran2
      external ran2
      do i = 1,10
         x(i) = ran2()
      enddo

      return 
      end


      subroutine check_momenta(pp,n)
C    check momentum conservation and on-shell relations
      implicit none
      double precision pp(0:3,*), mass
      integer n
      double precision tiny
      parameter (tiny=1d-5)
      double precision ptot(0:3), etot
      integer i, j
      double precision dot
      external dot
      include 'coupl.inc'
      mass = mdl_mt

      etot = pp(0,1) + pp(0,2)

      if (etot.lt.2*mass) then
          write(*,*) 'ERROR1', etot, mass
          call write_momenta(pp,n)
            call backtrace()
            stop
      endif

      ptot(:) = 0d0

      do i = 1, n
        if (i.le.2) ptot(:) = ptot(:) - pp(:,i)
        if (i.gt.2) ptot(:) = ptot(:) + pp(:,i)
      enddo

      do j = 0,3
        if (abs(ptot(j))/etot.gt.tiny) then
            write(*,*) 'ERROR2', j, ptot
            call write_momenta(pp,n)
            call backtrace()
            stop
        endif
      enddo

      ! check massless momenta
      do i = 1,n
        if (i.eq.3.or.i.eq.4) then
          if ((dot(pp(0,i),pp(0,i))-mass**2)/etot.gt.tiny) then
              write(*,*) 'ERROR3', i, dot(pp(0,i),pp(0,i)), etot
            call write_momenta(pp,n)
            call backtrace()
            stop
          endif

        else
          if (abs(dot(pp(0,i),pp(0,i))/etot**2).gt.tiny) then
              write(*,*) 'ERROR4', i, dot(pp(0,i),pp(0,i)), etot
            call write_momenta(pp,n)
            call backtrace()
            stop
          endif
        endif
      enddo 
      

      return
      end


      subroutine write_momenta(pp,n)
      implicit none
      double precision pp(0:3,*)
      integer n
      integer i

      do i = 1,n
        write(*,*) 'p(:,',i,')=(/', pp(0,i),',', pp(1,i),',',
     #                              pp(2,i),',',pp(3,i),'/)'
      enddo

      return
      end


      subroutine generate_momenta(shat, m, y1, y2, xi1, xi2, ph1, ph2, cth, phi,
     &                      p2, p1a, p1b, p0)
      implicit none
C starting from the kinematic variables, generates 4 sets of momenta:
C  - p2 with the double real emission
C  - p1a/b with a single real emission, from the first/second leg
C  - p0 without emissions
      double precision shat, m
      double precision y1, y2, xi1, xi2, ph1, ph2, phi, cth
      double precision p0(0:3,4), p1a(0:3,5), p1b(0:3,5), p2(0:3,6)
      double precision omega, sborn
      double precision prec(0:3)

      p0(:,:)=0d0
      p1a(:,:)=0d0
      p1b(:,:)=0d0
      p2(:,:)=0d0

!      write(*,*) 'shat =',shat
!      write(*,*) 'm=', m
!      write(*,*) 'y1=', y1
!      write(*,*) 'y2=', y2
!      write(*,*) 'xi1=', xi1
!      write(*,*) 'xi2=', xi2
!      write(*,*) 'ph1=', ph1
!      write(*,*) 'ph2=', ph2
!      write(*,*) 'cth=', cth
!      write(*,*) 'phi=', phi

      ! p0
      call generate_is(shat, p0(0,1))
      call generate_born_fs(shat,m,cth,phi,p0(0,3))
      ! p1a
      call generate_is(shat, p1a(0,1))
      call generate_born_fs(shat*(1d0-xi1),m,cth,phi,p1a(0,3))
      call generate_fks_momentum(shat,xi1,y1,ph1,1,p1a(0,5))
      call boost_momenta_born(p1a(0,3),p1a(0,5))
      ! p1b
      call generate_is(shat, p1b(0,1))
      call generate_born_fs(shat*(1d0-xi2),m,cth,phi,p1b(0,3))
      call generate_fks_momentum(shat,xi2,y2,ph2,2,p1b(0,5))
      call boost_momenta_born(p1b(0,3),p1b(0,5))
      ! p2
      omega = dsqrt(1d0-y1**2)*dsqrt(1d0-y2**2)*dcos(ph1-ph2)-y1*y2
      sborn = shat*(1d0-xi1-xi2+xi1*xi2*(1d0-omega)/2d0) 
      call generate_is(shat, p2(0,1))
      call generate_born_fs(sborn,m,cth,phi,p2(0,3))
      call generate_fks_momentum(shat,xi1,y1,ph1,1,p2(0,5))
      call generate_fks_momentum(shat,xi2,y2,ph2,2,p2(0,6))
      prec(:) = p2(:,5)+p2(:,6)
      call boost_momenta_born(p2(0,3),prec)

      return 
      end


      subroutine boost_momenta_born(pcm, prec)
      implicit none
C boost the born momenta in ther c.o.m so that they
C recoil against prec
      double precision pcm(0:3,2), prec(0:3)
      double precision pboost(0:3), minvb
      double precision sumdot
      external sumdot
      integer i

      double precision dot
      external dot

      ! the invariant mass of the born
      minvb = sumdot(pcm(0,1),pcm(0,2),1d0)
      ! the total momentum of the born system after the boost
      pboost(1:3) = -prec(1:3)
      pboost(0) = dsqrt(minvb+prec(1)**2+prec(2)**2+prec(3)**2)

      do i = 1,2
        call boostx(pcm(0,i),pboost,pcm(0,i))
      enddo

      return
      end


      subroutine generate_fks_momentum(shat,xi,y,ph,ileg,p)
      implicit none
C generate the FKS radiated particle, ie a massless particle
C with energy xi, angles y and phi, w.r.t. the initial leg i
      double precision shat,xi,y,ph
      integer ileg
      double precision p(0:3)
      double precision sth

      sth = dsqrt(max(1d0-y**2,0d0))

      p(:) = xi*dsqrt(shat)/2d0
      p(1) = p(1)*sth*dcos(ph)
      p(2) = p(2)*sth*dsin(ph)
      if (ileg.eq.1) then
          p(3) = p(3)*y
      else if (ileg.eq.2) then
          p(3) =-p(3)*y
      endif
      return
      end


      subroutine generate_is(shat, p)
      implicit none
C generate the center-of-mass initial-state momenta
C along the beam axis
      double precision shat, p(0:3,2) ! read only the first two momenta
      double precision sqrtso2


      sqrtso2 = sqrt(shat)/2d0
      p(0,1) = sqrtso2
      p(3,1) = sqrtso2
      p(0,2) = sqrtso2
      p(3,2) =-sqrtso2
      return
      end


      subroutine generate_born_fs(shat,m,cth,phi,p)
      implicit none
C generate the Born outgoing momenta, with angles ct,phi
C in their partonic c.o.m fram
      double precision shat, m, cth, phi, p(0:3,2) ! read only the first two momenta
      double precision pmod,sth

      sth = dsqrt(max(1d0-cth**2,0d0))

      pmod = sqrt(shat/4-m**2)
      p(0,1) = sqrt(shat)/2d0
      p(1,1) = pmod*sth*dcos(phi)
      p(2,1) = pmod*sth*dsin(phi)
      p(3,1) = pmod*cth
      p(0,2) = sqrt(shat)/2d0
      p(1,2) =-pmod*sth*dcos(phi)
      p(2,2) =-pmod*sth*dsin(phi)
      p(3,2) =-pmod*cth
      return
      end


      subroutine generate_kinematics(x, shat, thresh, icoll, y1, y2, omy1, omy2, xi1, xi2, ph1, ph2, phi, cth, jac)
      implicit none
C generates the kinematic variables (y_i,xi_i,ph_i, i=1,2) for each of
C the two collinear splittings
C Generates phi,cth, the angles in the 2->2 scattering
C   icoll:
C   1-> doubly resolved collinear emissions
C   2-> single resolved collinear emission (y1=1)
C   3-> single resolved collinear emission (y2=1)
C   4-> no resolved collinear emission (y1=y2=1)
C
C  jac includes the PS volumes, times flux (1/2shat) and converted to PB.
C  omy = 1-y (for better numerical accuracy)
      double precision x(10)
      double precision shat, thresh
      integer icoll
      double precision y1, y2, omy1, omy2, xi1, xi2, ph1, ph2, phi, cth, jac
      double precision omega, sborn
      double precision pi
      parameter (pi=3.14159265359d0)
      double precision conv
      parameter (conv=389379.66*1000)  !conv to picobarns
      ! born angles
      cth = x(1)*2d0-1d0
      jac = jac*2d0
      phi = x(2)*2d0*pi
      jac = jac*2d0*pi

      ! y, adaptive towards y->1
      if (icoll.ne.2.and.icoll.ne.4) then
        y1 = -2d0*x(3)**2+1d0
        omy1 = 2d0*x(3)**2
      else
        y1 = 1d0
        omy1 = 0d0
      endif
      jac = jac*4d0*x(3)
      if (icoll.ne.3.and.icoll.ne.4) then
        y2 = -2d0*x(4)**2+1d0
        omy2 = 2d0*x(4)**2
      else
        y2 = 1d0
        omy2 = 0d0
      endif
      jac = jac*4d0*x(4)
      ! phi, flat
      ph1 = x(5)*2d0*pi
      jac = jac*2d0*pi
      ph2 = x(6)*2d0*pi
      jac = jac*2d0*pi
      !write(*,*) 'FORCING y1', y1fix
      !y1 = y1fix
      !write(*,*) 'FORCING y2', y2fix
      !y2 = y2fix

      ! xi1/2 following the formulae on the note.
      ! randomize which one is generated first
      omega = sqrt(1d0-y1**2)*sqrt(1d0-y2**2)*dcos(ph1-ph2)-y1*y2

      if (x(7).lt.0.5d0) then
          xi1 = x(7)*2d0*(1d0-thresh)
          jac = jac*2d0*(1d0-thresh)
          xi2 = x(8)*2d0*(1d0-thresh-xi1)/(2d0-(1d0-omega)*xi1)
          jac = jac*2d0*(1d0-thresh-xi1)/(2d0-(1d0-omega)*xi1)
      else
          xi2 = (x(7)-0.5d0)*2d0*(1d0-thresh)
          jac = jac*2d0*(1d0-thresh)
          xi1 = x(8)*2d0*(1d0-thresh-xi2)/(2d0-(1d0-omega)*xi2)
          jac = jac*2d0*(1d0-thresh-xi2)/(2d0-(1d0-omega)*xi2)
      endif

      ! finally, turn jac into the proper phase-space volume
      ! this is the contribution from the two radiations
      jac = jac * (shat / 64d0 / pi**3)**2 * xi1 * xi2
      ! this is for the underlying born
      sborn = shat*(1d0-xi1-xi2+xi1*xi2*(1d0-omega)/2d0) 
      ! 1/32pi^2 p/sqrt[mt^2+p^2]
      jac = jac / 32d0 / pi**2 * sqrt(1d0-thresh*shat/sborn)
      ! to pb and flux
      jac = jac * conv / 2d0 /shat

      if (jac.ne.jac) jac = 0d0

      return
      end



      DOUBLE PRECISION FUNCTION SumDot(P1,P2,dsign)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3),dsign
c
c     Local
c      
      integer i
      double precision ptot(0:3)
c
c     External
c
      double precision dot
      external dot
c-----
c  Begin Code
c-----

      do i=0,3
         ptot(i)=p1(i)+dsign*p2(i)
      enddo
      SumDot = dot(ptot,ptot)
      RETURN
      END


      double precision function dot(p1,p2)
C****************************************************************************
C     4-Vector Dot product
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)

      if(dabs(dot).lt.1d-6)then ! solve numerical problem 
         dot=0d0
      endif

      end

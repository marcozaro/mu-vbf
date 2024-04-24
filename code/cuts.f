      logical function passcuts(p,pdgs,istatus)
      implicit none
      double precision p(0:3,6)
      integer pdgs(6), istatus(6)
      double precision pmumu(0:3)
      double precision dot
      external dot
      integer i
      logical foundmup, foundmum
      double precision rap2
      external rap2

      include 'input.inc'

      ! we ask that mu+ has rapidity larger than ymin
      ! and that mu- has rapidity smaller than -ymin
      ! mass is far away from the z (we pick mmin gev as threshold)
      passcuts = .true.
      foundmup = .false.
      foundmum = .false.
      pmumu(:) = 0d0
      do i = 1, 6
        if (istatus(i).eq.1.and.pdgs(i).eq.-13) then
         ! mu+ must have positive z component
         if (rap2(p(0,i)).lt.ymin) passcuts = .false.
         pmumu(:) = pmumu+p(:,i)
         foundmup = .true.
        endif
        if (istatus(i).eq.1.and.pdgs(i).eq.+13) then
         ! mu- must have negative z component
         if (rap2(p(0,i)).gt.-ymin) passcuts = .false.
         pmumu(:) = pmumu+p(:,i)
         foundmum = .true.
        endif

      enddo
      if (foundmup.and.foundmum) then
        if (dot(pmumu,pmumu).lt.mmin**2) passcuts = .false.
      endif

      if (foundmup.and.foundmum) then
        ! MZ to single out only the resolved contribution
        !if (p(1,5)**2+p(2,5)**2.lt.2500d0) passcuts = .false.
        !if (p(1,6)**2+p(2,6)**2.lt.2500d0) passcuts = .false.
        continue
      else if (foundmup.or.foundmum) then
        !if (p(1,5)**2+p(2,5)**2.lt.2500d0) passcuts = .false.
        continue
      endif

      return
      end



      DOUBLE PRECISION  FUNCTION rap2(p)
c************************************************************************
c     Returns rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision  p(0:3)
c
c     Local
c
      double precision pm
c
c     Global
c
c-----
c  Begin Code
c-----
c      pm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      pm = p(0)
      rap2 = .5d0*dlog((pm+p(3))/(pm-p(3)))
      end

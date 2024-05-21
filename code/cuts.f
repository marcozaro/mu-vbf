      logical function passcuts(p,pdgs,istatus)
      implicit none
      double precision p(0:3,6)
      integer pdgs(6), istatus(6)
      double precision pmumu(0:3)
      double precision dot
      external dot
      integer i, imup, imum
      logical foundmup, foundmum
      double precision rap2
      external rap2
      double precision pbeam(0:3,2)

      include 'input.inc'

      ! define the beam momenta
      pbeam(:,1) = (/ecm/2d0, 0d0, 0d0, ecm/2d0/)
      pbeam(:,2) = (/ecm/2d0, 0d0, 0d0,-ecm/2d0/)

      ! we ask that mu+ has rapidity larger than ymin
      ! and that mu- has rapidity smaller than -ymin
      ! mass is far away from the z (we pick mmin gev as threshold)

      imup = 0
      imum = 0
      passcuts = .true.
      foundmup = .false.
      foundmum = .false.
      pmumu(:) = 0d0
      do i = 1, 6
        if (istatus(i).eq.1.and.pdgs(i).eq.-13) then
         ! mu+ must have positive z component
         if (rap2(p(0,i)).lt.ymin) passcuts = .false.
         pmumu(:) = pmumu(:)+p(:,i)
         foundmup = .true.
         imup = i
        endif
        if (istatus(i).eq.1.and.pdgs(i).eq.+13) then
         ! mu- must have negative z component
         if (rap2(p(0,i)).gt.-ymin) passcuts = .false.
         pmumu(:) = pmumu(:)+p(:,i)
         foundmum = .true.
         imum = i
        endif
        if (istatus(i).eq.1.and.abs(pdgs(i)).eq.6) then
            if (ymaxtop.gt.0d0) then
                if (abs(rap2(p(0,i))).gt.ymaxtop) passcuts = .false.
            endif
        endif

      enddo

      ! MZ if mup/m have not been found, define them starting from the
      ! remnants
      if (.not. foundmup) then
          if (pdgs(1).ne.22) write(*,*) 'error mup', pdgs, foundmup
          foundmup = .true.
          pmumu(:) = pmumu(:) + pbeam(:,1) - p(:,1)
      endif
      if (.not. foundmum) then
          if (pdgs(2).ne.22) write(*,*) 'error mum', pdgs, foundmum
          foundmum = .true.
          pmumu(:) = pmumu(:) + pbeam(:,2) - p(:,2)
      endif


      if (foundmup.and.foundmum) then
        if (dot(pmumu,pmumu).lt.mllmin**2) passcuts = .false.
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

CC      if (imum.eq.0) then
CC          passcuts = .false.
CC          return
CC      endif
CC      if (foundmum) then
CC        if (p(1,imum)**2+p(2,imum)**2.lt.100d0**2) passcuts = .false.
CC      endif

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

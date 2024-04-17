      logical function passcuts(p,pdgs,istatus)
      implicit none
      double precision p(0:3,6)
      integer pdgs(6), istatus(6)
      double precision pmumu(0:3)
      double precision dot
      external dot
      integer i
      logical foundmup, foundmum
      ! to be changed/improved.
      !for the moment, we ask that the mu+(-) momentum (n 5(6)) 
      ! has z component > (<)0, and that the resulting invariant
      ! mass is far away from the z (we pick 200 gev as threshold)
      passcuts = .true.
      foundmup = .false.
      foundmum = .false.
      pmumu(:) = 0d0
      do i = 1, 6
        if (istatus(i).eq.1.and.pdgs(i).eq.-13) then
         ! mu+ must have positive z component
         if (p(3,i).lt.0d0) passcuts = .false.
         pmumu(:) = pmumu+p(:,i)
         foundmup = .true.
        endif
        if (istatus(i).eq.1.and.pdgs(i).eq.+13) then
         ! mu- must have negative z component
         if (p(3,i).gt.0d0) passcuts = .false.
         pmumu(:) = pmumu+p(:,i)
         foundmum = .true.
        endif

      enddo
      if (foundmup.and.foundmum) then
        if (dot(pmumu,pmumu).lt.200d0**2) passcuts = .false.
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

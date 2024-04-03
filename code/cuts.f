      logical function passcuts2(p)
      implicit none
      double precision p(0:3,6)
      double precision pmumu(0:3)
      double precision dot
      external dot
      ! to be changed/improved.
      !for the moment, we ask that the mu+(-) momentum (n 5(6)) 
      ! has z component > (<)0, and that the resulting invariant
      ! mass is far away from the z (we pick 200 gev as threshold)
      passcuts2 = .true.
      if (p(3,5).lt.0d0.or.p(3,6).gt.0) passcuts2 = .false.
      pmumu(:) = p(:,5) + p(:,6)
      if (dot(pmumu,pmumu).lt.200d0**2) passcuts2 = .false.

      ! MZ to single out only the resolved contribution
      if (p(1,5)**2+p(2,5)**2.lt.2500d0) passcuts2 = .false.
      if (p(1,6)**2+p(2,6)**2.lt.2500d0) passcuts2 = .false.

      return
      end

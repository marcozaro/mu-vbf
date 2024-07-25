      program check
      implicit none
      integer ip, iy, iph
      integer i, icoll, isoft

      double precision p(0:3,6,4)
      ! the last index of the momenta array:
      ! 1-> doubly resolved collinear emissions
      ! 2-> single resolved collinear emission (y1=1)
      ! 3-> single resolved collinear emission (y2=1)<-
      ! 4-> no resolved collinear emission (y1=y2=1)<-
      !!! note that xi are different in the various kinematics
      double precision jac2(4), jac1a(4), jac1b(4), jac0(4), me(4)
      double precision y1(4), y2(4), omy1(4), omy2(4), xi1(4), xi2(4), ph1(4), ph2(4), phi(4), cth(4)
      double precision p2(0:3,6,4), p1a(0:3,6,4), p1b(0:3,6,4),p0(0:3,6,4)
      double precision ans
      double precision x(12)
      double precision ran2
      external ran2
      double precision scoll, shat
      common /to_scoll/scoll
      double precision jac
      double precision thresh
      include 'coupl.inc'
      include 'input.inc'

      double precision tiny
      common/to_coll_cutoff/tiny

      common /to_shat/shat

      call setpara('Cards/param_card.dat')

      tiny = 1d-10


      scoll = 1d4**2
      shat = scoll
      thresh = 4*mdl_mt**2/shat
      isoft = 2 ! 1a
      isoft = 1 ! 1b


      do ip = 1,5
        write(*,*) 'NEWPOINT'
        do i = 1,12
          x(i) = ran2()
        enddo
        do iy = 1,11
          x(3) = 10d0**(-iy/2d0)! for y1
          do iph = 1,21
            x(5) = dble(iph)/21d0
            do icoll = 3, 4
            jac = 1d0
            jac1a(icoll) = jac


            call generate_kinematics(x, shat, thresh, icoll, 2, 
     &       y1(icoll), y2(icoll), omy1(icoll), omy2(icoll), xi1(icoll), xi2(icoll), 
     &       ph1(icoll), ph2(icoll), phi(icoll), cth(icoll),
     &       jac2(icoll), jac1a(icoll), jac1b(icoll), jac0(icoll))
            call generate_momenta(shat, mdl_mt, y1(icoll), y2(icoll), xi1(icoll), xi2(icoll), 
     &                      ph1(icoll), ph2(icoll), cth(icoll), phi(icoll),
     &                      p2(0,1,icoll), p1a(0,1,icoll), p1b(0,1,icoll), p0(0,1,icoll))
            enddo

            do icoll = 3, 4
            call compute_me_singlereal1a(p1a,y1(icoll),y2(icoll),omy1(icoll),omy2(icoll),xi1,
     &                             xi2,ph1(icoll),ph2(icoll),me(icoll))
            enddo
            write(*,*)  y1(3), omy1(3), ph1(3),  me(3), me(4)
          enddo
        enddo
      enddo

      return
      end


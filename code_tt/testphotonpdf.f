      program test
      implicit none
      double precision pi
      parameter (pi=3.14159265359d0)
      real*8 me
      data me /0.10566d0/

      double precision q2,x

      double precision gampdf, gampdf_emela, gampdf_ww
      double precision q2max, q2min
      integer i,j

      !gampdf = dble(gal(1)**2)/8d0/pi**2*

      call initfromgrid_lhaid(10005)

      q2 = me**2
      open(file = 'plot_gammapdf/gampdf_mmu.txt', unit=10, 
     #     status = 'unknown')
      do i = 1,100
      x = dble(i)/100d0
      gampdf = 0.007755544d0/2d0/pi*
     &           (2-2d0*x+x*x)/x*(dlog(q2/me**2)-2d0*dlog(x)-1)

      q2min= me**2*x**2/(1-x)
      q2max=q2
      if(q2min.lt.q2max) then 
             gampdf_ww =  0.007755544d0/2d0/pi*
     &           (2d0*me**2*x**2*(-1/q2min+1/q2max)+
     &           (2-2d0*x+x*x)*dlog(q2max/q2min))/x
      endif
      call elpdfq2(0, 22, x, 1d0-x, q2, 1d0, gampdf_emela)
      write(10,*)  x, gampdf, gampdf_ww, gampdf_emela
      enddo
      close(10)

      q2 = 50d0**2
      open(file = 'plot_gammapdf/gampdf_50.txt', unit=10, 
     #     status = 'unknown')
      do i = 1,100
      x = dble(i)/100d0
      gampdf = 0.007755544d0/2d0/pi*
     &           (2-2d0*x+x*x)/x*(dlog(q2/me**2)-2d0*dlog(x)-1)

      q2min= me**2*x**2/(1-x)
      q2max=q2
      if(q2min.lt.q2max) then 
             gampdf_ww =  0.007755544d0/2d0/pi*
     &           (2d0*me**2*x**2*(-1/q2min+1/q2max)+
     &           (2-2d0*x+x*x)*dlog(q2max/q2min))/x
      endif
      call elpdfq2(0, 22, x, 1d0-x, q2, 1d0, gampdf_emela)
      write(10,*)  x, gampdf, gampdf_ww, gampdf_emela
      enddo
      close(10)

      q2 = 1000d0**2
      open(file = 'plot_gammapdf/gampdf_1000.txt', unit=10, 
     #     status = 'unknown')
      do i = 1,100
      x = dble(i)/100d0
      gampdf = 0.007755544d0/2d0/pi*
     &           (2-2d0*x+x*x)/x*(dlog(q2/me**2)-2d0*dlog(x)-1)

      q2min= me**2*x**2/(1-x)
      q2max=q2
      if(q2min.lt.q2max) then 
             gampdf_ww =  0.007755544d0/2d0/pi*
     &           (2d0*me**2*x**2*(-1/q2min+1/q2max)+
     &           (2-2d0*x+x*x)*dlog(q2max/q2min))/x
      endif
      call elpdfq2(0, 22, x, 1d0-x, q2, 1d0, gampdf_emela)
      write(10,*)  x, gampdf, gampdf_ww, gampdf_emela
      enddo
      close(10)

      return
      end



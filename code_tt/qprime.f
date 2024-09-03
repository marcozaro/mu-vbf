      double precision function qprime(z,shat,mu2,delta)
      ! the Qprime term of the note, including alpha/2pi
      implicit none 
      double precision z,shat,mu2,delta
      double precision Kscheme
      external kscheme
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'coupl.inc'

      if (delta.gt.2d0.or.delta.le.0d0) then
          write(*,*)'Error, deltaI must be in (0,2]', delta
          stop 1
      endif
      qprime = (1d0+(1d0-z)**2) / z * (dlog(shat*delta/mu2/2d0)  
     $                                  +2d0*dlog(1d0-z)) + z - kscheme(z)
      !include alpha/2pi
      qprime = qprime * dble(gal(1))**2/8d0/pi**2
      return 
      end

      
      subroutine generate_qp_z(x, zmin, z, jac)
      ! generate z in the range [zmin, 1], using x as random number in [0,1]
      ! use an adaptive generation towards 1, to suppress the log(1-z) in Q
      implicit none
      double precision x, zmin, z, jac

      z = zmin + x * (1d0-zmin)
      z = 1d0 - x**2 * (1d0-zmin)
      jac = jac * 2*x*(1d0-zmin)
      return
      end


      double precision function Pgamu(z)
      implicit none
      double precision z
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'coupl.inc'

      Pgamu = (1+(1-z)**2)/z
      !include alpha/2pi
      Pgamu = Pgamu * dble(gal(1))**2/8d0/pi**2

      return
      end

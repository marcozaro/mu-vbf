      double precision function qprime(z,shat,mu2)
      ! the Qprime term of the note, including alpha/2pi
      implicit none 
      double precision z,shat,mu2
      double precision K
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'coupl.inc'
      include 'input.inc'

      k=0d0 ! put here change of scheme
      if (deltaI.gt.2d0) then
          write(*,*)'Error, deltaI must be in (0,2]', deltaI
          stop 1
      endif
      qprime = (1d0+(1d0-z)**2) / z * (dlog(shat*deltaI/mu2/2d0)  
     $                                  +2d0*dlog(1d0-z)) + z - k
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


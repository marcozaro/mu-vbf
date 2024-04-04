      double precision function qprime(z,shat,mu2)
      ! the Qprime term of the note, including alpha/2pi
      implicit none 
      double precision z,shat,mu2
      double precision K
      double precision pi
      parameter (pi=3.14159265359d0)
      include 'coupl.inc'

      k=0d0 ! put here change of scheme
      ! MZMZMZ the log(1-z) part may not converge
      qprime = (1d0+(1d0-z)**2) / z * (dlog(shat/mu2)+2d0*dlog(1d0-z)) + z - k
      !include alpha/2pi
      qprime = qprime * dble(gal(1))**2/8d0/pi
      return 
      end

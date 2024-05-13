      subroutine get_lum(ilum,rnd,scoll,smin,jac,lum,tau,ycm,x1,x2)
      ! ilum = 1 -> mu+ mu-
      ! ilum = 2 -> gamma mu-
      ! ilum = 3 -> mu+ gamma
      ! ilum = 4 -> gamma gamma
      implicit none
      integer ilum
      double precision rnd(2), scoll, smin, jac, lum, tau, ycm

      double precision x1, omx1, x2, omx2, jac_ee, mu2

      double precision tolerance
      parameter (tolerance=1e-3)
      double precision lx1, lx2

      double precision mupdf, gampdf, getscale
      external mupdf, gampdf, getscale

      ! generate the bjorken x's
      if (ilum.eq.1) then 
        ! mu mu scattering
        call generate_x_ee(rnd(1), smin/scoll, x1, omx1, jac_ee)
        jac = jac * jac_ee
        call generate_x_ee(rnd(2), smin/scoll/x1, x2, omx2, jac_ee)
        jac = jac * jac_ee
        mu2 = getscale(scoll, x1, x2)

        lum = mupdf(x1, omx1, mu2)
        lum = lum*mupdf(x2, omx2, mu2)

      else if (ilum.eq.2) then
        ! gamma mu scattering
        call generate_x_gam(rnd(1), smin/scoll, x1, omx1, jac_ee)
        jac = jac * jac_ee
        call generate_x_ee(rnd(2), smin/scoll/x1, x2, omx2, jac_ee)
        jac = jac * jac_ee
        mu2 = getscale(scoll, x1, x2)

        lum = gampdf(x1, omx1, mu2)
        lum = lum*mupdf(x2, omx2, mu2)

      else if (ilum.eq.3) then
        ! mu gamma scattering
        call generate_x_ee(rnd(1), smin/scoll, x1, omx1, jac_ee)
        jac = jac * jac_ee
        call generate_x_gam(rnd(2), smin/scoll/x1, x2, omx2, jac_ee)
        jac = jac * jac_ee
        mu2 = getscale(scoll, x1, x2)

        lum = mupdf(x1, omx1, mu2)
        lum = lum*gampdf(x2, omx2, mu2)

      else if (ilum.eq.4) then
        ! gamma gamma scattering
        call generate_x_gam(rnd(1), smin/scoll, x1, omx1, jac_ee)
        jac = jac * jac_ee
        call generate_x_gam(rnd(2), smin/scoll/x1, x2, omx2, jac_ee)
        jac = jac * jac_ee
        mu2 = getscale(scoll, x1, x2)

        lum = gampdf(x1, omx1, mu2)
        lum = lum*gampdf(x2, omx2, mu2)
      else
        write(*,*) 'ERROR: wrong ilum', ilum
        stop 1
      endif

      ! compute tau and y
      tau = x1*x2

      if (omx1.gt.tolerance) then
        lx1 = dlog(x1)
      else
        lx1 = -omx1-omx1**2/2d0-omx1**3/3d0-omx1**4/4d0-omx1**5/5d0
      endif
      ycm = 0.5d0*lx1

      if (omx2.gt.tolerance) then
        lx2 = dlog(x2)
      else
        lx2 = -omx2-omx2**2/2d0-omx2**3/3d0-omx2**4/4d0-omx2**5/5d0
      endif
      ycm = ycm-0.5d0*lx2

      return
      end


      double precision function getscale(scoll,x1,x2)
      implicit none
      double precision scoll, x1, x2
      include 'input.inc'

      if (imuf.eq.0) then
        getscale = scoll
      else if (imuf.eq.1) then
        getscale = scoll*x1*x2
      else if (imuf.eq.-1) then
        getscale = fixscale**2 
      endif

      return
      end


      subroutine generate_x_ee(rnd, xmin, x, omx, jac)
      implicit none
      ! generates the momentum fraction with importance
      !  sampling suitable for ee collisions
      ! rnd is generated uniformly in [0,1], 
      ! x is generated according to (1 -rnd)^-expo, starting
      ! from xmin
      ! jac is the corresponding jacobian
      ! omx is 1-x, stored to improve numerical accuracy
      double precision rnd, x, omx, jac, xmin
      double precision expo
      double precision get_ee_expo
      double precision tolerance
      parameter (tolerance=1.d-5)
      include 'input.inc'

      x = 1d0
      omx = 0d0
      jac = 1d0
      ! if we do not need to convolve with the muon PDF, just return here
      if (convolvemuon.eq.0) return

      expo = get_ee_expo()

      x = 1d0 - rnd ** (1d0/(1d0-expo))
      omx = rnd ** (1d0/(1d0-expo))
      if (x.ge.1d0) then
        if (x.lt.1d0+tolerance) then
          x=1d0
        else
          write(*,*) 'ERROR in generate_x_ee', rnd, x
          stop 1
        endif
      endif
      jac = 1d0/(1d0-expo) 
      ! then rescale it between xmin and 1
      x = x * (1d0 - xmin) + xmin
      omx = omx * (1d0 - xmin)
      jac = jac * (1d0 - xmin)**(1d0-expo)

      return 
      end


      subroutine generate_x_gam(rnd, xmin, x, omx, jac)
      implicit none
      ! generates the momentum fraction with importance
      !  sampling suitable for gamma e collisions
      ! rnd is generated uniformly in [0,1], 
      ! x is generated according to 1/rnd, starting
      ! from xmin
      ! jac is the corresponding jacobian
      ! omx is 1-x, stored to improve numerical accuracy
      double precision rnd, x, omx, jac, xmin
      double precision tolerance
      parameter (tolerance=1.d-5)


      x = exp(rnd*dlog(xmin))
      omx = 1d0-x 
      jac = abs(log(xmin))
      if (x.ge.1d0) then
        if (x.lt.1d0+tolerance) then
          x=1d0
        else
          write(*,*) 'ERROR in generate_x_ee', rnd, x
          stop 1
        endif
      endif

      return 
      end


      double precision function mupdf(x,omx,Q2)
      implicit none
      ! a wrapper to whatever implementation of the  muon PDF
      ! we are employing
      ! It includes a factor (1-x)^ get_ee_expo (included in the PS
      ! jacobian)
      double precision x, omx, q2
      double precision eepdf_tilde, eepdf_tilde_power, get_ee_expo
      external eepdf_tilde, eepdf_tilde_power, get_ee_expo
      double precision k_exp, ps_expo
      include 'input.inc'

      if (convolvemuon.eq.0) then
          mupdf = 1d0
          return
      endif

      mupdf = eepdf_tilde(x,Q2,1,13,13)
      k_exp = eepdf_tilde_power(Q2,1,13,13)
      ps_expo = get_ee_expo()

      if (k_exp.gt.ps_expo) then
          write(*,*) 'WARNING, e+e- exponent exceeding limit', k_exp, ps_expo
          stop 1
      endif
      mupdf = mupdf * omx**(-k_exp+ps_expo)
      return
      end


      double precision function gampdf(x,omx,Q2)
      implicit none
      ! a wrapper to whatever implementation of the photon PDF
      ! we are employing
      ! The function is multiplied by a factor x (included in the PS
      ! jacobian)
      ! For the moment, we use the crude WW formula
      include 'coupl.inc'
      double precision pi
      parameter (pi=3.14159265359d0)
      real*8 me
      data me /0.105658d0/
      double precision x, omx, q2
      double precision q2max,q2min

      q2min= me**2*x**2/(1-x)
      q2max=q2
      if(q2min.lt.q2max) then 
             gampdf = dble(gal(1)**2)/8d0/pi**2*
     &           (2d0*me**2*x**2*(-1/q2min+1/q2max)+
     &           (2-2d0*x+x*x)*dlog(q2max/q2min))
      endif

      return
      end


      double precision function get_ee_expo()
      ! return the exponent used in the
      ! importance-sampling transformation to sample
      ! the Bjorken x's
      implicit none
      integer idbeam
      double precision expo_e, expo_m
      parameter (expo_e=0.96d0)
      parameter (expo_m=0.975d0)
      get_ee_expo = expo_m
      return
      end


c     This function return the power of (1-x)
      real*8 function eepdf_tilde_power(Q2,n,partonid,beamid)
      implicit none
      real*8 me
      data me /0.105658d0/
      real*8 PI
      real*8 alphaem
c     In Gmu scheme
      data alphaem/0.007562397d0/
      real*8 beta,eta,Q2
      integer n,partonid,beamid
      real*8 k,b

      PI=4.D0*DATAN(1.D0)
      beta = alphaem/PI * (dlog(Q2/me/me)-1)
      eta = alphaem/PI * (dlog(Q2/me/me))
      b=-2.D0/3.D0

c     muon beam
      if (beamid .eq. 13) then
c     other partons are zero
        if (partonid .ne. 13) then
          k=0d0
        else
          if (n .eq. 1) then
            k=1d0-beta
          else
            k=0d0
          endif
        endif
      else if (beamid .eq. -13) then
        if (partonid .ne. -13) then
          k=0d0
        else
          if (n .eq. 1) then
            k=1d0-beta
          else
            k=0d0
          endif
        endif
      endif
      eepdf_tilde_power = k
      end


c     This function calculate the reduced structure function, with energy
c     fraction given by "x", at scale "Qsquare"
      real*8 function eepdf_tilde(x,Qsquare,n,partonid,beamid)
      implicit none
      real*8 x,Qsquare
      real*8 me
      data me /0.105658d0/
      real*8 PI
      real*8 alphaem
c     In Gmu scheme
      data alphaem/0.007562397d0/
      real*8 beta
      integer n,partonid,beamid
      real*8 isr_tilde_racoon
      real*8 res
      data res/0d0/

      PI=4.D0*DATAN(1.D0)
      
      beta = alphaem/PI * (dlog(Qsquare/me/me)-1)

c     electron beam
      if (beamid .eq. 13) then
c     other partons are zero
          if (partonid .ne. 13) then
             res = 0d0
          else
             if (n .eq. 1) then
                res = isr_tilde_racoon(x,beta)
             else
                res = 0d0
             endif
          endif
      else if (beamid .eq. -13) then
          if (partonid .ne. -13) then
              res = 0d0
          else
             if (n .eq. 1) then
                res = isr_tilde_racoon(x,beta)
             else
                res = 0d0
             endif
          endif
      endif
      eepdf_tilde = res
      end
      
c     https://arxiv.org/pdf/hep-ph/0302198.pdf, eq.(2.44)
c     note that beta_e in eq.(2.45) is twice our beta
c     so eq.(2.44) needs to be corrected by some factor of 2
      real*8 function isr_tilde_racoon(x,beta)
      implicit none
      real*8 x,beta
      real*8 res
      real*8 PI
      real*8 gE
      real*8 logx, logomx
      real*8 dlgam,DDILOG
      external dlgam,DDILOG
      PI=4.D0*DATAN(1.D0)
      gE=0.5772156649d0
      res=0d0
      if (x .lt. 0.9999999d0) then
         logx=dlog(x)
         logomx=dlog(1d0-x)
c     ----------------------------
c     order alpha
         res=-beta*(1d0+x)/2d0
c     order alpha^2
         res=res-(beta**2)/8d0*(
     c        (1d0+3d0*x*x)/(1d0-x)*logx
     c        +4d0*(1d0+x)*logomx+5d0+x)
c     order alpha^3
         res=res-(beta**3)/128d0*(
     c        (1d0+x)*(6d0*DDILOG(x)+12d0*(logomx**2)-3d0*PI**2)
     c        +(1d0/(1d0-x)) * (
     c        (3d0/2d0)*(1d0+8d0*x+3d0*x**2)*logx
     c        +6d0*(x+5d0)*(1d0-x)*logomx
     c        +12d0*(1d0+x**2)*logx*logomx
     c        -1d0/2d0*(1d0+7d0*x**2)*logx**2
     c        +1d0/4d0*(39d0-24d0*x-15d0*x**2)))
c     ----------------------------
         res=res*(1d0-x)**(1-beta)
      endif
      res=res+exp(beta*(-gE+3d0/4d0))/exp(dlgam(1d0+beta))*beta
      isr_tilde_racoon=res
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      function dlgam(xx)
c real logarithm of gamma function
      implicit real * 8 (a-h,o-z)
      real * 8 cof(6),stp,gt,g,cst
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     # -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data gt,g/5.5d0,5.0d0/
      data cst/4.081061466d0/
      x = xx - 1
      xpgt = x + gt
      xmp5  = x + .5d0
      s0 = 1
      do 1 j=1,6
        x = x + 1
        tmp = cof(j)/x
        s0  = s0 + tmp
  1     continue
      r10 = log(s0)
      dlgam = xmp5*(log(xpgt)-1) + r10 - cst
      return
      end      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION DDILOG(X)
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(0:19)
      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)
      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/
      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DDILOG=H
      RETURN
      END



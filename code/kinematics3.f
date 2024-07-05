c Given vegas x's, the partonic CM energy squared (shat), the mass of 
c the two massive particles (xm), and the corresponding threshold 
c (thresh=4*xm**2/shat), generates the phase-space variables (y1, y2,
c omy1=1-y1, omy2=1-y2, xi1, xi2, ph1, ph2, phi, cth) that correspond
c to one of the following configurations:
c   icoll=1  ->  generic y1 and y2
c   icoll=2  ->  y1=1, generic y2 
c   icoll=3  ->  generic y1, y2=1 
c   icoll=4  ->  y1=1, y2=1       
c In addition to generating the phase space variables, the routine returns
c the corresponding jacobians of the transformation x -> phase space 
c variables (that includes the GeV -> pb conversion factor), relevant to 
c the four body (jac4), three body (jac3a and jac3b), and two body (jac2) 
c phase spaces, TIMES the measure of the appropriate z variable(s) in the
c cases of the reduced kinematic configurations, namely: for icoll=2, 
c dz1 <-> jac3b; for icoll=3, dz2 <-> jac3a; and for icoll=4, dz1  <-> jac3b, 
c dz2 <-> jac3a, and dz1dz2 <-> jac2
c
c In addition, the variables xi1 and xi2 are also generated independently
c of those above for the cases of the three- and two-body phase spaces,
c and stored in the common block /cexternal1/. These will be used by
c generate_momenta3 to set up alternative three- and two-body kinematic
c configurations. There are two ways in which these alternative xi1 and xi2
c variables can be generated, controlled by nodz1dz2. When nodz1dz2=.false.
c one uses the parametrisation of sect.3 of muvbf.tex (rather than that
c of the appendix), but otherwise they behave the same as in the four-body
c case. When nodz1dz2=.true., they are set to values provided externally
c by the variables in the common block /cexternal3/ (which MUST respect
c the kinematics bounds -- no check is performed here). In such a case,
c the measure(s) associated with the reduced CM energies are REMOVED
c from the computations of the weights (dz2 for 3A, dz1 for 3B,
c and dz1dz2 for 2)'
c 
      subroutine generate_kinematics3(x,shat,thresh,icoll,isoft,
     #                                y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,
     #                                phi,cth,jac4,ejac3a,ejac3b,ejac2)
C     #                                phi,cth,jac4,jac3a,jac3b,jac2)
      implicit none
      double precision x(8)
      double precision shat,xm,thresh,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,
     # phi,cth,jac4,jac3a,jac3b,jac2,tmp,omega,tvar,lam,lamfun
      integer icoll,isoft
      double precision exi1b,exi2b,exi1a,exi2a,exi12,exi22,
     # ejac3b,ejac3a,ejac2
      common/cexternal1/exi1b,exi2b,exi1a,exi2a,exi12,exi22!,
C     # ejac3b,ejac3a,ejac2
      double precision eez2a,eez1b,eez12,eez22
      common/cexternal3/eez2a,eez1b,eez12,eez22
      double precision pi
      parameter (pi=3.14159265359d0)
      double precision conv
      parameter (conv=389379.66d0*1000)  !conv to picobarns
      double precision tiny
      parameter (tiny=1.d-4)
      logical nodz1dz2
      parameter (nodz1dz2=.false.)
      double precision mfin
      common /to_mfin/mfin
      xm = mfin
c
      if(abs(1-4*xm**2/(shat*thresh)).gt.1.d-10)then
        write(6,*)'Inconsistency in mass assignments:',xm,thresh
        stop
      endif
      jac3a = ejac3a
      jac3b = ejac3b
      jac2 = ejac2
      cth = x(1)*2d0-1d0
      phi = x(2)*2d0*pi
      jac4 = jac4*4d0*pi
      jac3a = jac3a*4d0*pi
      jac3b = jac3b*4d0*pi
      jac2 = jac2*4d0*pi
c
      if (icoll.eq.2.or.icoll.eq.4) then
        y1 = 1d0
        omy1 = 0d0
      else 
        y1 = -2d0*x(3)**2+1d0
        omy1 = 2d0*x(3)**2
      endif
      if (icoll.eq.3.or.icoll.eq.4) then
        y2 = 1d0
        omy2 = 0d0
      else
        y2 = -2d0*x(4)**2+1d0
        omy2 = 2d0*x(4)**2
      endif
      if(icoll.eq.3.or.icoll.eq.4)jac3a = jac3a*4d0*x(3)
      if(icoll.eq.2.or.icoll.eq.4)jac3b = jac3b*4d0*x(4)
      jac4 = jac4*4d0*x(3)
      jac4 = jac4*4d0*x(4)
c
      ph1 = x(5)*2d0*pi
      jac4 = jac4*2d0*pi
      ph2 = x(6)*2d0*pi
      jac4 = jac4*2d0*pi
      if (icoll.eq.2) then
        jac3b = jac3b*2d0*pi
        ph1 = 0.d0
      elseif (icoll.eq.3) then
        jac3a = jac3a*2d0*pi
        ph2 = 0.d0
      elseif (icoll.eq.4) then
        jac3a = jac3a*2d0*pi
        jac3b = jac3b*2d0*pi
        ph1 = 0.d0
        ph2 = 0.d0
      endif
c
      if(.not.nodz1dz2)then
        exi1b = x(7)*(1-thresh)
        exi2b = x(8)*2*(1-thresh-exi1b)/(2-(1+y2)*exi1b)
        ejac3b = jac3b*(1-thresh)*2*(1-thresh-exi1b)/(2-(1+y2)*exi1b)
        exi2a = x(7)*(1-thresh)
        exi1a = x(8)*2*(1-thresh-exi2a)/(2-(1+y1)*exi2a)
        ejac3a = jac3a*(1-thresh)*2*(1-thresh-exi2a)/(2-(1+y1)*exi2a)
        exi12 = x(7)*(1-thresh)
        exi22 = x(8)*(1-thresh-exi12)/(1-exi12)
        ejac2 = jac2*(1-thresh)*(1-thresh-exi12)/(1-exi12)
      else
        exi1b = eez1b
        exi2b = x(8)*2*(1-thresh-exi1b)/(2-(1+y2)*exi1b)
        ejac3b = jac3b*2*(1-thresh-exi1b)/(2-(1+y2)*exi1b)
        exi2a = eez2a
        exi1a = x(8)*2*(1-thresh-exi2a)/(2-(1+y1)*exi2a)
        ejac3a = jac3a*(1-thresh)*2*(1-thresh-exi2a)/(2-(1+y1)*exi2a)
        exi12 = eez12
        exi22 = eez22
        ejac2 = jac2
      endif
c
      omega=sqrt(1-y1**2)*sqrt(1-y2**2)*cos(ph1-ph2)-y1*y2
      tvar=x(7)
      tmp=(1-thresh)*(1-omega)*(1-tvar)*tvar
      if(tmp.le.tiny)then
        lamfun=(1-thresh)*(1+0.5d0*tmp+0.5d0*tmp**2+5/8.d0*tmp**3+
     #                     7/8.d0*tmp**4+21/16.d0*tmp**5)
      else
        lamfun=(1-thresh)*(1-sqrt(1-2*tmp))/tmp
      endif
      lam=x(8)*lamfun
      xi1=lam*tvar
      xi2=lam*(1-tvar)
      jac4 = jac4*lamfun*lam
      jac4 = jac4*conv
      if(icoll.eq.2)then
        jac3b = jac3b*lamfun*lam
      elseif(icoll.eq.3)then
        jac3a = jac3a*lamfun*lam
      elseif(icoll.eq.4)then
        jac3a = jac3a*lamfun*lam
        jac3b = jac3b*lamfun*lam
      endif
      jac3a = jac3a*conv
      jac3b = jac3b*conv
      jac2 = jac2*lamfun*lam
      jac2 = jac2*conv
c
      ejac3b = ejac3b*conv
      ejac3a = ejac3a*conv 
      ejac2 = ejac2*conv
c
      return
      end
      


c Given the partonic CM energy squared (shat), the mass of the two massive 
c particles (xm), and the FKS-type phase-space variables (y1...phi),
c generates the corresponding four momenta and phase-space weights.
c The conventions for the momenta are as follows:
c   i=0,1,2,3: i^th component, i=0 -> energy, i=3 -> longitudinal
c   j=1: incoming muon from the left
c   j=2: incoming muon from the right
c   j=3,4: massive particles
c   j=5: outgoing muon, possibly collinear to 1 
c   j=6: outgoing muon, possibly collinear to 2
c
c The momenta and the measures are relevant to the following multiplicities:
c   four body  <-> xmom4(i,j), meas4;
c   three body <-> xmom3a(i,j), meas3a; and xmom3b(i,j), meas3b
c   two body   <-> xmom2(i,j), meas2
c All of these quantities are computed for the kinematic configuration
c determined by whether y1 and y2 are equal or not equal to one.
c In particular, it is sensible to consider:
c   y1#1, y2#1 (icoll=1) -> four body (event)
c   y1=1, y2#1 (icoll=2) -> four body (cnt), three body 3B (event)
c   y1#1, y2=1 (icoll=3) -> four body (cnt), three body 3A (event)
c   y1=1, y2=1 (icoll=4) -> four body (cnt), three body 3A and 3B (cnt),
c                           two body (event).
c Regardless of which configuration is chosen, the momenta xmom4, xmom3a,
c xmom3b, and xmom2 are always given in the p1+p2 rest frame.
c In addition to the momenta, for the reduced configurations the routine
c also returns the corresponding cm energies squared (s3a,s3b,s2) -- this
c information is redundant, in that it could be computed from the momenta.
c
c Still in the p1+p2 rest frame, the momenta of the three- and two-body
c configurations are also given in the common block /cexternal2/ as they 
c emerge from setting xi1 and xi2 independently of what is done in the 
c four-body phase space; such a setting is handled by generate_kinematics3. 
c To inhibit this option, set mexternal=.false..
c
c In the common block /crsyst/ the routine also stores the momenta
c directly relevant to the collinear configuration selected: at variance
c with the former ones, these momenta are always in their rest frame.

c A check on momentum conservation is performed if mcheck=.true.
c The quantities relevant to the reduced system are computed if mboost=.true.
c
      subroutine generate_momenta3(shat,xm,y1,y2,xi1,xi2,ph1,ph2,
     #                             cth,phi,icoll,meas4,emeas3a,emeas3b,meas2,
     #                             es3a,es3b,es2,xmom4,exmom3a,exmom3b,
     #                             exmom2)
      implicit none
      double precision shat,m,y1,y2,xi1,xi2,ph1,ph2,cth,phi,meas4,
     # meas3a,meas3b,meas2,s3a,s3b,s2
      integer icoll
      double precision xmom4(0:3,1:6),xmom3a(0:3,1:6),xmom3b(0:3,1:6),
     # xmom2(0:3,1:6)
      double precision xm,s,sqs,tmp,prefact,q2,q2an,chybst,chybstmo,
     # shybst,weight,omega,xin(0:3),xout(0:3),vrec(1:3),rec(0:3),
     # xmomb(0:3,3:4),bvec(1:3)
      integer i,j,iseedsave,icount,ineg,idisc(2)
      double precision measred,sred,ey,ymom(0:3,1:6)
      common/crsyst/measred,sred,ey,ymom
      common/cissed/iseedsave
      double precision exi1b,exi2b,exi1a,exi2a,exi12,exi22
      common/cexternal1/exi1b,exi2b,exi1a,exi2a,exi12,exi22
      double precision emeas3a,emeas3b,emeas2,es3a,es3b,es2,
     # exmom3a(0:3,1:6),exmom3b(0:3,1:6),exmom2(0:3,1:6)
C      common/cexternal2/emeas3a,emeas3b,emeas2,es3a,es3b,es2,
C     # exmom3a,exmom3b,exmom2
      logical iscm,mcheck,mboost,mexternal
      parameter (mcheck=.true.)
      parameter (mboost=.true.)
      parameter (mexternal=.true.)
      double precision pi,vtiny
      parameter (pi=3.14159265359d0)
      parameter (vtiny=1.d-16)
c
      s=shat
      sqs=sqrt(s)
      prefact=s/(8*(2*pi)**3)
      meas4=meas4*prefact**2*xi1*xi2
      if(icoll.eq.1)then
        meas3a=0.d0
        meas3b=0.d0
        meas2=0.d0
        emeas3a=0.d0
        emeas3b=0.d0
        emeas2=0.d0
      elseif(icoll.eq.2)then
        meas3a=0.d0
        meas3b=meas3b*prefact*xi2
        meas2=0.d0
        emeas3a=0.d0
        emeas3b=emeas3b*prefact*exi2b
        emeas2=0.d0
      elseif(icoll.eq.3)then
        meas3a=meas3a*prefact*xi1
        meas3b=0.d0
        meas2=0.d0
        emeas3a=emeas3a*prefact*exi1a
        emeas3b=0.d0
        emeas2=0.d0
      elseif(icoll.eq.4)then
        meas3a=meas3a*prefact*xi1
        meas3b=meas3b*prefact*xi2
        meas2=meas2
        emeas3a=meas3a*prefact*exi1a
        emeas3b=meas3b*prefact*exi2b
        emeas2=emeas2
      endif
      xmom4(0,1)=sqs/2
      xmom4(1,1)=0.d0
      xmom4(2,1)=0.d0
      xmom4(3,1)=sqs/2
      xmom4(0,2)=sqs/2
      xmom4(1,2)=0.d0
      xmom4(2,2)=0.d0
      xmom4(3,2)=-sqs/2
      xmom4(0,5)=sqs/2*xi1
      xmom4(1,5)=xmom4(0,5)*sqrt(1-y1**2)*cos(ph1)
      xmom4(2,5)=xmom4(0,5)*sqrt(1-y1**2)*sin(ph1)
      xmom4(3,5)=xmom4(0,5)*y1
      xmom4(0,6)=sqs/2*xi2
      xmom4(1,6)=xmom4(0,6)*sqrt(1-y2**2)*cos(ph2)
      xmom4(2,6)=xmom4(0,6)*sqrt(1-y2**2)*sin(ph2)
      xmom4(3,6)=-xmom4(0,6)*y2
      tmp=0.d0
      do i=0,3
        rec(i)=xmom4(i,1)+xmom4(i,2)-(xmom4(i,5)+xmom4(i,6))
        rec(i)=-rec(i)
        if(i.gt.0)tmp=tmp+rec(i)**2
      enddo
      q2=rec(0)**2-tmp
      omega=sqrt(1-y1**2)*sqrt(1-y2**2)*cos(ph1-ph2)-y1*y2
      q2an=s+s/2.d0*xi1*xi2*(1-omega)-s*(xi1+xi2)
      if(abs(1-q2/q2an).gt.1.d-8)then
        write(6,*)'Error 1 in generate_momenta3',q2,q2an
        stop
      endif
      tmp=sqrt(tmp)
      do i=1,3
        vrec(i)=rec(i)/tmp
      enddo
      ey=sqrt( (rec(0)-tmp)/(rec(0)+tmp) )
      chybst=1/2.d0*(ey+1/ey)
      shybst=1/2.d0*(ey-1/ey)
      chybstmo=chybst-1
      call phsp2b(q2,xm,cth,phi,weight,xmomb)
      meas4=meas4*weight
      meas3a=meas3a*weight
      meas3b=meas3b*weight
      meas2=meas2*weight
      do j=3,4
        do i=0,3
          xin(i)=xmomb(i,j)
        enddo
        call boostwdir2(chybst,shybst,chybstmo,vrec,xin,xout)
        do i=0,3
          xmom4(i,j)=xout(i)
        enddo
      enddo
c
      if(icoll.eq.2)then
        do i=0,3
          xmom3b(i,1)=xmom4(i,1)-xmom4(i,5)
          xmom3b(i,2)=xmom4(i,2)
          xmom3b(i,3)=xmom4(i,3)
          xmom3b(i,4)=xmom4(i,4)
          xmom3b(i,5)=xmom4(i,6)
        enddo
        s3b=shat*(1-xi1)
      elseif(icoll.eq.3)then
        do i=0,3
          xmom3a(i,2)=xmom4(i,2)-xmom4(i,6)
          xmom3a(i,1)=xmom4(i,1)
          xmom3a(i,3)=xmom4(i,3)
          xmom3a(i,4)=xmom4(i,4)
          xmom3a(i,5)=xmom4(i,5)
        enddo
        s3a=shat*(1-xi2)
      elseif(icoll.eq.4)then
        do i=0,3
          xmom3a(i,2)=xmom4(i,2)-xmom4(i,6)
          xmom3a(i,1)=xmom4(i,1)
          xmom3a(i,3)=xmom4(i,3)
          xmom3a(i,4)=xmom4(i,4)
          xmom3a(i,5)=xmom4(i,5)
c
          xmom3b(i,1)=xmom4(i,1)-xmom4(i,5)
          xmom3b(i,2)=xmom4(i,2)
          xmom3b(i,3)=xmom4(i,3)
          xmom3b(i,4)=xmom4(i,4)
          xmom3b(i,5)=xmom4(i,6)
c
          xmom2(i,1)=xmom4(i,1)-xmom4(i,5)
          xmom2(i,2)=xmom4(i,2)-xmom4(i,6)
          xmom2(i,3)=xmom4(i,3)
          xmom2(i,4)=xmom4(i,4)
        enddo
        s2=shat*(1-xi1)*(1-xi2)
        s3a=shat*(1-xi2)
        s3b=shat*(1-xi1)
      endif
c
      if(mcheck)then
        iscm=.true.
        call checkmom(iseedsave,6,0,shat,xmom4,iscm)
        if((.not.iscm))then
          write(6,*)"CM error 1",iseedsave
          stop
        endif
        if(icoll.eq.3.or.icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,5,1,s3a,xmom3a,iscm)
        endif
        if(icoll.eq.2.or.icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,5,2,s3b,xmom3b,iscm)
        endif
        if(icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,4,3,s2,xmom2,iscm)
        endif
      endif
c
      if(.not.mexternal)goto 200
      es3a=shat*(1-exi2a)
      es3b=shat*(1-exi1b)
      es2=shat*(1-exi12)*(1-exi22)
      if(icoll.eq.2.or.icoll.eq.4)then
        exmom3b(0,1)=sqs*(1-exi1b)/2
        exmom3b(1,1)=0.d0
        exmom3b(2,1)=0.d0
        exmom3b(3,1)=sqs*(1-exi1b)/2
        exmom3b(0,2)=sqs/2
        exmom3b(1,2)=0.d0
        exmom3b(2,2)=0.d0
        exmom3b(3,2)=-sqs/2
        exmom3b(0,5)=sqs/2*exi2b
        exmom3b(1,5)=exmom3b(0,5)*sqrt(1-y2**2)*cos(ph2)
        exmom3b(2,5)=exmom3b(0,5)*sqrt(1-y2**2)*sin(ph2)
        exmom3b(3,5)=-exmom3b(0,5)*y2
      endif
      if(icoll.eq.3.or.icoll.eq.4)then
        exmom3a(0,1)=sqs/2
        exmom3a(1,1)=0.d0
        exmom3a(2,1)=0.d0
        exmom3a(3,1)=sqs/2
        exmom3a(0,2)=sqs*(1-exi2a)/2
        exmom3a(1,2)=0.d0
        exmom3a(2,2)=0.d0
        exmom3a(3,2)=-sqs*(1-exi2a)/2
        exmom3a(0,5)=sqs/2*exi1a
        exmom3a(1,5)=exmom3a(0,5)*sqrt(1-y1**2)*cos(ph1)
        exmom3a(2,5)=exmom3a(0,5)*sqrt(1-y1**2)*sin(ph1)
        exmom3a(3,5)=exmom3a(0,5)*y1
      endif
      if(icoll.eq.4)then
        exmom2(0,1)=sqs*(1-exi12)/2
        exmom2(1,1)=0.d0
        exmom2(2,1)=0.d0
        exmom2(3,1)=sqs*(1-exi12)/2
        exmom2(0,2)=sqs*(1-exi22)/2
        exmom2(1,2)=0.d0
        exmom2(2,2)=0.d0
        exmom2(3,2)=-sqs*(1-exi22)/2
      endif
c
      if(icoll.eq.2.or.icoll.eq.4)then
        tmp=0.d0
        do i=0,3
          rec(i)=exmom3b(i,1)+exmom3b(i,2)-exmom3b(i,5)
          rec(i)=-rec(i)
          if(i.gt.0)tmp=tmp+rec(i)**2
        enddo
        q2=rec(0)**2-tmp
        tmp=sqrt(tmp)
        do i=1,3
          vrec(i)=rec(i)/tmp
        enddo
        ey=sqrt( (rec(0)-tmp)/(rec(0)+tmp) )
        chybst=1/2.d0*(ey+1/ey)
        shybst=1/2.d0*(ey-1/ey)
        chybstmo=chybst-1
        call phsp2b(q2,xm,cth,phi,weight,xmomb)
        emeas3b=emeas3b*weight
        do j=3,4
          do i=0,3
            xin(i)=xmomb(i,j)
          enddo
          call boostwdir2(chybst,shybst,chybstmo,vrec,xin,xout)
          do i=0,3
            exmom3b(i,j)=xout(i)
          enddo
        enddo
      endif
c
      if(icoll.eq.3.or.icoll.eq.4)then
        tmp=0.d0
        do i=0,3
          rec(i)=exmom3a(i,1)+exmom3a(i,2)-exmom3a(i,5)
          rec(i)=-rec(i)
          if(i.gt.0)tmp=tmp+rec(i)**2
        enddo
        q2=rec(0)**2-tmp
        tmp=sqrt(tmp)
        do i=1,3
          vrec(i)=rec(i)/tmp
        enddo
        ey=sqrt( (rec(0)-tmp)/(rec(0)+tmp) )
        chybst=1/2.d0*(ey+1/ey)
        shybst=1/2.d0*(ey-1/ey)
        chybstmo=chybst-1
        call phsp2b(q2,xm,cth,phi,weight,xmomb)
        emeas3a=emeas3a*weight
        do j=3,4
          do i=0,3
            xin(i)=xmomb(i,j)
          enddo
          call boostwdir2(chybst,shybst,chybstmo,vrec,xin,xout)
          do i=0,3
            exmom3a(i,j)=xout(i)
          enddo
        enddo
      endif
c
      if(icoll.eq.4)then
        tmp=0.d0
        do i=0,3
          rec(i)=exmom2(i,1)+exmom2(i,2)
          rec(i)=-rec(i)
          if(i.gt.0)tmp=tmp+rec(i)**2
        enddo
        q2=rec(0)**2-tmp
        tmp=sqrt(tmp)
        do i=1,3
          vrec(i)=rec(i)/tmp
        enddo
        ey=sqrt( (rec(0)-tmp)/(rec(0)+tmp) )
        chybst=1/2.d0*(ey+1/ey)
        shybst=1/2.d0*(ey-1/ey)
        chybstmo=chybst-1
        call phsp2b(q2,xm,cth,phi,weight,xmomb)
        emeas2=emeas2*weight
        do j=3,4
          do i=0,3
            xin(i)=xmomb(i,j)
          enddo
          call boostwdir2(chybst,shybst,chybstmo,vrec,xin,xout)
          do i=0,3
            exmom2(i,j)=xout(i)
          enddo
        enddo
      endif
c
      if(mcheck)then
        if(icoll.eq.3.or.icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,5,11,es3a,exmom3a,iscm)
        endif
        if(icoll.eq.2.or.icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,5,12,es3b,exmom3b,iscm)
        endif
        if(icoll.eq.4)then
          iscm=.false.
          call checkmom(iseedsave,4,13,es2,exmom2,iscm)
        endif
      endif
c
 200  continue
      ey=-1.d16
      if(.not.mboost)goto 100
      bvec(1)=0
      bvec(2)=0
      bvec(3)=1
      if(icoll.eq.2)then
        ey=sqrt(1-xi1)
        idisc(1)=5
        idisc(2)=0
      elseif(icoll.eq.3)then
        ey=1/sqrt(1-xi2)
        idisc(1)=6
        idisc(2)=0
      elseif(icoll.eq.4)then
        ey=sqrt((1-xi1)/(1-xi2))
        idisc(1)=5
        idisc(2)=6
      else
        goto 100
      endif
      chybst=1/2.d0*(ey+1/ey)
      shybst=1/2.d0*(ey-1/ey)
      chybstmo=chybst-1
c
      icount=0
      ineg=1
      do j=1,6
        if(j.ne.idisc(ineg))then
          icount=icount+1
          do i=0,3
            ymom(i,icount)=xmom4(i,j)
          enddo
        else
          ineg=ineg+1
        endif
      enddo
      if(icoll.eq.2)then
        do i=0,3
          ymom(i,1)=ymom(i,1)-xmom4(i,5)
        enddo
        sred=shat*(1-xi1)
      elseif(icoll.eq.3)then
        do i=0,3
          ymom(i,2)=ymom(i,2)-xmom4(i,6)
        enddo
        sred=shat*(1-xi2)
      elseif(icoll.eq.4)then
        do i=0,3
          ymom(i,1)=ymom(i,1)-xmom4(i,5)
          ymom(i,2)=ymom(i,2)-xmom4(i,6)
        enddo
        sred=shat*(1-xi1)*(1-xi2)
      endif
      do j=1,icount
        do i=0,3
          xin(i)=ymom(i,j)
        enddo
        call boostwdir2(chybst,shybst,chybstmo,bvec,xin,xout)
        do i=0,3
          ymom(i,j)=xout(i)
        enddo
      enddo
      if(mcheck)then
        iscm=.true.
        call checkmom(iseedsave,icount,4,sred,ymom,iscm)
        if((.not.iscm))then
          write(6,*)"CM error 2",iseedsave
          stop
        endif
      endif
c
 100  continue
      measred = 0.d0
      if(icoll.eq.2)then
        measred=meas3b
      elseif(icoll.eq.3)then
        measred=meas3a
      elseif(icoll.eq.4)then
        measred=meas2
      endif
c
      return
      end



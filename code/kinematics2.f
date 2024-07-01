
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
      subroutine generate_kinematics2(x,shat,thresh,icoll,isoft,
     #                                y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,
     #                                phi,cth,jac4,jac3a,jac3b,jac2)
      implicit none
      double precision x(8)
      double precision shat,xm,thresh,y1,y2,omy1,omy2,xi1,xi2,ph1,ph2,
     # phi,cth,jac4,jac3a,jac3b,jac2,tmp,omega,tvar,lam,lamfun
      integer icoll,isoft
      double precision pi
      parameter (pi=3.14159265359d0)
      double precision conv
      parameter (conv=389379.66d0*1000)  !conv to picobarns
      double precision tiny
      parameter (tiny=1.d-4)
      double precision mfin
      common /to_mfin/mfin
      xm = mfin
c
      if(abs(1-4*xm**2/(shat*thresh)).gt.1.d-10)then
        write(6,*)'Inconsistency in mass assignments:',xm,thresh
        stop
      endif
      jac4 = 1.d0
      jac3a = 1.d0
      jac3b = 1.d0
      jac2 = 1.d0
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
      omega=sqrt(1-y1**2)*sqrt(1-y2**2)*cos(ph1-ph2)-y1*y2
      tvar=x(7)
      tmp=(1-thresh)*(1-omega)*(1-tvar)*tvar
      if(tmp.le.tiny)then
        lamfun=(1-thresh)*(1+0.5*tmp+0.5*tmp**2+5/8.d0*tmp**3+
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
      return
      end


      subroutine phsp2b(s,xm,cth,phi,weight,xmom)
      implicit none
      double precision s,xm,cth,phi,weight,xmom(0:3,3:4)
      double precision sqs,sth,pi,beta
      parameter (pi=3.14159265359d0)
c
      beta=sqrt(1-4*xm**2/s)
      weight=beta/(8*(2*pi)**2)
      sqs=sqrt(s)
      sth=sqrt(1-cth**2)
      xmom(0,3)=sqs/2.d0
      xmom(1,3)=sqs/2.d0*beta*sth*cos(phi)
      xmom(2,3)=sqs/2.d0*beta*sth*sin(phi)
      xmom(3,3)=sqs/2.d0*beta*cth
      xmom(0,4)=xmom(0,3)
      xmom(1,4)=-xmom(1,3)
      xmom(2,4)=-xmom(2,3)
      xmom(3,4)=-xmom(3,3)
      return
      end



      subroutine generate_momenta2(shat,xm,y1,y2,xi1,xi2,ph1,ph2,
     #                             cth,phi,icoll,meas4,meas3a,meas3b,meas2,
     #                             s3a,s3b,s2,xmom4,xmom3a,xmom3b,
     #                             xmom2)
      implicit none
      double precision shat,m,y1,y2,xi1,xi2,ph1,ph2,cth,phi,meas4,
     # meas3a,meas3b,meas2,s3a,s3b,s2
      double precision xmom4(0:3,1:6),xmom3a(0:3,1:6),xmom3b(0:3,1:6),
     # xmom2(0:3,1:6)
      double precision xm,s,sqs,tmp,prefact,q2,q2an,chybst,chybstmo,
     # shybst,weight,omega,xin(0:3),xout(0:3),vrec(1:3),rec(0:3),
     # xmomb(0:3,3:4),bvec(1:3)
      integer i,j,iseedsave,icoll,icount,ineg,idisc(2)
      double precision measred,sred,ey,ymom(0:3,1:6)
      common/crsyst/measred,sred,ey,ymom
      common/cissed/iseedsave
      logical iscm,mcheck,mboost
      parameter (mcheck=.true.)
      parameter (mboost=.true.)
      double precision pi,vtiny
      parameter (pi=3.14159265359d0)
      parameter (vtiny=1.d-16)
c
      xmom4(:,:) = 0d0
      xmom3a(:,:) = 0d0
      xmom3b(:,:) = 0d0
      xmom2(:,:) = 0d0
      s=shat
      sqs=sqrt(s)
      prefact=s/(8*(2*pi)**3)
      meas4=meas4*prefact**2*xi1*xi2
      if(icoll.eq.1)then
        meas3a=0.d0
        meas3b=0.d0
        meas2=0.d0
      elseif(icoll.eq.2)then
        meas3a=0.d0
        meas3b=meas3b*prefact*xi2
        meas2=0.d0
      elseif(icoll.eq.3)then
        meas3a=meas3a*prefact*xi1
        meas3b=0.d0
        meas2=0.d0
      elseif(icoll.eq.4)then
        meas3a=meas3a*prefact*xi1
        meas3b=meas3b*prefact*xi2
        meas2=meas2
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
        write(6,*)'Error 1 in generate_momenta2',q2,q2an
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



      subroutine checkmom(iseedsave,maxmom,icall,s,xmom,iscm)
      implicit none
      real*8 s,snew,xmom(0:3,1:6),check(0:3),cm(0:3)
      integer i,j,maxmom,iseedsave,icall
      logical iscm
c
      snew=xmom(0,1)*xmom(0,2)-xmom(1,1)*xmom(1,2)-
     #     xmom(2,1)*xmom(2,2)-xmom(3,1)*xmom(3,2)
      snew=2*snew
      if(abs(s/snew-1).gt.1.d-8)then
        write(6,*)"wrong energy:",s,snew,iseedsave,icall
        stop
      endif
      do i=0,3
        check(i)=0.d0
        cm(i)=0.d0
        do j=1,maxmom
          if(j.le.2)then
            check(i)=check(i)+xmom(i,j)
            cm(i)=cm(i)+xmom(i,j)
          else
            check(i)=check(i)-xmom(i,j)
          endif
        enddo
      enddo
      do i=0,3
        if(abs(check(i)).gt.1.d-8)then
          write(6,*)"momentum is not conserved:",check(i),i,iseedsave
          stop
        endif
      enddo
      iscm=(iscm.and.abs(cm(0)**2/snew-1).le.1.d-8)
      do i=1,3
        iscm=(iscm.and.abs(cm(i)/sqrt(snew)).le.1.d-8)
      enddo
c
      return
      end

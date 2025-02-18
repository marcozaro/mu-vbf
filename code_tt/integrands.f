      double precision function integrand(x, vegas_wgt)
      implicit none
      double precision x(12), vegas_wgt
      include 'coupl.inc'
      double precision  mmin
      common /to_mmin/mmin
      double precision mfin
      common /to_mfin/mfin
      double precision integrand_mumu,integrand_muga,integrand_gamu,integrand_gaga
      external integrand_mumu,integrand_muga,integrand_gamu,integrand_gaga
      logical fill_histos
      common /to_fill_histos/fill_histos

      include 'mfin.inc'
      mmin = 2d0*mfin

      integrand = 0d0

      ! mu-mu in initial state
      integrand = integrand + integrand_mumu(x,vegas_wgt) 
      ! gam-gam in initial state
      integrand = integrand + integrand_gaga(x,vegas_wgt) 
      ! mu-gam in initial state
      integrand = integrand + integrand_muga(x,vegas_wgt) 
      ! gam-mu in initial state
      integrand = integrand + integrand_gamu(x,vegas_wgt) 

      if (fill_histos) call HwU_add_points()

      return
      end



      double precision function integrand_mumu(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE MUON-MUON CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, xbk(2), omxbk(2)
      integer ilum
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, zq(2), jac_pdf_save

      double precision compute_subtracted_me_2, compute_subtracted_me_1b,
     $ compute_subtracted_me_1a, compute_subtracted_me_0, qprime, getscale
      external compute_subtracted_me_2, compute_subtracted_me_1b,
     $ compute_subtracted_me_1a, compute_subtracted_me_0, qprime, getscale
      double precision mu2
      double precision qq, pploglog
      double precision pgamu
      external pgamu

      logical mumu_doublereal
      parameter (mumu_doublereal=.true.)

      integer orders_tag ! 0->LO,1->NLO,2->NNLO
      common/to_orderstag/orders_tag

      double precision delta_used
      common/to_delta_used/delta_used
      double precision tiny_used
      common/to_coll_cutoff/tiny_used
      include 'input.inc'

      integrand_mumu = 0d0
      !
      ! generate the mu mu luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      ilum = 1
      call get_xbk(ilum,x(9:10),scoll,mmin**2,jac_pdf,tau,ycm,xbk,omxbk)

      ! THE DOUBLE-REAL CONTRIBUTION FOR THE MUON PAIR
      if (.not.mumu_doublereal) goto 10

      orders_tag = 2
      delta_used = deltaI
      tiny_used = tinycoll
      integrand_mumu = integrand_mumu + compute_subtracted_me_2(x,vegas_wgt,xbk,omxbk,tau,ycm,jac_pdf,ilum)

 10   continue

      ! THE CONVOLUTION OF M_GAM MU WITH Q'(Z1)
      jac_pdf_save = jac_pdf
      zq(2) = -1d0
      call generate_qp_z(x(11),tau_min/tau,zq(1),jac_pdf)

      orders_tag = 2
      delta_used = deltaIb
      tiny_used = tinycoll*2/(2*zq(1)+tinycoll*(1-zq(1)))
      !lum*qprime(z1,scoll*tau,mu2,deltaI)
      integrand_mumu = integrand_mumu +
     $ compute_subtracted_me_1b(x,vegas_wgt,xbk,omxbk,zq,
     $               tau*zq(1),ycm+0.5*dlog(zq(1)),jac_pdf,ilum)

      ! THE CONVOLUTION OF M MU GAM WITH Q'(Z2)
      jac_pdf = jac_pdf_save
      zq(1) = -1d0
      call generate_qp_z(x(11),tau_min/tau,zq(2),jac_pdf)

      orders_tag = 2
      delta_used = deltaIb
      tiny_used = tinycoll*2/(2*zq(2)+tinycoll*(1-zq(2)))
      integrand_mumu = integrand_mumu +
     $ compute_subtracted_me_1a(x,vegas_wgt,xbk,omxbk,zq,
     $               tau*zq(2),ycm-0.5*dlog(zq(2)),jac_pdf,ilum)

      ! THE CONVOLUTION OF M GAM GAM WITH Q'(Z1) Q'(Z2)
      jac_pdf = jac_pdf_save
      call generate_qp_z(x(11),tau_min/tau,zq(1),jac_pdf)
      call generate_qp_z(x(12),tau_min/tau/zq(1),zq(2),jac_pdf)
      orders_tag = 2
      !! caareful here with the scale that enters in the logs inside q'
      ! it is the only non-trivial place. In Q'(z_i), the com energy
      ! that enters is z_1 * z_2 * tau * shat / z_i

      !qq = qprime(z1,scoll*tau*z2,mu2,deltaIb)*qprime(z2,scoll*tau*z1,mu2,deltaIb)
      !pploglog = Pgamu(z1)*Pgamu(z2)*dlog(z1*deltaIb/deltaI)*dlog(z2*deltaIb/deltaI)
        !lum*(qq-pploglog)

      integrand_mumu = integrand_mumu + 
     $ compute_subtracted_me_0(x,vegas_wgt,xbk,omxbk,zq,
     $  tau*zq(1)*zq(2),ycm+0.5*dlog(zq(1))-0.5*dlog(zq(2)),jac_pdf,ilum)

      return
      end


      double precision function integrand_gaga(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE GAMMA-GAMMA CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, xbk(2), omxbk(2), zq(2)
      integer ilum
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin

      double precision compute_subtracted_me_0
      external compute_subtracted_me_0

      logical gaga_born
      parameter (gaga_born=.true.)

      integer orders_tag ! 0->LO,1->NLO,2->NNLO
      common/to_orderstag/orders_tag
       
      integrand_gaga = 0d0

      ! generate the gamma-gamma luminosity
      jac_pdf = 1d0
      ilum = 4
      call get_xbk(ilum,x(9:10),scoll,mmin**2,jac_pdf,tau,ycm,xbk,omxbk)

      ! THE BORN CONTRIBUTION FOR THE PHOTON PAIR
      if (.not.gaga_born) goto 10
      orders_tag = 0
      zq(:)=-1d0
      integrand_gaga = integrand_gaga + compute_subtracted_me_0(x,vegas_wgt,xbk,omxbk,zq,tau,ycm,jac_pdf,ilum)

 10   continue

      return
      end



      double precision function integrand_muga(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE MUON-GAMMA CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, xbk(2), omxbk(2)
      integer ilum
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, zq(2)

      logical muga_singlereal
      parameter (muga_singlereal=.true.)

      double precision compute_subtracted_me_1a, compute_subtracted_me_0, qprime, getscale
      external compute_subtracted_me_1a, compute_subtracted_me_0, qprime, getscale
      double precision mu2

      integer orders_tag ! 0->LO,1->NLO,2->NNLO
      common/to_orderstag/orders_tag
      include 'input.inc'
      double precision delta_used
      common/to_delta_used/delta_used
      double precision tiny_used
      common/to_coll_cutoff/tiny_used

      integrand_muga = 0d0
      !
      ! generate the mu gam luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      ilum = 3
      call get_xbk(ilum,x(9:10),scoll,mmin**2,jac_pdf,tau,ycm,xbk,omxbk)

      ! THE SINGLE-REAL CONTRIBUTION 
      if (.not.muga_singlereal) goto 10
      orders_tag = 1
      delta_used = deltaI
      tiny_used = tinycoll
      zq(:) = -1d0
      integrand_muga = integrand_muga +
     $ compute_subtracted_me_1a(x,vegas_wgt,xbk,omxbk,zq,
     $               tau,ycm,jac_pdf,ilum)

 10   continue

      ! THE CONVOLUTION OF M_GAM GAM WITH Q'(Z1)
      call generate_qp_z(x(11),tau_min/tau,zq(1),jac_pdf)

      !!mu2 = getscale(scoll, x1bk, x2bk)
      !lum*qprime(z1,scoll*tau,mu2,deltaI)

      orders_tag = 1
      integrand_muga = integrand_muga +
     $ compute_subtracted_me_0(x,vegas_wgt,xbk,omxbk,zq,
     $               tau*zq(1),ycm+0.5*dlog(zq(1)),jac_pdf,ilum)

      return
      end


      double precision function integrand_gamu(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE GAMMA-MUON CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, xbk(2), omxbk(2)
      integer ilum
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, zq(2)

      logical gamu_singlereal
      parameter (gamu_singlereal=.true.)

      double precision compute_subtracted_me_1b, compute_subtracted_me_0, qprime, getscale
      external compute_subtracted_me_1b, compute_subtracted_me_0, qprime, getscale
      double precision mu2

      integer orders_tag ! 0->LO,1->NLO,2->NNLO
      common/to_orderstag/orders_tag
      include 'input.inc'
      double precision delta_used
      common/to_delta_used/delta_used
      double precision tiny_used
      common/to_coll_cutoff/tiny_used

      integrand_gamu = 0d0
      !
      ! generate the gam mu luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      ilum = 2
      call get_xbk(ilum,x(9:10),scoll,mmin**2,jac_pdf,tau,ycm,xbk,omxbk)

      ! THE SINGLE-REAL CONTRIBUTION 
      if (.not.gamu_singlereal) goto 10

      orders_tag = 1
      delta_used = deltaI
      tiny_used = tinycoll
      zq(:) = -1d0
      integrand_gamu = integrand_gamu +
     $ compute_subtracted_me_1b(x,vegas_wgt,xbk,omxbk,zq,
     $               tau,ycm,jac_pdf,ilum)

 10   continue

      ! THE CONVOLUTION OF M_GAM GAM WITH Q'(Z2)
      !We use jac0, since we convolve with the born-like matrix element
      call generate_qp_z(x(11),tau_min/tau,zq(2),jac_pdf)

      !mu2 = getscale(scoll, x1bk, x2bk)
      !lum*qprime(z2,scoll*tau,mu2,deltaI)
      orders_tag = 1
      integrand_gamu = integrand_gamu +
     $ compute_subtracted_me_0(x,vegas_wgt,xbk,omxbk,zq,
     $               tau*zq(2),ycm-0.5*dlog(zq(2)),jac_pdf,ilum)

      return
      end


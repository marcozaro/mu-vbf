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

      integer npoints 
      data npoints/0/
      npoints = npoints+1

      mfin = mdl_mt
      mmin = 2d0*mfin
      write(*,*) 'XX', x, vegas_wgt

      integrand = 0d0

      write(*,*) 'zero', integrand
      ! mu-mu in initial state
      !integrand = integrand + integrand_mumu(x,vegas_wgt) 
      write(*,*) 'mumu', integrand
      ! gam-gam in initial state
      integrand = integrand + integrand_gaga(x,vegas_wgt) 
      write(*,*) 'gaga', integrand
      ! mu-gam in initial state
      integrand = integrand + integrand_muga(x,vegas_wgt) 
      write(*,*) 'muga', integrand
      ! gam-mu in initial state
      integrand = integrand + integrand_gamu(x,vegas_wgt) 
      write(*,*) 'gamu', integrand

      if (fill_histos) call HwU_add_points()
      if (npoints.eq.10) stop

      return
      end



      double precision function integrand_mumu(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE MUON-MUON CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, x1bk, x2bk
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, z1, z2, jac_pdf_save

      double precision compute_subtracted_me_2, compute_subtracted_me_1b,
     $ compute_subtracted_me_1a, compute_subtracted_me_0_qq, qprime, getscale
      external compute_subtracted_me_2, compute_subtracted_me_1b,
     $ compute_subtracted_me_1a, compute_subtracted_me_0_qq, qprime, getscale
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
      include 'input.inc'

      integrand_mumu = 0d0
      !
      ! generate the mu mu luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      call get_lum(1,x(9:10),scoll,mmin**2,jac_pdf,lum,tau,ycm,x1bk,x2bk)
      jac_pdf_save = jac_pdf

      ! THE DOUBLE-REAL CONTRIBUTION FOR THE MUON PAIR
      if (.not.mumu_doublereal) goto 10

      orders_tag = 2
      delta_used = deltaI
      integrand_mumu = integrand_mumu + compute_subtracted_me_2(x,vegas_wgt,lum,tau,ycm,jac_pdf)

 10   continue

      ! THE CONVOLUTION OF M_GAM MU WITH Q'(Z1)
      integrand_mumu = integrand_mumu +
     $ compute_subtracted_me_1b(x,vegas_wgt,lum,
     $               tau,ycm,jac_pdf,1)

      ! THE CONVOLUTION OF M MU GAM WITH Q'(Z2)
      integrand_mumu = integrand_mumu +
     $ compute_subtracted_me_1a(x,vegas_wgt,lum,
     $               tau,ycm,jac_pdf,2)

      ! THE CONVOLUTION OF M GAM GAM WITH Q'(Z1) Q'(Z2)
      integrand_mumu = integrand_mumu + 
     $ compute_subtracted_me_0_qq(x,vegas_wgt,lum,
     $               tau,ycm,jac_pdf)

      return
      end


      double precision function integrand_gaga(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE GAMMA-GAMMA CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, x1bk, x2bk
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
      call get_lum(4,x(9:10),scoll,mmin**2,jac_pdf,lum,tau,ycm,x1bk,x2bk)

      ! THE BORN CONTRIBUTION FOR THE PHOTON PAIR
      if (.not.gaga_born) goto 10
      orders_tag = 0
      integrand_gaga = integrand_gaga + compute_subtracted_me_0(x,vegas_wgt,lum,tau,ycm,jac_pdf)

 10   continue

      return
      end



      double precision function integrand_muga(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE MUON-GAMMA CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, x1bk,x2bk
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, z1, z2

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

      integrand_muga = 0d0
      !
      ! generate the mu gam luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      call get_lum(3,x(9:10),scoll,mmin**2,jac_pdf,lum,tau,ycm,x1bk,x2bk)

      ! THE SINGLE-REAL CONTRIBUTION 
      if (.not.muga_singlereal) goto 10
      orders_tag = 1
      delta_used = deltaI
      integrand_muga = integrand_muga +
     $ compute_subtracted_me_1a(x,vegas_wgt,lum,
     $               tau,ycm,jac_pdf,0)
      write(*,*) 'SR', integrand_muga

 10   continue

      ! THE CONVOLUTION OF M_GAM GAM WITH Q'(Z1)
      call generate_qp_z(x(11),tau_min/tau,z1,jac_pdf)

      mu2 = getscale(scoll, x1bk*x2bk*scoll)

      orders_tag = 1
       write(*,*) 'Z1',z1,x(11), tau_min/tau
      integrand_muga = integrand_muga +
     $ compute_subtracted_me_0(x,vegas_wgt,lum*qprime(z1,scoll*tau,mu2),
     $               tau*z1,ycm+0.5*dlog(z1),jac_pdf)
      write(*,*) 'QQ', integrand_muga

      return
      end


      double precision function integrand_gamu(x,vegas_wgt)
      implicit none
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! THE GAMMA-MUON CONTRIBUTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision x(12), vegas_wgt
      double precision jac_pdf, lum, tau, ycm, x1bk, x2bk
      double precision scoll
      common /to_scoll/scoll
      double precision mmin
      common /to_mmin/mmin
      double precision tau_min, z1, z2

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

      integrand_gamu = 0d0
      !
      ! generate the gam mu luminosity
      jac_pdf = 1d0
      tau_min = mmin**2/scoll
      call get_lum(2,x(9:10),scoll,mmin**2,jac_pdf,lum,tau,ycm,x1bk,x2bk)

      ! THE SINGLE-REAL CONTRIBUTION 
      if (.not.gamu_singlereal) goto 10

      orders_tag = 1
      delta_used = deltaI
      integrand_gamu = integrand_gamu +
     $ compute_subtracted_me_1b(x,vegas_wgt,lum,
     $               tau,ycm,jac_pdf,0)

 10   continue

      ! THE CONVOLUTION OF M_GAM GAM WITH Q'(Z2)
      !We use jac0, since we convolve with the born-like matrix element
      call generate_qp_z(x(11),tau_min/tau,z2,jac_pdf)

      mu2 = getscale(scoll, x1bk*x2bk*scoll)

      orders_tag = 1
      integrand_gamu = integrand_gamu +
     $ compute_subtracted_me_0(x,vegas_wgt,lum*qprime(z2,scoll*tau,mu2),
     $               tau*z2,ycm-0.5*dlog(z2),jac_pdf)

      return
      end


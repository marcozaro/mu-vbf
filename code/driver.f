
      double precision scoll
      common /to_scoll/scoll
      double precision x(12)
      integer i,j
      logical check
      parameter(check=.true.)

      double precision ans
      include 'coupl.inc'

      double precision integrand
      external integrand
      double precision integral,error,prob
      integer nprn
      logical fill_histos
      common /to_fill_histos/fill_histos
      logical sameflav_diags
      common /to_sameflav/sameflav_diags
      double precision tiny
      common/to_coll_cutoff/tiny
      include 'input.inc'

      sameflav_diags = sameflav.gt.0

      tiny = tinycoll

      scoll = ecm**2

      call setpara('Cards/param_card.dat')
      call printout()
      call print_run()
      !MZ leave this for now, to keep the RN sequence
      call fill_vegas_x(x)

      nprn = 0
      ! fill histogram only in the refine phase
      fill_histos = .false.
      call vegas01(12,integrand,0,40000,
     1        10,nprn,integral,error,prob)

      ! for the analysis
      fill_histos = .true.
      call set_error_estimation(1)
      call analysis_begin(1,"central value")

      call vegas01(12,integrand,1,400000,
     1        4,nprn,integral,error,prob)
      call analysis_end(1d0)

      return
      end



      subroutine fill_vegas_x(x)
C     fill the vegas x.
      implicit none
      double precision x(12)
      integer i
      double precision ran2
      external ran2
      do i = 1,12
         x(i) = ran2()
      enddo

      return 
      end


      subroutine print_run()
      implicit none
      include "input.inc"

      include "printout.inc"

      return
      end

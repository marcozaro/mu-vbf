C ************************************************************
C Source for the library implementing a bias function that 
C populates the large pt tale of the leading jet. 
C
C The two options of this subroutine, that can be set in
C the run card are:
C    > (double precision) ptj_bias_target_ptj : target ptj value
C    > (double precision) ptj_bias_enhancement_power : exponent
C
C Schematically, the functional form of the enhancement is
C    bias_wgt = [ptj(evt)/mean_ptj]^enhancement_power
C ************************************************************
C
C The following lines are read by MG5aMC to set what are the 
C relevant parameters for this bias module.
C
C  parameters = {'ptj_bias_target_ptj': 1000.0,
C               'ptj_bias_enhancement_power': 4.0}
C

      subroutine bias_wgt(p, original_weight, bias_weight)
          implicit none
C
C Parameters
C
          include '../../maxparticles.inc'          
          include '../../nexternal.inc'

C
C Arguments
C
          double precision p(0:3,nexternal)
          double precision original_weight, bias_weight
C
C local variables
C
          integer i
          double precision ptj(nexternal)
          double precision max_ptj
c
c local variables defined in the run_card
c
          double precision ptj_bias_target_ptj
          double precision ptj_bias_enhancement_power
C
C Global variables
C
C
C Mandatory common block to be defined in bias modules
C
          double precision stored_bias_weight
          data stored_bias_weight/1.0d0/          
          logical impact_xsec, requires_full_event_info
C         We only want to bias distributions, but not impact the xsec. 
          data impact_xsec/.False./
C         Of course this module does not require the full event
C         information (color, resonances, helicities, etc..)
          data requires_full_event_info/.False./ 
          common/bias/stored_bias_weight,impact_xsec,
     &                requires_full_event_info

          double precision minvtwo

C
C Accessingt the details of the event
C
          logical is_a_j(nexternal),is_a_l(nexternal),
     &            is_a_b(nexternal),is_a_a(nexternal),
     &            is_a_onium(nexternal),is_a_nu(nexternal),
     &            is_heavy(nexternal),do_cuts(nexternal)
          common/to_specisa/is_a_j,is_a_a,is_a_l,is_a_b,is_a_nu,
     &                      is_heavy,is_a_onium,do_cuts

C
C    Setup the value of the parameters from the run_card    
C
      include '../bias.inc'

C --------------------
C BEGIN IMPLEMENTATION
C --------------------
          
c          do i=1,nexternal
c            ptj(i)=-1.0d0
c            if (is_a_j(i)) then
c              ptj(i)=sqrt(p(1,i)**2+p(2,i)**2)
c            endif
c         enddo

c          max_ptj=-1.0d0
c          do i=1,nexternal
c            max_ptj = max(max_ptj,ptj(i))
c          enddo
c          if (max_ptj.lt.0.0d0) then
c            bias_weight = 1.0d0
c            return
c          endif

c          bias_weight = (max_ptj/ptj_bias_target_ptj)
c     &                                      **ptj_bias_enhancement_power

           minvtwo=(p(0,3)+p(0,4))**2 - (p(1,3)+p(1,4))**2 - (p(2,3)+p(2,4))**2 - (p(3,3)+p(3,4))**2

c          bias_weight =max(1d0,(p(1,3)**2+p(2,3)**2)**3/(1000d0)**6)
         
           bias_weight =max(1d0,(minvtwo)**(1.7d0)/(1000d0)**6)

          return

      end subroutine bias_wgt

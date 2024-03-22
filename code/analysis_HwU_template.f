c
c This file contains the default histograms for fixed order runs: it
c only plots the total rate as an example. It can be used as a template
c to make distributions for other observables.
c
c This uses the HwU package and generates histograms in the HwU/GnuPlot
c format. This format is human-readable. After running, the histograms
c can be found in the Events/run_XX/ directory.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_begin(nwgt,weights_info)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called once at the start of each run. Here the
c histograms should be declared. 
c
c Declare the histograms using 'HwU_book'.
c     o) The first argument is an integer that labels the histogram. In
c     the analysis_end and analysis_fill subroutines this label is used
c     to keep track of the histogram. The label should be a number
c     starting at 1 and be increased for each plot.
c     o) The second argument is a string that will apear above the
c     histogram. Do not use brackets "(" or ")" inside this string.
c     o) The third, forth and fifth arguments are the number of bis, the
c     lower edge of the first bin and the upper edge of the last
c     bin.
c     o) When including scale and/or PDF uncertainties, declare a
c     histogram for each weight, and compute the uncertainties from the
c     final set of histograms
c
      implicit none
c When including scale and/or PDF uncertainties the total number of
c weights considered is nwgt
      integer nwgt
c In the weights_info, there is an text string that explains what each
c weight will mean. The size of this array of strings is equal to nwgt.
      character*(*) weights_info(*)
c Initialize the histogramming package (HwU). Pass the number of
c weights and the information on the weights:
      call HwU_inithist(nwgt,weights_info)
c declare (i.e. book) the histograms
      call HwU_book(1,'total rate      ', 5,0.5d0,5.5d0)
      call HwU_book(2,'total rate real ', 5,0.5d0,5.5d0)
      call HwU_book(3,'mtt ',   50,0d0,1000d0)
      call HwU_book(4,'thetatt ', 50,-1.5707963268d0,1.5707963268d0)
      call HwU_book(5,'pttt ', 50,0d0,1000d0)
      call HwU_book(6,'ptmup', 50,0d0,1000d0)
      call HwU_book(7,'ptmum', 50,0d0,1000d0)
      call HwU_book(8,'pttop', 50,0d0,1000d0)
      call HwU_book(9,'ptatop', 50,0d0,1000d0)
      call HwU_book(10,'mmumu ',   50,0d0,1000d0)
      call HwU_book(11,'thetamup', 50,-1.5707963268d0,1.5707963268d0)
      call HwU_book(12,'thetamum', 50,-1.5707963268d0,1.5707963268d0)
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_end(dummy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called once at the end of the run. Here the
c histograms are written to disk. Note that this is done for each
c integration channel separately. There is an external script that will
c read the HwU data files in each of the integration channels and
c combines them by summing all the bins in a final single HwU data file
c to be put in the Events/run_XX directory, together with a gnuplot
c file to convert them to a postscript histogram file.
      implicit none
      double precision dummy
      call HwU_write_file
      return                
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine analysis_fill(p,istatus,ipdg,wgts,ibody)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine is called for each n-body and (n+1)-body configuration
c that passes the generation cuts. Here the histrograms are filled.
      implicit none
c This includes the 'nexternal' parameter that labels the number of
c particles in the (n+1)-body process
c This is an array which is '-1' for initial state and '1' for final
c state particles
      integer istatus(6)
c This is an array with (simplified) PDG codes for the particles. Note
c that channels that are combined (i.e. they have the same matrix
c elements) are given only 1 set of PDG codes. This means, e.g., that
c when using a 5-flavour scheme calculation (massless b quark), no
c b-tagging can be applied.
      integer iPDG(6)
c The array of the momenta and masses of the initial and final state
c particles in the lab frame. The format is "E, px, py, pz, mass", while
c the second dimension loops over the particles in the process. Note
c that these are the (n+1)-body particles; for the n-body there is one
c momenta equal to all zeros (this is not necessarily the last particle
c in the list). If one uses IR-safe obserables only, there should be no
c difficulty in using this.
      double precision p(0:3,6)
c The weight of the current phase-space point is wgts(1). If scale
c and/or PDF uncertainties are included through reweighting, the rest of
c the array contains the list of weights in the same order as described
c by the weigths_info strings in analysis_begin
      double precision wgts(*)
c The ibody variable is:
c     ibody=1 : (n+1)-body contribution
c     ibody=2 : n-body contribution (excluding the Born)
c     ibody=3 : Born contribution
c The histograms need to be filled for all these contribution to get the
c physics NLO results. (Note that the adaptive phase-space integration
c is optimized using the sum of the contributions, therefore plotting
c them separately might lead to larger than expected statistical
c fluctuations).
      integer ibody
c local variable
      double precision var
      double precision p_tt(0:3), p_mm(0:3)
      double precision pt_tt, th_tt, m_tt
      double precision th_mp, th_mm, m_mm
      double precision pt_mp, pt_mm, pt_t, pt_tx
c
c Fill the histograms here using a call to the HwU_fill()
c subroutine. The first argument is the histogram label, the second is
c the numerical value of the variable to plot for the current
c phase-space point and the final argument is the weight of the current
c phase-space point.
      var=1d0
      p_tt(:) = p(:,3) + p(:,4)
      p_mm(:) = p(:,5) + p(:,6)
      pt_tt = dsqrt(p_tt(1)**2 + p_tt(2)**2)
      th_tt = datan(pt_tt / p_tt(3))
      m_tt = p_tt(0)**2 - p_tt(1)**2 - p_tt(2)**2 - p_tt(3)**2
      m_mm = p_mm(0)**2 - p_mm(1)**2 - p_mm(2)**2 - p_mm(3)**2
      m_tt = dsqrt(m_tt)
      m_mm = dsqrt(m_mm)

      pt_t  = dsqrt(p(1,3)**2 + p(2,3)**2)
      pt_tx = dsqrt(p(1,4)**2 + p(2,4)**2)
      pt_mp = dsqrt(p(1,5)**2 + p(2,5)**2)
      pt_mm = dsqrt(p(1,6)**2 + p(2,6)**2)
      th_mp = datan(pt_mp / p(3,5))
      th_mm = datan(pt_mm / p(3,6))
c always fill the total rate
      call HwU_fill(1,var,wgts)
c only fill the total rate for the Born when ibody=3
      if (ibody.eq.3) call HwU_fill(2,var,wgts)

      call HwU_fill(3,m_tt,wgts)
      call HwU_fill(4,th_tt,wgts)
      call HwU_fill(5,pt_tt,wgts)

      call HwU_fill(6,pt_mp,wgts)
      call HwU_fill(7,pt_mm,wgts)
      call HwU_fill(8,pt_t ,wgts)
      call HwU_fill(9,pt_tx,wgts)
      call HwU_fill(10,m_mm,wgts)
      call HwU_fill(11,th_mp,wgts)
      call HwU_fill(12,th_mm,wgts)
      return
      end

      integer max_jet_number
      parameter(max_jet_number=20)
      integer njet,jet(max_jet_number)
      double precision pjet(0:3,max_jet_number)
      double precision ptjet(max_jet_number),etajet(max_jet_number),yjet(max_jet_number)
      common/cfastjeti/njet,jet
      common/cfastjetr/pjet,ptjet,etajet,yjet
      REAL*8 JETALGO,JETRADIUS,PTJ,ETAJ ! jet cuts
      REAL*8 PTL,ETAL,DRLL,DRLL_SF,MLL,MLL_SF ! lepton cuts
      COMMON /TO_NEW_AUTO_CUTS/
     &     JETALGO,JETRADIUS,PTJ,ETAJ,
     &     PTL,ETAL,DRLL,DRLL_SF,MLL,MLL_SF
      LOGICAL GAMMA_IS_J
      COMMON /TO_GAMMA_IS_J/GAMMA_IS_J
      REAL*8 PTGMIN,R0GAMMA,XN,EPSGAMMA,ETAGAMMA  ! photon cuts
      LOGICAL ISOEM                           ! photon cuts
      COMMON /TO_ISOGAMMA_CUTS/PTGMIN,R0GAMMA,XN,EPSGAMMA,ETAGAMMA
      COMMON /TO_ISOGAMMA_EM/ISOEM
      integer maxjetflavor
      COMMON /to_min_max_cuts/maxjetflavor
      DOUBLE PRECISION RPHRECO, ETAPHRECO
      LOGICAL LEPPHRECO, QUARKPHRECO
      COMMON /TO_PHRECO/ RPHRECO, ETAPHRECO, LEPPHRECO, QUARKPHRECO
C
C     SPECIAL CUTS (FOR MERGING)
C
      real*8 xmtc,xqcut
      common /to_specxpt/xmtc,xqcut

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        MDL_YB = 0.000000D+00

        MDL_CKM31 = 0.000000D+00

        MDL_CKM32 = 0.000000D+00

        MDL_MB__EXP__4 = 0.000000D+00

        MDL_MB__EXP__2 = 0.000000D+00

        MDL_CKM11 = 1.000000D+00

        MDL_SCKM33 = 1.000000D+00

        MDL_SCKM11 = 1.000000D+00

        MDL_VEP = 0.000000D+00

        MDL_CKM33 = 1.000000D+00

        MDL_CMSPARAM = 0.000000D+00

        MDL_MC__EXP__2 = 0.000000D+00

        MDL_SCKM33__EXP__2 = 1.000000D+00

        MDL_QL__EXP__2 = 1.000000D+00

        MDL_CONJG__CKM22 = 1.000000D+00

        MDL_SCKM22 = 1.000000D+00

        MDL_CKM22 = 1.000000D+00

        MDL_SCKM22__EXP__2 = 1.000000D+00

        MDL_MD__EXP__2 = 0.000000D+00

        MDL_MU__EXP__2 = 0.000000D+00

        MDL_MM__EXP__2 = 0.000000D+00

        MDL_MTA__EXP__2 = 0.000000D+00

        MDL_QL__EXP__4 = 1.000000D+00

        MDL_CONJG__CKM11 = 1.000000D+00

        MDL_LHV = 1.000000D+00

        MDL_CONJG__CKM31 = 0.000000D+00

        MDL_CONJG__CKM32 = 0.000000D+00

        MDL_CONJG__CKM33 = 1.000000D+00

        MDL_NCOL = 3.000000D+00

        MDL_CA = 3.000000D+00

        MDL_TF = 5.000000D-01

        MDL_CF = (4.000000D+00/3.000000D+00)

        MDL_COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

        MDL_MZ__EXP__2 = MDL_MZ**2

        MDL_MW__EXP__2 = MDL_MW**2

        MDL_SQRT__2 = SQRT(DCMPLX(2.000000D+00))

        MDL_SW2 = 1.000000D+00-MDL_MW__EXP__2/MDL_MZ__EXP__2

        MDL_CW = SQRT(DCMPLX(1.000000D+00-MDL_SW2))

        MDL_SQRT__SW2 = SQRT(DCMPLX(MDL_SW2))

        MDL_SW = MDL_SQRT__SW2

        MDL_MH__EXP__2 = MDL_MH**2

        MDL_I3D = (-1.000000D+00/2.000000D+00)

        MDL_I3U = (1.000000D+00/2.000000D+00)

        MDL_I3L = (-1.000000D+00/2.000000D+00)

        MDL_I3V = (1.000000D+00/2.000000D+00)

        MDL_QD = (-1.000000D+00/3.000000D+00)

        MDL_QU = (2.000000D+00/3.000000D+00)

        MDL_QL = (-1.000000D+00)

        MDL_NCOL__EXP__2 = MDL_NCOL**2

        MDL_SW__EXP__2 = MDL_SW**2

        MDL_MT__EXP__2 = MDL_MT**2

        MDL_CW__EXP__2 = MDL_CW**2

        MDL_MT__EXP__4 = MDL_MT**4

        MDL_SW__EXP__4 = MDL_SW**4

        MDL_MW__EXP__4 = MDL_MW**4

        MDL_MT__EXP__3 = MDL_MT**3

        MDL_MH__EXP__4 = MDL_MH**4

        MDL_MZ__EXP__4 = MDL_MZ**4

        MDL_MT__EXP__6 = MDL_MT**6

        MDL_MW__EXP__6 = MDL_MW**6

        MDL_MT__EXP__8 = MDL_MT**8

        MDL_MW__EXP__8 = MDL_MW**8

        MDL_MT__EXP__10 = MDL_MT**10

        MDL_MW__EXP__10 = MDL_MW**10

        MDL_MZ__EXP__6 = MDL_MZ**6

        MDL_SW__EXP__3 = MDL_SW**3

        MDL_MZ__EXP__3 = MDL_MZ**3

        MDL_CW__EXP__4 = MDL_CW**4

        MDL_CW__EXP__6 = MDL_CW**6

        MDL_MH__EXP__8 = MDL_MH**8

        MDL_MH__EXP__6 = MDL_MH**6

        MDL_CW__EXP__3 = MDL_CW**3

        MDL_MW__EXP__3 = MDL_MW**3

        MDL_QD__EXP__3 = MDL_QD**3

        MDL_QD__EXP__2 = MDL_QD**2

        MDL_I3D__EXP__2 = MDL_I3D**2

        MDL_I3D__EXP__3 = MDL_I3D**3

        MDL_QU__EXP__3 = MDL_QU**3

        MDL_QU__EXP__2 = MDL_QU**2

        MDL_I3U__EXP__2 = MDL_I3U**2

        MDL_I3U__EXP__3 = MDL_I3U**3

        MDL_QL__EXP__3 = MDL_QL**3

        MDL_I3L__EXP__2 = MDL_I3L**2

        MDL_I3L__EXP__3 = MDL_I3L**3

        MDL_I3V__EXP__3 = MDL_I3V**3

        MDL_SW__EXP__6 = MDL_SW**6

        MDL_SW__EXP__8 = MDL_SW**8

        MDL_QU__EXP__4 = MDL_QU**4

        MDL_QD__EXP__4 = MDL_QD**4

        MDL_I3L__EXP__4 = MDL_I3L**4

        MDL_I3U__EXP__4 = MDL_I3U**4

        MDL_I3D__EXP__4 = MDL_I3D**4

        MDL_I3V__EXP__4 = MDL_I3V**4

        MDL_I3V__EXP__2 = MDL_I3V**2

        MDL_AEW = 1.000000D+00/AEWM1

        MDL_GF = -MDL_AEW*MDL_MZ__EXP__2*PI/(MDL_SQRT__2*MDL_MW__EXP__2
     $*(MDL_MW__EXP__2-MDL_MZ__EXP__2))

        MDL_SQRT__AEW = SQRT(DCMPLX(MDL_AEW))

        MDL_EE = 2.000000D+00*MDL_SQRT__AEW*SQRT(DCMPLX(PI))

        MDL_G1 = MDL_EE/MDL_CW

        MDL_GW = MDL_EE/MDL_SW

        MDL_VEV = (2.000000D+00*MDL_MW*MDL_SW)/MDL_EE

        MDL_VEV__EXP__2 = MDL_VEV**2

        MDL_LAM = MDL_MH__EXP__2/(2.000000D+00*MDL_VEV__EXP__2)

        MDL_YT = (MDL_YMT*MDL_SQRT__2)/MDL_VEV

        MDL_MUH = SQRT(DCMPLX(MDL_LAM*MDL_VEV__EXP__2))

        MDL_I233 = MDL_YT*MDL_CONJG__CKM33

        MDL_I333 = MDL_CKM33*MDL_YT

        MDL_AXIALZUP = (3.000000D+00/2.000000D+00)*(-(MDL_EE*MDL_SW)/
     $   (6.000000D+00*MDL_CW))-(1.000000D+00/2.000000D+00)*(
     $   (MDL_CW*MDL_EE)/(2.000000D+00*MDL_SW))

        MDL_AXIALZDOWN = (-1.000000D+00/2.000000D+00)*(-(MDL_CW*MDL_EE)
     $   /(2.000000D+00*MDL_SW))+(-3.000000D+00/2.000000D+00)*(-
     $   (MDL_EE*MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_VECTORZUP = (1.000000D+00/2.000000D+00)*((MDL_CW*MDL_EE)/
     $   (2.000000D+00*MDL_SW))+(5.000000D+00/2.000000D+00)*(-
     $   (MDL_EE*MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_VECTORZDOWN = (1.000000D+00/2.000000D+00)*(-(MDL_CW*MDL_EE)
     $   /(2.000000D+00*MDL_SW))+(-1.000000D+00/2.000000D+00)*(-
     $   (MDL_EE*MDL_SW)/(6.000000D+00*MDL_CW))

        MDL_AXIALG0UP3 = -MDL_YT/(MDL_SQRT__2)

        MDL_VECTORHUP3 = -MDL_COMPLEXI*MDL_YT/(MDL_SQRT__2)

        MDL_VECTORAUP = (2.000000D+00*MDL_EE)/3.000000D+00

        MDL_VECTORADOWN = -(MDL_EE)/3.000000D+00

        MDL_VECTORWMDXU = (1.000000D+00/2.000000D+00)*((MDL_EE)/
     $   (MDL_SW*MDL_SQRT__2))

        MDL_AXIALWMDXU = (-1.000000D+00/2.000000D+00)*((MDL_EE)/
     $   (MDL_SW*MDL_SQRT__2))

        MDL_VECTORWPUXD = (1.000000D+00/2.000000D+00)*((MDL_EE)/
     $   (MDL_SW*MDL_SQRT__2))

        MDL_AXIALWPUXD = -(1.000000D+00/2.000000D+00)*((MDL_EE)/
     $   (MDL_SW*MDL_SQRT__2))

        MDL_VECTORGPUX3D3 = ((-1.000000D+00)/2.000000D+00)
     $   *MDL_CONJG__CKM33*(MDL_YT-MDL_YB)

        MDL_AXIALGPUX3D3 = ((-1.000000D+00)/2.000000D+00)
     $   *MDL_CONJG__CKM33*(MDL_YT+MDL_YB)

        MDL_VECTORGMDX3U3 = (1.000000D+00/2.000000D+00)*MDL_CKM33*
     $   (MDL_YT-MDL_YB)

        MDL_AXIALGMDX3U3 = (-1.000000D+00/2.000000D+00)*MDL_CKM33*
     $   (MDL_YT+MDL_YB)

        MDL_EE__EXP__2 = MDL_EE**2

        MDL_EE__EXP__3 = MDL_EE**3

        MDL_EE__EXP__4 = MDL_EE**4

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      AS = G**2/4/PI

      MDL_SQRT__AS = SQRT(DCMPLX(AS))

      MDL_G__EXP__4 = G**4

      MDL_G__EXP__2 = G**2

      MDL_G_UVC_1EPS_ = ((MDL_G__EXP__2)/(2.000000D+00*4.800000D
     $ +01*PI**2))*4.000000D+00*MDL_TF

      MDL_GWCFT_UV_T_1EPS_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $ *4.000000D+00*MDL_TF))

      MDL_TWCFT_UV_1EPS_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))
     $ *3.000000D+00*MDL_CF))

      MDL_G__EXP__3 = G**3

      MDL_MU_R__EXP__2 = MU_R**2

      MDL_GWCFT_UV_T_FIN_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(((MDL_G__EXP__2)/(2.000000D+00*4.800000D+01*PI**2))
     $ *4.000000D+00*MDL_TF*REGLOG(DCMPLX(MDL_MT__EXP__2
     $ /MDL_MU_R__EXP__2))))

      MDL_TWCFT_UV_FIN_ = COND(DCMPLX(MDL_MT),DCMPLX(0.000000D+00)
     $ ,DCMPLX(-((MDL_G__EXP__2)/(2.000000D+00*1.600000D+01*PI**2))
     $ *MDL_CF*(4.000000D+00-3.000000D+00*REGLOG(DCMPLX(MDL_MT__EXP__2
     $ /MDL_MU_R__EXP__2)))))

C     
C     Definition of the EW coupling used in the write out of aqed
C     
      GAL(1) = 3.5449077018110318D0 / DSQRT(AEWM1)
      GAL(2) = 1D0


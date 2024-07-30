ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C     
      IF(READLHA) THEN

        MP__G = 2 * SQRT(MP__AS*MP__PI)  ! for the first init

        MP__MDL_YB = 0.000000E+00_16

        MP__MDL_CKM31 = 0.000000E+00_16

        MP__MDL_CKM32 = 0.000000E+00_16

        MP__MDL_MB__EXP__4 = 0.000000E+00_16

        MP__MDL_MB__EXP__2 = 0.000000E+00_16

        MP__MDL_CKM11 = 1.000000E+00_16

        MP__MDL_SCKM33 = 1.000000E+00_16

        MP__MDL_SCKM11 = 1.000000E+00_16

        MP__MDL_VEP = 0.000000E+00_16

        MP__MDL_CKM33 = 1.000000E+00_16

        MP__MDL_CMSPARAM = 0.000000E+00_16

        MP__MDL_MC__EXP__2 = 0.000000E+00_16

        MP__MDL_SCKM33__EXP__2 = 1.000000E+00_16

        MP__MDL_QL__EXP__2 = 1.000000E+00_16

        MP__MDL_CONJG__CKM22 = 1.000000E+00_16

        MP__MDL_SCKM22 = 1.000000E+00_16

        MP__MDL_CKM22 = 1.000000E+00_16

        MP__MDL_SCKM22__EXP__2 = 1.000000E+00_16

        MP__MDL_MD__EXP__2 = 0.000000E+00_16

        MP__MDL_MU__EXP__2 = 0.000000E+00_16

        MP__MDL_MM__EXP__2 = 0.000000E+00_16

        MP__MDL_MTA__EXP__2 = 0.000000E+00_16

        MP__MDL_QL__EXP__4 = 1.000000E+00_16

        MP__MDL_CONJG__CKM11 = 1.000000E+00_16

        MP__MDL_LHV = 1.000000E+00_16

        MP__MDL_CONJG__CKM31 = 0.000000E+00_16

        MP__MDL_CONJG__CKM32 = 0.000000E+00_16

        MP__MDL_CONJG__CKM33 = 1.000000E+00_16

        MP__MDL_NCOL = 3.000000E+00_16

        MP__MDL_CA = 3.000000E+00_16

        MP__MDL_TF = 5.000000E-01_16

        MP__MDL_CF = (4.000000E+00_16/3.000000E+00_16)

        MP__MDL_COMPLEXI = CMPLX(0.000000E+00_16,1.000000E+00_16
     $   ,KIND=16)

        MP__MDL_MZ__EXP__2 = MP__MDL_MZ**2

        MP__MDL_MW__EXP__2 = MP__MDL_MW**2

        MP__MDL_SQRT__2 = SQRT(CMPLX((2.000000E+00_16),KIND=16))

        MP__MDL_SW2 = 1.000000E+00_16-MP__MDL_MW__EXP__2
     $   /MP__MDL_MZ__EXP__2

        MP__MDL_CW = SQRT(CMPLX((1.000000E+00_16-MP__MDL_SW2),KIND=16))

        MP__MDL_SQRT__SW2 = SQRT(CMPLX((MP__MDL_SW2),KIND=16))

        MP__MDL_SW = MP__MDL_SQRT__SW2

        MP__MDL_MH__EXP__2 = MP__MDL_MH**2

        MP__MDL_I3D = (-1.000000E+00_16/2.000000E+00_16)

        MP__MDL_I3U = (1.000000E+00_16/2.000000E+00_16)

        MP__MDL_I3L = (-1.000000E+00_16/2.000000E+00_16)

        MP__MDL_I3V = (1.000000E+00_16/2.000000E+00_16)

        MP__MDL_QD = (-1.000000E+00_16/3.000000E+00_16)

        MP__MDL_QU = (2.000000E+00_16/3.000000E+00_16)

        MP__MDL_QL = (-1.000000E+00_16)

        MP__MDL_NCOL__EXP__2 = MP__MDL_NCOL**2

        MP__MDL_SW__EXP__2 = MP__MDL_SW**2

        MP__MDL_MT__EXP__2 = MP__MDL_MT**2

        MP__MDL_CW__EXP__2 = MP__MDL_CW**2

        MP__MDL_MT__EXP__4 = MP__MDL_MT**4

        MP__MDL_SW__EXP__4 = MP__MDL_SW**4

        MP__MDL_MW__EXP__4 = MP__MDL_MW**4

        MP__MDL_MT__EXP__3 = MP__MDL_MT**3

        MP__MDL_MH__EXP__4 = MP__MDL_MH**4

        MP__MDL_MZ__EXP__4 = MP__MDL_MZ**4

        MP__MDL_MT__EXP__6 = MP__MDL_MT**6

        MP__MDL_MW__EXP__6 = MP__MDL_MW**6

        MP__MDL_MT__EXP__8 = MP__MDL_MT**8

        MP__MDL_MW__EXP__8 = MP__MDL_MW**8

        MP__MDL_MT__EXP__10 = MP__MDL_MT**10

        MP__MDL_MW__EXP__10 = MP__MDL_MW**10

        MP__MDL_MZ__EXP__6 = MP__MDL_MZ**6

        MP__MDL_SW__EXP__3 = MP__MDL_SW**3

        MP__MDL_MZ__EXP__3 = MP__MDL_MZ**3

        MP__MDL_CW__EXP__4 = MP__MDL_CW**4

        MP__MDL_CW__EXP__6 = MP__MDL_CW**6

        MP__MDL_MH__EXP__8 = MP__MDL_MH**8

        MP__MDL_MH__EXP__6 = MP__MDL_MH**6

        MP__MDL_CW__EXP__3 = MP__MDL_CW**3

        MP__MDL_MW__EXP__3 = MP__MDL_MW**3

        MP__MDL_QD__EXP__3 = MP__MDL_QD**3

        MP__MDL_QD__EXP__2 = MP__MDL_QD**2

        MP__MDL_I3D__EXP__2 = MP__MDL_I3D**2

        MP__MDL_I3D__EXP__3 = MP__MDL_I3D**3

        MP__MDL_QU__EXP__3 = MP__MDL_QU**3

        MP__MDL_QU__EXP__2 = MP__MDL_QU**2

        MP__MDL_I3U__EXP__2 = MP__MDL_I3U**2

        MP__MDL_I3U__EXP__3 = MP__MDL_I3U**3

        MP__MDL_QL__EXP__3 = MP__MDL_QL**3

        MP__MDL_I3L__EXP__2 = MP__MDL_I3L**2

        MP__MDL_I3L__EXP__3 = MP__MDL_I3L**3

        MP__MDL_I3V__EXP__3 = MP__MDL_I3V**3

        MP__MDL_SW__EXP__6 = MP__MDL_SW**6

        MP__MDL_SW__EXP__8 = MP__MDL_SW**8

        MP__MDL_QU__EXP__4 = MP__MDL_QU**4

        MP__MDL_QD__EXP__4 = MP__MDL_QD**4

        MP__MDL_I3L__EXP__4 = MP__MDL_I3L**4

        MP__MDL_I3U__EXP__4 = MP__MDL_I3U**4

        MP__MDL_I3D__EXP__4 = MP__MDL_I3D**4

        MP__MDL_I3V__EXP__4 = MP__MDL_I3V**4

        MP__MDL_I3V__EXP__2 = MP__MDL_I3V**2

        MP__MDL_AEW = 1.000000E+00_16/MP__AEWM1

        MP__MDL_GF = -MP__MDL_AEW*MP__MDL_MZ__EXP__2*MP__PI/
     $   (MP__MDL_SQRT__2*MP__MDL_MW__EXP__2*(MP__MDL_MW__EXP__2
     $   -MP__MDL_MZ__EXP__2))

        MP__MDL_SQRT__AEW = SQRT(CMPLX((MP__MDL_AEW),KIND=16))

        MP__MDL_EE = 2.000000E+00_16*MP__MDL_SQRT__AEW*SQRT(CMPLX(
     $   (MP__PI),KIND=16))

        MP__MDL_G1 = MP__MDL_EE/MP__MDL_CW

        MP__MDL_GW = MP__MDL_EE/MP__MDL_SW

        MP__MDL_VEV = (2.000000E+00_16*MP__MDL_MW*MP__MDL_SW)
     $   /MP__MDL_EE

        MP__MDL_VEV__EXP__2 = MP__MDL_VEV**2

        MP__MDL_LAM = MP__MDL_MH__EXP__2/(2.000000E
     $   +00_16*MP__MDL_VEV__EXP__2)

        MP__MDL_YT = (MP__MDL_YMT*MP__MDL_SQRT__2)/MP__MDL_VEV

        MP__MDL_MUH = SQRT(CMPLX((MP__MDL_LAM*MP__MDL_VEV__EXP__2)
     $   ,KIND=16))

        MP__MDL_I233 = MP__MDL_YT*MP__MDL_CONJG__CKM33

        MP__MDL_I333 = MP__MDL_CKM33*MP__MDL_YT

        MP__MDL_AXIALZUP = (3.000000E+00_16/2.000000E+00_16)*(-
     $   (MP__MDL_EE*MP__MDL_SW)/(6.000000E+00_16*MP__MDL_CW))-
     $   (1.000000E+00_16/2.000000E+00_16)*((MP__MDL_CW*MP__MDL_EE)/
     $   (2.000000E+00_16*MP__MDL_SW))

        MP__MDL_AXIALZDOWN = (-1.000000E+00_16/2.000000E+00_16)*(-
     $   (MP__MDL_CW*MP__MDL_EE)/(2.000000E+00_16*MP__MDL_SW))+(
     $   -3.000000E+00_16/2.000000E+00_16)*(-(MP__MDL_EE*MP__MDL_SW)/
     $   (6.000000E+00_16*MP__MDL_CW))

        MP__MDL_VECTORZUP = (1.000000E+00_16/2.000000E+00_16)*(
     $   (MP__MDL_CW*MP__MDL_EE)/(2.000000E+00_16*MP__MDL_SW))+
     $   (5.000000E+00_16/2.000000E+00_16)*(-(MP__MDL_EE*MP__MDL_SW)/
     $   (6.000000E+00_16*MP__MDL_CW))

        MP__MDL_VECTORZDOWN = (1.000000E+00_16/2.000000E+00_16)*(-
     $   (MP__MDL_CW*MP__MDL_EE)/(2.000000E+00_16*MP__MDL_SW))+(
     $   -1.000000E+00_16/2.000000E+00_16)*(-(MP__MDL_EE*MP__MDL_SW)/
     $   (6.000000E+00_16*MP__MDL_CW))

        MP__MDL_AXIALG0UP3 = -MP__MDL_YT/(MP__MDL_SQRT__2)

        MP__MDL_VECTORHUP3 = -MP__MDL_COMPLEXI*MP__MDL_YT/
     $   (MP__MDL_SQRT__2)

        MP__MDL_VECTORAUP = (2.000000E+00_16*MP__MDL_EE)/3.000000E
     $   +00_16

        MP__MDL_VECTORADOWN = -(MP__MDL_EE)/3.000000E+00_16

        MP__MDL_VECTORWMDXU = (1.000000E+00_16/2.000000E+00_16)*(
     $   (MP__MDL_EE)/(MP__MDL_SW*MP__MDL_SQRT__2))

        MP__MDL_AXIALWMDXU = (-1.000000E+00_16/2.000000E+00_16)*(
     $   (MP__MDL_EE)/(MP__MDL_SW*MP__MDL_SQRT__2))

        MP__MDL_VECTORWPUXD = (1.000000E+00_16/2.000000E+00_16)*(
     $   (MP__MDL_EE)/(MP__MDL_SW*MP__MDL_SQRT__2))

        MP__MDL_AXIALWPUXD = -(1.000000E+00_16/2.000000E+00_16)*(
     $   (MP__MDL_EE)/(MP__MDL_SW*MP__MDL_SQRT__2))

        MP__MDL_VECTORGPUX3D3 = ((-1.000000E+00_16)/2.000000E+00_16)
     $   *MP__MDL_CONJG__CKM33*(MP__MDL_YT-MP__MDL_YB)

        MP__MDL_AXIALGPUX3D3 = ((-1.000000E+00_16)/2.000000E+00_16)
     $   *MP__MDL_CONJG__CKM33*(MP__MDL_YT+MP__MDL_YB)

        MP__MDL_VECTORGMDX3U3 = (1.000000E+00_16/2.000000E+00_16)
     $   *MP__MDL_CKM33*(MP__MDL_YT-MP__MDL_YB)

        MP__MDL_AXIALGMDX3U3 = (-1.000000E+00_16/2.000000E+00_16)
     $   *MP__MDL_CKM33*(MP__MDL_YT+MP__MDL_YB)

        MP__MDL_EE__EXP__2 = MP__MDL_EE**2

        MP__MDL_EE__EXP__3 = MP__MDL_EE**3

        MP__MDL_EE__EXP__4 = MP__MDL_EE**4

      ENDIF
C     
C     Parameters that should be recomputed at an event by even basis.
C     
      MP__AS = MP__G**2/4/MP__PI

      MP__MDL_SQRT__AS = SQRT(CMPLX((MP__AS),KIND=16))

      MP__MDL_G__EXP__4 = MP__G**4

      MP__MDL_G__EXP__2 = MP__G**2

      MP__MDL_G_UVC_1EPS_ = ((MP__MDL_G__EXP__2)/(2.000000E
     $ +00_16*4.800000E+01_16*MP__PI**2))*4.000000E+00_16*MP__MDL_TF

      MP__MDL_GWCFT_UV_T_1EPS_ = MP_COND(CMPLX(MP__MDL_MT,KIND=16)
     $ ,CMPLX(0.000000E+00_16,KIND=16),CMPLX(-((MP__MDL_G__EXP__2)/
     $ (2.000000E+00_16*4.800000E+01_16*MP__PI**2))*4.000000E
     $ +00_16*MP__MDL_TF,KIND=16))

      MP__MDL_TWCFT_UV_1EPS_ = MP_COND(CMPLX(MP__MDL_MT,KIND=16),CMPLX
     $ (0.000000E+00_16,KIND=16),CMPLX(-((MP__MDL_G__EXP__2)/
     $ (2.000000E+00_16*1.600000E+01_16*MP__PI**2))*3.000000E
     $ +00_16*MP__MDL_CF,KIND=16))

      MP__MDL_G__EXP__3 = MP__G**3

      MP__MDL_MU_R__EXP__2 = MP__MU_R**2

      MP__MDL_GWCFT_UV_T_FIN_ = MP_COND(CMPLX(MP__MDL_MT,KIND=16)
     $ ,CMPLX(0.000000E+00_16,KIND=16),CMPLX(((MP__MDL_G__EXP__2)/
     $ (2.000000E+00_16*4.800000E+01_16*MP__PI**2))*4.000000E
     $ +00_16*MP__MDL_TF*MP_REGLOG(CMPLX((MP__MDL_MT__EXP__2
     $ /MP__MDL_MU_R__EXP__2),KIND=16)),KIND=16))

      MP__MDL_TWCFT_UV_FIN_ = MP_COND(CMPLX(MP__MDL_MT,KIND=16),CMPLX
     $ (0.000000E+00_16,KIND=16),CMPLX(-((MP__MDL_G__EXP__2)/
     $ (2.000000E+00_16*1.600000E+01_16*MP__PI**2))*MP__MDL_CF*
     $ (4.000000E+00_16-3.000000E+00_16*MP_REGLOG(CMPLX(
     $ (MP__MDL_MT__EXP__2/MP__MDL_MU_R__EXP__2),KIND=16))),KIND=16))

C     
C     Definition of the EW coupling used in the write out of aqed
C     
      MP__GAL(1) = 2 * SQRT(MP__PI/MP__AEWM1)
      MP__GAL(2) = 1D0


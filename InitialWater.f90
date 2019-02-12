!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MODULE InitialWaterSystem
  IMPLICIT NONE
  
  ! Physical System Parameters
  REAL(8), PARAMETER :: Lx = 10         ! Physical Domain Length in X
  REAL(8), PARAMETER :: Ly = 10         ! Physical Domain Length in Y

  ! Computational System Parameters
  INTEGER, PARAMETER :: nM = 256         ! Number of Collocation Pts in x
  INTEGER, PARAMETER :: nN = 256         ! Number of Collocation Pts in y

  ! Derived Parameters
  INTEGER, PARAMETER :: mM = nM - 1      ! Maximum Mode Number in X
  INTEGER, PARAMETER :: mN = nN - 1      ! Maximum Mode Number in Y  
  
  REAL(8), PARAMETER :: dx = Lx / nM     ! Physical x-spacing
  REAL(8), PARAMETER :: dy = Ly / nN     ! Physical y-spacing
  
  ! Physical System Variables
  REAL(8),    DIMENSION(0:mM,0:mN) :: u, v, h         ! x,y velocity, height
  COMPLEX(8), DIMENSION(0:mM,0:mN) :: Apv, Agp, Agm   ! Modal flow components
  
  ! Intermediary FFTW Variables
  COMPLEX(8), DIMENSION(0:nM/2, 0:mN) :: uF, vF, hF         ! FFT of u, v, h
  COMPLEX(8), DIMENSION(0:mM, 0:mN)   :: ApvF, AgpF, AgmF   ! FFT of modal flow
  
  ! FFTW Control Variables
  INTEGER(8) :: ApvPlan, AgpPlan, AgmPlan         ! FFTW:  A -> AF
  INTEGER(8) :: uFPlan, vFPlan, hFPlan            ! FFTW: uF -> u
  
  ! NetCDF Control Variables
  INTEGER :: ncFileID, &
             xDimID, yDimID, &
             xVarID, yVarID, uVarID, vVarID, hVarID, &
             lxVarID, lyVarID, dxVarID, dyVarID, &
             nMVarID, nNVarID
  CHARACTER (len = *), PARAMETER :: fileName = "InitialData.nc"
  
  ! Numerical Parameters
  REAL(4),  PARAMETER :: FPI = 3.14159
  REAL(8),  PARAMETER :: DPI = 3.14159265358979D+0
  ! REAL(16), PARAMETER :: QPI = 3.14159265358979323846264338327950Q+0
  
END MODULE InitialWaterSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM InitialWater
!
! Program Description:
!   Constructs the (u,v,h) shallow water initial state, given a placement of
!   PV and gravity wave distributions
!
! --Marshall Ward
!
! Log:
!  2 November 2006 - Program begun

  ! Modules
  USE typeSizes
  USE netcdf
  USE InitialWaterSystem

  IMPLICIT NONE
  INCLUDE 'fftw3.f'
    
  ! Computational Control Variables
  INTEGER :: i, j              ! Spatial/modal timestepping index
  INTEGER :: t                 ! Timestepping index

  REAL(8) :: kmode, lmode, Km, omega
  
  COMPLEX(8) :: upvF, vpvF, hpvF
  COMPLEX(8) :: ugpF, vgpF, hgpF
  COMPLEX(8) :: ugmF, vgmF, hgmF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT*, "Begin InitialWater"
  
  ! (1) Create FFTW Plans
  ! * This clobbers any data in the arrays
  CALL dfftw_plan_dft_2d(ApvPlan, nM, nN, Apv, ApvF, FFTW_FORWARD, & 
                         FFTW_MEASURE)
  CALL dfftw_plan_dft_2d(AgpPlan, nM, nN, Agp, AgpF, FFTW_FORWARD, & 
                         FFTW_MEASURE)
  CALL dfftw_plan_dft_2d(AgmPlan, nM, nN, Agm, AgmF, FFTW_FORWARD, & 
                         FFTW_MEASURE)

  CALL dfftw_plan_dft_c2r_2d(uFPlan, nM, nN, uF, u, FFTW_BACKWARD, & 
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vFPlan, nM, nN, vF, v, FFTW_BACKWARD, & 
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hFPlan, nM, nN, hF, h, FFTW_BACKWARD, & 
                             FFTW_MEASURE)
  
  PRINT*, "FFTW Initialized"
  
  ! (1) Initialize Data File
  ! (1a) Create Data File
  CALL check(nf90_create(path = trim(fileName), &
                         cmode = or(nf90_clobber, nf90_64bit_offset), &
                         ncid=ncFileID))
             
  ! (1b) Define Dimensions
  CALL check(nf90_def_dim(ncid = ncFileID, name = "x", len = nM, &
                          dimID = xDimID))
  
  CALL check(nf90_def_dim(ncid = ncFileID, name = "y", len = nN, &
                          dimID = yDimID))
  
  ! (1c) Define Variables
  CALL check(nf90_def_var(ncid = ncFileID, name = "x", xtype = nf90_double, &
                          dimids = xDimID, varID = xVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "y", xtype = nf90_double, &
                          dimids = yDimID, varID = yVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "u", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID /), &
                          varID = uVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "v", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID /), &
                          varID = vVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "h", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID /), &
                          varID = hVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "Lx", xtype = nf90_double, &
                          varID = LxVarID))

  CALL check(nf90_def_var(ncid = ncFileID, name = "Ly", xtype = nf90_double, &
                          varID = LyVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "dx", xtype = nf90_double, &
                          varID = dxVarID))

  CALL check(nf90_def_var(ncid = ncFileID, name = "dy", xtype = nf90_double, &
                          varID = dyVarID))

  CALL check(nf90_def_var(ncid = ncFileID, name = "nM", xtype = nf90_int, &
                          varID = nMVarID))
  
  CALL check(nf90_def_var(ncid = ncFileID, name = "nN", xtype = nf90_int, &
                          varID = nNVarID))
  
  CALL check(nf90_enddef(ncFileID))
  
  PRINT*, "Output File Created"
  
  ! (1d) Fill in Grid Values (i.e. Dimension Grid)
  DO i = 1, nM
    CALL check(nf90_put_var(ncFileID, xVarID, (i-1)*dx, start = (/ i /) ))
  END DO

  DO i = 1, nN
    CALL check(nf90_put_var(ncFileID, yVarID, (i-1)*dy, start = (/ i /) ))
  END DO

  ! (1d') Fill in Parameters
  CALL check(nf90_put_var(ncFileID, LxVarID, Lx))
  CALL check(nf90_put_var(ncFileID, LyVarID, Ly))
  
  CALL check(nf90_put_var(ncFileID, dxVarID, dx))
  CALL check(nf90_put_var(ncFileID, dyVarID, dy))

  CALL check(nf90_put_var(ncFileID, nMVarID, nM))
  CALL check(nf90_put_var(ncFileID, nNVarID, nN))

  PRINT*, "PV and Gravity Mode Placement"
  
  ! Spatially define the PV and gravity wave distrubutions  
  ! Actually I think I will force the waves in the model
  DO i = 0, mM
    DO j = 0, mN
!      IF((i*dx - Lx/2)**2 + (j*dy - Ly/2)**2 .LT. 150**2) THEN
!        Apv(i,j) = COS(2*DPI*(i*dx*8/Lx - j*dy*8/Ly))
!      ELSE
!        Apv(i,j) = 0.0
!      END IF

!      Apv(i,j) = EXP(-((i*dx-Lx/2)**2 + (j*dy-Ly/2)**2)/(2)**2)
!      Agp(i,j) = 0.125*EXP((0,1)*(2*DPI)*32*i*dx/Lx);
!      Agm(i,j) = CONJG(Agp(i,j))

!      IF((i*dx - Lx/8)**2 .KT. 3**2) THEN
!        Agp(i,j) = 0.25*COS(2*DPI*i*dx*8/Lx)
!        Agm(i,j) = 0.25*COS(2*DPI*i*dx*8/Lx)
!      ELSE
!        Agp(i,j) = 0;
!        Agm(i,j) = 0;
!      END IF
      
!      Apv(i,j) = COS(2*DPI*(i*dx*8*(SQRT(2.)-2)/Lx - j*dy*8*SQRT(2.)/Ly)/2) &
!                + COS(2*DPI*(i*dx*8*(SQRT(2.)-2)/Lx + j*dy*8*SQRT(2.)/Ly)/2)

      Apv(i,j) = COS(2*DPI*(-i*dx*15/Lx + j*dy*26/Ly)) &
                + COS(2*DPI*(-i*dx*15/Lx - j*dy*26/Ly))
      Agp(i,j) = 0.25*EXP((0,1)*(2*DPI)*i*dx*30/Lx)
      Agm(i,j) = CONJG(Agp(i,j))
    END DO
  END DO
  
  ! Transform modal data to spectral form
  CALL dfftw_execute(ApvPlan)
  CALL dfftw_execute(AgpPlan)
  CALL dfftw_execute(AgmPlan)
  
  ! Construct modal (u,v,h) flow
  DO i = 0, nM/2
    DO j = 0, mN
      kmode = 2*DPI * i / Lx
      lmode = 2*DPI * (j - nN*(j/(nN/2 + 1))) / Ly
      
      Km = SQRT(kmode**2 + lmode**2)
      omega = SQRT(1 + Km**2)
      
      upvF =  (0,1)*lmode*ApvF(i,j) / omega
      vpvF = -(0,1)*kmode*ApvF(i,j) / omega
      hpvF = -ApvF(i,j) / omega
      
      ugpF = ( omega*kmode + (0,1)*lmode) * AgpF(i,j) &
                   / (SQRT(2d0)*omega*Km)
      vgpF = (-(0,1)*kmode + omega*lmode) * AgpF(i,j) &
                   / (SQRT(2d0)*omega*Km)
      hgpF = Km * AgpF(i,j) / (SQRT(2d0)*omega)
      
      ugmF = (-omega*kmode + (0,1)*lmode) * AgmF(i,j) &
                   / (SQRT(2d0)*omega*Km)
      vgmF = (-(0,1)*kmode - omega*lmode) * AgmF(i,j) &
                   / (SQRT(2d0)*omega*Km)
      hgmF = Km * AgmF(i,j) / (SQRT(2d0)*omega)
      
      uF(i,j) = (upvF + ugpF + ugmF)/(nM*nN)
      vF(i,j) = (vpvF + vgpF + vgmF)/(nM*nN)
      hF(i,j) = (hpvF + hgpF + hgmF)/(nM*nN)
    END DO
  END DO
    
  !** Redo k = l = 0 (ha ha ha) **!   
  upvF =  0d0
  vpvF =  0d0
  hpvF = -1d0

  ugpF = AgpF(0,0) / SQRT(2d0)
  vgpF = -(0,1)*AgpF(0,0) / SQRT(2d0)
  hgpF = 0d0

  ugmF = -Agm(0,0) / SQRT(2d0)
  vgmF = -(0,1)*Agm(0,0) / SQRT(2d0)
  hgmF = 0d0
  
  uF(0,0) = (upvF + ugpF + ugmF)/(nM*nN)
  vF(0,0) = (vpvF + vgpF + vgmF)/(nM*nN)
  hF(0,0) = (hpvF + hgpF + hgmF)/(nM*nN)
  
  CALL dfftw_execute(uFPlan)
  CALL dfftw_execute(vFPlan)
  CALL dfftw_execute(hFPlan)  

  !** Write (u,v,h) to "InitialData.nc" **
  CALL check(nf90_put_var(ncFileID, uVarID, u, &
                          start = (/ 1, 1 /), count = (/ nM, nN /) ))
  CALL check(nf90_put_var(ncFileID, vVarID, v, &
                          start = (/ 1, 1 /), count = (/ nM, nN /) ))
  CALL check(nf90_put_var(ncFileID, hVarID, h, &
                          start = (/ 1, 1 /), count = (/ nM, nN /) ))
  
  ! ** CLEANUP **
  ! ** Deallocate FFTW Plans **
  CALL dfftw_destroy_plan(ApvPlan)
  CALL dfftw_destroy_plan(AgpPlan)
  CALL dfftw_destroy_plan(AgmPlan)
  
  CALL dfftw_destroy_plan(uFPlan)
  CALL dfftw_destroy_plan(vFPlan)
  CALL dfftw_destroy_plan(hFPlan)
  
  ! Close NetCDF Data File
  CALL check(nf90_close(ncFileID))
  
!! ***RETURN TO SHELL***
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
CONTAINS

  ! Internal subroutine - checks error status after each netcdf, prints out 
  ! text message each time an error code is returned. 
  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status
    
    IF(status /= nf90_noerr) THEN 
      PRINT *, trim(nf90_strerror(status))
    END IF
  END SUBROUTINE check  

!!!!!!!!!!!!!!!

END PROGRAM InitialWater

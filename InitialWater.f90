!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MODULE InitialWaterSystem
  IMPLICIT NONE
  
  ! Numerical Parameters
  REAL(4),  PARAMETER :: FPI = 3.14159
  REAL(8),  PARAMETER :: DPI = 3.14159265358979D+0
  ! REAL(16), PARAMETER :: QPI = 3.14159265358979323846264338327950Q+0
  
  ! My Special Variables
  REAL(8), PARAMETER :: Jn = 0.01   ! Control parameter for jet spread
  REAL(8), PARAMETER :: N0 = 2      ! Mode number for singlemode
  
  ! Physical System Parameters
  REAL(8), PARAMETER :: Lx = 2*DPI*1    ! Physical Domain Length in X
  REAL(8), PARAMETER :: Ly = 2*DPI*1    ! Physical Domain Length in Y

  ! Computational System Parameters
  INTEGER, PARAMETER :: nM = 128         ! Number of Collocation Pts in x
  INTEGER, PARAMETER :: nN = 128         ! Number of Collocation Pts in y

  ! Derived Parameters
  INTEGER, PARAMETER :: mM = nM - 1      ! Maximum Mode Number in X
  INTEGER, PARAMETER :: mN = nN - 1      ! Maximum Mode Number in Y  
` 
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

CONTAINS
  ! Initial PV Mode Equation
  REAL(8) FUNCTION PVEqn(x, y)
    REAL(8), INTENT(IN) :: x, y
!    PVEqn = 0.0
!    PVEqn = -5.0*DTANH((x-Lx/2))*EXP(-Jn*(x-Lx/2)**2)
!    PVEqn = COS(2*DPI*(x*(0-N0)/Lx + y*(N0-0)/Ly))
  END FUNCTION PVEqn

  ! Initial +Gravity Wave Equation
  COMPLEX(8) FUNCTION GpEqn(x, y)
    REAL(8), INTENT(IN) :: x, y
!    GpEqn = 0.0
    GpEqn = 0.5*EXP((0,1)*(2*DPI)*x*N0/Lx)
  END FUNCTION GpEqn

END MODULE InitialWaterSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM InitialWater
!
! Program Description:
!   Constructs the (u,v,h) shallow water initial state, given a placement of
!   PV and gravity wave distributions
!
! --Marshall Ward

  ! Modules
  USE InitialWaterSystem

  IMPLICIT NONE
  INCLUDE 'fftw3.f'
    
  ! Computational Control Variables
  INTEGER :: i, j              ! Spatial/modal timestepping index
  INTEGER :: t                 ! Timestepping index

  REAL(8) :: kmode, lmode, Km, omega
  REAL(8), DIMENSION(0:mM) :: xGrid
  REAL(8), DIMENSION(0:mN) :: yGrid
  
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
  OPEN(unit = 20, file = 'InitialData.raw', status = 'REPLACE', &
       form = 'UNFORMATTED')
  
  PRINT *, "Output File Created"
  
  ! Write controlling parameters to "InitialData.raw"
  WRITE(20) Lx; WRITE(20) Ly; WRITE(20) dx; WRITE(20) dy
  WRITE(20) nM; WRITE(20) nN; WRITE(20) mM; WRITE(20) mN

  ! Fill in Grid Values
  DO i = 0, mM
    xGrid(i) = i*dx
  END DO

  DO i = 0, mN
    yGrid(i) = i*dy
  END DO
  
  WRITE(20) xGrid; WRITE(20) yGrid

  PRINT*, "PV and Gravity Mode Placement"
  
  ! Spatially define the PV and gravity wave distrubutions  
  ! Actually I think I will force the waves in the model
  DO i = 0, mM
    DO j = 0, mN
    
!      Apv(i,j) = COS(2*DPI*(i*dx*(0-1)/Lx + j*dy*(1-0)/Ly))
!      Agp(i,j) = 0.25*EXP((0,1)*(2*DPI)*i*dx*1/Lx)
!      Agm(i,j) = CONJG(Agp(i,j))

      Apv(i,j) = PVEqn(i*dx, j*dy)
      Agp(i,j) = GpEqn(i*dx, j*dy)
      Agm(i,j) = CONJG(GpEqn(i*dx, j*dy))

    END DO
  END DO
  
  ! Transform modal data to spectral form
  CALL dfftw_execute(ApvPlan)
  CALL dfftw_execute(AgpPlan)
  CALL dfftw_execute(AgmPlan)
  
  ! Construct modal (u,v,h) flow
  DO j = 0, mN
    DO i = 0, nM/2
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
  hpvF = -ApvF(0,0)

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

  !** Write (u,v,h) to "InitialData.raw" **
  WRITE(20) u; WRITE(20) v; WRITE(20) h
  
  ! ** CLEANUP **
  ! ** Deallocate FFTW Plans **
  CALL dfftw_destroy_plan(ApvPlan)
  CALL dfftw_destroy_plan(AgpPlan)
  CALL dfftw_destroy_plan(AgmPlan)
  
  CALL dfftw_destroy_plan(uFPlan)
  CALL dfftw_destroy_plan(vFPlan)
  CALL dfftw_destroy_plan(hFPlan)
  
  ! Close Data File
  CLOSE(20)
  
!! ***RETURN TO SHELL***
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END PROGRAM InitialWater

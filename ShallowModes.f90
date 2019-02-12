!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MODULE ShallowModesSystem
  IMPLICIT NONE
  
  ! Spatial System Parameters (from input file)
  REAL(8) :: Lx, Ly                      ! Physical Domain Length in x and y
  REAL(8) :: dx, dy                      ! Physical x- and y-spacing
  REAL(8) :: dt                          ! Timestep size in the computation
  INTEGER :: nM, nN                      ! Num of Collocation Pts in x and y
  INTEGER :: mM, mN                      ! Maximum Mode Number in x and y  
  INTEGER :: nT                          ! Number of Timesteps
  REAL(8) :: epsilon                     ! Rossby Number
  INTEGER :: dataRate                    ! nT/DataRate = number of stored steps
  INTEGER :: buffSize                    ! Number of stored steps in buffer

  REAL(8) :: Tmax

  ! Physical System Variables
  REAL(8), POINTER :: u(:,:), v(:,:), h(:,:)    ! velocity and height
  REAL(8), POINTER :: dataBuffer(:,:,:,:)  
  
  ! Grids
  REAL(8), POINTER :: xGrid(:), yGrid(:), tGrid(:)  ! x,y,t stencils

  ! Spectal System Variables
  COMPLEX(8), POINTER :: uF(:,:),   vF(:,:),   hF(:,:)    ! velocity and height
  COMPLEX(8), POINTER :: Apv(:,:), Agp(:,:), Agm(:,:)  ! modal components

  ! FFTW Control Variables
  INTEGER(8) :: uPlan, vPlan, hPlan            ! FFTW:  u -> uF

  ! Numerical Parameters
  REAL(4),  PARAMETER :: FPI = 3.14159
  REAL(8),  PARAMETER :: DPI = 3.14159265358979D+0
  !REAL(16), PARAMETER :: QPI = 3.14159265358979323846264338327950Q+0
  
END MODULE ShallowModesSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ShallowModes
!
! Program Description:
!   Decomposes a flow into its PV containing flow and the two gravity modes
!
! --Marshall Ward
!

  ! Modules
  USE ShallowModesSystem

  IMPLICIT NONE
  INCLUDE 'fftw3.f'
  
  ! Computational Control Variables
  INTEGER :: i, j, k           ! Spatial/modal timestepping index
  INTEGER :: t                 ! Timestepping index
  REAL(8) :: kmode, lmode      ! Wavenumbers
  REAL(8) :: Km, omega         ! useful junk
  
  INTEGER :: bufferIndex
    
!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT*, "MODAL DECOMPOSITION OF SHALLOW WATER DATA"
  PRINT*, "*****************************************"
  
  !****************************************************************************

  !**** Open Input File ****!
  OPEN(unit = 20, file = 'ShallowData.raw', status = 'OLD', &
       form = 'UNFORMATTED')
  
  !**** Read parameter data ****
  READ(20) Lx; READ(20) Ly; READ(20) dx; READ(20) dy; READ(20) dt
  READ(20) nM; READ(20) nN; READ(20) nT; READ(20) mM; READ(20) mN
  READ(20) epsilon; READ(20) dataRate; READ(20) buffSize

  Tmax = dt*nT
  
  !**** Allocate variable memory
  ALLOCATE(u(0:mM, 0:mN)); ALLOCATE(uF(0:nM/2, 0:mN))
  ALLOCATE(v(0:mM, 0:mN)); ALLOCATE(vF(0:nM/2, 0:mN))
  ALLOCATE(h(0:mM, 0:mN)); ALLOCATE(hF(0:nM/2, 0:mN))

  ALLOCATE(Apv(0:mM, 0:mN));
  ALLOCATE(Agp(0:mM, 0:mN));
  ALLOCATE(Agm(0:mM, 0:mN));
  
  ALLOCATE(dataBuffer(0:mM, 0:mN, 1:buffSize, 3))
  
  ALLOCATE(xGrid(0:mM)); ALLOCATE(yGrid(0:mN)); ALLOCATE(tGrid(0:nT/dataRate))
  
  READ(20) xGrid; READ(20) yGrid; READ(20) tGrid
  
  PRINT*, "FFTW Initialization"
  
  ! (1) Create FFTW Plans
  ! * This clobbers any data in the arrays
  CALL dfftw_plan_dft_r2c_2d(uPlan, nM, nN, u, uF, FFTW_FORWARD, FFTW_MEASURE)
  CALL dfftw_plan_dft_r2c_2d(vPlan, nM, nN, v, vF, FFTW_FORWARD, FFTW_MEASURE)
  CALL dfftw_plan_dft_r2c_2d(hPlan, nM, nN, h, hF, FFTW_FORWARD, FFTW_MEASURE)

  !****************************************************************************
  
  PRINT*, "Output File Initialization"
  
  ! (1) Initialize Data File
  ! (1a) Create Data File
  OPEN(unit = 30, file = 'ShallowModes.raw', status = 'REPLACE', &
        form = 'UNFORMATTED')
  
  ! (1b) Save parameters
  WRITE(30) Lx; WRITE(30) Ly; WRITE(30) dx; WRITE(30) dy; WRITE(30) dt
  WRITE(30) nM; WRITE(30) nN; WRITE(30) nT; WRITE(30) mM; WRITE(30) mN
  WRITE(30) epsilon; WRITE(30) dataRate; WRITE(30) buffSize
  
  WRITE(30) xGrid; WRITE(30) yGrid; WRITE(30) tGrid
  
  ! Write a k and l grid? It would be more convenient
    
  !***************************************************************************!

  PRINT*, "Computing Modal Decomposition"

  ! Need to do t=0 separately because of my dumdum file format

  DO t = 0, (nT/dataRate)
    IF(t == 0) THEN
      READ(20) u; READ(20) v; READ(20) h
    ELSE
      bufferIndex = MOD(t-1,buffSize)
    
      IF(MOD(t-1, buffSize) == 0) THEN
        READ(20) dataBuffer
      END IF

      DO j = 0, nN
        DO i = 0, nM
          u(i,j) = dataBuffer(i,j,bufferIndex, 1)
        END DO
      END DO
      
      DO j = 0, nN
        DO i = 0, nM
          v(i,j) = dataBuffer(i,j,bufferIndex, 2)
        END DO
      END DO
      
      DO j = 0, nN
        DO i = 0, nM
          h(i,j) = dataBuffer(i,j,bufferIndex, 3)
        END DO
      END DO  
    END IF
    
    !** Compute DFT **!
    CALL dfftw_execute(uPlan)
    CALL dfftw_execute(vPlan)
    CALL dfftw_execute(hPlan)
    
    !** Compute mode amplitudes **!
    DO j = 0, mN
      DO i = 0, nM/2
        kmode = 2*DPI * i / Lx
        lmode = 2*DPI * (j - nN*(j/(nN/2 + 1))) / Ly
      
        Km = SQRT(kmode**2 + lmode**2)
        omega = SQRT(1 + Km**2)

        Apv(i,j) = (-(0,1)*lmode*uF(i,j) + (0,1)*kmode*vF(i,j) - hF(i,j)) &
                  / omega / (nM*nN)
        
        Agp(i,j) = ((omega*kmode - (0,1)*lmode)/Km*uF(i,j) &
                  + ((0,1)*kmode + omega*lmode)/Km*vF(i,j) &
                  + Km*hF(i,j))/(SQRT(2d0)*omega) / (nM*nN)
               
        Agm(i,j) = ((-omega*kmode - (0,1)*lmode)/Km*uF(i,j) &
                  + ((0,1)*kmode - omega*lmode)/Km*vF(i,j) &
                  + Km*hF(i,j))/(SQRT(2d0)*omega) / (nM*nN)
                       
      END DO
    END DO
    
    !** Redo k = l = 0 (ha ha ha) **!
    Agp(0,0) = ( uF(0,0) + (0,1)*vF(0,0)) / SQRT(2d0) / (nM*nN) 
    Agm(0,0) = (-uF(0,0) + (0,1)*vF(0,0)) / SQRT(2d0) / (nM*nN)
    
    WRITE(30) Apv; WRITE(30) Agp; WRITE(30) Agm;
    
    ! End timestep
  END DO  
  
  !****************************************************************************
  ! CLEANUP
  
  ! Deallocate plans
  CALL dfftw_destroy_plan(uPlan)
  CALL dfftw_destroy_plan(vPlan)
  CALL dfftw_destroy_plan(hPlan)

  CLOSE(20); CLOSE(30)
  
!! ***RETURN TO SHELL***
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END PROGRAM ShallowModes

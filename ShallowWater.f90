!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MODULE ShallowWaterSystem
  IMPLICIT NONE

  ! Spatial System Parameters (Stored in 'InitialData.nc')
  REAL(8) :: Lx, Ly                       ! Physical Domain Length in x and y
  INTEGER :: nM, nN                       ! Num of Collocation Pts in x and y
  !   (Derived quantities)
  REAL(8) :: dx, dy                       ! Physical x- and y-spacing
  INTEGER :: mM, mN                       ! Maximum Mode Number in x and y  

  ! Control Parameters
  LOGICAL, PARAMETER :: dealias = .TRUE.
  
  ! Not Implemented Yet
  LOGICAL, PARAMETER :: shockfilter = .FALSE.

  ! Storage Paramters
  INTEGER, PARAMETER :: buffSize = 1e2
  INTEGER, PARAMETER :: dataRate = 1e1

  ! Temporal System Parameters
  REAL(8), PARAMETER :: dt = 1D-2         ! Integration Timestep
  INTEGER, PARAMETER :: nT = 1e3          ! Number of Timesteps
  REAL(8), PARAMETER :: Tmax = dt*nT      ! Integration Time

  ! Physical Scaling Parameters
  REAL(8), PARAMETER :: epsilon = 0.1     ! Rossby Number

  ! Grids
  REAL(8), POINTER :: xGrid(:), yGrid(:)  ! Space Grid
  REAL(8), POINTER :: kGrid(:), lGrid(:)  ! Mode Grid
  REAL(8), POINTER :: tGrid(:)            ! Time Grid

  ! Physical System Variables
  REAL(8), POINTER :: u(:,:), ux(:,:), uy(:,:), Fu(:,:)     ! x-velocity
  REAL(8), POINTER :: v(:,:), vx(:,:), vy(:,:), Fv(:,:)     ! y-velocity
  REAL(8), POINTER :: h(:,:), hx(:,:), hy(:,:), Fh(:,:)     ! Surface Ht

  COMPLEX(8), POINTER :: Fq(:,:)    ! PV Mode Forcing
  COMPLEX(8), POINTER :: Fp(:,:)    ! Gravity Mode Forcing (NB: Fm = (Fp)*)

  ! Intermediary FFTW Variables
  COMPLEX(8), POINTER :: uF(:,:), uxF(:,:), uyF(:,:), FuF(:,:)   ! FFT u data
  COMPLEX(8), POINTER :: vF(:,:), vxF(:,:), vyF(:,:), FvF(:,:)   ! FFT v data
  COMPLEX(8), POINTER :: hF(:,:), hxF(:,:), hyF(:,:), FhF(:,:)   ! FFT h data

  COMPLEX(8), POINTER :: FqF(:,:), FpF(:,:) ! FFT Modal Forcing

  ! FFTW Control Variables
  INTEGER(8) :: uPlan, uxPlan, uyPlan, FuPlan        ! FFTW u plans
  INTEGER(8) :: vPlan, vxPlan, vyPlan, FvPlan        ! FFTW v plans
  INTEGER(8) :: hPlan, hxPlan, hyPlan, FhPlan        ! FFTW h plans
  INTEGER(8) :: FqPlan, FpPlan                       ! FFTW Forcing Plans

  INTEGER(8) :: uFilterPlan, vFilterPlan, hFilterPlan

  ! Numerical Parameters
  REAL(4),  PARAMETER :: FPI = 3.14159
  REAL(8),  PARAMETER :: DPI = 3.14159265358979D+0
  ! REAL(16), PARAMETER :: QPI = 3.14159265358979323846264338327950Q+0

CONTAINS
  ! PV Mode Forcing Equation
  REAL(8) FUNCTION FqEqn(x, y, t)
    REAL(8), INTENT(IN) :: x, y, t
    FqEqn = 0.0
  END FUNCTION FqEqn

  ! +Gravity Wave Forcing Equation
  COMPLEX(8) FUNCTION FpEqn(x, y, t)
    REAL(8), INTENT(IN) :: x, y, t
!    FpEqn = EXP((0.0,1.0)*(5*x-sqrt(1+5.0)**2)*t)
!    FpEqn = EXP(-1**2*(x-5.0)**2)*EXP((0.0,1.0)*(1.6473*x - 1.9271*t))
    FpEqn = 0.0
  END FUNCTION FpEqn
  
  ! Jump Discontinuity Filter Function
  ! This is not a smart idea
  ! - It produces nicer results, not more accurate results
  ! - I need a physically motivated smoothing, a closure model of some sort
  REAL(8) FUNCTION FilterEqn(x)
    REAL(8), INTENT(IN) :: x
  ! [Raised Cosine]
    FilterEqn = 0.5*(1+DCOS(x))
  END FUNCTION FilterEqn
  
END MODULE ShallowWaterSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ShallowWaterModel
!
! Program Description:
!   Pseudospectral shallow water model
!   Periodic Rectangular Boundary Conditions
!   Free Parameters:
!     Rossby Number (epsilon)
!     Domain Size (Lx, Ly) [Specified in InitialWater]
!     Initial Conditions [Specified in InitialWater]
!
! --Marshall Ward
!

  ! Modules
  USE ShallowWaterSystem

  IMPLICIT NONE
  INCLUDE 'fftw3.f'
    
  ! Computational Control Variables
  INTEGER :: i, j              ! Spatial/modal timestepping index
  INTEGER :: t                 ! Timestepping index
  
  ! Physical Control Variables
  REAL(8), POINTER :: uLast(:,:), uTemp(:,:), uF2(:,:)   ! Prior and current u
  REAL(8), POINTER :: vLast(:,:), vTemp(:,:), vF2(:,:)   ! Prior and current v
  REAL(8), POINTER :: hLast(:,:), hTemp(:,:), hF2(:,:)   ! Prior and current h
  REAL(8) :: energy
  
  ! I/O Write Buffer
  REAL(8), POINTER :: dataBuffer(:,:,:,:)

  INTEGER :: tFileWrite, tFileBase, tBufferWrite, tBufferIndex

!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT*, "Begin ShallowWater"
  
  !** Load Input File **!
  OPEN(unit = 20, file = 'InitialData.raw', status = 'OLD', &
       form = 'UNFORMATTED')
  
  !**** Get the data size variables
  READ(20) Lx; READ(20) Ly; READ(20) dx; READ(20) dy
  READ(20) nM; READ(20) nN; READ(20) mM; READ(20) mN

  !**** Allocate variable memory
  ALLOCATE( u(0:mM, 0:mN)); ALLOCATE( uF(0:nM/2, 0:mN))
  ALLOCATE( v(0:mM, 0:mN)); ALLOCATE( vF(0:nM/2, 0:mN))
  ALLOCATE( h(0:mM, 0:mN)); ALLOCATE( hF(0:nM/2, 0:mN))

  ALLOCATE(ux(0:mM, 0:mN)); ALLOCATE(uxF(0:nM/2, 0:mN))
  ALLOCATE(vx(0:mM, 0:mN)); ALLOCATE(vxF(0:nM/2, 0:mN))
  ALLOCATE(hx(0:mM, 0:mN)); ALLOCATE(hxF(0:nM/2, 0:mN))

  ALLOCATE(uy(0:mM, 0:mN)); ALLOCATE(uyF(0:nM/2, 0:mN))
  ALLOCATE(vy(0:mM, 0:mN)); ALLOCATE(vyF(0:nM/2, 0:mN))
  ALLOCATE(hy(0:mM, 0:mN)); ALLOCATE(hyF(0:nM/2, 0:mN))

  ALLOCATE(Fu(0:mM, 0:mN)); ALLOCATE(FuF(0:nM/2, 0:mN))
  ALLOCATE(Fv(0:mM, 0:mN)); ALLOCATE(FvF(0:nM/2, 0:mN))
  ALLOCATE(Fh(0:mM, 0:mN)); ALLOCATE(FhF(0:nM/2, 0:mN))

  ALLOCATE(Fq(0:mM, 0:mN)); ALLOCATE(FqF(0:mM, 0:mN))
  ALLOCATE(Fp(0:mM, 0:mN)); ALLOCATE(FpF(0:mM, 0:mN))

  ALLOCATE(uLast(0:mM, 0:mN))
  ALLOCATE(uTemp(0:mM, 0:mN))
  ALLOCATE(  uF2(0:mM, 0:mN))

  ALLOCATE(vLast(0:mM, 0:mN))
  ALLOCATE(vTemp(0:mM, 0:mN))
  ALLOCATE(  vF2(0:mM, 0:mN))

  ALLOCATE(hLast(0:mM, 0:mN))
  ALLOCATE(hTemp(0:mM, 0:mN))
  ALLOCATE(  hF2(0:mM, 0:mN))

  ALLOCATE(dataBuffer(0:mM, 0:mN, 0:buffSize-1, 3))
  
  ALLOCATE(xGrid(0:mM)); ALLOCATE(yGrid(0:mN))
  ALLOCATE(kGrid(0:mM)); ALLOCATE(lGrid(0:mN))
  ALLOCATE(tGrid(0:(nT/dataRate)))
  
  PRINT*, "Memory Allocated"
  
  !** Create FFTW Plans **!
  CALL dfftw_plan_dft_r2c_2d(uPlan,  nM, nN, u,   uF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(uxPlan, nM, nN, uxF, ux, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(uyPlan, nM, nN, uyF, uy, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  CALL dfftw_plan_dft_r2c_2d(vPlan,  nM, nN, v,   vF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vxPlan, nM, nN, vxF, vx, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vyPlan, nM, nN, vyF, vy, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  CALL dfftw_plan_dft_r2c_2d(hPlan,  nM, nN, h,   hF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hxPlan, nM, nN, hxF, hx, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hyPlan, nM, nN, hyF, hy, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  !** Forcing FFTW Plans **!
  CALL dfftw_plan_dft_2d(FqPlan, nM, nN, Fq, FqF, FFTW_FORWARD, &
                         FFTW_MEASURE)
  CALL dfftw_plan_dft_2d(FpPlan, nM, nN, Fp, FpF, FFTW_FORWARD, &
                         FFTW_MEASURE)

  CALL dfftw_plan_dft_c2r_2d(FuPlan, nM, nN, FuF, Fu, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(FvPlan, nM, nN, FvF, Fv, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(FhPlan, nM, nN, FhF, Fh, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  !** Dealias FFTW Plans **!
  CALL dfftw_plan_dft_c2r_2d(uFilterPlan, nM, nN, uF, u, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vFilterPlan, nM, nN, vF, v, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hFilterPlan, nM, nN, hF, h, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  
  PRINT*, "FFTW Initialized"
  
  ! Get (x,y) grid
  READ(20) xGrid; READ(20) yGrid

  ! Get (u0,v0,h0)
  READ(20) u; READ(20) v; READ(20) h
  
  ! Close 'InitialData.raw'
  CLOSE(20)

  ! Initialize Output File
  OPEN(unit = 30, file = 'ShallowData.raw', status = 'REPLACE', &
        form = 'UNFORMATTED')
  OPEN(unit = 31, file = 'LastData.raw', status = 'REPLACE', &
       form = 'UNFORMATTED')
    
  PRINT*, "Output Files Created"
  
  ! Create k grid
  DO i = 0, mM
    kGrid(i) = (2.0*DPI * (i - nM*(i/(nM/2 + 1)))) / Lx   ! (Check this)
  END DO

  ! Create l grid
  DO j = 0, mN
    lGrid(j) = (2.0*DPI * (j - nN*(j/(nN/2 + 1)))) / Ly   ! (Check this)
  END DO

  ! Create t grid
  DO t = 0, (nT/DataRate)
    tGrid(t) = t*dt*DataRate
  END DO

  ! Save parameters to 'ShallowData.raw'
  WRITE(30) Lx; WRITE(30) Ly; WRITE(30) dx; WRITE(30) dy; WRITE(30) dt
  WRITE(30) nM; WRITE(30) nN; WRITE(30) nT; WRITE(30) mM; WRITE(30) mN
  WRITE(30) epsilon; WRITE(30) dataRate; WRITE(30) buffSize
  
  WRITE(30) xGrid; WRITE(30) yGrid; WRITE(30) tGrid

  ! Save parameters to 'InitialData.raw'
  WRITE(31) Lx; WRITE(31) Ly; WRITE(31) dx; WRITE(31) dy
  WRITE(31) nM; WRITE(31) nN; WRITE(31) mM; WRITE(31) mN

  WRITE(31) xGrid; WRITE(31) yGrid
  
  IF(dealias) THEN
    CALL DealiasFilter
  END IF
  
  !** Calculate dumdum energy
  energy = 0
  DO j = 0, mN
    DO i = 0, mM
      energy = energy + 0.5*h(i,j)**2 &
              + 0.5*(1 + epsilon*h(i,j))*(u(i,j)**2 + v(i,j)**2)
    END DO
  END DO
  energy = energy/(nM*nN)
  
  WRITE(30) u; WRITE(30) v; WRITE(30) h
  
  PRINT*, "Energy:", energy
  
  DO t = 1, nT
    ! NOTE: At timestep t, we are computing u(t) from u(t-1) and u(t-2)
    !  (except at t=1)

    !--------------------------------------
    ! Calculate the Forcing along the grid
    !  - Determine (q,p,m) forcing
    !  - FFT
    !  - Calculate (u,v,h) forcing
    !  - IFFT to spatial forcing
    !--------------------------------------
    DO j = 0, mN
      DO i = 0, mM
        Fq(i,j) = FqEqn(i*dx, j*dy, t*dt)
        Fp(i,j) = FpEqn(i*dx, j*dy, t*dt)
      END DO
    END DO

    CALL GetForcing

!
!   Use this code if you want to force (u,v,h) instead of (q,p,m)

!    DO j = 0, mN
!      DO i = 0, mM
!        Fu(i,j) = 0.0
!        Fv(i,j) = 0.0
!        Fh(i,j) = DEXP(-(i*dx-5.0)**2) * DSIN(10*t*dt)
!      END DO
!    END DO

    IF(t == 1) THEN
      !------------------------------------
      ! [[RK2 Timestep]]
      !------------------------------------
      ! Incoming: u = u(0)
      !------------------------------------

      !------------------------------------
      ! uLast <- u(0)
      !------------------------------------
      DO j = 0, mN
        DO i = 0, mM
          uLast(i,j) = u(i,j)
          vLast(i,j) = v(i,j)
          hLast(i,j) = h(i,j)
        END DO
      END DO
      
      !------------------------------------
      ! uTemp <- F(u(0))
      !------------------------------------      
      CALL GetDerivatives
      
      DO j = 0, mN
        DO i = 0, mM
          uTemp(i,j) =  vLast(i,j) - hx(i,j) &
                      - epsilon*(uLast(i,j)*ux(i,j) + vLast(i,j)*uy(i,j)) &
                      + Fu(i,j)
          vTemp(i,j) = -uLast(i,j) - hy(i,j) &
                      - epsilon*(uLast(i,j)*vx(i,j) + vLast(i,j)*vy(i,j)) &
                      + Fv(i,j)
          hTemp(i,j) = -(ux(i,j) + vy(i,j)) &
                      - epsilon*(uLast(i,j)*hx(i,j) + vLast(i,j)*hy(i,j) &
                                 + hLast(i,j)*(ux(i,j) + vy(i,j))) &
                      + Fh(i,j)
        END DO
      END DO
  
      !------------------------------------
      !  u <- u(0) + dt*F(u(0)) == u(midpt)
      !------------------------------------
      DO j = 0, mN
        DO i = 0, mM
          u(i,j) = uLast(i,j) + dt*uTemp(i,j)
          v(i,j) = vLast(i,j) + dt*vTemp(i,j)
          h(i,j) = hLast(i,j) + dt*hTemp(i,j)
        END DO
      END DO
      
      IF(dealias) THEN
        CALL DealiasFilter
      END IF
  
      !------------------------------------
      ! uF2 <- F(u')
      !   u <- u(1) == u(0) + dt*(1/2)*(F(u(0)) + F(u(midpt)))
      !------------------------------------
      CALL GetDerivatives
  
      energy = 0
      DO j = 0, mN
        DO i = 0, mM
          uF2(i,j) =  v(i,j) - hx(i,j) &
                    - epsilon*(u(i,j)*ux(i,j) + v(i,j)*uy(i,j)) &
                    + Fu(i,j)
          vF2(i,j) = -u(i,j) - hy(i,j) &
                    - epsilon*(u(i,j)*vx(i,j) + v(i,j)*vy(i,j)) &
                    + Fv(i,j)
          hF2(i,j) = -(ux(i,j) + vy(i,j)) &
                    - epsilon*(u(i,j)*hx(i,j) + v(i,j)*hy(i,j) &
                               + h(i,j)*(ux(i,j) + vy(i,j))) &
                    + Fh(i,j)
                      
          u(i,j) = uLast(i,j) + 0.5*dt*(uTemp(i,j) + uF2(i,j))
          v(i,j) = vLast(i,j) + 0.5*dt*(vTemp(i,j) + vF2(i,j))
          h(i,j) = hLast(i,j) + 0.5*dt*(hTemp(i,j) + hF2(i,j))

          energy = energy + 0.5*h(i,j)**2 &
                  + 0.5*(1 + epsilon*h(i,j))*(u(i,j)**2 + v(i,j)**2)
        END DO
      END DO
      energy = energy/(nM*nN)

      IF(dealias) THEN
        CALL DealiasFilter
      END IF
      
      !------------------------------------
      ! We leave timestep t=0 with the following:
      !   uLast = u(0)
      !   uTemp = F(u(0))
      !       u = u(1)
      !------------------------------------
      
    ELSE
      !------------------------------------
      ! [[Leapfrog Trapezoidal Timestep]]
      !------------------------------------
      ! We should be coming in with these values at t:
      !      u = u(t-1)
      !  uLast = u(t-2)
      !------------------------------------
      
      !------------------------------------
      ! uTemp <- u(t-1)
      !------------------------------------
      DO j = 0, mN
        DO i = 0, mM
          uTemp(i,j) = u(i,j)
          vTemp(i,j) = v(i,j)
          hTemp(i,j) = h(i,j)
        END DO
      END DO
    
      !------------------------------------
      !  uF2 = F(u(t-1))
      !------------------------------------
      CALL GetDerivatives
    
      energy = 0
      DO j = 0, mN
        DO i = 0, mM
          uF2(i,j) =  vTemp(i,j) - hx(i,j) &
                    - epsilon*(uTemp(i,j)*ux(i,j) + vTemp(i,j)*uy(i,j)) &
                    + Fu(i,j)
                  
          vF2(i,j) = -uTemp(i,j) - hy(i,j) &
                    - epsilon*(uTemp(i,j)*vx(i,j) + vTemp(i,j)*vy(i,j)) &
                    + Fv(i,j)
          
          hF2(i,j) = -(ux(i,j) + vy(i,j)) &
                    - epsilon*(hTemp(i,j)*(ux(i,j) + vy(i,j))       &
                               + uTemp(i,j)*hx(i,j) + vTemp(i,j)*hy(i,j)) &
                    + Fh(i,j)
        END DO
      END DO
      
      !------------------------------------
      ! u(t) = u(t-2) + 2*dt*F(u(t-1))
      !------------------------------------
      DO j = 0, mN
        DO i = 0, mM
          u(i,j) = uLast(i,j) + 2.0*dt*uF2(i,j)
          v(i,j) = vLast(i,j) + 2.0*dt*vF2(i,j)
          h(i,j) = hLast(i,j) + 2.0*dt*hF2(i,j)
        END DO
      END DO
      
      IF(dealias) THEN
        CALL DealiasFilter
      END IF
          
      !------------------------------------
      ! uLast <- u(t-1)
      !------------------------------------
      DO i = 0, mM
        DO j = 0, mN
          uLast(i,j) = uTemp(i,j)
          vLast(i,j) = vTemp(i,j)
          hLast(i,j) = hTemp(i,j)
        END DO
      END DO
                
      !------------------------------------
      ! uTemp <- F(u(t-1))
      !------------------------------------
      CALL GetDerivatives
            
      energy = 0
      DO j = 0, mN
        DO i = 0, mM
          uTemp(i,j) =  v(i,j) - hx(i,j) &
                      - epsilon*(u(i,j)*ux(i,j) + v(i,j)*uy(i,j)) &
                      + Fu(i,j)
                  
          vTemp(i,j) = -u(i,j) - hy(i,j) &
                      - epsilon*(u(i,j)*vx(i,j) + v(i,j)*vy(i,j)) &
                      + Fv(i,j)
          
          hTemp(i,j) = -(ux(i,j) + vy(i,j)) &
                      - epsilon*(h(i,j)*(ux(i,j) + vy(i,j)) &
                                 + u(i,j)*hx(i,j) + v(i,j)*hy(i,j)) &
                      + Fh(i,j)
        END DO
      END DO
      
      !------------------------------------
      !     u <- u(t-1) + (1/2)*dt*( F(u(t-1)) + F(u(*)) )
      !
      ! [Again I use u(*), and again I hope FFTW preserves u(*)]
      !------------------------------------
      DO j = 0, mN
        DO i = 0, mM
          u(i,j) = uLast(i,j) + 0.5*dt*(uTemp(i,j) + uF2(i,j))
          v(i,j) = vLast(i,j) + 0.5*dt*(vTemp(i,j) + vF2(i,j))
          h(i,j) = hLast(i,j) + 0.5*dt*(hTemp(i,j) + hF2(i,j))

          energy = energy + 0.5*h(i,j)**2 &
                  + 0.5*(1 + epsilon*h(i,j))*(u(i,j)**2 + v(i,j)**2)
        END DO
      END DO
      energy = energy/(nM*nN)
      
      IF(dealias) THEN
        CALL DealiasFilter
      END IF
      
      !------------------------------------
      ! Leave t = n with these values:
      !     u <- u(n)
      ! uLast <- u(n-1)
      !------------------------------------

    ! End LFT Timestep
    END IF
    
    
    !------------------------------------
    ! Output Results
    !------------------------------------
    ! There may be some bugs in here
    !*********************************************************************!
    !*** This won't work if BuffSize and DataRate don't divide into nT ***!
    !*********************************************************************!
    
    tFileWrite = MOD(t, BuffSize*DataRate)
    tFileBase = (t-1)/(DataRate*BuffSize)
    
    tBufferWrite = MOD(t, DataRate)
    tBufferIndex = MOD((t-1)/DataRate, BuffSize)
    
    ! (There may be a faster way to do this in F90)
    IF(tBufferWrite == 0) THEN
      DO j = 0, mN
        DO i = 0, mM
          dataBuffer(i,j,tBufferIndex,1) = u(i,j)
          dataBuffer(i,j,tBufferIndex,2) = v(i,j)
          dataBuffer(i,j,tBufferIndex,3) = h(i,j)
        END DO
      END DO
    END IF
    
    IF(tFileWrite == 0) THEN
      PRINT*, "Writing at t =", t
      WRITE(30) dataBuffer
    END IF
    
    IF(MOD(100*t, 10*nT) == 0) THEN
      PRINT*, "Percentage Completed:", 100*t/nT
      PRINT*, "Energy:", energy
    END IF
    
    ! Need a better blowup check here
    IF(energy > 1e3) THEN
      PRINT*, "BLOWUP"
      PRINT*, "Energy:", energy
      PRINT*, "t:", t
      EXIT
    END IF
    
    ! End timestep
  END DO
  
  ! Write final state to 'LastData.raw'
  WRITE(31) u; WRITE(31) v; WRITE(31) h
  
  !---------------
  ! CLEANUP
  !---------------
  
  ! Deallocate FFTW Plans
  CALL dfftw_destroy_plan(uPlan)
  CALL dfftw_destroy_plan(vPlan)
  CALL dfftw_destroy_plan(hPlan)

  CALL dfftw_destroy_plan(uxPlan)
  CALL dfftw_destroy_plan(vxPlan)
  CALL dfftw_destroy_plan(hxPlan)

  CALL dfftw_destroy_plan(uyPlan)
  CALL dfftw_destroy_plan(vyPlan)
  CALL dfftw_destroy_plan(hyPlan)

  CALL dfftw_destroy_plan(uFilterPlan)
  CALL dfftw_destroy_plan(vFilterPlan)
  CALL dfftw_destroy_plan(hFilterPlan)

  CALL dfftw_destroy_plan(FqPlan)
  CALL dfftw_destroy_plan(FpPlan)
  CALL dfftw_destroy_plan(FuPlan)
  CALL dfftw_destroy_plan(FvPlan)
  CALL dfftw_destroy_plan(FhPlan)
  
  ! Close 'ShallowData.raw'
  CLOSE(30)
  CLOSE(31)
  
!! ***RETURN TO SHELL***
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
CONTAINS

! ** SUBROUTINE GetDerivatives **
! Calculate derivatives by collocation
!  * Preserves values of u, v, h and clobbers initial contents of ux, vx,...
  SUBROUTINE GetDerivatives
    USE ShallowWaterSystem
    
    ! Indexing variables
    INTEGER :: i, j
    
    ! (1) Get u, v, h in spectral form
    CALL dfftw_execute(uPlan)
    CALL dfftw_execute(vPlan)
    CALL dfftw_execute(hPlan)    
    
    ! (2) Calculate spectral derivatives
    !  * I preadjust for the unnormalized FFT weights here
    DO j = 0, mN
      DO i = 0, nM/2
        uxF(i,j) = (0.0, 1.0) * kGrid(i) * uF(i,j) / (nM*nN)
        uyF(i,j) = (0.0, 1.0) * lGrid(j) * uF(i,j) / (nM*nN)
        
        vxF(i,j) = (0.0, 1.0) * kGrid(i) * vF(i,j) / (nM*nN)
        vyF(i,j) = (0.0, 1.0) * lGrid(j) * vF(i,j) / (nM*nN)

        hxF(i,j) = (0.0, 1.0) * kGrid(i) * hF(i,j) / (nM*nN)
        hyF(i,j) = (0.0, 1.0) * lGrid(j) * hF(i,j) / (nM*nN)
      END DO
    END DO
    
    ! (3) Convert derivatives back to spatial form
    CALL dfftw_execute(uxPlan)
    CALL dfftw_execute(uyPlan)

    CALL dfftw_execute(vxPlan)
    CALL dfftw_execute(vyPlan)

    CALL dfftw_execute(hxPlan)
    CALL dfftw_execute(hyPlan)
    
  END SUBROUTINE
  
!!!!!!!!!!!!!!!

  SUBROUTINE GetForcing
    USE ShallowWaterSystem

    ! Indexing variables
    INTEGER :: i, j

    ! Physical variables
    REAL(8) :: omega, Km
    COMPLEX(8) :: FupvF, FvpvF, FhpvF
    COMPLEX(8) :: FugpF, FvgpF, FhgpF
    COMPLEX(8) :: FugmF, FvgmF, FhgmF

    ! - Convert (Fq, Fp, Fm) to the spectral basis
    CALL dfftw_execute(FqPlan)
    CALL dfftw_execute(FpPlan)

    ! - Calculate (Fu, Fv, Fh)
    !  * I preadjust for the unnormalized FFT weights here
    DO j = 0, mN
      DO i = 0, nM/2
        Km = SQRT(kGrid(i)**2 + lGrid(j)**2)
        omega = SQRT(1 + Km**2)

        FupvF =  (0,1)*lGrid(j)*FqF(i,j) / omega
        FvpvF = -(0,1)*kGrid(i)*FqF(i,j) / omega
        FhpvF = -FqF(i,j) / omega

        FugpF = ( omega*kGrid(i) + (0,1)*lGrid(j)) * FpF(i,j) &
                   / (SQRT(2d0)*omega*Km)
        FvgpF = (-(0,1)*kGrid(i) + omega*lGrid(j)) * FpF(i,j) &
                   / (SQRT(2d0)*omega*Km)
        FhgpF = Km * FpF(i,j) / (SQRT(2d0)*omega)

        FugmF = (-omega*kGrid(i) + (0,1)*lGrid(j)) * CONJG(FpF(nM-i,nN-j)) &
                   / (SQRT(2d0)*omega*Km)
        FvgmF = (-(0,1)*kGrid(i) - omega*lGrid(j)) * CONJG(FpF(nM-i,nN-j)) &
                   / (SQRT(2d0)*omega*Km)
        FhgmF = Km * CONJG(FpF(nM-i,nN-j)) / (SQRT(2d0)*omega)

        FuF(i,j) = (FupvF + FugpF + FugmF)/(nM*nN)
        FvF(i,j) = (FvpvF + FvgpF + FvgmF)/(nM*nN)
        FhF(i,j) = (FhpvF + FhgpF + FhgmF)/(nM*nN)
      END DO
    END DO
    
    !** Redo (k,l) = (0,0) **!   
    FupvF =  0d0
    FvpvF =  0d0
    FhpvF = -FqF(0,0)

    FugpF = FpF(0,0) / SQRT(2d0)
    FvgpF = -(0,1)*FpF(0,0) / SQRT(2d0)
    FhgpF = 0d0

    FugmF = -(-CONJG(FpF(0,0))) / SQRT(2d0)
    FvgmF = -(0,1)*(-CONJG(FpF(0,0))) / SQRT(2d0)
    FhgmF = 0d0
  
    FuF(0,0) = (FupvF + FugpF + FugmF)/(nM*nN)
    FvF(0,0) = (FvpvF + FvgpF + FvgmF)/(nM*nN)
    FhF(0,0) = (FhpvF + FhgpF + FhgmF)/(nM*nN)

    ! - Convert (Fu, Fv, Fh) to the spatial basis
    CALL dfftw_execute(FuPlan)
    CALL dfftw_execute(FvPlan)
    CALL dfftw_execute(FhPlan)

  END SUBROUTINE GetForcing   

!!!!!!!!!!!!!!!

  SUBROUTINE DealiasFilter
    USE ShallowWaterSystem
    
    INTEGER :: i, j
    
    CALL dfftw_execute(uPlan)
    CALL dfftw_execute(vPlan)
    CALL dfftw_execute(hPlan)
    
    DO j = 0, mN
      DO i = 0, nM/2
        IF(i .GE. (nM/3 - 1) .OR. (j .GE. (nN/3 - 1) .AND. j .LE. (2*nN/3 + 1))) THEN
          uF(i,j) = 0
          vF(i,j) = 0
          hF(i,j) = 0
        ELSE
          uF(i,j) = uF(i,j)/(nM*nN)
          vF(i,j) = vF(i,j)/(nM*nN)
          hF(i,j) = hF(i,j)/(nM*nN)
        END IF
      END DO
    END DO
    
    CALL dfftw_execute(uFilterPlan)
    CALL dfftw_execute(vFilterPlan)
    CALL dfftw_execute(hFilterPlan)
    
  END SUBROUTINE DealiasFilter

!!!!!!!!!!!!!!!

END PROGRAM ShallowWaterModel

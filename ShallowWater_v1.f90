!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
MODULE ShallowWaterSystem
  IMPLICIT NONE
  
  ! Spatial System Parameters (Stored in "InitialData.nc")
  REAL(8) :: Lx, Ly                       ! Physical Domain Length in x and y
  INTEGER :: nM, nN                       ! Num of Collocation Pts in x and y
  !   (Derived quantities)
  REAL(8) :: dx, dy                       ! Physical x- and y-spacing
  INTEGER :: mM, mN                       ! Maximum Mode Number in x and y  

  ! Control Parameters
  LOGICAL, PARAMETER :: dealias = .TRUE.
  
  ! Storage Paramters
  INTEGER, PARAMETER :: BuffSize = 1e2
  INTEGER, PARAMETER :: DataRate = 1e2
  
  ! Temporal System Parameters
  REAL(8), PARAMETER :: dt = 1D-3         ! Integration Timestep
  INTEGER, PARAMETER :: nT = 4e4          ! Number of Timesteps
  REAL(8), PARAMETER :: Tmax = dt*nT      ! Integration Time

  ! Physical Scaling Parameters
  ! (I'm going to phase these out, they are redundant)
  REAL(8), PARAMETER :: epsilon = 0.1     ! Rossby Number
  
  ! Numerical Parameters
  ! This is redundant, but still present; dealiasing achieves the same effect
  REAL(8), PARAMETER :: visc = 0e-5      ! Numerical Viscosity
  INTEGER, PARAMETER :: vN = 1           ! Viscosity Order, (nabla)^2n
  
  ! Physical System Variables
  REAL(8), POINTER :: u(:,:), ux(:,:), uy(:,:), Du(:,:)     ! x-velocity
  REAL(8), POINTER :: v(:,:), vx(:,:), vy(:,:), Dv(:,:)     ! y-velocity
  REAL(8), POINTER :: h(:,:), hx(:,:), hy(:,:), Dh(:,:)     ! Surface Ht
  
  ! Intermediary FFTW Variables
  COMPLEX(8), POINTER :: uF(:,:), uxF(:,:), uyF(:,:), DuF(:,:) ! Fourier u data
  COMPLEX(8), POINTER :: vF(:,:), vxF(:,:), vyF(:,:), DvF(:,:) ! Fourier v data
  COMPLEX(8), POINTER :: hF(:,:), hxF(:,:), hyF(:,:), DhF(:,:) ! Fourier h data
  
  ! FFTW Control Variables
  INTEGER(8) :: uPlan, uxPlan, uyPlan, DuPlan        ! FFTW u plans
  INTEGER(8) :: vPlan, vxPlan, vyPlan, DvPlan        ! FFTW v plans
  INTEGER(8) :: hPlan, hxPlan, hyPlan, DhPlan        ! FFTW h plans
  
  INTEGER(8) :: uFilterPlan, vFilterPlan, hFilterPlan

  ! NetCDF Control Variables
  INTEGER :: ncInputID, ncOutputID, &
             xDimID, yDimID, tDimID, &
             uVInputID, vVInputID, hVInputID, &
             LxVInputID, LyVInputID, nMVInputID, nNVInputID, &
             xVarID, yVarID, tVarID, uVarID, vVarID, hVarID, &
             epsilonVarID, dataRateVarID, &
             LxVarID, LyVarID, dxVarID, dyVarID, dtVarID, &
             nMVarID, nNVarID, ntVarID, mMVarID, mNVarID, &
             energyVarID

  CHARACTER (len = *), PARAMETER :: inputName  = "InitialData.nc"
  CHARACTER (len = *), PARAMETER :: outputName = "ShallowData.nc"
  
  ! Numerical Parameters
  REAL(4),  PARAMETER :: FPI = 3.14159
  REAL(8),  PARAMETER :: DPI = 3.14159265358979D+0
  ! REAL(16), PARAMETER :: QPI = 3.14159265358979323846264338327950Q+0
  
END MODULE ShallowWaterSystem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ShallowWaterModel
!
! Program Description:
!   Weak attempt at a pseudospectral shallow water model, with a minimum of
!   cheats.
!
! --Marshall Ward
!
! Log:
! 26 January   2006 - Begin Project
! 27 February  2006 - Began implementation of FFTW
!  9 March     2006 - Finished first draft
!    May       2006 - Really trucking along here...
! 18 September 2006 - Semi-functional, exploring different timestepping methods
!  2 October   2006 - Nonlinear instability
!    December  2006 - Producing physically meaningful results

  ! Modules
  USE typeSizes
  USE netcdf
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
  REAL(8), POINTER :: uBuff(:,:,:), vBuff(:,:,:), hBuff(:,:,:), eBuff(:)
  INTEGER :: tFileWrite, tFileBase, tBufferWrite, tBufferIndex
  
!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  PRINT*, "Begin ShallowWater"
  
  !** Load Input File **!
  CALL check(nf90_open(trim(inputName), nf90_nowrite, ncInputID))
  
  !**** Get variable IDs
  CALL check(nf90_inq_varid(ncInputID, "u", uVInputID))
  CALL check(nf90_inq_varid(ncInputID, "v", vVInputID))
  CALL check(nf90_inq_varid(ncInputID, "h", hVInputID))

  CALL check(nf90_inq_varid(ncInputID, "Lx", LxVInputID))
  CALL check(nf90_inq_varid(ncInputID, "Ly", LyVInputID))
  
  CALL check(nf90_inq_varid(ncInputID, "nM", nMVInputID))
  CALL check(nf90_inq_varid(ncInputID, "nN", nNVInputID))

  !**** Get the data size variables
  CALL check(nf90_get_var(ncInputID, LxVInputID, Lx))
  CALL check(nf90_get_var(ncInputID, LyVInputID, Ly))
  CALL check(nf90_get_var(ncInputID, nMVInputID, nM))
  CALL check(nf90_get_var(ncInputID, nNVInputID, nN))
  
  !**** Compute derived variables
  dx = Lx / nM;  dy = Ly / nN
  mM = nM - 1;   mN = nN - 1
  
  !**** Allocate variable memory
  ALLOCATE( u(0:mM, 0:mN)); ALLOCATE( uF(0:nM/2, 0:mN));
  ALLOCATE( v(0:mM, 0:mN)); ALLOCATE( vF(0:nM/2, 0:mN));
  ALLOCATE( h(0:mM, 0:mN)); ALLOCATE( hF(0:nM/2, 0:mN));

  ALLOCATE(ux(0:mM, 0:mN)); ALLOCATE(uxF(0:nM/2, 0:mN));
  ALLOCATE(vx(0:mM, 0:mN)); ALLOCATE(vxF(0:nM/2, 0:mN));
  ALLOCATE(hx(0:mM, 0:mN)); ALLOCATE(hxF(0:nM/2, 0:mN));

  ALLOCATE(uy(0:mM, 0:mN)); ALLOCATE(uyF(0:nM/2, 0:mN));
  ALLOCATE(vy(0:mM, 0:mN)); ALLOCATE(vyF(0:nM/2, 0:mN));
  ALLOCATE(hy(0:mM, 0:mN)); ALLOCATE(hyF(0:nM/2, 0:mN));

  ALLOCATE(Du(0:mM, 0:mN)); ALLOCATE(DuF(0:nM/2, 0:mN));
  ALLOCATE(Dv(0:mM, 0:mN)); ALLOCATE(DvF(0:nM/2, 0:mN));
  ALLOCATE(Dh(0:mM, 0:mN)); ALLOCATE(DhF(0:nM/2, 0:mN));
  
  ALLOCATE(uLast(0:mM, 0:mN));
  ALLOCATE(uTemp(0:mM, 0:mN));
  ALLOCATE(  uF2(0:mM, 0:mN));

  ALLOCATE(vLast(0:mM, 0:mN));
  ALLOCATE(vTemp(0:mM, 0:mN));
  ALLOCATE(  vF2(0:mM, 0:mN));

  ALLOCATE(hLast(0:mM, 0:mN));
  ALLOCATE(hTemp(0:mM, 0:mN));
  ALLOCATE(  hF2(0:mM, 0:mN));

  ALLOCATE(uBuff(0:mM, 0:mN, 0:BuffSize-1));
  ALLOCATE(vBuff(0:mM, 0:mN, 0:BuffSize-1));
  ALLOCATE(hBuff(0:mM, 0:mN, 0:BuffSize-1));
  ALLOCATE(eBuff(0:BuffSize-1));
  
  PRINT*, "Memory Allocated"
  
  !** Create FFTW Plans **!
  CALL dfftw_plan_dft_r2c_2d(uPlan,  nM, nN, u,   uF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(uxPlan, nM, nN, uxF, ux, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(uyPlan, nM, nN, uyF, uy, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(DuPlan, nM, nN, DuF, Du, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  CALL dfftw_plan_dft_r2c_2d(vPlan,  nM, nN, v,   vF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vxPlan, nM, nN, vxF, vx, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vyPlan, nM, nN, vyF, vy, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(DvPlan, nM, nN, DvF, Dv, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  CALL dfftw_plan_dft_r2c_2d(hPlan,  nM, nN, h,   hF, FFTW_FORWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hxPlan, nM, nN, hxF, hx, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hyPlan, nM, nN, hyF, hy, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(DhPlan, nM, nN, DhF, Dh, FFTW_BACKWARD, &
                             FFTW_MEASURE)

  !** Dealias FFTW Plans **!
  CALL dfftw_plan_dft_c2r_2d(uFilterPlan, nM, nN, uF, u, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(vFilterPlan, nM, nN, vF, v, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  CALL dfftw_plan_dft_c2r_2d(hFilterPlan, nM, nN, hF, h, FFTW_BACKWARD, &
                             FFTW_MEASURE)
  
  PRINT*, "FFTW Initialized"

  !**** Fill (u,v,h) variables
  CALL check(nf90_get_var(ncInputID, uVInputID, u, (/ 1, 1 /), (/ nM, nN /) ))
  CALL check(nf90_get_var(ncInputID, vVInputID, v, (/ 1, 1 /), (/ nM, nN /) ))
  CALL check(nf90_get_var(ncInputID, hVInputID, h, (/ 1, 1 /), (/ nM, nN /) ))
  
  ! (1) Initialize Output File
  ! (1a) Create Output File
  CALL check(nf90_create(trim(outputName), &
                         or(nf90_clobber, nf90_64bit_offset), &
  !*(alt cmode)*         cmode = nf90_clobber, &
                         ncOutputID))
             
  ! (1b) Define Dimensions
  CALL check(nf90_def_dim(ncOutputID, "x", nM, xDimID))
  CALL check(nf90_def_dim(ncOutputID, "y", nN, yDimID))  
  CALL check(nf90_def_dim(ncOutputID, "t", (nT/DataRate)+1, tDimID))
  
  ! (1c) Define Variables
  CALL check(nf90_def_var(ncOutputID, "x", nf90_double, xDimID, xVarID))  
  CALL check(nf90_def_var(ncOutputID, "y", nf90_double, yDimID, yVarID))
  CALL check(nf90_def_var(ncOutputID, "t", nf90_double, tDimID, tVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "u", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID, tDimID /), &
                          varID = uVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "v", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID, tDimID /), &
                          varID = vVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "h", xtype = nf90_double, &
                          dimids = (/ xDimID, yDimID, tDimID /), &
                          varID = hVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "energy", &
                          xtype = nf90_double, &
                          dimids = (/ tDimID /), &
                          varID = energyVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "epsilon", &
                          xtype = nf90_double, &
                          varID = epsilonVarID))

  CALL check(nf90_def_var(ncid = ncOutputID, name = "Lx", xtype = nf90_double, &
                          varID = lxVarID))

  CALL check(nf90_def_var(ncid = ncOutputID, name = "Ly", xtype = nf90_double, &
                          varID = lyVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "dx", xtype = nf90_double, &
                          varID = dxVarID))

  CALL check(nf90_def_var(ncid = ncOutputID, name = "dy", xtype = nf90_double, &
                          varID = dyVarID))

  CALL check(nf90_def_var(ncid = ncOutputID, name = "dt", xtype = nf90_double, &
                          varID = dtVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "nM", xtype = nf90_int, &
                          varID = nMVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "nN", xtype = nf90_int, &
                          varID = nNVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "nT", xtype = nf90_int, &
                          varID = nTVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "mM", xtype = nf90_int, &
                          varID = mMVarID))
  
  CALL check(nf90_def_var(ncid = ncOutputID, name = "mN", xtype = nf90_int, &
                          varID = mNVarID))
  
  CALL check(nf90_def_var(ncOutputID, "dataRate", nf90_int, dataRateVarID))

  CALL check(nf90_enddef(ncOutputID))
  
  PRINT*, "Output File Created"
  
  ! (1d) Fill in Grid Values (i.e. Dimension Grid)
  DO i = 1, nM
    CALL check(nf90_put_var(ncOutputID, xVarID, (i-1)*dx, start = (/ i /) ))
  END DO

  DO i = 1, nN
    CALL check(nf90_put_var(ncOutputID, yVarID, (i-1)*dy, start = (/ i /) ))
  END DO

  DO i = 1, (nT/DataRate)+1
    CALL check(nf90_put_var(ncOutputID, tVarID, (i-1)*dt*DataRate, &
                            start = (/ i /) ))
  END DO

  ! (1d') Fill in Parameters
  CALL check(nf90_put_var(ncOutputID, epsilonVarID, epsilon))
  CALL check(nf90_put_var(ncOutputID, dataRateVarID, dataRate))
  
  CALL check(nf90_put_var(ncOutputID, lxVarID, Lx))
  CALL check(nf90_put_var(ncOutputID, lyVarID, Ly))
  
  CALL check(nf90_put_var(ncOutputID, dxVarID, dx))
  CALL check(nf90_put_var(ncOutputID, dyVarID, dy))
  CALL check(nf90_put_var(ncOutputID, dtVarID, dt))

  CALL check(nf90_put_var(ncOutputID, nMVarID, nM))
  CALL check(nf90_put_var(ncOutputID, nNVarID, nN))
  CALL check(nf90_put_var(ncOutputID, ntVarID, nT))

  CALL check(nf90_put_var(ncOutputID, mMVarID, mM))
  CALL check(nf90_put_var(ncOutputID, mNVarID, mN))

  IF(dealias) THEN
    CALL DealiasFilter
  END IF
  
  !** Calculate energy
  energy = 0
  DO j = 0, mN
    DO i = 0, mM
      energy = energy + 0.5*h(i,j)**2 &
              + 0.5*(1 + epsilon*h(i,j))*(u(i,j)**2 + v(i,j)**2)
    END DO
  END DO
  energy = energy/(nM*nN)
  
  CALL check(nf90_put_var(ncOutputID, uVarID, u, &
                            start = (/ 1, 1, 1 /),   &
                            count = (/ nM, nN, 1 /) ))
  CALL check(nf90_put_var(ncOutputID, vVarID, v, &
                            start = (/ 1, 1, 1 /),   &
                            count = (/ nM, nN, 1 /) ))
  CALL check(nf90_put_var(ncOutputID, hVarID, h, &
                            start = (/ 1, 1, 1 /),   &
                            count = (/ nM, nN, 1 /) ))
  CALL check(nf90_put_var(ncOutputID, energyVarID, energy, &
                            start = (/ 1 /) ))  
  
  PRINT*, "Energy:", energy
  
  DO t = 1, nT
    ! NOTE: At timestep t, we are computing u(t) from u(t-1) and u(t-2)

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
                      + Du(i,j)
          vTemp(i,j) = -uLast(i,j) - hy(i,j) &
                      - epsilon*(uLast(i,j)*vx(i,j) + vLast(i,j)*vy(i,j)) &
                      + Dv(i,j)
          hTemp(i,j) = -(ux(i,j) + vy(i,j)) &
                      - epsilon*(uLast(i,j)*hx(i,j) + vLast(i,j)*hy(i,j) &
                              + hLast(i,j)*(ux(i,j) + vy(i,j))) + Dh(i,j)
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
      !
      ! [verify that FFTW in GetDerivatives doesn't clobber u]
      !------------------------------------
      CALL GetDerivatives
  
      energy = 0
      DO j = 0, mN
        DO i = 0, mM
          uF2(i,j) =  v(i,j) - hx(i,j) &
                    - epsilon*(u(i,j)*ux(i,j) + v(i,j)*uy(i,j)) + Du(i,j)
          vF2(i,j) = -u(i,j) - hy(i,j) &
                    - epsilon*(u(i,j)*vx(i,j) + v(i,j)*vy(i,j)) + Dv(i,j)
          hF2(i,j) = -(ux(i,j) + vy(i,j)) &
                    - epsilon*(u(i,j)*hx(i,j) + v(i,j)*hy(i,j) &
                          + h(i,j)*(ux(i,j) + vy(i,j))) + Dh(i,j)
                      
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
                    + Du(i,j)
                  
          vF2(i,j) = -uTemp(i,j) - hy(i,j) &
                    - epsilon*(uTemp(i,j)*vx(i,j) + vTemp(i,j)*vy(i,j)) &
                    + Dv(i,j)
          
          hF2(i,j) = -(ux(i,j) + vy(i,j)) &
                    - epsilon*(hTemp(i,j)*(ux(i,j) + vy(i,j))       &
                               + uTemp(i,j)*hx(i,j) + vTemp(i,j)*hy(i,j)) &
                    + Dh(i,j)
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
                      + Du(i,j)
                  
          vTemp(i,j) = -u(i,j) - hy(i,j) &
                      - epsilon*(u(i,j)*vx(i,j) + v(i,j)*vy(i,j)) &
                      + Dv(i,j)
          
          hTemp(i,j) = -(ux(i,j) + vy(i,j)) &
                      - epsilon*(h(i,j)*(ux(i,j) + vy(i,j)) &
                                 + u(i,j)*hx(i,j) + v(i,j)*hy(i,j)) &
                      + Dh(i,j)
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
    
    !*****************************************************************!
    !*** This won't work if BuffSize and DataRate don't go into nT ***!
    !*****************************************************************!
    
    tFileWrite = MOD(t, BuffSize*DataRate)
    tFileBase = (t-1)/(DataRate*BuffSize)
    
    tBufferWrite = MOD(t, DataRate)
    tBufferIndex = MOD((t-1)/DataRate, BuffSize)
    
    ! (There may be a faster way to do this in F90)
    IF(tBufferWrite == 0) THEN
      DO j = 0, mN
        DO i = 0, mM
          uBuff(i,j,tBufferIndex) = u(i,j)
          vBuff(i,j,tBufferIndex) = v(i,j)
          hBuff(i,j,tBufferIndex) = h(i,j)
        END DO
      END DO
      eBuff(tBufferIndex) = energy
    END IF
    
    IF(tFileWrite == 0) THEN
      PRINT*, "Writing at t =", t
      CALL check(nf90_put_var(ncOutputID, uVarID, uBuff, &
                              start = (/ 1, 1, tFileBase*BuffSize+2 /), &
                              count = (/ nM, nN, BuffSize /) ))
      CALL check(nf90_put_var(ncOutputID, vVarID, vBuff, &
                              start = (/ 1, 1, tFileBase*BuffSize+2 /), &
                              count = (/ nM, nN, BuffSize /) ))
      CALL check(nf90_put_var(ncOutputID, hVarID, hBuff, &
                              start = (/ 1, 1, tFileBase*BuffSize+2 /), &
                              count = (/ nM, nN, BuffSize /) ))
      CALL check(nf90_put_var(ncOutputID, energyVarID, eBuff, &
                              start = (/ tFileBase*BuffSize+2 /), &
                              count = (/ BuffSize /) ))
    END IF
    
    IF(MOD(100*t, 10*nT) == 0) THEN
      PRINT*, "Percentage Completed:", 100*t/nT
      PRINT*, "Energy:", energy
    END IF
    
    IF(energy > 1e3) THEN
      PRINT*, "BLOWUP"
      PRINT*, "Energy:", energy
      PRINT*, "t:", t
      EXIT
    END IF
    
    ! End timestep
  END DO
  
  ! CLEANUP
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

  CALL dfftw_destroy_plan(DuPlan)
  CALL dfftw_destroy_plan(DvPlan)
  CALL dfftw_destroy_plan(DhPlan)

  CALL dfftw_destroy_plan(uFilterPlan)
  CALL dfftw_destroy_plan(vFilterPlan)
  CALL dfftw_destroy_plan(hFilterPlan)
  
  ! Close NetCDF Data File
  CALL check(nf90_close(ncOutputID))
  
!! ***RETURN TO SHELL***
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
CONTAINS

! ** SUBROUTINE GetDerivatives **
! Calculate derivatives by collocation
!  * Preserves values of u, v, h and clobbers initial contents of ux, vx, etc
  SUBROUTINE GetDerivatives
    USE ShallowWaterSystem
    
    !Physical Control Variables
    REAL(8) :: kmode             ! Physical wavenumber for x (k = 2pi m / Lx)
    REAL(8) :: lmode             ! Physical wavenumber for y (l = 2pi n / Ly)
        
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
        kmode = 2.0*DPI * i / Lx
        lmode = (2.0*DPI * (j - nN*(j/(nN/2 + 1)))) / Ly   ! (Check this)
      
        uxF(i,j) = (0.0, 1.0) * kmode * uF(i,j) / (nM*nN)
        uyF(i,j) = (0.0, 1.0) * lmode * uF(i,j) / (nM*nN)
        DuF(i,j) = -visc * (kmode**2 + lmode**2)**vN * uF(i,j) / (nM*nN)
        
        vxF(i,j) = (0.0, 1.0) * kmode * vF(i,j) / (nM*nN)
        vyF(i,j) = (0.0, 1.0) * lmode * vF(i,j) / (nM*nN)
        DvF(i,j) = -visc * (kmode**2 + lmode**2)**vN * vF(i,j) / (nM*nN)

        hxF(i,j) = (0.0, 1.0) * kmode * hF(i,j) / (nM*nN)
        hyF(i,j) = (0.0, 1.0) * lmode * hF(i,j) / (nM*nN)
        DhF(i,j) = -visc * (kmode**2 + lmode**2)**vN * hF(i,j) / (nM*nN)
      END DO
    END DO
    
    ! (3) Convert derivatives back to spatial form
    CALL dfftw_execute(uxPlan)
    CALL dfftw_execute(uyPlan)
    CALL dfftw_execute(DuPlan)

    CALL dfftw_execute(vxPlan)
    CALL dfftw_execute(vyPlan)
    CALL dfftw_execute(DvPlan)

    CALL dfftw_execute(hxPlan)
    CALL dfftw_execute(hyPlan)
    CALL dfftw_execute(DhPlan)
    
  END SUBROUTINE
  
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
!        IF(i .GE. nM/4 .OR. (j .GE. nN/4 .AND. j .LE. 3*nN/4)) THEN
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


  ! Internal subroutine - checks error status after each netcdf, prints out 
  ! text message each time an error code is returned. 
  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status
    
    IF(status /= nf90_noerr) THEN 
      PRINT *, trim(nf90_strerror(status))
    END IF
  END SUBROUTINE check  

!!!!!!!!!!!!!!!


END PROGRAM ShallowWaterModel

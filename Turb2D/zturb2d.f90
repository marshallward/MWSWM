!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! zTurb2D (zturb2d.f90) : Two dimensional Turbulence Simulator
!
! Spectral evolution of 2D vorticity
!
! Using Singleton's F77 FFT code. There is probably a better alternative,
!  i.e. one that assumes real data in spatial form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM Turbulence2D
  IMPLICIT NONE

  REAL, PARAMETER :: pi = 3.1415926535897932384626433832795029

!!!!!!!!!!!!!!!!!
!               !
!  DECLARATION  !
!               !
!!!!!!!!!!!!!!!!!

  LOGICAL, PARAMETER :: resume = .FALSE.
  ! .TRUE.  =  use .last file
  ! .FALSE. =  initialize with new data

  LOGICAL, PARAMETER :: alias = .TRUE.
  ! .TRUE.  = use 2/3 dealiasing
  ! .FALSE. = don't dealias

  INTEGER, PARAMETER :: Kmax = 127       ! 2x-1 the max k-mode / max x-node
  INTEGER, PARAMETER :: Lmax = 127       ! 2x-1 the max l-mode / max y-node

  INTEGER, PARAMETER :: NK = Kmax + 1    ! Number of k-modes
  INTEGER, PARAMETER :: NL = Lmax + 1    ! Number of l-modes

  INTEGER, PARAMETER :: Emax = MIN(NK/2,NL/2)  ! Maximum spectral energy mode

  REAL   , PARAMETER :: h = 1E-3          ! Timestep size
  INTEGER, PARAMETER :: Tmax = 1*1e3     ! Number of timesteps

!  REAL   , PARAMETER :: Re_T = 1e9       ! Temporal Reynolds number ((L^2/T) / viscosity)
!  REAL   , PARAMETER :: Re_L = Re_T      ! Spatial Reynolds number  (L^2 zeta / viscosity)

  REAL   , PARAMETER :: nu = 1e-6         ! nu = 1/Re_T

  COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: zeta       ! Spectral vorticity
  COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: lastZeta   ! prior zeta
  COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: tempZeta   ! temporary storage of zeta

  COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: Jcbn       ! Jacobian
  COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: lastJcbn   ! previous Jacobian

  INTEGER, DIMENSION(0:Kmax, 0:Lmax) :: iso    ! (k,l) -> (K) referencer

  INTEGER, DIMENSION(0:Kmax) :: mK       ! k mode associated with i
  INTEGER, DIMENSION(0:Lmax) :: mL       ! l mode associated with j

  REAL, DIMENSION(0:Emax) :: ESpect      ! Isotropic spectral energy density

  REAL :: InitialEnstrophy

  INTEGER :: i, j                        ! Summation indices
  REAL :: k, l                           ! a = k(i), b = l(j)
  INTEGER :: t                           ! Timestep index


!!!!!!!!!!!!!!!!!!!!
!                  !
!  INITIALIZATION  !
!                  !
!!!!!!!!!!!!!!!!!!!!

  ! Construct the isotropic referencer
  CALL BuildRefs(iso, mK, mL)

  ! Either initalize the vorticity or load up a prior vorticity state
  IF(resume) THEN
    CALL LoadVorticity(zeta)
  ELSE
    CALL InitVorticity(zeta)
  END IF

  ! Create the output files

  OPEN(16, FILE='zeta.dat', STATUS='REPLACE', ACCESS='SEQUENTIAL', &
           ACTION='WRITE')

  OPEN(18, FILE='spect.dat', STATUS='REPLACE', ACCESS='SEQUENTIAL', &
           ACTION='WRITE')

  t = 0

  CALL BuildSpectrum(zeta, iso, ESpect)

  CALL OutputVorticity(zeta)
  CALL OutputSpectrum(ESpect)

  PRINT*, "t = ", h*t, "Energy = ", Energy(zeta), "Enstrophy = ", Enstrophy(zeta)

  InitialEnstrophy = Enstrophy(zeta)

!!!!!!!!!!!!!!!
!             !
!  EVOLUTION  !
!             !
!!!!!!!!!!!!!!!

  !*** We assume Re_L = Re_T ***!

  ! Euler step to the second iteration

  t = 1

  CALL GetJacobian(zeta, Jcbn)
!  CALL Jac2(zeta, Jcbn)

  IF(Tmax > 0) THEN

  DO i = 0, Kmax
    DO j = 0, Lmax

      k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
      l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

      zeta(i,j) = zeta(i,j) - h*(Jcbn(i,j) + nu*(k**2 + l**2)**2*zeta(i,j))

      lastJcbn(i,j) = Jcbn(i,j)

    END DO
  END DO

  CALL BuildSpectrum(zeta, iso, ESpect)

  IF(MODULO(t,NINT(0.1/h)) == 0) THEN
    CALL OutputVorticity(zeta)
    CALL OutputSpectrum(ESpect)
  END IF

  PRINT*, "t = ", t*h, "Energy = ", Energy(zeta), "Enstrophy = ", Enstrophy(zeta)

  END IF

  DO t = 2, Tmax

    CALL GetJacobian(zeta,Jcbn)
!    CALL Jac2(zeta, Jcbn)

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        zeta(i,j) = ((1 - 0.5*h*nu*(k**2 + l**2)**2) * zeta(i,j)  &
                             - h*(1.5*Jcbn(i,j) - 0.5*lastJcbn(i,j)))  &
                         / (1 + 0.5*h*nu*(k**2 + l**2)**2)

        lastJcbn(i,j) = Jcbn(i,j)

      END DO
    END DO

    ! Lame abort check for blowup
    IF(ABS(Enstrophy(zeta)) > 100*InitialEnstrophy) THEN
      PRINT*, "Warning! Possible instability. Fuck Crank-Nicolson."
      EXIT
    END IF

    CALL BuildSpectrum(zeta, iso, ESpect)

    IF(MODULO(t,NINT(0.1/h)) == 0) THEN
      CALL OutputVorticity(zeta)
      CALL OutputSpectrum(ESpect)
    END IF

    PRINT*, "t = ", t*h, "Energy = ", Energy(zeta), "Enstrophy = ", Enstrophy(zeta)

  END DO

! Write 'zeta.last'
  CALL BackupVorticity(zeta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!  ****************  END MAIN PROGRAM  *****************  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! BuildRefs: The following referencers are built:
!
! Isotropic referencer:
!    Constructs the (k,l) -> (K = sqrt(k**2 + l**2)) referencer
!
! k/l mode referencers:
!    Contructs the mode referencer (i=0..N-1 -> k=0..N/2,-N/2+1..-1)
!
! INPUT : iso, mK, mL
! OUTPUT: iso, mK, mL
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BuildRefs(iso, mK, mL)

  INTEGER, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: iso
  INTEGER, DIMENSION(0:Kmax) :: mK
  INTEGER, DIMENSION(0:Lmax) :: mL

  INTEGER :: i, j

  DO i = 0, Kmax
    mK(i) = MODULO(i+(NK/2-1),NK) - (NK/2-1)
  END DO

  DO j = 0, Lmax
    mL(j) = MODULO(j+(NL/2-1),NL) - (NL/2-1)
  END DO

  DO i = 0, Kmax
    DO j = 0, Lmax
      iso(i,j) = NINT(SQRT(REAL(mK(i)**2 + mL(j)**2)))
    END DO
  END DO

  END SUBROUTINE BuildRefs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! LoadVorticity: Loads the vorticity from the restart file (zeta.last). This
!                allows us to resume the fluid evolution from a prior state.
!
! INPUT : Vorticity (zeta)
! OUTPUT: Vorticity (zeta)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE LoadVorticity(zeta)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta

    INTEGER :: i, j

    OPEN(15, FILE='zeta.start', STATUS='OLD', ACCESS='SEQUENTIAL', &
             ACTION='READ')

    DO i = 0, Kmax
      DO j = 0, Lmax
        READ(15, *) zeta(i,j)
      END DO
    END DO

  END SUBROUTINE LoadVorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! InitVorticity: Calculates an "interesting" initial vorticity contribution.
!
! INPUT : Vorticity (zeta)
! OUTPUT: Vorticity (zeta)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE InitVorticity(zeta)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: zetaRe, zetaIm

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: h1, h2
    REAL :: K2

    INTEGER :: i, j
    REAL :: x, y, k, l

    CALL RANDOM_NUMBER(h1)
    CALL RANDOM_NUMBER(h2)

  ! Psi is given in spatial form, then converted to spectral form

  DO i = 0, Kmax
    DO j = 0, Lmax

      x = MODULO(i+(NK/2-1),NK) - (NK/2-1)
      y = MODULO(j+(NL/2-1),NL) - (NL/2-1)

!       IF(j .EQ. Kmax/4)   zeta(i,j) = 10
!       IF(j .EQ. 3*Kmax/4) zeta(i,j) = -10

!      zeta(i,j) = 10*(2*h2(i,j)-1)

!      zeta(i,j) = 5 * EXP(-REAL((x/6.)**2 + (y)**2)/3.**2)

!      zeta(i,j) = 1 * EXP(-REAL((x-1)**2 + (y)**2)/5.**2) &
!                 - 1 * EXP(-REAL((x+1)**2 + (y)**2)/5.**2)

      zeta(i,j) = 10*x*y*EXP(-REAL((x**2 + y**2)/5.**2))

      zetaRe(i,j) = REAL(zeta(i,j))
      zetaIm(i,j) = AIMAG(zeta(i,j))

    END DO
  END DO

  CALL FFT(zetaRe, zetaIm, NK*NL, NK, NK, -1)
  CALL FFT(zetaRe, zetaIm, NK*NL, NL, NK*NL, -1)

  DO i = 0, Kmax
    DO j = 0, Lmax

      zeta(i,j) = (zetaRe(i,j) + (0,1)*zetaIm(i,j))/(NK*NL)
!      zeta(i,j) = 0

    END DO
  END DO

!!! *** zeta = cos(x) + cos(y) + cos(x+y) + cos(x-y) ***

!    zeta(10,0) = 1. ; zeta(Kmax-9,0) = 1.
!    zeta(0,10) = 1. ; zeta(0,Lmax-9) = 1.
!    zeta(10,10) = 1. ; zeta(Kmax-9,Lmax-9) = 1.
!    zeta(10,Lmax-9) = 1. ; zeta(Kmax-9,10) = 1.

!    zeta(1,0) = -1/2. ; zeta(Kmax-0,0) = -1/2.
!    zeta(0,1) = -1/2. ; zeta(0,Lmax-0) = -1/2.
!    zeta(1,1) = -1. ; zeta(Kmax-0,Lmax-0) = -1.
!    zeta(1,Lmax-0) = -1. ; zeta(Kmax-0,1) = -1.

!    zeta(2,0) = -1/2.*(0.707,0.707) ; zeta(Kmax-1,0) = -1/2.*(0.707,-0.707)
!    zeta(0,2) = -1/2.*(0.707,0.707) ; zeta(0,Lmax-1) = -1/2.*(0.707,-0.707)
!    zeta(2,2) = -1. ; zeta(Kmax-1,Lmax-1) = -1.
!    zeta(2,Lmax-1) = -1. ; zeta(Kmax-1,2) = -1.

  END SUBROUTINE InitVorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Energy: Calculates the energy of the system
!         (or at least proportional to it)
!
! INPUT: Vorticity (zeta)
! OUTPUT: Energy of system (Energy)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL FUNCTION Energy(zeta)
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: zeta

    INTEGER :: i, j
    REAL :: k, l

    Energy = 0

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        IF(k == 0 .AND. l == 0) CYCLE

        Energy = Energy + zeta(i,j)*zeta(MODULO(NK-i,NK), MODULO(NL-j,NL)) &
                              / (k**2 + l**2)
      END DO
    END DO

  END FUNCTION Energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Enstrophy: Calculates the enstrophy of the system
!            (or at least something proportional to it)
!
! INPUT: Vorticity (zeta)
! OUTPUT: Enstrophy of system (Enstrophy)
!
! STATUS: So far consistent
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL FUNCTION Enstrophy(zeta)
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: zeta

    INTEGER :: i, j

    Enstrophy = 0

    DO i = 0, Kmax
      DO j = 0, Lmax

        Enstrophy = Enstrophy + zeta(i,j)*zeta(MODULO(NK-i,NK), MODULO(NL-j,NL))

      END DO
    END DO

  END FUNCTION Enstrophy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! BuildSpectrum: Calculates the energy in each wavenumber magnitude
!                We ignore energy without a complete angle sweep
!                (Except the asymmetric mode, which is included for some
!                 stupid reason...)
!
! INPUT: Amplitudes (zeta), referencer (iso)
! OUTPUT: Energy spectrum (ESpect)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BuildSpectrum(zeta, iso, ESpect)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: zeta
    INTEGER, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: iso

    REAL, DIMENSION(0:Emax), INTENT(OUT) :: ESpect

    INTEGER :: i, j
    REAL :: k, l

! Initialize the spectrum
    DO i = 0, Emax
      ESpect(i) = 0
    END DO

! Compute the modes
    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2 - 1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2 - 1)

        IF(k == 0 .AND. l == 0) CYCLE 
        IF((k**2 + l**2) > Emax**2) CYCLE

        ESpect(iso(i,j)) = ESpect(iso(i,j)) &
             + zeta(i,j)*zeta(MODULO(NK-i, NK), MODULO(NL-j, NL)) &
                 / (k**2 + l**2)

      END DO
    END DO

  END SUBROUTINE BuildSpectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! OutputVorticity: Places the spatial vorticity in 'zeta.dat'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE OutputVorticity(zeta)
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: zeta

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: zetaRe, zetaIm

    INTEGER :: i, j

    DO i = 0, Kmax
      DO j = 0, Lmax

        zetaRe(i,j) = REAL(zeta(i,j))
        zetaIm(i,j) = AIMAG(zeta(i,j))

      END DO
    END DO

    CALL FFT(zetaRe, zetaIm, NK*NL, NK, NK, 1)
    CALL FFT(zetaRe, zetaIm, NK*NL, NL, NK*NL, 1)

    DO i = 0, Kmax
      DO j = 0, Lmax
        WRITE(16, "(E $)") zetaRe(MODULO(i + NK/2, NK), MODULO(j + NL/2, NL))
      END DO
      WRITE(16,*) ""
    END DO

  END SUBROUTINE OutputVorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! OutputSpectrum: Output the energy spectrum into 'spect.dat'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE OutputSpectrum(ESpect)
    REAL, DIMENSION(0:Emax), INTENT(IN) :: ESpect

    INTEGER :: i

    DO i = 0, Emax
      WRITE(18, *) ESpect(i)
    END DO

  END SUBROUTINE OutputSpectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! GetJacobian: Calculates the spectral coefficients of the Jacobian
!              by a Galerkin-type method
!
! INPUT: Vorticity (zeta)
! OUTPUT: Jacobian/Inertial Interaction (Jcbn)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GetJacobian(zeta, Jcbn)
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(OUT) :: Jcbn

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: PsiX, PsiY, Vort, J1, J2

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: PsiX_Re, PsiX_Im, PsiY_Re, PsiY_Im, &
                                       Vort_Re, Vort_Im

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: J1_Re, J1_Im, J2_Re, J2_Im

    INTEGER :: i, j
    REAL :: k, l

!!!! DeAlias
    IF(alias) THEN
      DO i = 0, Kmax
        DO j = 0, Lmax

          k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
          l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

          IF((k**2 + l**2) > (0.666*Emax)**2) zeta(i,j) = 0

        END DO
      END DO
    END IF
!!!!

    ! STEP 1: Calculate the F's in spectral form

    PsiX(0,0) = 0.
    PsiY(0,0) = 0.
    Vort(0,0) = 0.

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        IF(k == 0 .AND. l == 0) CYCLE

        PsiX(i,j) = -(0,1) * k * zeta(i,j) / (k**2 + l**2)
        PsiY(i,j) = -(0,1) * l * zeta(i,j) / (k**2 + l**2)
        Vort(i,j) = zeta(i,j)

      END DO
    END DO

    ! STEP 2: Transform F's to spatial form

      ! This is how ignorant people split a complex matrix into its parts
      ! You need to learn F90 pointers, Marshall

    DO i = 0, Kmax
      DO  j = 0, Lmax

        PsiX_Re(i,j) = REAL(PsiX(i,j))
        PsiX_Im(i,j) = AIMAG(PsiX(i,j))

        PsiY_Re(i,j) = REAL(PsiY(i,j))
        PsiY_Im(i,j) = AIMAG(PsiY(i,j))

        Vort_Re(i,j) = REAL(Vort(i,j))
        Vort_Im(i,j) = AIMAG(Vort(i,j))

      END DO
    END DO

    CALL FFT(PsiX_Re, PsiX_Im, NK*NL, NK, NK, 1)
    CALL FFT(PsiX_Re, PsiX_Im, NK*NL, NL, NK*NL, 1)

    CALL FFT(PsiY_Re, PsiY_Im, NK*NL, NK, NK, 1)
    CALL FFT(PsiY_Re, PsiY_Im, NK*NL, NL, NK*NL, 1)

    CALL FFT(Vort_Re, Vort_Im, NK*NL, NK, NK, 1)
    CALL FFT(Vort_Re, Vort_Im, NK*NL, NL, NK*NL, 1)

    DO i = 0, Kmax
      DO j = 0, Lmax

        PsiX(i,j) = PsiX_Re(i,j) + (0,1)*PsiX_Im(i,j)
        PsiY(i,j) = PsiY_Re(i,j) + (0,1)*PsiY_Im(i,j)
        Vort(i,j) = Vort_Re(i,j) + (0,1)*Vort_Im(i,j)

!!! Help nature along?
!        PsiX(i,j) = PsiX_Re(i,j)
!        PsiY(i,j) = PsiY_Re(i,j)
!        Vort(i,j) = Vort_Re(i,j)

        J1(i,j) = PsiX(i,j) * Vort(i,j)
        J2(i,j) = PsiY(i,j) * Vort(i,j)

      END DO
    END DO

    ! STEP 4: Convert the Jacobian to spectral form

    DO i = 0, Kmax
      DO j = 0, Lmax

        J1_Re(i,j) = REAL(J1(i,j))
        J1_Im(i,j) = AIMAG(J1(i,j))

        J2_Re(i,j) = REAL(J2(i,j))
        J2_Im(i,j) = AIMAG(J2(i,j))

      END DO
    END DO

    CALL FFT(J1_Re, J1_Im, NK*NL, NK, NK, -1)
    CALL FFT(J1_Re, J1_Im, NK*NL, NL, NK*NL, -1)

    CALL FFT(J2_Re, J2_Im, NK*NL, NK, NK, -1)
    CALL FFT(J2_Re, J2_Im, NK*NL, NL, NK*NL, -1)

    DO i = 0, Kmax
      DO j = 0, Lmax
        IF(i == 0 .AND. j == 0) CYCLE

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        J1(i,j) = (J1_Re(i,j) + (0,1)*J1_Im(i,j))/(NK*NL)
        J2(i,j) = (J2_Re(i,j) + (0,1)*J2_Im(i,j))/(NK*NL)

        Jcbn(i,j) = (0,1)*l*J1(i,j) - (0,1)*k*J2(i,j)

      END DO
    END DO

    Jcbn(0,0) = 0

!!!! DeAlias
    IF(alias) THEN
      DO i = 0, Kmax
        DO j = 0, Lmax

          k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
          l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

          IF((k**2 + l**2) > (0.666*Emax)**2) Jcbn(i,j) = 0

        END DO
      END DO
    END IF
!!!!

  END SUBROUTINE GetJacobian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is an explicit convolution. It is very slow but is more predictable
! and has less roundoff error.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE Jac2(zeta, Jcbn)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta
    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(OUT) :: Jcbn

    INTEGER :: i, j, m, n
    REAL :: ii, jj, mm, nn

!!!! DeAlias

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        IF((k**2 + l**2) > (0.666*Emax)**2) zeta(i,j) = 0

      END DO
    END DO

!!!!

    DO i = 0, Kmax
      DO j = 0, Lmax

        ii = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        jj = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        Jcbn(i,j) = 0

        DO m = 0, Kmax
          DO n = 0, Lmax

            mm = MODULO(m+(NK/2-1),NK) - (NK/2-1)
            nn = MODULO(n+(NL/2-1),NL) - (NL/2-1)

            IF(mm == 0 .AND. nn == 0) CYCLE

            Jcbn(i,j) = Jcbn(i,j) + zeta(m,n) &
                            * zeta(MODULO(NK-m+i,NK),MODULO(NL-n+j, NL)) &
                            * (mm*jj - ii*nn) / (mm**2 + nn**2)

          END DO
        END DO

      END DO
    END DO

!!!! Dealias

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2-1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2-1)

        IF((k**2 + l**2) > (0.666*Emax)**2) Jcbn(i,j) = 0

      END DO
    END DO

!!!!

  END SUBROUTINE Jac2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Dealias(zeta)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta

    INTEGER :: i, j

    DO i = 0, Kmax
      DO j = 0, Lmax

        k = MODULO(i+(NK/2-1),NK) - (NK/2 - 1)
        l = MODULO(j+(NL/2-1),NL) - (NL/2 - 1)

        !** Dealias the Beast! SAAATAN!
        IF((k**2 + l**2) > (0.666*Emax)**2) zeta(i,j) = 0

      END DO
    END DO

  END SUBROUTINE Dealias

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE BackupVorticity(zeta)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(IN) :: zeta

    INTEGER :: i, j

    OPEN(19, FILE='zeta.last', STATUS='REPLACE', ACCESS='SEQUENTIAL', &
             ACTION='WRITE')

    DO i = 0, Kmax
      DO j = 0, Lmax
        WRITE(19,*) zeta(i,j)
      END DO
    END DO

  END SUBROUTINE BackupVorticity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Change(zeta)

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax), INTENT(INOUT) :: zeta

    COMPLEX, DIMENSION(0:Kmax, 0:Lmax) :: temp

    REAL, DIMENSION(0:Kmax, 0:Lmax) :: zre, zim

    INTEGER :: i, j
    REAL :: x, y

    DO i = 0, Kmax
      DO  j = 0, Lmax

        zre(i,j) = REAL(zeta(i,j))
        zim(i,j) = AIMAG(zeta(i,j))

      END DO
    END DO

    CALL FFT(zre, zim, NK*NL, NK, NK, 1)
    CALL FFT(zre, zim, NK*NL, NL, NK*NL, 1)

    DO i = 0, Kmax
      DO j = 0, Lmax

        temp(i,j) = zre(i,j) + (0,1)*zim(i,j)

      END DO
    END DO

! change temp
   DO i = 0, Kmax
     DO j = 0, Lmax

       x = MODULO(i+(NK/2-1),NK) - (NK/2 - 1)
       y = MODULO(j+(NL/2-1),NL) - (NL/2 - 1)

       temp(i,j) = temp(i,j) - 1.*EXP(-(x**2 + y**2)/5.**2)

     END DO
  END DO

    DO i = 0, Kmax
      DO j = 0, Lmax

        zre(i,j) = REAL(temp(i,j))
        zim(i,j) = AIMAG(temp(i,j))

      END DO
    END DO

    CALL FFT(zre, zim, NK*NL, NK, NK, -1)
    CALL FFT(zre, zim, NK*NL, NL, NK*NL, -1)

    DO i = 0, Kmax
      DO j = 0, Lmax

        zeta(i,j) = (zre(i,j) + (0,1)*zim(i,j))/(NK*NL)

      END DO
    END DO

  END SUBROUTINE Change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END PROGRAM Turbulence2D

! FUNCTION TO COMPUTE SPHERICAL BESSEL FUNCTION J_N(X)
FUNCTION SPH_BESSEL_J(N, X) RESULT(JN)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION :: JN
    INTEGER, PARAMETER :: NMAX = 50  ! Choose a sufficiently large upper bound
    DOUBLE PRECISION :: J(NMAX+1)
    INTEGER :: I

    ! Special cases
    IF (X == 0.0D0) THEN
        IF (N == 0) THEN
            JN = 1.0D0
        ELSE
            JN = 0.0D0
        END IF
        RETURN
    END IF

    IF (N <= 2) THEN
        IF (N == 0) JN = DSIN(X)/X
        IF (N == 1) JN = -(DCOS(X)/X) + DSIN(X)/X**2
        RETURN
    ENDIF

    ! Initial conditions: large n approximation
    J(NMAX) = 0.0D0
    J(NMAX-1) = 1.0D-10  ! Arbitrary small number for stable backward recursion

    ! Backward recursion
    DO I = NMAX-2, 0, -1
        J(I) = (2.0D0*I + 3.0D0) / X * J(I+1) - J(I+2)
    END DO

    ! Normalization using j_0(x) = sin(x)/x
    JN = J(N) * (SIN(X) / X) / J(0)
END FUNCTION SPH_BESSEL_J


! FUNCTION TO COMPUTE SPHERICAL NEUMANN FUNCTION Y_N(X)
FUNCTION SPH_BESSEL_Y(N, X) RESULT(YN)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION :: YN, Y0, Y1, YPREV, YPREV2
    INTEGER :: I

    ! Special case: singularity at x = 0
    IF (X == 0.0D0) THEN
        YN = -1.0D99  ! Large negative number to indicate singularity
        RETURN
    END IF

    ! Base cases using explicit formulas
    Y0 = -COS(X) / X
    Y1 = (-COS(X) / X**2) - (SIN(X) / X)

    IF (N == 0) THEN
        YN = Y0
        RETURN
    ELSE IF (N == 1) THEN
        YN = Y1
        RETURN
    END IF

    ! Backward recursion for stability
    YPREV2 = Y0
    YPREV = Y1
    DO I = 2, N
        YN = ((2.0D0 * I - 1.0D0) / X) * YPREV - YPREV2
        YPREV2 = YPREV
        YPREV = YN
    END DO
END FUNCTION SPH_BESSEL_Y

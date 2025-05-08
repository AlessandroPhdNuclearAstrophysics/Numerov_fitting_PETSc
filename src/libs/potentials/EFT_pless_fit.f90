SUBROUTINE EFT_PLESS_PW_FIT(R, L, S, J, TZ, V, CF, RF, NCOEFF)
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: HTC = 197.32697D0
  DOUBLE PRECISION, PARAMETER :: PI= 4.D0*DATAN(1.D0)
  
  INTEGER, INTENT(IN) :: L, S, J, TZ, NCOEFF
  DOUBLE PRECISION, INTENT(IN) :: R
  DOUBLE PRECISION, INTENT(IN) :: CF(NCOEFF), RF
  DOUBLE PRECISION, INTENT(OUT):: V(2,2)

  LOGICAL, SAVE :: FIRST_CALL = .TRUE., NEW_CALL
  INTEGER, SAVE :: T, LS, TTZ, I2(2,2), LOLD=-1, SOLD=-1, JOLD=-1, TZOLD=-10
  DOUBLE PRECISION, SAVE :: S12(2,2)
  DOUBLE PRECISION :: VR(NCOEFF)

  NEW_CALL = L.NE.LOLD .OR. S.NE.SOLD .OR. J.NE.JOLD .OR. TZ.NE.TZOLD
  IF (FIRST_CALL.OR.NEW_CALL) THEN
    FIRST_CALL = .FALSE.
    
    T = MOD(MOD(L+S, 2) + 1, 2)
    I2(1,1) = 1
    I2(1,2) = 0
    I2(2,1) = 0
    I2(2,2) = 1
    
    LS = (J*(J+1) - L*(L+1) - S*(S+1))/2

    IF (TZ.EQ.0) THEN
      TTZ = -4
    ELSE
      TTZ = 2
    ENDIF
    
    S12 = 0.D0
    IF (L.EQ.(J+1)) THEN
      S12(1,1) = -2*(J+2)/(2*J+1)
    ELSEIF (L.EQ.J) THEN
      S12(1,1) = 2.d0
    ELSEIF (L.EQ.(J-1)) THEN
      S12(1,1) = -2*(J-1)/(2*J+1)
      S12(1,2) = 6*DSQRT(J*(J+1.D0))/(2*J+1)
      S12(2,1) = 6*DSQRT(J*(J+1.D0))/(2*J+1)
      S12(2,2) = -2*(J+2)/(2*J+1)
    ELSE
      WRITE(*,*) "ERROR: EFT_PLESS_PW_FIT, L, S, J = ", L, S, J
      STOP
    ENDIF
    LOLD = L
    SOLD = S
    JOLD = J
    TZOLD = TZ
    ! WRITE(*,*) "S12 = ", S12(1,1), S12(1,2), S12(2,1), S12(2,2)
    ! WRITE(*,*) "LS = ", LS
    ! WRITE(*,*) "TTZ = ", TTZ
    ! WRITE(*,*) "T = ", T
    ! WRITE(*,*) "CF = ", CF
  ENDIF
  
  
  V = 0.D0
  IF ( T.EQ.0 .AND. S.EQ.0 ) THEN
    V(1,1) = CF(1)*(-4*R**2 + 6*RF**2)/(DEXP(R**2/RF**2)*PI**1.5D0*RF**7)
  ENDIF
  IF ( T.EQ.1 .AND. S.EQ.1 ) THEN
    VR(1) = (-4*R**2 + 6*RF**2)/(RF**4)
    VR(2) =-(4*R**2)/(RF**4)
    VR(3) = 2/RF**2
    VR(4) = 1.D0
    V =  CF(1)*I2*      VR(1) &
        +CF(2)*S12*     VR(2) &
        +CF(3)*LS*I2*   VR(3) &
        +CF(4)*TTZ*I2*  VR(4)
    ! V =  CF(1)*(-4*R**2 + 6*RF**2)/RF**4*I2 &
    !     -CF(2)*S12*(4*R**2)/RF**4 &
    !     +CF(3)*LS*2/RF**2*I2 &
    !     +CF(4)*TTZ*I2
    V = V/(DEXP(R**2/RF**2)*PI**1.5D0*RF**3)
  ENDIF
  V = V*HTC

  RETURN
END SUBROUTINE EFT_PLESS_PW_FIT
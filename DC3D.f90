!----------------------------------------------------------------------
subroutine  DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4, &
                  UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
!********************************************************************
!*****                                                          *****
!*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****
!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
!*****                         CODED BY  Y.OKADA ... SEP.1991   *****
!*****                         REVISED     NOV.1991, MAY.2002   *****
!--------------------------------------------------------------------
!-----                     converted to Fortran90 (free form)   -----
!-----                                T. Miyashita, Jan. 2020   -----
!--------------------------------------------------------------------
!*****                                                          *****
!********************************************************************
!
!***** INPUT
!*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
!*****   X,Y,Z : COORDINATE OF OBSERVING POINT
!*****   DEPTH : SOURCE DEPTH
!*****   DIP   : DIP-ANGLE (DEGREE)
!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
!*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2
!*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3
!*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4
!
!***** OUTPUT
!*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /
!*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )
!*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /
!*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )
!*****   UXZ,UYZ,UZZ : Z-DERIVATIVE
!*****   IRET        : return CODE
!*****               :   =0....NORMAL
!*****               :   =1....SINGULAR
!*****               :   =2....POSITIVE Z WAS GIVEN
    implicit none

    ! arguments
    real*8, intent(in) :: ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4
    real*8, intent(out) :: UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
    integer, intent(out) :: IRET

    ! local variables
    real*8 :: ALP(5)
    real*8 :: U(12),DUA(12),DUB(12),DUC(12)
    real*8 :: DU
    real*8 :: SD,CD,R
    real*8 :: XX,YY,ZZ,DD
    ! parameters
    real*8, parameter :: EPS=1.0d-6
    ! loop counter
    integer :: I

    ! initialization
    U(1:12) = 0.0d0
    DUA(1:12) = 0.0d0
    DUB(1:12) = 0.0d0
    DUC(1:12) = 0.0d0
    call varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)

!### CAUTION ### if X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO
    XX=X
    YY=Y
    ZZ=Z
    DD=DEPTH+Z
    if(dabs(XX)<EPS) XX=0.0d0
    if(dabs(YY)<EPS) YY=0.0d0
    if(dabs(DD)<EPS) DD=0.0d0
!-----
    IRET=0
    if(Z>0.0d0) then
      IRET=2
      return
    endif
!-----
    R=dsqrt(XX*XX+YY*YY+DD*DD)
    if (R==0.0d0) then !IN CASE OF SINGULAR (R=0)
      IRET=1
      return
    endif

!**********************************************************************
!*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    *****
!**********************************************************************
    call assign_alpha(ALPHA,ALP)
    call DCCON0(DIP,SD,CD)
!======================================
!=====  REAL-SOURCE CONTRIBUTION  =====
!======================================
!-----
    call UA0(XX,YY,DD,POT1,POT2,POT3,POT4,SD,CD,ALP,DUA)
!-----
    do 222 I=1,12
      if(I<10) U(I)=U(I)-DUA(I)
      if(I>=10) U(I)=U(I)+DUA(I)
222 enddo
!=======================================
!=====  IMAGE-SOURCE CONTRIBUTION  =====
!=======================================
    DD=DEPTH-Z
    call UA0(XX,YY,DD,POT1,POT2,POT3,POT4,SD,CD,ALP,DUA)
    call UB0(XX,YY,DD,ZZ,POT1,POT2,POT3,POT4,SD,CD,ALP,DUB)
    call UC0(XX,YY,DD,ZZ,POT1,POT2,POT3,POT4,SD,CD,ALP,DUC)
!-----
    do 333 I=1,12
      DU=DUA(I)+DUB(I)+ZZ*DUC(I)
      if(I>=10) DU=DU+DUC(I-9)
      U(I)=U(I)+DU
333 enddo
!=====
    call varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)
    return
end subroutine DC3D0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UA0(X,Y,D,POT1,POT2,POT3,POT4,SD,CD,ALP,U)
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****
!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM
!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
    implicit none

    ! arguments
    real*8, intent(in) :: X,Y,D
    real*8, intent(in) :: POT1,POT2,POT3,POT4
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)

    ! local variables
    real*8 :: DU(12)
    real*8 :: P,Q,S,T
    real*8 :: R, R2, R3, R5
    real*8 :: S2D,C2D
    real*8 :: QR,QRX
    real*8 :: XY,X2,Y2,D2,A3,A5,B3,C3
    real*8 :: UY,VY,WY,UZ,VZ,WZ
    ! parameters
    real*8, parameter :: PI2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F3=3.0d0, F5=5.0d0
    ! loop counter
    integer :: I

    U(1:12)=0.0d0
    DU(1:12)=0.0d0

    call DCCON1(X,Y,D,SD,CD,R,P,Q,S,T)
    S2D=2.0d0*SD*CD
    C2D=CD*CD-SD*SD

    XY=X*Y
    X2=X*X
    Y2=Y*Y
    D2=D*D

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    QR=F3*Q/R5      ! 3q/(R^5)
    QRX=F5*QR*X/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*X2/R2
    A5=F1-F5*X2/R2
    B3=F1-F3*Y2/R2
    C3=F1-F3*D2/R2
!-----
    ! -- Tab.4 and 5, pp.1027-1028 --
    UY=SD-F5*Y*Q/R2   ! U  in Tab.4
    UZ=CD+F5*D*Q/R2   ! U' in Tab.5
    VY=S -F5*Y*P*Q/R2 ! V  in Tab.4
    VZ=T +F5*D*P*Q/R2 ! V' in Tab.5
    WY=UY+SD          ! W  in Tab.4
    WZ=UZ+CD          ! W' in Tab.5
!-----
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(POT1/=F0) then
      DU( 1)= ALP(1)*Q/R3    +ALP(2)*X2*QR
      DU( 2)= ALP(1)*X/R3*SD +ALP(2)*XY*QR
      DU( 3)=-ALP(1)*X/R3*CD +ALP(2)*X*D*QR
      DU( 4)= X*QR*(-ALP(1) +ALP(2)*(F1+A5) )
      DU( 5)= ALP(1)*A3/R3*SD +ALP(2)*Y*QR*A5
      DU( 6)=-ALP(1)*A3/R3*CD +ALP(2)*D*QR*A5
      DU( 7)= ALP(1)*(SD/R3-Y*QR) +ALP(2)*F3*X2/R5*UY
      DU( 8)= F3*X/R5*(-ALP(1)*Y*SD +ALP(2)*(Y*UY+Q) )
      DU( 9)= F3*X/R5*( ALP(1)*Y*CD +ALP(2)*D*UY )
      DU(10)= ALP(1)*(CD/R3+D*QR) +ALP(2)*F3*X2/R5*UZ
      DU(11)= F3*X/R5*( ALP(1)*D*SD +ALP(2)*Y*UZ )
      DU(12)= F3*X/R5*(-ALP(1)*D*CD +ALP(2)*(D*UZ-Q) )
      do 222 I=1,12
        U(I)=U(I)+POT1/PI2*DU(I)
 222  enddo
    endif
!===================================
!=====  DIP-SLIP CONTRIBUTION  =====
!===================================
    if(POT2/=F0) then
      DU( 1)=              ALP(2)*X*P*QR
      DU( 2)= ALP(1)*S/R3 +ALP(2)*Y*P*QR
      DU( 3)=-ALP(1)*T/R3 +ALP(2)*D*P*QR
      DU( 4)=              ALP(2)*P*QR*A5
      DU( 5)=-ALP(1)*F3*X*S/R5 -ALP(2)*Y*P*QRX
      DU( 6)= ALP(1)*F3*X*T/R5 -ALP(2)*D*P*QRX
      DU( 7)=                            ALP(2)*F3*X/R5*VY
      DU( 8)= ALP(1)*(S2D/R3-F3*Y*S/R5) +ALP(2)*(F3*Y/R5*VY+P*QR)
      DU( 9)=-ALP(1)*(C2D/R3-F3*Y*T/R5) +ALP(2)*F3*D/R5*VY
      DU(10)=                            ALP(2)*F3*X/R5*VZ
      DU(11)= ALP(1)*(C2D/R3+F3*D*S/R5) +ALP(2)*F3*Y/R5*VZ
      DU(12)= ALP(1)*(S2D/R3-F3*D*T/R5) +ALP(2)*(F3*D/R5*VZ-P*QR)
      do 333 I=1,12
        U(I)=U(I)+POT2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(POT3/=F0) then
      DU( 1)= ALP(1)*X/R3      -ALP(2)*X*Q*QR
      DU( 2)= ALP(1)*T/R3      -ALP(2)*Y*Q*QR
      DU( 3)= ALP(1)*S/R3      -ALP(2)*D*Q*QR
      DU( 4)= ALP(1)*A3/R3     -ALP(2)*Q*QR*A5
      DU( 5)=-ALP(1)*F3*X*T/R5 +ALP(2)*Y*Q*QRX
      DU( 6)=-ALP(1)*F3*X*S/R5 +ALP(2)*D*Q*QRX
      DU( 7)=-ALP(1)*F3*XY/R5           -ALP(2)*X*QR*WY
      DU( 8)= ALP(1)*(C2D/R3-F3*Y*T/R5) -ALP(2)*(Y*WY+Q)*QR
      DU( 9)= ALP(1)*(S2D/R3-F3*Y*S/R5) -ALP(2)*D*QR*WY
      DU(10)= ALP(1)*F3*X*D/R5          -ALP(2)*X*QR*WZ
      DU(11)=-ALP(1)*(S2D/R3-F3*D*T/R5) -ALP(2)*Y*QR*WZ
      DU(12)= ALP(1)*(C2D/R3+F3*D*S/R5) -ALP(2)*(D*WZ-Q)*QR
      do 444 I=1,12
        U(I)=U(I)+POT3/PI2*DU(I)
 444  enddo
    endif
!=========================================
!=====  INFLATE SOURCE CONTRIBUTION  =====
!=========================================
    if(POT4/=F0) then
      DU( 1)=-ALP(1)*X/R3
      DU( 2)=-ALP(1)*Y/R3
      DU( 3)=-ALP(1)*D/R3
      DU( 4)=-ALP(1)*A3/R3
      DU( 5)= ALP(1)*F3*XY/R5
      DU( 6)= ALP(1)*F3*X*D/R5
      DU( 7)= DU(5)
      DU( 8)=-ALP(1)*B3/R3
      DU( 9)= ALP(1)*F3*Y*D/R5
      DU(10)=-DU(6)
      DU(11)=-DU(9)
      DU(12)= ALP(1)*C3/R3
      do 555 I=1,12
        U(I)=U(I)+POT4/PI2*DU(I)
 555  enddo
    endif
    return
end subroutine  UA0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UB0(X,Y,D,Z,POT1,POT2,POT3,POT4,SD,CD,ALP,U)
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM
!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
    implicit none

    ! arguments
    real*8, intent(in) :: X,Y,D,Z
    real*8, intent(in) :: POT1,POT2,POT3,POT4
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)


    ! local variables
    real*8 :: DU(12)
    real*8 :: P,Q,S,T
    real*8 :: R, R2, R3, R5
    real*8 :: QR,QRX
    real*8 :: XY,X2,Y2,D2,A3,A5,B3,C3
    real*8 :: UY,VY,WY,UZ,VZ,WZ
    real*8 :: C,RD
    real*8 :: D12,D32,D33,D53,D54
    real*8 :: FI1,FI2,FI3,FI4,FI5
    real*8 :: FJ1,FJ2,FJ3,FJ4
    real*8 :: FK1,FK2,FK3
    ! parameters
    real*8, parameter :: PI2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0, F4=4.0d0, F5=5.0d0, F8=8.0d0, F9=9.0d0
    ! loop counter
    integer :: I

    DU(1:12)=0.0d0
    U(1:12)=0.0d0

    call DCCON1(X,Y,D,SD,CD,R,P,Q,S,T)
    XY=X*Y
    X2=X*X
    Y2=Y*Y
    D2=D*D

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    QR=F3*Q/R5      ! 3q/(R^5)
    QRX=F5*QR*X/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*X2/R2
    A5=F1-F5*X2/R2
    B3=F1-F3*Y2/R2
    C3=F1-F3*D2/R2
!-----
    ! -- Tab.4 and 5, pp.1027-1028 --
    UY=SD-F5*Y*Q/R2   ! U  in Tab.4
    UZ=CD+F5*D*Q/R2   ! U' in Tab.5
    VY=S -F5*Y*P*Q/R2 ! V  in Tab.4
    VZ=T +F5*D*P*Q/R2 ! V' in Tab.5
    WY=UY+SD          ! W  in Tab. 4
    WZ=UZ+CD          ! W' in Tab.5
!-----
    C=D+Z
    RD=R+D
    D12=F1/(R*RD*RD)
    D32=D12*(F2*R+D)/R2
    D33=D12*(F3*R+D)/(R2*RD)
    D53=D12*(F8*R2+F9*R*D+F3*D2)/(R2*R2*RD)
    D54=D12*(F5*R2+F4*R*D+D2)/R3*D12
!-----
    FI1= Y*(D12-X2*D33)
    FI2= X*(D12-Y2*D33)
    FI3= X/R3-FI2
    FI4=-XY*D32
    FI5= F1/(R*RD)-X2*D32
    FJ1=-F3*XY*(D33-X2*D54)
    FJ2= F1/R3-F3*D12+F3*X2*Y2*D54
    FJ3= A3/R3-FJ2
    FJ4=-F3*XY/R5-FJ1
    FK1=-Y*(D32-X2*D53)
    FK2=-X*(D32-Y2*D53)
    FK3=-F3*X*D/R5-FK2
!-----
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(POT1/=F0) then
      DU( 1)=-X2*QR  -ALP(3)*FI1*SD
      DU( 2)=-XY*QR  -ALP(3)*FI2*SD
      DU( 3)=-C*X*QR -ALP(3)*FI4*SD
      DU( 4)=-X*QR*(F1+A5) -ALP(3)*FJ1*SD
      DU( 5)=-Y*QR*A5      -ALP(3)*FJ2*SD
      DU( 6)=-C*QR*A5      -ALP(3)*FK1*SD
      DU( 7)=-F3*X2/R5*UY      -ALP(3)*FJ2*SD
      DU( 8)=-F3*XY/R5*UY-X*QR -ALP(3)*FJ4*SD
      DU( 9)=-F3*C*X/R5*UY     -ALP(3)*FK2*SD
      DU(10)=-F3*X2/R5*UZ  +ALP(3)*FK1*SD
      DU(11)=-F3*XY/R5*UZ  +ALP(3)*FK2*SD
      DU(12)= F3*X/R5*(-C*UZ +ALP(3)*Y*SD)
      do 222 I=1,12
        U(I)=U(I)+POT1/PI2*DU(I)
 222  enddo
    endif
!===================================
!=====  DIP-SLIP CONTRIBUTION  =====
!===================================
    if(POT2/=F0) then
      DU( 1)=-X*P*QR          +ALP(3)*FI3*SD*CD
      DU( 2)=-Y*P*QR          +ALP(3)*FI1*SD*CD
      DU( 3)=-C*P*QR          +ALP(3)*FI5*SD*CD
      DU( 4)=-P*QR*A5         +ALP(3)*FJ3*SD*CD
      DU( 5)= Y*P*QRX         +ALP(3)*FJ1*SD*CD
      DU( 6)= C*P*QRX         +ALP(3)*FK3*SD*CD
      DU( 7)=-F3*X/R5*VY      +ALP(3)*FJ1*SD*CD
      DU( 8)=-F3*Y/R5*VY-P*QR +ALP(3)*FJ2*SD*CD
      DU( 9)=-F3*C/R5*VY      +ALP(3)*FK1*SD*CD
      DU(10)=-F3*X/R5*VZ      -ALP(3)*FK3*SD*CD
      DU(11)=-F3*Y/R5*VZ      -ALP(3)*FK1*SD*CD
      DU(12)=-F3*C/R5*VZ      +ALP(3)*A3/R3*SD*CD
      do 333 I=1,12
        U(I)=U(I)+POT2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(POT3/=F0) then
      DU( 1)= X*Q*QR      -ALP(3)*FI3*SD*SD
      DU( 2)= Y*Q*QR      -ALP(3)*FI1*SD*SD
      DU( 3)= C*Q*QR      -ALP(3)*FI5*SD*SD
      DU( 4)= Q*QR*A5     -ALP(3)*FJ3*SD*SD
      DU( 5)=-Y*Q*QRX     -ALP(3)*FJ1*SD*SD
      DU( 6)=-C*Q*QRX     -ALP(3)*FK3*SD*SD
      DU( 7)= X*QR*WY     -ALP(3)*FJ1*SD*SD
      DU( 8)= QR*(Y*WY+Q) -ALP(3)*FJ2*SD*SD
      DU( 9)= C*QR*WY     -ALP(3)*FK1*SD*SD
      DU(10)= X*QR*WZ     +ALP(3)*FK3*SD*SD
      DU(11)= Y*QR*WZ     +ALP(3)*FK1*SD*SD
      DU(12)= C*QR*WZ     -ALP(3)*A3/R3*SD*SD
      do 444 I=1,12
        U(I)=U(I)+POT3/PI2*DU(I)
 444  enddo
    endif
!=========================================
!=====  INFLATE SOURCE CONTRIBUTION  =====
!=========================================
    if(POT4/=F0) then
      DU( 1)= ALP(3)*X/R3
      DU( 2)= ALP(3)*Y/R3
      DU( 3)= ALP(3)*D/R3
      DU( 4)= ALP(3)*A3/R3
      DU( 5)=-ALP(3)*F3*XY/R5
      DU( 6)=-ALP(3)*F3*X*D/R5
      DU( 7)= DU(5)
      DU( 8)= ALP(3)*B3/R3
      DU( 9)=-ALP(3)*F3*Y*D/R5
      DU(10)=-DU(6)
      DU(11)=-DU(9)
      DU(12)=-ALP(3)*C3/R3
      do 555 I=1,12
        U(I)=U(I)+POT4/PI2*DU(I)
 555  enddo
    endif
    return
end subroutine UB0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UC0(X,Y,D,Z,POT1,POT2,POT3,POT4,SD,CD,ALP,U)
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
!*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM
!*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
    implicit none

    ! arguments
    real*8, intent(in) :: X,Y,D,Z
    real*8, intent(in) :: POT1,POT2,POT3,POT4
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)

    ! local variables
    real*8 :: DU(12)
    real*8 :: P,Q,S,T
    real*8 :: R, R2, R3, R5, R7
    real*8 :: S2D,C2D
    real*8 :: Q2,QR,QRX,QR5,QR7,DR5
    real*8 :: XY,X2,Y2,D2,C
    real*8 :: A3,A5,A7,B3,B5,B7,C3,C5,C7,D7
    ! parameters
    real*8, parameter :: PI2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0, F5=5.0d0, F7=7.0d0, F10=1.0d1, F15=1.5d1
    ! loop counter
    integer :: I

    U(1:12)=0.0d0
    DU(1:12)=0.0d0

    call DCCON1(X,Y,D,SD,CD,R,P,Q,S,T)
    XY=X*Y
    X2=X*X
    Y2=Y*Y
    D2=D*D

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    QR=F3*Q/R5      ! 3q/(R^5)
    QRX=F5*QR*X/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*X2/R2
    A5=F1-F5*X2/R2
    B3=F1-F3*Y2/R2
    C3=F1-F3*D2/R2
!-----
    S2D=2.0d0*SD*CD
    C2D=CD*CD-SD*SD
!-----
    C=D+Z
    Q2=Q*Q
    R7=R5*R2
    A7=F1-F7*X2/R2
    B5=F1-F5*Y2/R2
    B7=F1-F7*Y2/R2
    C5=F1-F5*D2/R2
    C7=F1-F7*D2/R2
    D7=F2-F7*Q2/R2
    QR5=F5*Q/R2
    QR7=F7*Q/R2
    DR5=F5*D/R2
!-----
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(POT1/=F0) then
      DU( 1)=-ALP(4)*A3/R3*CD  +ALP(5)*C*QR*A5
      DU( 2)= F3*X/R5*( ALP(4)*Y*CD +ALP(5)*C*(SD-Y*QR5) )
      DU( 3)= F3*X/R5*(-ALP(4)*Y*SD +ALP(5)*C*(CD+D*QR5) )
      DU( 4)= ALP(4)*F3*X/R5*(F2+A5)*CD   -ALP(5)*C*QRX*(F2+A7)
      DU( 5)= F3/R5*( ALP(4)*Y*A5*CD +ALP(5)*C*(A5*SD-Y*QR5*A7) )
      DU( 6)= F3/R5*(-ALP(4)*Y*A5*SD +ALP(5)*C*(A5*CD+D*QR5*A7) )
      DU( 7)= DU(5)
      DU( 8)= F3*X/R5*( ALP(4)*B5*CD -ALP(5)*F5*C/R2*(F2*Y*SD+Q*B7) )
      DU( 9)= F3*X/R5*(-ALP(4)*B5*SD +ALP(5)*F5*C/R2*(D*B7*SD-Y*C7*CD) )
      DU(10)= F3/R5*   (-ALP(4)*D*A5*CD +ALP(5)*C*(A5*CD+D*QR5*A7) )
      DU(11)= F15*X/R7*( ALP(4)*Y*D*CD  +ALP(5)*C*(D*B7*SD-Y*C7*CD) )
      DU(12)= F15*X/R7*(-ALP(4)*Y*D*SD  +ALP(5)*C*(F2*D*CD-Q*C7) )
      do 222 I=1,12
        U(I)=U(I)+POT1/PI2*DU(I)
 222  enddo
    endif
!===================================
!=====  DIP-SLIP CONTRIBUTION  =====
!===================================
    if(POT2/=F0) then
      DU( 1)= ALP(4)*F3*X*T/R5          -ALP(5)*C*P*QRX
      DU( 2)=-ALP(4)/R3*(C2D-F3*Y*T/R2) +ALP(5)*F3*C/R5*(S-Y*P*QR5)
      DU( 3)=-ALP(4)*A3/R3*SD*CD        +ALP(5)*F3*C/R5*(T+D*P*QR5)
      DU( 4)= ALP(4)*F3*T/R5*A5              -ALP(5)*F5*C*P*QR/R2*A7
      DU( 5)= F3*X/R5*(ALP(4)*(C2D-F5*Y*T/R2)-ALP(5)*F5*C/R2*(S-Y*P*QR7))
      DU( 6)= F3*X/R5*(ALP(4)*(F2+A5)*SD*CD  -ALP(5)*F5*C/R2*(T+D*P*QR7))
      DU( 7)= DU(5)
      DU( 8)= F3/R5*( ALP(4)*(F2*Y*C2D+T*B5) &
                     +ALP(5)*C*(S2D-F10*Y*S/R2-P*QR5*B7))
      DU( 9)= F3/R5*(ALP(4)*Y*A5*SD*CD-ALP(5)*C*((F3+A5)*C2D+Y*P*DR5*QR7))
      DU(10)= F3*X/R5*(-ALP(4)*(S2D-T*DR5) -ALP(5)*F5*C/R2*(T+D*P*QR7))
      DU(11)= F3/R5*(-ALP(4)*(D*B5*C2D+Y*C5*S2D) &
                     -ALP(5)*C*((F3+A5)*C2D+Y*P*DR5*QR7))
      DU(12)= F3/R5*(-ALP(4)*D*A5*SD*CD-ALP(5)*C*(S2D-F10*D*T/R2+P*QR5*C7))
      do 333 I=1,12
        U(I)=U(I)+POT2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(POT3/=F0) then
      DU( 1)= F3*X/R5*(-ALP(4)*S +ALP(5)*(C*Q*QR5-Z))
      DU( 2)= ALP(4)/R3*(S2D-F3*Y*S/R2)+ALP(5)*F3/R5*(C*(T-Y+Y*Q*QR5)-Y*Z)
      DU( 3)=-ALP(4)/R3*(F1-A3*SD*SD)  -ALP(5)*F3/R5*(C*(S-D+D*Q*QR5)-D*Z)
      DU( 4)=-ALP(4)*F3*S/R5*A5 +ALP(5)*(C*QR*QR5*A7-F3*Z/R5*A5)
      DU( 5)= F3*X/R5*(-ALP(4)*(S2D-F5*Y*S/R2) &
                       -ALP(5)*F5/R2*(C*(T-Y+Y*Q*QR7)-Y*Z))
      DU( 6)= F3*X/R5*( ALP(4)*(F1-(F2+A5)*SD*SD) &
                       +ALP(5)*F5/R2*(C*(S-D+D*Q*QR7)-D*Z))
      DU( 7)= DU(5)
      DU( 8)= F3/R5*(-ALP(4)*(F2*Y*S2D+S*B5) &
                     -ALP(5)*(C*(F2*SD*SD+F10*Y*(T-Y)/R2-Q*QR5*B7)+Z*B5))
      DU( 9)= F3/R5*( ALP(4)*Y*(F1-A5*SD*SD) &
                     +ALP(5)*(C*(F3+A5)*S2D-Y*DR5*(C*D7+Z)))
      DU(10)= F3*X/R5*(-ALP(4)*(C2D+S*DR5) &
                       +ALP(5)*(F5*C/R2*(S-D+D*Q*QR7)-F1-Z*DR5))
      DU(11)= F3/R5*( ALP(4)*(D*B5*S2D-Y*C5*C2D) &
                     +ALP(5)*(C*((F3+A5)*S2D-Y*DR5*D7)-Y*(F1+Z*DR5)))
      DU(12)= F3/R5*(-ALP(4)*D*(F1-A5*SD*SD) &
                     -ALP(5)*(C*(C2D+F10*D*(S-D)/R2-Q*QR5*C7)+Z*(F1+C5)))
      do 444 I=1,12
        U(I)=U(I)+POT3/PI2*DU(I)
 444  enddo
    endif
!=========================================
!=====  INFLATE SOURCE CONTRIBUTION  =====
!=========================================
    if(POT4/=F0) then
      DU( 1)= ALP(4)*F3*X*D/R5
      DU( 2)= ALP(4)*F3*Y*D/R5
      DU( 3)= ALP(4)*C3/R3
      DU( 4)= ALP(4)*F3*D/R5*A5
      DU( 5)=-ALP(4)*F15*XY*D/R7
      DU( 6)=-ALP(4)*F3*X/R5*C5
      DU( 7)= DU(5)
      DU( 8)= ALP(4)*F3*D/R5*B5
      DU( 9)=-ALP(4)*F3*Y/R5*C5
      DU(10)= DU(6)
      DU(11)= DU(9)
      DU(12)= ALP(4)*F3*D/R5*(F2+C5)
      do 555 I=1,12
        U(I)=U(I)+POT4/PI2*DU(I)
 555  enddo
    endif
    return
end subroutine UC0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine DC3D(ALPHA,X,Y,Z,DEPTH,DIP, &
                AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, &
                UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
!********************************************************************
!*****                                                          *****
!*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****
!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
!*****              CODED BY  Y.OKADA ... SEP.1991              *****
!*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****
!*****                          JUL.1993, MAY.2002              *****
!--------------------------------------------------------------------
!-----                     converted to Fortran90 (free form)   -----
!-----                                T. Miyashita, Jan. 2020   -----
!--------------------------------------------------------------------
!********************************************************************
!
!***** INPUT
!*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
!*****   X,Y,Z : COORDINATE OF OBSERVING POINT
!*****   DEPTH : DEPTH OF REFERENCE POINT
!*****   DIP   : DIP-ANGLE (DEGREE)
!*****   AL1,AL2   : FAULT LENGTH RANGE
!*****   AW1,AW2   : FAULT WIDTH RANGE
!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
!
!***** OUTPUT
!*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)
!*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /
!*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )
!*****   UXZ,UYZ,UZZ : Z-DERIVATIVE
!*****   IRET        : return CODE
!*****               :   =0....NORMAL
!*****               :   =1....SINGULAR
!*****               :   =2....POSITIVE Z WAS GIVEN
    implicit none

    ! arguments
    real*8, intent(in) :: ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
    real*8, intent(out) :: UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
    integer, intent(out) :: IRET

    ! local variables
    real*8 :: SD,CD
    real*8 :: XI(2),ET(2)
    integer :: KXI(2),KET(2)
    real*8 :: U(12),DU(12),DUA(12),DUB(12),DUC(12)
    real*8 :: DD1, DD2, DD3
    real*8 :: ZZ, D, P, Q, R12, R21, R22
    real*8 :: ALP(5)
    ! parameters
    real*8, parameter :: F0 = 0.0d0
    real*8, parameter :: EPS = 1.0d-6
    ! loop counters
    integer :: I,J,K

    ! initialization
    U(1:12) = 0.0d0
    DU(1:12) = 0.0d0
    DUA(1:12) = 0.0d0
    DUB(1:12) = 0.0d0
    DUC(1:12) = 0.0d0
    call varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)

!-----
      IRET=0
    if(Z>0.0d0) then
      IRET=2
      return
    endif
!-----
    ZZ=Z
    DD1=DISL1
    DD2=DISL2
    DD3=DISL3
!### CAUTION ### if XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0
    XI(1)=X-AL1
    XI(2)=X-AL2
    if(dabs(XI(1))<EPS) XI(1)=F0
    if(dabs(XI(2))<EPS) XI(2)=F0
!======================================
!=====  REAL-SOURCE CONTRIBUTION  =====
!======================================
    call assign_alpha(ALPHA,ALP)
    call DCCON0(DIP,SD,CD)
!-----
    D=DEPTH+Z
    P=Y*CD+D*SD
    Q=Y*SD-D*CD
    ET(1)=P-AW1
    ET(2)=P-AW2
    if(dabs(Q)<EPS)  Q=F0
    if(dabs(ET(1))<EPS) ET(1)=F0
    if(dabs(ET(2))<EPS) ET(2)=F0
!--------------------------------
!----- REJECT SINGULAR CASE -----
!--------------------------------
!----- ON FAULT EDGE
    if(Q==F0 .and. &
      (    (XI(1)*XI(2)<=F0 .and. ET(1)*ET(2)==F0) &
       .or.(ET(1)*ET(2)<=F0 .and. XI(1)*XI(2)==F0) )) then
      IRET=1
      return
    endif
!----- ON NEGATIVE EXTENSION OF FAULT EDGE
    KXI(1)=0
    KXI(2)=0
    KET(1)=0
    KET(2)=0
    R12=dsqrt(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)
    R21=dsqrt(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)
    R22=dsqrt(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)
    if(XI(1)<F0 .and. R21+XI(2)<EPS) KXI(1)=1
    if(XI(1)<F0 .and. R22+XI(2)<EPS) KXI(2)=1
    if(ET(1)<F0 .and. R12+ET(2)<EPS) KET(1)=1
    if(ET(1)<F0 .and. R22+ET(2)<EPS) KET(2)=1
!=====
    do 223 K=1,2
    do 222 J=1,2
      call UA(XI(J),ET(K),Q,DD1,DD2,DD3,KXI(K),KET(J),SD,CD,ALP,DUA)
!-----
      do 220 I=1,10,3
        DU(I)  =-DUA(I)
        DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD
        DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD
        if(I<10) cycle
        DU(I)  =-DU(I)
        DU(I+1)=-DU(I+1)
        DU(I+2)=-DU(I+2)
220   enddo
      do 221 I=1,12
        if(J+K/=3) U(I)=U(I)+DU(I)
        if(J+K==3) U(I)=U(I)-DU(I)
221   enddo
!-----
222 enddo
223 enddo
!=======================================
!=====  IMAGE-SOURCE CONTRIBUTION  =====
!=======================================
    D=DEPTH-Z
    P=Y*CD+D*SD
    Q=Y*SD-D*CD
    ET(1)=P-AW1
    ET(2)=P-AW2
    if(dabs(Q)<EPS)  Q=F0
    if(dabs(ET(1))<EPS) ET(1)=F0
    if(dabs(ET(2))<EPS) ET(2)=F0
!--------------------------------
!----- REJECT SINGULAR CASE -----
!--------------------------------
!----- ON FAULT EDGE
    if(Q==F0 .and. &
      (    (XI(1)*XI(2)<=F0 .and. ET(1)*ET(2)==F0) &
       .or.(ET(1)*ET(2)<=F0 .and. XI(1)*XI(2)==F0) )) then
      IRET=1
      call varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)
      return
    endif
!----- ON NEGATIVE EXTENSION OF FAULT EDGE
    KXI(1)=0
    KXI(2)=0
    KET(1)=0
    KET(2)=0
    R12=dsqrt(XI(1)*XI(1)+ET(2)*ET(2)+Q*Q)
    R21=dsqrt(XI(2)*XI(2)+ET(1)*ET(1)+Q*Q)
    R22=dsqrt(XI(2)*XI(2)+ET(2)*ET(2)+Q*Q)
    if(XI(1)<F0 .and. R21+XI(2)<EPS) KXI(1)=1
    if(XI(1)<F0 .and. R22+XI(2)<EPS) KXI(2)=1
    if(ET(1)<F0 .and. R12+ET(2)<EPS) KET(1)=1
    if(ET(1)<F0 .and. R22+ET(2)<EPS) KET(2)=1
!=====
    do 334 K=1,2
    do 333 J=1,2
      call UA(XI(J),ET(K),Q   ,DD1,DD2,DD3,KXI(K),KET(J),SD,CD,ALP,DUA)
      call UB(XI(J),ET(K),Q   ,DD1,DD2,DD3,KXI(K),KET(J),SD,CD,ALP,DUB)
      call UC(XI(J),ET(K),Q,ZZ,DD1,DD2,DD3,KXI(K),KET(J),SD,CD,ALP,DUC)
!-----
      do 330 I=1,10,3
        DU(I)=DUA(I)+DUB(I)+Z*DUC(I)
        DU(I+1)=(DUA(I+1)+DUB(I+1)+Z*DUC(I+1))*CD &
               -(DUA(I+2)+DUB(I+2)+Z*DUC(I+2))*SD
        DU(I+2)=(DUA(I+1)+DUB(I+1)-Z*DUC(I+1))*SD &
               +(DUA(I+2)+DUB(I+2)-Z*DUC(I+2))*CD
        if(I<10) cycle
        DU(10)=DU(10)+DUC(1)
        DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD
        DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD
330   enddo
      do 331 I=1,12
        if(J+K/=3) U(I)=U(I)+DU(I)
        if(J+K==3) U(I)=U(I)-DU(I)
331   enddo
!-----
333 enddo
334 enddo
!=====
    call varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)
    return
end subroutine DC3D
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UA(XI,ET,Q,DISL1,DISL2,DISL3,KXI,KET,SD,CD,ALP,U)
!
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****
!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
!*****   KXI,KET     : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
!*****   SD,CD       : SIN, COS OF DIP-ANGLE
!*****   ALP         : COEFFICIENTS RELATED TO ALPHA VALUE
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
!-----
    implicit none
    ! arguments
    real*8, intent(in) :: XI,ET,Q,DISL1,DISL2,DISL3
    integer, intent(in) :: KXI,KET
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)
    ! local variables
    real*8 :: DU(12)
    real*8 :: R,Ytilde,Dtilde
    real*8 :: THETA,ALX,ALE
    real*8 :: X11,X32,Y11,Y32
    real*8 :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    real*8 :: XY,QX,QY
    real*8 :: XI2,Q2
    real*8 :: R3
    ! parameters
    real*8, parameter :: F0=0.0d0, F2=2.0d0
    real*8, parameter :: PI2=6.283185307179586d0
    ! loop counter
    integer :: I

    U(1:12) = 0.0d0
    DU(1:12) = 0.0d0

    call DCCON2(XI,ET,Q,SD,CD,R,Ytilde,Dtilde)
    call math_singularity(XI,ET,Q,R,KXI,KET,THETA,ALX,ALE)
    call eq14_xy(XI,ET,KXI,KET,R,X11,X32,Y11,Y32)
    call auxvar_derivative_yz(XI,Q,SD,CD,R,Ytilde,Dtilde,X11,X32,Y32, &
                              EY,EZ,FY,FZ,GY,GZ,HY,HZ)

    R3=R**3.0d0
    XI2=XI*XI
    Q2=Q*Q

    XY=XI*Y11
    QX=Q *X11
    QY=Q *Y11
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(DISL1/=F0) then
      DU( 1)= THETA/F2     +ALP(2)*XI*QY
      DU( 2)=               ALP(2)*Q/R
      DU( 3)= ALP(1)*ALE   -ALP(2)*Q*QY
      DU( 4)=-ALP(1)*QY    -ALP(2)*XI2*Q*Y32
      DU( 5)=              -ALP(2)*XI*Q/R3
      DU( 6)= ALP(1)*XY    +ALP(2)*XI*Q2*Y32
      DU( 7)= ALP(1)*XY*SD +ALP(2)*XI*FY+Dtilde/F2*X11
      DU( 8)=               ALP(2)*EY
      DU( 9)= ALP(1)*(CD/R+QY*SD) -ALP(2)*Q*FY
      DU(10)= ALP(1)*XY*CD        +ALP(2)*XI*FZ+Ytilde/F2*X11
      DU(11)=               ALP(2)*EZ
      DU(12)=-ALP(1)*(SD/R-QY*CD) -ALP(2)*Q*FZ
      do 222 I=1,12
        U(I)=U(I)+DISL1/PI2*DU(I)
 222  enddo
    endif
!======================================
!=====    DIP-SLIP CONTRIBUTION   =====
!======================================
    if(DISL2/=F0) then
      DU( 1)=             ALP(2)*Q/R
      DU( 2)= THETA/F2   +ALP(2)*ET*QX
      DU( 3)= ALP(1)*ALX -ALP(2)*Q*QX
      DU( 4)=            -ALP(2)*XI*Q/R3
      DU( 5)= -QY/F2     -ALP(2)*ET*Q/R3
      DU( 6)= ALP(1)/R   +ALP(2)*Q2/R3
      DU( 7)=             ALP(2)*EY
      DU( 8)= ALP(1)*Dtilde*X11+XY/F2*SD +ALP(2)*ET*GY
      DU( 9)= ALP(1)*Ytilde*X11          -ALP(2)*Q*GY
      DU(10)=             ALP(2)*EZ
      DU(11)= ALP(1)*Ytilde*X11+XY/F2*CD +ALP(2)*ET*GZ
      DU(12)=-ALP(1)*Dtilde*X11          -ALP(2)*Q*GZ
      do 333 I=1,12
        U(I)=U(I)+DISL2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(DISL3/=F0) then
      DU( 1)=-ALP(1)*ALE -ALP(2)*Q*QY
      DU( 2)=-ALP(1)*ALX -ALP(2)*Q*QX
      DU( 3)= THETA/F2   -ALP(2)*(ET*QX+XI*QY)
      DU( 4)=-ALP(1)*XY  +ALP(2)*XI*Q2*Y32
      DU( 5)=-ALP(1)/R   +ALP(2)*Q2/R3
      DU( 6)=-ALP(1)*QY  -ALP(2)*Q*Q2*Y32
      DU( 7)=-ALP(1)*(CD/R+QY*SD)  -ALP(2)*Q*FY
      DU( 8)=-ALP(1)*Ytilde*X11         -ALP(2)*Q*GY
      DU( 9)= ALP(1)*(Dtilde*X11+XY*SD) +ALP(2)*Q*HY
      DU(10)= ALP(1)*(SD/R-QY*CD)  -ALP(2)*Q*FZ
      DU(11)= ALP(1)*Dtilde*X11         -ALP(2)*Q*GZ
      DU(12)= ALP(1)*(Ytilde*X11+XY*CD) +ALP(2)*Q*HZ
      do 444 I=1,12
        U(I)=U(I)+DISL3/PI2*DU(I)
 444  enddo
    endif
    return
end subroutine UA
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UB(XI,ET,Q,DISL1,DISL2,DISL3,KXI,KET,SD,CD,ALP,U)
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****
!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
!*****   KXI,KET     : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
!*****   SD,CD       : SIN, COS OF DIP-ANGLE
!*****   ALP         : COEFFICIENTS RELATED TO ALPHA VALUE
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
!-----
    implicit none
    ! arguments
    real*8, intent(in) :: XI,ET,Q,DISL1,DISL2,DISL3
    integer, intent(in) :: KXI,KET
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)
    ! local variables
    real*8 :: DU(12)
    real*8 :: R,Ytilde,Dtilde
    real*8 :: THETA,ALE
    real*8 :: X11,X32,Y11,Y32
    real*8 :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    real*8 :: X
    real*8 :: XY,QX,QY
    real*8 :: XI2,Q2
    real*8 :: R3
    real*8 :: RD,RD2
    real*8 :: D11
    real*8 :: AI1,AI2,AI3,AI4
    real*8 :: AJ1,AJ2,AJ3,AJ4,AJ5,AJ6
    real*8 :: AK1,AK2,AK3,AK4
    ! parameters
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0
    real*8, parameter :: PI2=6.283185307179586d0
    ! loop counter
    integer :: I
    ! unused (dummy) vars
    real*8 :: dum ! ALX

    U(1:12)=0.0d0
    DU(1:12)=0.0d0

    call DCCON2(XI,ET,Q,SD,CD,R,Ytilde,Dtilde)
    call math_singularity(XI,ET,Q,R,KXI,KET,THETA,dum,ALE)
    call eq14_xy(XI,ET,KXI,KET,R,X11,X32,Y11,Y32)
    call auxvar_derivative_yz(XI,Q,SD,CD,R,Ytilde,Dtilde,X11,X32,Y32, &
                              EY,EZ,FY,FZ,GY,GZ,HY,HZ)

    R3=R**3.0d0
    XI2=XI*XI
    Q2=Q*Q

    RD=R+Dtilde
    D11=F1/(R*RD)
    AJ2=XI*Ytilde/RD*D11
    AJ5=-(Dtilde+Ytilde*Ytilde/RD)*D11
    if(CD/=F0) then
      ! p.1034 (ii)
      if(XI==F0) then
        AI4=F0
      else
        X=dsqrt(XI2+Q2)
        AI4=F1/(CD*CD)*( XI/RD*SD*CD &
           +F2*datan((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)) )
      endif
      AI3=(Ytilde*CD/RD-ALE+SD*dlog(RD))/(CD*CD)
      AK1=XI*(D11-Y11*SD)/CD
      AK3=(Q*Y11-Ytilde*D11)/CD
      AJ3=(AK1-AJ2*SD)/CD
      AJ6=(AK3-AJ5*SD)/CD
    else
      RD2=RD*RD
      AI3=(ET/RD+Ytilde*Q/RD2-ALE)/F2
      AI4=XI*Ytilde/RD2/F2
      AK1=XI*Q/RD*D11
      AK3=SD/RD*(XI2*D11-F1)
      AJ3=-XI/RD2*(Q2*D11-F1/F2)
      AJ6=-Ytilde/RD2*(XI2*D11-F1/F2)
    endif
!-----
    XY=XI*Y11
    AI1=-XI/RD*CD-AI4*SD
    AI2= dlog(RD)+AI3*SD
    AK2= F1/R+AK3*SD
    AK4= XY*CD-AK1*SD
    AJ1= AJ5*CD-AJ6*SD
    AJ4=-XY-AJ2*CD+AJ3*SD
!=====
    QX=Q*X11
    QY=Q*Y11
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(DISL1/=F0) then
      DU( 1)=-XI*QY-THETA -ALP(3)*AI1*SD
      DU( 2)=-Q/R         +ALP(3)*Ytilde/RD*SD
      DU( 3)= Q*QY        -ALP(3)*AI2*SD
      DU( 4)= XI2*Q*Y32   -ALP(3)*AJ1*SD
      DU( 5)= XI*Q/R3     -ALP(3)*AJ2*SD
      DU( 6)=-XI*Q2*Y32   -ALP(3)*AJ3*SD
      DU( 7)=-XI*FY-Dtilde*X11 +ALP(3)*(XY+AJ4)*SD
      DU( 8)=-EY          +ALP(3)*(F1/R+AJ5)*SD
      DU( 9)= Q*FY        -ALP(3)*(QY-AJ6)*SD
      DU(10)=-XI*FZ-Ytilde*X11 +ALP(3)*AK1*SD
      DU(11)=-EZ          +ALP(3)*Ytilde*D11*SD
      DU(12)= Q*FZ        +ALP(3)*AK2*SD
      do 222 I=1,12
        U(I)=U(I)+DISL1/PI2*DU(I)
 222  enddo
    endif
!======================================
!=====    DIP-SLIP CONTRIBUTION   =====
!======================================
    if(DISL2/=F0) then
      DU( 1)=-Q/R         +ALP(3)*AI3*SD*CD
      DU( 2)=-ET*QX-THETA -ALP(3)*XI/RD*SD*CD
      DU( 3)= Q*QX        +ALP(3)*AI4*SD*CD
      DU( 4)= XI*Q/R3     +ALP(3)*AJ4*SD*CD
      DU( 5)= ET*Q/R3+QY  +ALP(3)*AJ5*SD*CD
      DU( 6)=-Q2/R3       +ALP(3)*AJ6*SD*CD
      DU( 7)=-EY          +ALP(3)*AJ1*SD*CD
      DU( 8)=-ET*GY-XY*SD +ALP(3)*AJ2*SD*CD
      DU( 9)= Q*GY        +ALP(3)*AJ3*SD*CD
      DU(10)=-EZ          -ALP(3)*AK3*SD*CD
      DU(11)=-ET*GZ-XY*CD -ALP(3)*XI*D11*SD*CD
      DU(12)= Q*GZ        -ALP(3)*AK4*SD*CD
      do 333 I=1,12
        U(I)=U(I)+DISL2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(DISL3/=F0) then
      DU( 1)= Q*QY              -ALP(3)*AI3*SD*SD
      DU( 2)= Q*QX              +ALP(3)*XI/RD*SD*SD
      DU( 3)= ET*QX+XI*QY-THETA -ALP(3)*AI4*SD*SD
      DU( 4)=-XI*Q2*Y32 -ALP(3)*AJ4*SD*SD
      DU( 5)=-Q2/R3     -ALP(3)*AJ5*SD*SD
      DU( 6)= Q*Q2*Y32  -ALP(3)*AJ6*SD*SD
      DU( 7)= Q*FY -ALP(3)*AJ1*SD*SD
      DU( 8)= Q*GY -ALP(3)*AJ2*SD*SD
      DU( 9)=-Q*HY -ALP(3)*AJ3*SD*SD
      DU(10)= Q*FZ +ALP(3)*AK3*SD*SD
      DU(11)= Q*GZ +ALP(3)*XI*D11*SD*SD
      DU(12)=-Q*HZ +ALP(3)*AK4*SD*SD
      do 444 I=1,12
        U(I)=U(I)+DISL3/PI2*DU(I)
 444  enddo
    endif
    return
end subroutine UB
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  UC(XI,ET,Q,Z,DISL1,DISL2,DISL3,KXI,KET,SD,CD,ALP,U)
!********************************************************************
!*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****
!*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****
!********************************************************************
!
!***** INPUT
!*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM
!*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS
!*****   KXI,KET     : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
!*****   SD,CD       : SIN, COS OF DIP-ANGLE
!*****   ALP         : COEFFICIENTS RELATED TO ALPHA VALUE
!***** OUTPUT
!*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES
    implicit none
    ! arguments
    real*8, intent(in) :: XI,ET,Z,Q,DISL1,DISL2,DISL3
    integer, intent(in) :: KXI,KET
    real*8, intent(in) :: SD,CD
    real*8, intent(in) :: ALP(5)
    real*8, intent(out) :: U(12)
    ! local variables
    real*8 :: DU(12)
    real*8 :: R,Ytilde,Dtilde
    real*8 :: X11,Y11,X32,Y32
    real*8 :: XY,QX,QY
    real*8 :: XI2,ET2,Q2
    real*8 :: R2,R3,R5
    real*8 :: C
    real*8 :: X53,Y53,Z32,Z53,H
    real*8 :: Y0,Z0
    real*8 :: PPY,PPZ
    real*8 :: QQ,QQY,QQZ,QR,CQX,CDR,YY0
    ! parameters
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0
    real*8, parameter :: PI2=6.283185307179586d0
    ! loop counter
    integer :: I

    U(1:12) = 0.0d0
    DU(1:12) = 0.0d0

    call DCCON2(XI,ET,Q,SD,CD,R,Ytilde,Dtilde)
    call eq14_xy(XI,ET,KXI,KET,R,X11,X32,Y11,Y32)

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    XI2=XI*XI
    ET2=ET*ET
    Q2=Q*Q
!-----
    C=Dtilde+Z
    X53=(8.D0*R2+9.D0*R*XI+F3*XI2)*X11*X11*X11/R2
    Y53=(8.D0*R2+9.D0*R*ET+F3*ET2)*Y11*Y11*Y11/R2
    H=Q*CD-Z
    Z32=SD/R3-H*Y32
    Z53=F3*SD/R5-H*Y53
    Y0=Y11-XI2*Y32
    Z0=Z32-XI2*Z53
    PPY=CD/R3+Q*Y32*SD
    PPZ=SD/R3-Q*Y32*CD
    QQ=Z*Y32+Z32+Z0
    QQY=F3*C*Dtilde/R5-QQ*SD
    QQZ=F3*C*Ytilde/R5-QQ*CD+Q*Y32
    XY=XI*Y11
    QX=Q*X11
    QY=Q*Y11
    QR=F3*Q/R5
    CQX=C*Q*X53
    CDR=(C+Dtilde)/R3
    YY0=Ytilde/R3-Y0*CD
!=====
!======================================
!=====  STRIKE-SLIP CONTRIBUTION  =====
!======================================
    if(DISL1/=F0) then
      DU( 1)= ALP(4)*XY*CD           -ALP(5)*XI*Q*Z32
      DU( 2)= ALP(4)*(CD/R+F2*QY*SD) -ALP(5)*C*Q/R3
      DU( 3)= ALP(4)*QY*CD           -ALP(5)*(C*ET/R3-Z*Y11+XI2*Z32)
      DU( 4)= ALP(4)*Y0*CD                  -ALP(5)*Q*Z0
      DU( 5)=-ALP(4)*XI*(CD/R3+F2*Q*Y32*SD) +ALP(5)*C*XI*QR
      DU( 6)=-ALP(4)*XI*Q*Y32*CD  +ALP(5)*XI*(F3*C*ET/R5-QQ)
      DU( 7)=-ALP(4)*XI*PPY*CD    -ALP(5)*XI*QQY
      DU( 8)= ALP(4)*F2*(Dtilde/R3-Y0*SD)*SD-Ytilde/R3*CD &
             -ALP(5)*(CDR*SD-ET/R3-C*Ytilde*QR)
      DU( 9)=-ALP(4)*Q/R3+YY0*SD  +ALP(5)*(CDR*CD+C*Dtilde*QR-(Y0*CD+Q*Z0)*SD)
      DU(10)= ALP(4)*XI*PPZ*CD    -ALP(5)*XI*QQZ
      DU(11)= ALP(4)*F2*(Ytilde/R3-Y0*CD)*SD+Dtilde/R3*CD -ALP(5)*(CDR*CD+C*Dtilde*QR)
      DU(12)=        YY0*CD    -ALP(5)*(CDR*SD-C*Ytilde*QR-Y0*SD*SD+Q*Z0*CD)
      do 222 I=1,12
        U(I)=U(I)+DISL1/PI2*DU(I)
 222  enddo
    endif
!======================================
!=====    DIP-SLIP CONTRIBUTION   =====
!======================================
    if(DISL2/=F0) then
      DU( 1)= ALP(4)*CD/R -QY*SD      -ALP(5)*C*Q/R3
      DU( 2)= ALP(4)*Ytilde*X11       -ALP(5)*C*ET*Q*X32
      DU( 3)=       -Dtilde*X11-XY*SD -ALP(5)*C*(X11-Q2*X32)
      DU( 4)=-ALP(4)*XI/R3*CD         +ALP(5)*C*XI*QR +XI*Q*Y32*SD
      DU( 5)=-ALP(4)*Ytilde/R3        +ALP(5)*C*ET*QR
      DU( 6)=        Dtilde/R3-Y0*SD  +ALP(5)*C/R3*(F1-F3*Q2/R2)
      DU( 7)=-ALP(4)*ET/R3+Y0*SD*SD   -ALP(5)*(CDR*SD-C*Ytilde*QR)
      DU( 8)= ALP(4)*(X11-Ytilde*Ytilde*X32) -ALP(5)*C*((Dtilde+F2*Q*CD)*X32-Ytilde*ET*Q*X53)
      DU( 9)=   XI*PPY*SD+Ytilde*Dtilde*X32  +ALP(5)*C*((Ytilde+F2*Q*SD)*X32-Ytilde*Q2*X53)
      DU(10)=   -Q/R3+Y0*SD*CD               -ALP(5)*(CDR*CD+C*Dtilde*QR)
      DU(11)= ALP(4)*Ytilde*Dtilde*X32       -ALP(5)*C*((Ytilde-F2*Q*SD)*X32+Dtilde*ET*Q*X53)
      DU(12)=-XI*PPZ*SD+X11-Dtilde*Dtilde*X32-ALP(5)*C*((Dtilde-F2*Q*CD)*X32-Dtilde*Q2*X53)
      do 333 I=1,12
        U(I)=U(I)+DISL2/PI2*DU(I)
 333  enddo
    endif
!========================================
!=====  TENSILE-FAULT CONTRIBUTION  =====
!========================================
    if(DISL3/=F0) then
      DU( 1)=-ALP(4)*(SD/R+QY*CD)        -ALP(5)*(Z*Y11-Q2*Z32)
      DU( 2)= ALP(4)*F2*XY*SD+Dtilde*X11 -ALP(5)*C*(X11-Q2*X32)
      DU( 3)= ALP(4)*(Ytilde*X11+XY*CD)  +ALP(5)*Q*(C*ET*X32+XI*Z32)
      DU( 4)= ALP(4)*XI/R3*SD+XI*Q*Y32*CD+ALP(5)*XI*(F3*C*ET/R5-F2*Z32-Z0)
      DU( 5)= ALP(4)*F2*Y0*SD-Dtilde/R3  +ALP(5)*C/R3*(F1-F3*Q2/R2)
      DU( 6)=-ALP(4)*YY0                 -ALP(5)*(C*ET*QR-Q*Z0)
      DU( 7)= ALP(4)*(Q/R3+Y0*SD*CD)     +ALP(5)*(Z/R3*CD+C*Dtilde*QR-Q*Z0*SD)
      DU( 8)=-ALP(4)*F2*XI*PPY*SD-Ytilde*Dtilde*X32 &
             +ALP(5)*C*((Ytilde+F2*Q*SD)*X32-Ytilde*Q2*X53)
      DU( 9)=-ALP(4)*(XI*PPY*CD-X11+Ytilde*Ytilde*X32) &
             +ALP(5)*(C*((Dtilde+F2*Q*CD)*X32-Ytilde*ET*Q*X53)+XI*QQY)
      DU(10)=  -ET/R3+Y0*CD*CD           -ALP(5)*(Z/R3*SD-C*Ytilde*QR-Y0*SD*SD+Q*Z0*CD)
      DU(11)= ALP(4)*F2*XI*PPZ*SD-X11+Dtilde*Dtilde*X32 &
             -ALP(5)*C*((Dtilde-F2*Q*CD)*X32-Dtilde*Q2*X53)
      DU(12)= ALP(4)*(XI*PPZ*CD+Ytilde*Dtilde*X32) &
             +ALP(5)*(C*((Ytilde-F2*Q*SD)*X32+Dtilde*ET*Q*X53)+XI*QQZ)
      do 444 I=1,12
        U(I)=U(I)+DISL3/PI2*DU(I)
 444  enddo
    endif
    return
end subroutine UC
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  DCCON0(DIP,SD,CD)
!*******************************************************************
!*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****
!*******************************************************************
!
!***** INPUT
!*****   DIP   : DIP-ANGLE (DEGREE)
!***** OUTPUT
!*****   SD,CD : SIN, COS OF DIP-ANGLE
    implicit none
    real*8, intent(in) :: DIP
    real*8, intent(out) :: SD,CD

    real*8, parameter :: F0=0.0d0, F1=1.0d0
    real*8, parameter :: PI2=6.283185307179586d0
    real*8, parameter :: EPS=1.0d-6
    real*8, parameter :: P18=PI2/360.D0
!-----
    SD=DSIN(DIP*P18)
    CD=DCOS(DIP*P18)
!### CAUTION ### if COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO
    if(dabs(CD)<EPS) then
      CD=F0
      if(SD>F0) SD= F1
      if(SD<F0) SD=-F1
    endif
    return
end subroutine DCCON0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  DCCON1(X,Y,D,SD,CD,R,P,Q,S,T)
!**********************************************************************
!*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    *****
!**********************************************************************
!***** INPUT
!*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM
!*****   SD,CD : SIN, COS OF DIP-ANGLE
    implicit none
    real*8, intent(in) :: X,Y,D
    real*8, intent(in) :: SD,CD
    real*8, intent(out) :: R
    real*8, intent(out) :: P,Q,S,T
!-----
    P=Y*CD+D*SD
    Q=Y*SD-D*CD
    S=P*SD+Q*CD
    T=P*CD-Q*SD
    R=dsqrt(X*X+Y*Y+D*D)
    return
end subroutine DCCON1
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  DCCON2(XI,ET,Q,SD,CD,R,Ytilde,Dtilde)
!**********************************************************************
!*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   *****
!**********************************************************************
!***** INPUT
!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
!*****   SD,CD   : SIN, COS OF DIP-ANGLE
    implicit none
    real*8, intent(in) :: XI,ET,Q
    real*8, intent(in) :: SD,CD
    real*8, intent(out) :: R,Ytilde,Dtilde
!-----
    R = dsqrt(XI*XI+ET*ET+Q*Q)
    Ytilde = ET*CD+Q*SD
    Dtilde = ET*SD-Q*CD
    return
end subroutine DCCON2
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine math_singularity(XI,ET,Q,R,KXI,KET,THETA,ALX,ALE)
!***** INPUT
!*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM
!*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
    implicit none
    real*8, intent(in) :: XI,ET,Q,R
    integer, intent(in) :: KXI,KET
    real*8, intent(out) :: THETA,ALX,ALE

    ! p.1034 (i)
    if(Q==0.0d0) then
      THETA=0.0d0
    else
      THETA=datan(XI*ET/(Q*R))
    endif
    ! p.1034 (iii)
    if(KXI==1) then
      ALX=-dlog(R-XI)
    else
      ALX=dlog(R+XI)
    endif
    ! p.1034 (iv)
    if(KET==1) then
      ALE=-dlog(R-ET)
    else
      ALE=dlog(R+ET)
    endif

    return
end subroutine math_singularity
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine eq14_xy(XI,ET,KXI,KET,R,X11,X32,Y11,Y32)
!***** INPUT
!*****   XI,ET   : STATION COORDINATES IN FAULT SYSTEM
!*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY
    implicit none
    real*8, intent(in) :: XI,ET
    integer, intent(in) :: KXI,KET
    real*8, intent(in) :: R
    real*8, intent(out) :: X11,X32,Y11,Y32

    X11 = 0.0d0
    X32 = 0.0d0
    Y11 = 0.0d0
    Y32 = 0.0d0
!-----
    if(KXI/=1) then
      X11=1.0d0/(R*(R+XI))        ! Eq.(14), p.1026
      X32=(2.0d0*R+XI)*X11*X11/R  ! Eq.(14), p.1026
    endif
!-----
    if(KET/=1) then
      Y11=1.0d0/(R*(R+ET))        ! Eq.(14), p.1026
      Y32=(2.0d0*R+ET)*Y11*Y11/R  ! Eq.(14), p.1026
    endif
!-----
    return
end subroutine eq14_xy
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine auxvar_derivative_yz(XI,Q,SD,CD,R,Ytilde,Dtilde,X11,X32,Y32, &
                                EY,EZ,FY,FZ,GY,GZ,HY,HZ)
!***** INPUT
!*****   XI,Q  : STATION COORDINATES IN FAULT SYSTEM
!*****   SD,CD : SIN, COS OF DIP-ANGLE
    implicit none
    ! arguments
    real*8, intent(in) :: XI,Q,SD,CD,R,Ytilde,Dtilde
    real*8, intent(in) :: X11,X32,Y32
    real*8, intent(out) :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    ! local variables
    real*8 :: R3
    R3 = R**3.0d0

    EY=SD/R-Ytilde*Q/R3           ! E  in Tab.8, p.1032
    EZ=CD/R+Dtilde*Q/R3           ! E' in Tab.9, p.1033
    FY=Dtilde/R3+XI*XI*Y32*SD     ! F  in Tab.8, p.1032
    FZ=Ytilde/R3+XI*XI*Y32*CD     ! F' in Tab.9, p.1033
    GY=2.0d0*X11*SD-Ytilde*Q*X32  ! G  in Tab.8, p.1032
    GZ=2.0d0*X11*CD+Dtilde*Q*X32  ! G' in Tab.9, p.1033
    HY=Dtilde*Q*X32+XI*Q*Y32*SD   ! H  in Tab.8, p.1032
    HZ=Ytilde*Q*X32+XI*Q*Y32*CD   ! H' in Tab.9, p.1033
    return
end subroutine auxvar_derivative_yz
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine assign_alpha(ALPHA,ALP)
!***** INPUT
!*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)
!***** OUTPUT
!*****   ALP   : COEFFICIENTS RELATED TO ALPHA VALUE
    implicit none
    real*8, intent(in) :: ALPHA
    real*8, intent(out) :: ALP(5)
    ALP(1)=(1.0d0-ALPHA)/2.0d0
    ALP(2)= ALPHA/2.0d0
    ALP(3)=(1.0d0-ALPHA)/ALPHA
    ALP(4)=1.0d0-ALPHA
    ALP(5)=ALPHA
    return
end subroutine assign_alpha
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine varout(U,UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ)
    implicit none
    real*8, intent(in) :: U(12)
    real*8, intent(out) :: UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
    UX=U(1)
    UY=U(2)
    UZ=U(3)
    UXX=U(4)
    UYX=U(5)
    UZX=U(6)
    UXY=U(7)
    UYY=U(8)
    UZY=U(9)
    UXZ=U(10)
    UYZ=U(11)
    UZZ=U(12)
    return
end subroutine varout
!----------------------------------------------------------------------

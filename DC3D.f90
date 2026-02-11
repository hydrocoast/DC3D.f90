!----------------------------------------------------------------------
subroutine  DC3D0(alpha,x,y,z,depth,dip,pot1,pot2,pot3,pot4, &
                  ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
!********************************************************************
!*****                                                          *****
!*****   displacement and strain at depth                       *****
!*****   due to a buried point source in a semi-infinite medium *****
!*****                         coded by  Y.Okada ... Sep.1991   *****
!*****                         revised     Nov.1991, May.2002   *****
!--------------------------------------------------------------------
!-----                     converted to Fortran90 (free form)   -----
!-----                                T. Miyashita, Jan. 2020   -----
!--------------------------------------------------------------------
!*****                                                          *****
!********************************************************************
!
!***** input
!*****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
!*****   x,y,z : coordinate of observing point
!*****   depth : source depth
!*****   dip   : dip-angle (degree)
!*****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
!*****       potency=(  moment of double-couple  )/myu     for pot1,2
!*****       potency=(intensity of isotropic part)/lambda  for pot3
!*****       potency=(intensity of linear dipole )/myu     for pot4
!
!***** output
!*****   ux, uy, uz  : displacement ( unit=(unit of potency) /
!*****               :                     (unit of x,y,z,depth)**2  )
!*****   uxx,uyx,uzx : x-derivative ( unit= unit of potency) /
!*****   uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth)**3  )
!*****   uxz,uyz,uzz : z-derivative
!*****   iret        : return code
!*****               :   =0....normal
!*****               :   =1....singular
!*****               :   =2....positive z was given
    implicit none

    ! arguments
    real*8, intent(in) :: alpha,x,y,z,depth,dip,pot1,pot2,pot3,pot4
    real*8, intent(out) :: ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer, intent(out) :: iret

    ! local variables
    real*8 :: alp(5)
    real*8 :: u(12),duA(12),duB(12),duC(12)
    real*8 :: du
    real*8 :: sd,cd,R
    real*8 :: xx,yy,zz,dd
    ! parameters
    real*8, parameter :: eps=1.0d-6
    ! loop counter
    integer :: i

    ! initialization
    u(1:12) = 0.0d0
    duA(1:12) = 0.0d0
    duB(1:12) = 0.0d0
    duC(1:12) = 0.0d0
    call varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)

!### caution ### if x,y,d are sufficiently small, they are set to zero
    xx=x
    yy=y
    zz=z
    dd=depth+z
    if(dabs(xx)<eps) xx=0.0d0
    if(dabs(yy)<eps) yy=0.0d0
    if(dabs(dd)<eps) dd=0.0d0
!-----
    iret=0
    if(z>0.0d0) then
      iret=2
      return
    endif
!-----
    R=dsqrt(xx*xx+yy*yy+dd*dd)
    if (R==0.0d0) then !in case of singular (R=0)
      iret=1
      return
    endif

!**********************************************************************
!*****   calculate station geometry constants for point source    *****
!**********************************************************************
    call assign_alpha(alpha,alp)
    call dccon0(dip,sd,cd)
!======================================
!=====  real-source contribution  =====
!======================================
!-----
    call uA0(xx,yy,dd,pot1,pot2,pot3,pot4,sd,cd,alp,duA)
!-----
    do 222 i=1,12
      if(i<10) u(i)=u(i)-duA(i)
      if(i>=10) u(i)=u(i)+duA(i)
222 enddo
!=======================================
!=====  image-source contribution  =====
!=======================================
    dd=depth-z
    call uA0(xx,yy,dd,pot1,pot2,pot3,pot4,sd,cd,alp,duA)
    call uB0(xx,yy,dd,zz,pot1,pot2,pot3,pot4,sd,cd,alp,duB)
    call uC0(xx,yy,dd,zz,pot1,pot2,pot3,pot4,sd,cd,alp,duC)
!-----
    do 333 i=1,12
      du=duA(i)+duB(i)+zz*duC(i)
      if(i>=10) du=du+duC(i-9)
      u(i)=u(i)+du
333 enddo
!=====
    call varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)
    return
end subroutine DC3D0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uA0(x,y,d,pot1,pot2,pot3,pot4,sd,cd,alp,u)
!********************************************************************
!*****    displacement and strain at depth (part-A)             *****
!*****    due to buried point source in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   x,y,d : station coordinates in fault system
!*****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
!***** output
!*****   u(12) : displacement and their derivatives
    implicit none

    ! arguments
    real*8, intent(in) :: x,y,d
    real*8, intent(in) :: pot1,pot2,pot3,pot4
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)

    ! local variables
    real*8 :: du(12)
    real*8 :: p,q,s,t
    real*8 :: R, R2, R3, R5
    real*8 :: s2d,c2d
    real*8 :: qR,qRx
    real*8 :: xy,x2,y2,d2,A3,A5,B3,C3
    real*8 :: uy,vy,wy,uz,vz,wz
    ! parameters
    real*8, parameter :: pi2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F3=3.0d0, F5=5.0d0
    ! loop counter
    integer :: i

    u(1:12)=0.0d0
    du(1:12)=0.0d0

    call dccon1(x,y,d,sd,cd,R,p,q,s,t)
    s2d=2.0d0*sd*cd
    c2d=cd*cd-sd*sd

    xy=x*y
    x2=x*x
    y2=y*y
    d2=d*d

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    qR=F3*q/R5      ! 3q/(R^5)
    qRx=F5*qR*x/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*x2/R2
    A5=F1-F5*x2/R2
    B3=F1-F3*y2/R2
    C3=F1-F3*d2/R2
!-----
    ! -- Tab.4 and 5, pp.1027-1028 --
    uy=sd-F5*y*q/R2   ! U  in Tab.4
    uz=cd+F5*d*q/R2   ! U' in Tab.5
    vy=s -F5*y*p*q/R2 ! V  in Tab.4
    vz=t +F5*d*p*q/R2 ! V' in Tab.5
    vz=uy+sd          ! W  in Tab.4
    wz=uz+cd          ! W' in Tab.5
!-----
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(pot1/=F0) then
      du( 1)= alp(1)*q/R3    +alp(2)*x2*qR
      du( 2)= alp(1)*x/R3*sd +alp(2)*xy*qR
      du( 3)=-alp(1)*x/R3*cd +alp(2)*x*d*qR
      du( 4)= x*qR*(-alp(1) +alp(2)*(F1+A5) )
      du( 5)= alp(1)*A3/R3*sd +alp(2)*y*qR*A5
      du( 6)=-alp(1)*A3/R3*cd +alp(2)*d*qR*A5
      du( 7)= alp(1)*(sd/R3-y*qR) +alp(2)*F3*x2/R5*uy
      du( 8)= F3*x/R5*(-alp(1)*y*sd +alp(2)*(y*uy+q) )
      du( 9)= F3*x/R5*( alp(1)*y*cd +alp(2)*d*uy )
      du(10)= alp(1)*(cd/R3+d*qR) +alp(2)*F3*x2/R5*uz
      du(11)= F3*x/R5*( alp(1)*d*sd +alp(2)*y*uz )
      du(12)= F3*x/R5*(-alp(1)*d*cd +alp(2)*(d*uz-q) )
      do 222 i=1,12
        u(i)=u(i)+pot1/pi2*du(i)
 222  enddo
    endif
!===================================
!=====  dip-slip contribution  =====
!===================================
    if(pot2/=F0) then
      du( 1)=              alp(2)*x*p*qR
      du( 2)= alp(1)*s/R3 +alp(2)*y*p*qR
      du( 3)=-alp(1)*t/R3 +alp(2)*d*p*qR
      du( 4)=              alp(2)*p*qR*A5
      du( 5)=-alp(1)*F3*x*s/R5 -alp(2)*y*p*qRx
      du( 6)= alp(1)*F3*x*t/R5 -alp(2)*d*p*qRx
      du( 7)=                            alp(2)*F3*x/R5*vy
      du( 8)= alp(1)*(s2d/R3-F3*y*s/R5) +alp(2)*(F3*y/R5*vy+p*qR)
      du( 9)=-alp(1)*(c2d/R3-F3*y*t/R5) +alp(2)*F3*d/R5*vy
      du(10)=                            alp(2)*F3*x/R5*vz
      du(11)= alp(1)*(c2d/R3+F3*d*s/R5) +alp(2)*F3*y/R5*vz
      du(12)= alp(1)*(s2d/R3-F3*d*t/R5) +alp(2)*(F3*d/R5*vz-p*qR)
      do 333 i=1,12
        u(i)=u(i)+pot2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(pot3/=F0) then
      du( 1)= alp(1)*x/R3      -alp(2)*x*q*qR
      du( 2)= alp(1)*t/R3      -alp(2)*y*q*qR
      du( 3)= alp(1)*s/R3      -alp(2)*d*q*qR
      du( 4)= alp(1)*A3/R3     -alp(2)*q*qR*A5
      du( 5)=-alp(1)*F3*x*t/R5 +alp(2)*y*q*qRx
      du( 6)=-alp(1)*F3*x*s/R5 +alp(2)*d*q*qRx
      du( 7)=-alp(1)*F3*xy/R5           -alp(2)*x*qR*wy
      du( 8)= alp(1)*(c2d/R3-F3*y*t/R5) -alp(2)*(y*wy+q)*qR
      du( 9)= alp(1)*(s2d/R3-F3*y*s/R5) -alp(2)*d*qR*wy
      du(10)= alp(1)*F3*x*d/R5          -alp(2)*x*qR*wz
      du(11)=-alp(1)*(s2d/R3-F3*d*t/R5) -alp(2)*y*qR*wz
      du(12)= alp(1)*(c2d/R3+F3*d*s/R5) -alp(2)*(d*wz-q)*qR
      do 444 i=1,12
        u(i)=u(i)+pot3/pi2*du(i)
 444  enddo
    endif
!=========================================
!=====  inflate source contribution  =====
!=========================================
    if(pot4/=F0) then
      du( 1)=-alp(1)*x/R3
      du( 2)=-alp(1)*y/R3
      du( 3)=-alp(1)*d/R3
      du( 4)=-alp(1)*A3/R3
      du( 5)= alp(1)*F3*xy/R5
      du( 6)= alp(1)*F3*x*d/R5
      du( 7)= du(5)
      du( 8)=-alp(1)*B3/R3
      du( 9)= alp(1)*F3*y*d/R5
      du(10)=-du(6)
      du(11)=-du(9)
      du(12)= alp(1)*C3/R3
      do 555 i=1,12
        u(i)=u(i)+pot4/pi2*du(i)
 555  enddo
    endif
    return
end subroutine  uA0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uB0(x,y,d,z,pot1,pot2,pot3,pot4,sd,cd,alp,u)
!********************************************************************
!*****    displacement and strain at depth (part-B)             *****
!*****    due to buried point source in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   x,y,d,z : station coordinates in fault system
!*****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
!***** output
!*****   u(12) : displacement and their derivatives
    implicit none

    ! arguments
    real*8, intent(in) :: x,y,d,z
    real*8, intent(in) :: pot1,pot2,pot3,pot4
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)


    ! local variables
    real*8 :: du(12)
    real*8 :: p,q,s,t
    real*8 :: R, R2, R3, R5
    real*8 :: qR,qRx
    real*8 :: xy,x2,y2,d2,A3,A5,B3,C3
    real*8 :: uy,vy,wy,uz,vz,wz
    real*8 :: C,Rd
    real*8 :: D12,D32,D33,D53,D54
    real*8 :: FI1,FI2,FI3,FI4,FI5
    real*8 :: FJ1,FJ2,FJ3,FJ4
    real*8 :: FK1,FK2,FK3
    ! parameters
    real*8, parameter :: pi2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0, F4=4.0d0, F5=5.0d0, F8=8.0d0, F9=9.0d0
    ! loop counter
    integer :: i

    du(1:12)=0.0d0
    u(1:12)=0.0d0

    call dccon1(x,y,d,sd,cd,R,p,q,s,t)
    xy=x*y
    x2=x*x
    y2=y*y
    d2=d*d

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    qR=F3*q/R5      ! 3q/(R^5)
    qRx=F5*qR*x/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*x2/R2
    A5=F1-F5*x2/R2
    B3=F1-F3*y2/R2
    C3=F1-F3*d2/R2
!-----
    ! -- Tab.4 and 5, pp.1027-1028 --
    uy=sd-F5*y*q/R2   ! u  in Tab.4
    uz=cd+F5*d*q/R2   ! u' in Tab.5
    vy=s -F5*y*p*q/R2 ! V  in Tab.4
    vz=t +F5*d*p*q/R2 ! V' in Tab.5
    wy=uy+sd          ! W  in Tab. 4
    wz=uz+cd          ! W' in Tab.5
!-----
    C=d+z
    Rd=R+d
    D12=F1/(R*Rd*Rd)
    D32=D12*(F2*R+d)/R2
    D33=D12*(F3*R+d)/(R2*Rd)
    D53=D12*(F8*R2+F9*R*d+F3*d2)/(R2*R2*Rd)
    D54=D12*(F5*R2+F4*R*d+d2)/R3*D12
!-----
    FI1= y*(D12-x2*D33)
    FI2= x*(D12-y2*D33)
    FI3= x/R3-FI2
    FI4=-xy*D32
    FI5= F1/(R*Rd)-x2*D32
    FJ1=-F3*xy*(D33-x2*D54)
    FJ2= F1/R3-F3*D12+F3*x2*y2*D54
    FJ3= A3/R3-FJ2
    FJ4=-F3*xy/R5-FJ1
    FK1=-y*(D32-x2*D53)
    FK2=-x*(D32-y2*D53)
    FK3=-F3*x*d/R5-FK2
!-----
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(pot1/=F0) then
      du( 1)=-x2*qR  -alp(3)*FI1*sd
      du( 2)=-xy*qR  -alp(3)*FI2*sd
      du( 3)=-C*x*qR -alp(3)*FI4*sd
      du( 4)=-x*qR*(F1+A5) -alp(3)*FJ1*sd
      du( 5)=-y*qR*A5      -alp(3)*FJ2*sd
      du( 6)=-C*qR*A5      -alp(3)*FK1*sd
      du( 7)=-F3*x2/R5*uy      -alp(3)*FJ2*sd
      du( 8)=-F3*xy/R5*uy-x*qR -alp(3)*FJ4*sd
      du( 9)=-F3*C*x/R5*uy     -alp(3)*FK2*sd
      du(10)=-F3*x2/R5*uz  +alp(3)*FK1*sd
      du(11)=-F3*xy/R5*uz  +alp(3)*FK2*sd
      du(12)= F3*x/R5*(-C*uz +alp(3)*y*sd)
      do 222 i=1,12
        u(i)=u(i)+pot1/pi2*du(i)
 222  enddo
    endif
!===================================
!=====  dip-slip contribution  =====
!===================================
    if(pot2/=F0) then
      du( 1)=-x*p*qR          +alp(3)*FI3*sd*cd
      du( 2)=-y*p*qR          +alp(3)*FI1*sd*cd
      du( 3)=-C*p*qR          +alp(3)*FI5*sd*cd
      du( 4)=-p*qR*A5         +alp(3)*FJ3*sd*cd
      du( 5)= y*p*qRx         +alp(3)*FJ1*sd*cd
      du( 6)= C*p*qRx         +alp(3)*FK3*sd*cd
      du( 7)=-F3*x/R5*vy      +alp(3)*FJ1*sd*cd
      du( 8)=-F3*y/R5*vy-p*qR +alp(3)*FJ2*sd*cd
      du( 9)=-F3*C/R5*vy      +alp(3)*FK1*sd*cd
      du(10)=-F3*x/R5*vz      -alp(3)*FK3*sd*cd
      du(11)=-F3*y/R5*vz      -alp(3)*FK1*sd*cd
      du(12)=-F3*C/R5*vz      +alp(3)*A3/R3*sd*cd
      do 333 i=1,12
        u(i)=u(i)+pot2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(pot3/=F0) then
      du( 1)= x*q*qR      -alp(3)*FI3*sd*sd
      du( 2)= y*q*qR      -alp(3)*FI1*sd*sd
      du( 3)= C*q*qR      -alp(3)*FI5*sd*sd
      du( 4)= q*qR*A5     -alp(3)*FJ3*sd*sd
      du( 5)=-y*q*qRx     -alp(3)*FJ1*sd*sd
      du( 6)=-C*q*qRx     -alp(3)*FK3*sd*sd
      du( 7)= x*qR*wy     -alp(3)*FJ1*sd*sd
      du( 8)= qR*(y*wy+q) -alp(3)*FJ2*sd*sd
      du( 9)= C*qR*wy     -alp(3)*FK1*sd*sd
      du(10)= x*qR*wz     +alp(3)*FK3*sd*sd
      du(11)= y*qR*wz     +alp(3)*FK1*sd*sd
      du(12)= C*qR*wz     -alp(3)*A3/R3*sd*sd
      do 444 i=1,12
        u(i)=u(i)+pot3/pi2*du(i)
 444  enddo
    endif
!=========================================
!=====  inflate source contribution  =====
!=========================================
    if(pot4/=F0) then
      du( 1)= alp(3)*x/R3
      du( 2)= alp(3)*y/R3
      du( 3)= alp(3)*d/R3
      du( 4)= alp(3)*A3/R3
      du( 5)=-alp(3)*F3*xy/R5
      du( 6)=-alp(3)*F3*x*d/R5
      du( 7)= du(5)
      du( 8)= alp(3)*B3/R3
      du( 9)=-alp(3)*F3*y*d/R5
      du(10)=-du(6)
      du(11)=-du(9)
      du(12)=-alp(3)*C3/R3
      do 555 i=1,12
        u(i)=u(i)+pot4/pi2*du(i)
 555  enddo
    endif
    return
end subroutine uB0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uC0(x,y,d,z,pot1,pot2,pot3,pot4,sd,cd,alp,u)
!********************************************************************
!*****    displacement and strain at depth (part-B)             *****
!*****    due to buried point source in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   x,y,d,z : station coordinates in fault system
!*****   pot1-pot4 : strike-, dip-, tensile- and inflate-potency
!***** output
!*****   u(12) : displacement and their derivatives
    implicit none

    ! arguments
    real*8, intent(in) :: x,y,d,z
    real*8, intent(in) :: pot1,pot2,pot3,pot4
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)

    ! local variables
    real*8 :: du(12)
    real*8 :: p,q,s,t
    real*8 :: R, R2, R3, R5, R7
    real*8 :: s2d,c2d
    real*8 :: q2,qR,qRx,QR5,QR7,DR5
    real*8 :: xy,x2,y2,d2,C
    real*8 :: A3,A5,A7,B3,B5,B7,C3,C5,C7,D7
    ! parameters
    real*8, parameter :: pi2=6.283185307179586D0
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0, F5=5.0d0, F7=7.0d0, F10=1.0d1, F15=1.5d1
    ! loop counter
    integer :: i

    u(1:12)=0.0d0
    du(1:12)=0.0d0

    call dccon1(x,y,d,sd,cd,R,p,q,s,t)
    xy=x*y
    x2=x*x
    y2=y*y
    d2=d*d

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    qR=F3*q/R5      ! 3q/(R^5)
    qRx=F5*qR*x/R2  ! 15qx/(R^7)

! -- Eq.(8), p.1024 --
    A3=F1-F3*x2/R2
    A5=F1-F5*x2/R2
    B3=F1-F3*y2/R2
    C3=F1-F3*d2/R2
!-----
    s2d=2.0d0*sd*cd
    c2d=cd*cd-sd*sd
!-----
    C=d+z
    q2=q*q
    R7=R5*R2
    A7=F1-F7*x2/R2
    B5=F1-F5*y2/R2
    B7=F1-F7*y2/R2
    C5=F1-F5*d2/R2
    C7=F1-F7*d2/R2
    D7=F2-F7*q2/R2
    QR5=F5*q/R2
    QR7=F7*q/R2
    DR5=F5*d/R2
!-----
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(pot1/=F0) then
      du( 1)=-alp(4)*A3/R3*cd  +alp(5)*C*qR*A5
      du( 2)= F3*x/R5*( alp(4)*y*cd +alp(5)*C*(sd-y*QR5) )
      du( 3)= F3*x/R5*(-alp(4)*y*sd +alp(5)*C*(cd+d*QR5) )
      du( 4)= alp(4)*F3*x/R5*(F2+A5)*cd   -alp(5)*C*qRx*(F2+A7)
      du( 5)= F3/R5*( alp(4)*y*A5*cd +alp(5)*C*(A5*sd-y*QR5*A7) )
      du( 6)= F3/R5*(-alp(4)*y*A5*sd +alp(5)*C*(A5*cd+d*QR5*A7) )
      du( 7)= du(5)
      du( 8)= F3*x/R5*( alp(4)*B5*cd -alp(5)*F5*C/R2*(F2*y*sd+q*B7) )
      du( 9)= F3*x/R5*(-alp(4)*B5*sd +alp(5)*F5*C/R2*(d*B7*sd-y*C7*cd) )
      du(10)= F3/R5*   (-alp(4)*d*A5*cd +alp(5)*C*(A5*cd+d*QR5*A7) )
      du(11)= F15*x/R7*( alp(4)*y*d*cd  +alp(5)*C*(d*B7*sd-y*C7*cd) )
      du(12)= F15*x/R7*(-alp(4)*y*d*sd  +alp(5)*C*(F2*d*cd-q*C7) )
      do 222 i=1,12
        u(i)=u(i)+pot1/pi2*du(i)
 222  enddo
    endif
!===================================
!=====  dip-slip contribution  =====
!===================================
    if(pot2/=F0) then
      du( 1)= alp(4)*F3*x*t/R5          -alp(5)*C*p*qRx
      du( 2)=-alp(4)/R3*(c2d-F3*y*t/R2) +alp(5)*F3*C/R5*(s-y*p*QR5)
      du( 3)=-alp(4)*A3/R3*sd*cd        +alp(5)*F3*C/R5*(t+d*p*QR5)
      du( 4)= alp(4)*F3*t/R5*A5              -alp(5)*F5*C*p*qR/R2*A7
      du( 5)= F3*x/R5*(alp(4)*(c2d-F5*y*t/R2)-alp(5)*F5*C/R2*(s-y*p*QR7))
      du( 6)= F3*x/R5*(alp(4)*(F2+A5)*sd*cd  -alp(5)*F5*C/R2*(t+d*p*QR7))
      du( 7)= du(5)
      du( 8)= F3/R5*( alp(4)*(F2*y*c2d+t*B5) &
                     +alp(5)*C*(s2d-F10*y*s/R2-p*QR5*B7))
      du( 9)= F3/R5*(alp(4)*y*A5*sd*cd-alp(5)*C*((F3+A5)*c2d+y*p*DR5*QR7))
      du(10)= F3*x/R5*(-alp(4)*(s2d-t*DR5) -alp(5)*F5*C/R2*(t+d*p*QR7))
      du(11)= F3/R5*(-alp(4)*(d*B5*c2d+y*C5*s2d) &
                     -alp(5)*C*((F3+A5)*c2d+y*p*DR5*QR7))
      du(12)= F3/R5*(-alp(4)*d*A5*sd*cd-alp(5)*C*(s2d-F10*d*t/R2+p*QR5*C7))
      do 333 i=1,12
        u(i)=u(i)+pot2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(pot3/=F0) then
      du( 1)= F3*x/R5*(-alp(4)*s +alp(5)*(C*q*QR5-z))
      du( 2)= alp(4)/R3*(s2d-F3*y*s/R2)+alp(5)*F3/R5*(C*(t-y+y*q*QR5)-y*z)
      du( 3)=-alp(4)/R3*(F1-A3*sd*sd)  -alp(5)*F3/R5*(C*(s-d+d*q*QR5)-d*z)
      du( 4)=-alp(4)*F3*s/R5*A5 +alp(5)*(C*qR*QR5*A7-F3*z/R5*A5)
      du( 5)= F3*x/R5*(-alp(4)*(s2d-F5*y*s/R2) &
                       -alp(5)*F5/R2*(C*(t-y+y*q*QR7)-y*z))
      du( 6)= F3*x/R5*( alp(4)*(F1-(F2+A5)*sd*sd) &
                       +alp(5)*F5/R2*(C*(s-d+d*q*QR7)-d*z))
      du( 7)= du(5)
      du( 8)= F3/R5*(-alp(4)*(F2*y*s2d+s*B5) &
                     -alp(5)*(C*(F2*sd*sd+F10*y*(t-y)/R2-q*QR5*B7)+z*B5))
      du( 9)= F3/R5*( alp(4)*y*(F1-A5*sd*sd) &
                     +alp(5)*(C*(F3+A5)*s2d-y*DR5*(C*D7+z)))
      du(10)= F3*x/R5*(-alp(4)*(c2d+s*DR5) &
                       +alp(5)*(F5*C/R2*(s-d+d*q*QR7)-F1-z*DR5))
      du(11)= F3/R5*( alp(4)*(d*B5*s2d-y*C5*c2d) &
                     +alp(5)*(C*((F3+A5)*s2d-y*DR5*D7)-y*(F1+z*DR5)))
      du(12)= F3/R5*(-alp(4)*d*(F1-A5*sd*sd) &
                     -alp(5)*(C*(c2d+F10*d*(s-d)/R2-q*QR5*C7)+z*(F1+C5)))
      do 444 i=1,12
        u(i)=u(i)+pot3/pi2*du(i)
 444  enddo
    endif
!=========================================
!=====  inflate source contribution  =====
!=========================================
    if(pot4/=F0) then
      du( 1)= alp(4)*F3*x*d/R5
      du( 2)= alp(4)*F3*y*d/R5
      du( 3)= alp(4)*C3/R3
      du( 4)= alp(4)*F3*d/R5*A5
      du( 5)=-alp(4)*F15*xy*d/R7
      du( 6)=-alp(4)*F3*x/R5*C5
      du( 7)= du(5)
      du( 8)= alp(4)*F3*d/R5*B5
      du( 9)=-alp(4)*F3*y/R5*C5
      du(10)= du(6)
      du(11)= du(9)
      du(12)= alp(4)*F3*d/R5*(F2+C5)
      do 555 i=1,12
        u(i)=u(i)+pot4/pi2*du(i)
 555  enddo
    endif
    return
end subroutine uC0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine DC3D(alpha,x,y,z,depth,dip, &
                AL1,AL2,AW1,AW2,disl1,disl2,disl3, &
                ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
!********************************************************************
!*****                                                          *****
!*****    displacement and strain at depth                      *****
!*****    due to buried finite fault in a semi-infinite medium  *****
!*****              CODED BY  y.OKADA ... SEP.1991              *****
!*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****
!*****                          JUL.1993, MAY.2002              *****
!--------------------------------------------------------------------
!-----                     converted to Fortran90 (free form)   -----
!-----                                t. Miyashita, Jan. 2020   -----
!--------------------------------------------------------------------
!********************************************************************
!
!***** input
!*****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
!*****   x,y,z : coordinate of observing point
!*****   depth : depth of reference point
!*****   dip   : dip-angle (degree)
!*****   AL1,AL2   : fault length range
!*****   AW1,AW2   : fault width range
!*****   disl1-disl3 : strike-, dip-, tensile-dislocations
!
!***** output
!*****   ux, uy, uz  : displacement ( unit=(unit of dislocation) )
!*****   uxx,uyx,uzx : x-derivative ( unit=(unit of dislocation) /
!*****   uxy,uyy,uzy : y-derivative        (unit of x,y,z,depth,AL,AW) )
!*****   uxz,uyz,uzz : z-derivative
!*****   iret        : return code
!*****               :   =0....normal
!*****               :   =1....singular
!*****               :   =2....positive z was given
    implicit none

    ! arguments
    real*8, intent(in) :: alpha,x,y,z,depth,dip,AL1,AL2,AW1,AW2,disl1,disl2,disl3
    real*8, intent(out) :: ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer, intent(out) :: iret

    ! local variables
    real*8 :: sd,cd
    real*8 :: xi(2),et(2)
    integer :: kxi(2),ket(2)
    real*8 :: u(12),du(12),duA(12),duB(12),duC(12)
    real*8 :: dd1, dd2, dd3
    real*8 :: zz, d, p, q, R12, R21, R22
    real*8 :: alp(5)
    ! parameters
    real*8, parameter :: F0 = 0.0d0
    real*8, parameter :: eps = 1.0d-6
    ! loop counters
    integer :: i,j,k

    ! initialization
    u(1:12) = 0.0d0
    du(1:12) = 0.0d0
    duA(1:12) = 0.0d0
    duB(1:12) = 0.0d0
    duC(1:12) = 0.0d0
    call varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)

!-----
      iret=0
    if(z>0.0d0) then
      iret=2
      return
    endif
!-----
    zz=z
    dd1=disl1
    dd2=disl2
    dd3=disl3
!### caution ### if xi,et,q are sufficiently small, they are set to zero
    xi(1)=x-AL1
    xi(2)=x-AL2
    if(dabs(xi(1))<eps) xi(1)=F0
    if(dabs(xi(2))<eps) xi(2)=F0
!======================================
!=====  real-source contribution  =====
!======================================
    call assign_alpha(alpha,alp)
    call dccon0(dip,sd,cd)
!-----
    d=depth+z
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-AW1
    et(2)=p-AW2
    if(dabs(q)<eps)  q=F0
    if(dabs(et(1))<eps) et(1)=F0
    if(dabs(et(2))<eps) et(2)=F0
!--------------------------------
!----- reject singular case -----
!--------------------------------
!----- on fault edge
    if(q==F0 .and. &
      (    (xi(1)*xi(2)<=F0 .and. et(1)*et(2)==F0) &
       .or.(et(1)*et(2)<=F0 .and. xi(1)*xi(2)==F0) )) then
      iret=1
      return
    endif
!----- on negative extension of fault edge
    kxi(1)=0
    kxi(2)=0
    ket(1)=0
    ket(2)=0
    R12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    R21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    R22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
    if(xi(1)<F0 .and. R21+xi(2)<eps) kxi(1)=1
    if(xi(1)<F0 .and. R22+xi(2)<eps) kxi(2)=1
    if(et(1)<F0 .and. R12+et(2)<eps) ket(1)=1
    if(et(1)<F0 .and. R22+et(2)<eps) ket(2)=1
!=====
    do 223 k=1,2
    do 222 j=1,2
      call uA(xi(j),et(k),q,dd1,dd2,dd3,kxi(k),ket(j),sd,cd,alp,duA)
!-----
      do 220 i=1,10,3
        du(i)  =-duA(i)
        du(i+1)=-duA(i+1)*cd+duA(i+2)*sd
        du(i+2)=-duA(i+1)*sd-duA(i+2)*cd
        if(i<10) cycle
        du(i)  =-du(i)
        du(i+1)=-du(i+1)
        du(i+2)=-du(i+2)
220   enddo
      do 221 i=1,12
        if(j+k/=3) u(i)=u(i)+du(i)
        if(j+k==3) u(i)=u(i)-du(i)
221   enddo
!-----
222 enddo
223 enddo
!=======================================
!=====  image-source contribution  =====
!=======================================
    d=depth-z
    p=y*cd+d*sd
    q=y*sd-d*cd
    et(1)=p-AW1
    et(2)=p-AW2
    if(dabs(q)<eps)  q=F0
    if(dabs(et(1))<eps) et(1)=F0
    if(dabs(et(2))<eps) et(2)=F0
!--------------------------------
!----- reject singular case -----
!--------------------------------
!----- on fault edge
    if(q==F0 .and. &
      (    (xi(1)*xi(2)<=F0 .and. et(1)*et(2)==F0) &
       .or.(et(1)*et(2)<=F0 .and. xi(1)*xi(2)==F0) )) then
      iret=1
      call varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)
      return
    endif
!----- on negative extension of fault edge
    kxi(1)=0
    kxi(2)=0
    ket(1)=0
    ket(2)=0
    R12=dsqrt(xi(1)*xi(1)+et(2)*et(2)+q*q)
    R21=dsqrt(xi(2)*xi(2)+et(1)*et(1)+q*q)
    R22=dsqrt(xi(2)*xi(2)+et(2)*et(2)+q*q)
    if(xi(1)<F0 .and. R21+xi(2)<eps) kxi(1)=1
    if(xi(1)<F0 .and. R22+xi(2)<eps) kxi(2)=1
    if(et(1)<F0 .and. R12+et(2)<eps) ket(1)=1
    if(et(1)<F0 .and. R22+et(2)<eps) ket(2)=1
!=====
    do 334 k=1,2
    do 333 j=1,2
      call uA(xi(j),et(k),q   ,dd1,dd2,dd3,kxi(k),ket(j),sd,cd,alp,duA)
      call uB(xi(j),et(k),q   ,dd1,dd2,dd3,kxi(k),ket(j),sd,cd,alp,duB)
      call uC(xi(j),et(k),q,zz,dd1,dd2,dd3,kxi(k),ket(j),sd,cd,alp,duC)
!-----
      do 330 i=1,10,3
        du(i)=duA(i)+duB(i)+z*duC(i)
        du(i+1)=(duA(i+1)+duB(i+1)+z*duC(i+1))*cd &
               -(duA(i+2)+duB(i+2)+z*duC(i+2))*sd
        du(i+2)=(duA(i+1)+duB(i+1)-z*duC(i+1))*sd &
               +(duA(i+2)+duB(i+2)-z*duC(i+2))*cd
        if(i<10) cycle
        du(10)=du(10)+duC(1)
        du(11)=du(11)+duC(2)*cd-duC(3)*sd
        du(12)=du(12)-duC(2)*sd-duC(3)*cd
330   enddo
      do 331 i=1,12
        if(j+k/=3) u(i)=u(i)+du(i)
        if(j+k==3) u(i)=u(i)-du(i)
331   enddo
!-----
333 enddo
334 enddo
!=====
    call varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)
    return
end subroutine DC3D
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uA(xi,et,q,disl1,disl2,disl3,kxi,ket,sd,cd,alp,u)
!
!********************************************************************
!*****    displacement and strain at depth (part-A)             *****
!*****    due to buried finite fault in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   xi,et,q : station coordinates in fault system
!*****   disl1-disl3 : strike-, dip-, tensile-dislocations
!*****   kxi,ket     : kxi=1, ket=1 means R+xi<eps, R+et<eps, respectively
!*****   sd,cd       : sin, cos of dip-angle
!*****   alp         : coefficients related to alpha value
!***** output
!*****   u(12) : displacement and their derivatives
!-----
    implicit none
    ! arguments
    real*8, intent(in) :: xi,et,q,disl1,disl2,disl3
    integer, intent(in) :: kxi,ket
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)
    ! local variables
    real*8 :: du(12)
    real*8 :: R,ytilde,dtilde
    real*8 :: THETA,ALX,ALE
    real*8 :: X11,X32,Y11,Y32
    real*8 :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    real*8 :: xy,qX,qY
    real*8 :: xi2,q2
    real*8 :: R3
    ! parameters
    real*8, parameter :: F0=0.0d0, F2=2.0d0
    real*8, parameter :: pi2=6.283185307179586d0
    ! loop counter
    integer :: i

    u(1:12) = 0.0d0
    du(1:12) = 0.0d0

    call dccon2(xi,et,q,sd,cd,R,ytilde,dtilde)
    call math_singularity(xi,et,q,R,kxi,ket,THETA,ALX,ALE)
    call eq14_xy(xi,et,kxi,ket,R,X11,X32,Y11,Y32)
    call auxvar_derivative_yz(xi,q,sd,cd,R,ytilde,dtilde,X11,X32,Y32, &
                              EY,EZ,FY,FZ,GY,GZ,HY,HZ)

    R3=R**3.0d0
    xi2=xi*xi
    q2=q*q

    xy=xi*Y11
    qX=q *X11
    qY=q *Y11
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(disl1/=F0) then
      du( 1)= THETA/F2     +alp(2)*xi*qY
      du( 2)=               alp(2)*q/R
      du( 3)= alp(1)*ALE   -alp(2)*q*qY
      du( 4)=-alp(1)*qY    -alp(2)*xi2*q*Y32
      du( 5)=              -alp(2)*xi*q/R3
      du( 6)= alp(1)*xy    +alp(2)*xi*q2*Y32
      du( 7)= alp(1)*xy*sd +alp(2)*xi*FY+dtilde/F2*X11
      du( 8)=               alp(2)*EY
      du( 9)= alp(1)*(cd/R+qY*sd) -alp(2)*q*FY
      du(10)= alp(1)*xy*cd        +alp(2)*xi*FZ+ytilde/F2*X11
      du(11)=               alp(2)*EZ
      du(12)=-alp(1)*(sd/R-qY*cd) -alp(2)*q*FZ
      do 222 i=1,12
        u(i)=u(i)+disl1/pi2*du(i)
 222  enddo
    endif
!======================================
!=====    dip-slip contribution   =====
!======================================
    if(disl2/=F0) then
      du( 1)=             alp(2)*q/R
      du( 2)= THETA/F2   +alp(2)*et*qX
      du( 3)= alp(1)*ALX -alp(2)*q*qX
      du( 4)=            -alp(2)*xi*q/R3
      du( 5)= -qY/F2     -alp(2)*et*q/R3
      du( 6)= alp(1)/R   +alp(2)*q2/R3
      du( 7)=             alp(2)*EY
      du( 8)= alp(1)*dtilde*X11+xy/F2*sd +alp(2)*et*GY
      du( 9)= alp(1)*ytilde*X11          -alp(2)*q*GY
      du(10)=             alp(2)*EZ
      du(11)= alp(1)*ytilde*X11+xy/F2*cd +alp(2)*et*GZ
      du(12)=-alp(1)*dtilde*X11          -alp(2)*q*GZ
      do 333 i=1,12
        u(i)=u(i)+disl2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(disl3/=F0) then
      du( 1)=-alp(1)*ALE -alp(2)*q*qY
      du( 2)=-alp(1)*ALX -alp(2)*q*qX
      du( 3)= THETA/F2   -alp(2)*(et*qX+xi*qY)
      du( 4)=-alp(1)*xy  +alp(2)*xi*q2*Y32
      du( 5)=-alp(1)/R   +alp(2)*q2/R3
      du( 6)=-alp(1)*qY  -alp(2)*q*q2*Y32
      du( 7)=-alp(1)*(cd/R+qY*sd)  -alp(2)*q*FY
      du( 8)=-alp(1)*ytilde*X11         -alp(2)*q*GY
      du( 9)= alp(1)*(dtilde*X11+xy*sd) +alp(2)*q*HY
      du(10)= alp(1)*(sd/R-qY*cd)  -alp(2)*q*FZ
      du(11)= alp(1)*dtilde*X11         -alp(2)*q*GZ
      du(12)= alp(1)*(ytilde*X11+xy*cd) +alp(2)*q*HZ
      do 444 i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
 444  enddo
    endif
    return
end subroutine uA
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uB(xi,et,q,disl1,disl2,disl3,kxi,ket,sd,cd,alp,u)
!********************************************************************
!*****    displacement and strain at depth (part-B)             *****
!*****    due to buried finite fault in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   xi,et,q : station coordinates in fault system
!*****   disl1-disl3 : strike-, dip-, tensile-dislocations
!*****   kxi,ket     : kxi=1, ket=1 means R+xi<eps, R+et<eps, respectively
!*****   sd,cd       : sin, cos of dip-angle
!*****   alp         : coefficients related to alpha value
!***** output
!*****   u(12) : displacement and their derivatives
!-----
    implicit none
    ! arguments
    real*8, intent(in) :: xi,et,q,disl1,disl2,disl3
    integer, intent(in) :: kxi,ket
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)
    ! local variables
    real*8 :: du(12)
    real*8 :: R,ytilde,dtilde
    real*8 :: THETA,ALE
    real*8 :: X11,X32,Y11,Y32
    real*8 :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    real*8 :: x
    real*8 :: xy,qX,qY
    real*8 :: xi2,q2
    real*8 :: R3
    real*8 :: Rd,Rd2
    real*8 :: D11
    real*8 :: AI1,AI2,AI3,AI4
    real*8 :: AJ1,AJ2,AJ3,AJ4,AJ5,AJ6
    real*8 :: AK1,AK2,AK3,AK4
    ! parameters
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0
    real*8, parameter :: pi2=6.283185307179586d0
    ! loop counter
    integer :: i
    ! unused (dummy) vars
    real*8 :: dum ! ALX

    u(1:12)=0.0d0
    du(1:12)=0.0d0

    call dccon2(xi,et,q,sd,cd,R,ytilde,dtilde)
    call math_singularity(xi,et,q,R,kxi,ket,THETA,dum,ALE)
    call eq14_xy(xi,et,kxi,ket,R,X11,X32,Y11,Y32)
    call auxvar_derivative_yz(xi,q,sd,cd,R,ytilde,dtilde,X11,X32,Y32, &
                              EY,EZ,FY,FZ,GY,GZ,HY,HZ)

    R3=R**3.0d0
    xi2=xi*xi
    q2=q*q

    Rd=R+dtilde
    D11=F1/(R*Rd)
    AJ2=xi*ytilde/Rd*D11
    AJ5=-(dtilde+ytilde*ytilde/Rd)*D11
    if(cd/=F0) then
      ! p.1034 (ii)
      if(xi==F0) then
        AI4=F0
      else
        x=dsqrt(xi2+q2)
        AI4=F1/(cd*cd)*( xi/Rd*sd*cd &
           +F2*datan((et*(x+q*cd)+x*(R+x)*sd)/(xi*(R+x)*cd)) )
      endif
      AI3=(ytilde*cd/Rd-ALE+sd*dlog(Rd))/(cd*cd)
      AK1=xi*(D11-Y11*sd)/cd
      AK3=(q*Y11-ytilde*D11)/cd
      AJ3=(AK1-AJ2*sd)/cd
      AJ6=(AK3-AJ5*sd)/cd
    else
      Rd2=Rd*Rd
      AI3=(et/Rd+ytilde*q/Rd2-ALE)/F2
      AI4=xi*ytilde/Rd2/F2
      AK1=xi*q/Rd*D11
      AK3=sd/Rd*(xi2*D11-F1)
      AJ3=-xi/Rd2*(q2*D11-F1/F2)
      AJ6=-ytilde/Rd2*(xi2*D11-F1/F2)
    endif
!-----
    xy=xi*Y11
    AI1=-xi/Rd*cd-AI4*sd
    AI2= dlog(Rd)+AI3*sd
    AK2= F1/R+AK3*sd
    AK4= xy*cd-AK1*sd
    AJ1= AJ5*cd-AJ6*sd
    AJ4=-xy-AJ2*cd+AJ3*sd
!=====
    qX=q*X11
    qY=q*Y11
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(disl1/=F0) then
      du( 1)=-xi*qY-THETA -alp(3)*AI1*sd
      du( 2)=-q/R         +alp(3)*ytilde/Rd*sd
      du( 3)= q*qY        -alp(3)*AI2*sd
      du( 4)= xi2*q*Y32   -alp(3)*AJ1*sd
      du( 5)= xi*q/R3     -alp(3)*AJ2*sd
      du( 6)=-xi*q2*Y32   -alp(3)*AJ3*sd
      du( 7)=-xi*FY-dtilde*X11 +alp(3)*(xy+AJ4)*sd
      du( 8)=-EY          +alp(3)*(F1/R+AJ5)*sd
      du( 9)= q*FY        -alp(3)*(qY-AJ6)*sd
      du(10)=-xi*FZ-ytilde*X11 +alp(3)*AK1*sd
      du(11)=-EZ          +alp(3)*ytilde*D11*sd
      du(12)= q*FZ        +alp(3)*AK2*sd
      do 222 i=1,12
        u(i)=u(i)+disl1/pi2*du(i)
 222  enddo
    endif
!======================================
!=====    dip-slip contribution   =====
!======================================
    if(disl2/=F0) then
      du( 1)=-q/R         +alp(3)*AI3*sd*cd
      du( 2)=-et*qX-THETA -alp(3)*xi/Rd*sd*cd
      du( 3)= q*qX        +alp(3)*AI4*sd*cd
      du( 4)= xi*q/R3     +alp(3)*AJ4*sd*cd
      du( 5)= et*q/R3+qY  +alp(3)*AJ5*sd*cd
      du( 6)=-q2/R3       +alp(3)*AJ6*sd*cd
      du( 7)=-EY          +alp(3)*AJ1*sd*cd
      du( 8)=-et*GY-xy*sd +alp(3)*AJ2*sd*cd
      du( 9)= q*GY        +alp(3)*AJ3*sd*cd
      du(10)=-EZ          -alp(3)*AK3*sd*cd
      du(11)=-et*GZ-xy*cd -alp(3)*xi*D11*sd*cd
      du(12)= q*GZ        -alp(3)*AK4*sd*cd
      do 333 i=1,12
        u(i)=u(i)+disl2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(disl3/=F0) then
      du( 1)= q*qY              -alp(3)*AI3*sd*sd
      du( 2)= q*qX              +alp(3)*xi/Rd*sd*sd
      du( 3)= et*qX+xi*qY-THETA -alp(3)*AI4*sd*sd
      du( 4)=-xi*q2*Y32 -alp(3)*AJ4*sd*sd
      du( 5)=-q2/R3     -alp(3)*AJ5*sd*sd
      du( 6)= q*q2*Y32  -alp(3)*AJ6*sd*sd
      du( 7)= q*FY -alp(3)*AJ1*sd*sd
      du( 8)= q*GY -alp(3)*AJ2*sd*sd
      du( 9)=-q*HY -alp(3)*AJ3*sd*sd
      du(10)= q*FZ +alp(3)*AK3*sd*sd
      du(11)= q*GZ +alp(3)*xi*D11*sd*sd
      du(12)=-q*HZ +alp(3)*AK4*sd*sd
      do 444 i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
 444  enddo
    endif
    return
end subroutine uB
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  uC(xi,et,q,z,disl1,disl2,disl3,kxi,ket,sd,cd,alp,u)
!********************************************************************
!*****    displacement and strain at depth (part-C)             *****
!*****    due to buried finite fault in a semi-infinite medium  *****
!********************************************************************
!
!***** input
!*****   xi,et,q,z   : station coordinates in fault system
!*****   disl1-disl3 : strike-, dip-, tensile-dislocations
!*****   kxi,ket     : kxi=1, ket=1 means R+xi<eps, R+et<eps, respectively
!*****   sd,cd       : sin, cos of dip-angle
!*****   alp         : coefficients related to alpha value
!***** output
!*****   u(12) : displacement and their derivatives
    implicit none
    ! arguments
    real*8, intent(in) :: xi,et,z,q,disl1,disl2,disl3
    integer, intent(in) :: kxi,ket
    real*8, intent(in) :: sd,cd
    real*8, intent(in) :: alp(5)
    real*8, intent(out) :: u(12)
    ! local variables
    real*8 :: du(12)
    real*8 :: R,ytilde,dtilde
    real*8 :: X11,Y11,X32,Y32
    real*8 :: xy,qX,qY
    real*8 :: xi2,ET2,q2
    real*8 :: R2,R3,R5
    real*8 :: C
    real*8 :: X53,Y53,Z32,Z53,H
    real*8 :: Y0,Z0
    real*8 :: PPY,PPZ
    real*8 :: QQ,QQY,QQZ,qR,CqX,CDR,YY0
    ! parameters
    real*8, parameter :: F0=0.0d0, F1=1.0d0, F2=2.0d0, F3=3.0d0
    real*8, parameter :: pi2=6.283185307179586d0
    ! loop counter
    integer :: i

    u(1:12) = 0.0d0
    du(1:12) = 0.0d0

    call dccon2(xi,et,q,sd,cd,R,ytilde,dtilde)
    call eq14_xy(xi,et,kxi,ket,R,X11,X32,Y11,Y32)

    R2=R**2.0d0
    R3=R**3.0d0
    R5=R**5.0d0
    xi2=xi*xi
    ET2=et*et
    q2=q*q
!-----
    C=dtilde+z
    X53=(8.d0*R2+9.d0*R*xi+F3*xi2)*X11*X11*X11/R2
    Y53=(8.d0*R2+9.d0*R*et+F3*ET2)*Y11*Y11*Y11/R2
    H=q*cd-z
    Z32=sd/R3-H*Y32
    Z53=F3*sd/R5-H*Y53
    Y0=Y11-xi2*Y32
    Z0=Z32-xi2*Z53
    PPY=cd/R3+q*Y32*sd
    PPZ=sd/R3-q*Y32*cd
    QQ=z*Y32+Z32+Z0
    QQY=F3*C*dtilde/R5-QQ*sd
    QQZ=F3*C*ytilde/R5-QQ*cd+q*Y32
    xy=xi*Y11
    qX=q*X11
    qY=q*Y11
    qR=F3*q/R5
    CqX=C*q*X53
    CDR=(C+dtilde)/R3
    YY0=ytilde/R3-Y0*cd
!=====
!======================================
!=====  strike-slip contribution  =====
!======================================
    if(disl1/=F0) then
      du( 1)= alp(4)*xy*cd           -alp(5)*xi*q*Z32
      du( 2)= alp(4)*(cd/R+F2*qY*sd) -alp(5)*C*q/R3
      du( 3)= alp(4)*qY*cd           -alp(5)*(C*et/R3-z*Y11+xi2*Z32)
      du( 4)= alp(4)*Y0*cd                  -alp(5)*q*Z0
      du( 5)=-alp(4)*xi*(cd/R3+F2*q*Y32*sd) +alp(5)*C*xi*qR
      du( 6)=-alp(4)*xi*q*Y32*cd  +alp(5)*xi*(F3*C*et/R5-QQ)
      du( 7)=-alp(4)*xi*PPY*cd    -alp(5)*xi*QQY
      du( 8)= alp(4)*F2*(dtilde/R3-Y0*sd)*sd-ytilde/R3*cd &
             -alp(5)*(CDR*sd-et/R3-C*ytilde*qR)
      du( 9)=-alp(4)*q/R3+YY0*sd  +alp(5)*(CDR*cd+C*dtilde*qR-(Y0*cd+q*Z0)*sd)
      du(10)= alp(4)*xi*PPZ*cd    -alp(5)*xi*QQZ
      du(11)= alp(4)*F2*(ytilde/R3-Y0*cd)*sd+dtilde/R3*cd -alp(5)*(CDR*cd+C*dtilde*qR)
      du(12)=        YY0*cd    -alp(5)*(CDR*sd-C*ytilde*qR-Y0*sd*sd+q*Z0*cd)
      do 222 i=1,12
        u(i)=u(i)+disl1/pi2*du(i)
 222  enddo
    endif
!======================================
!=====    dip-slip contribution   =====
!======================================
    if(disl2/=F0) then
      du( 1)= alp(4)*cd/R -qY*sd      -alp(5)*C*q/R3
      du( 2)= alp(4)*ytilde*X11       -alp(5)*C*et*q*X32
      du( 3)=       -dtilde*X11-xy*sd -alp(5)*C*(X11-q2*X32)
      du( 4)=-alp(4)*xi/R3*cd         +alp(5)*C*xi*qR +xi*q*Y32*sd
      du( 5)=-alp(4)*ytilde/R3        +alp(5)*C*et*qR
      du( 6)=        dtilde/R3-Y0*sd  +alp(5)*C/R3*(F1-F3*q2/R2)
      du( 7)=-alp(4)*et/R3+Y0*sd*sd   -alp(5)*(CDR*sd-C*ytilde*qR)
      du( 8)= alp(4)*(X11-ytilde*ytilde*X32) -alp(5)*C*((dtilde+F2*q*cd)*X32-ytilde*et*q*X53)
      du( 9)=   xi*PPY*sd+ytilde*dtilde*X32  +alp(5)*C*((ytilde+F2*q*sd)*X32-ytilde*q2*X53)
      du(10)=   -q/R3+Y0*sd*cd               -alp(5)*(CDR*cd+C*dtilde*qR)
      du(11)= alp(4)*ytilde*dtilde*X32       -alp(5)*C*((ytilde-F2*q*sd)*X32+dtilde*et*q*X53)
      du(12)=-xi*PPZ*sd+X11-dtilde*dtilde*X32-alp(5)*C*((dtilde-F2*q*cd)*X32-dtilde*q2*X53)
      do 333 i=1,12
        u(i)=u(i)+disl2/pi2*du(i)
 333  enddo
    endif
!========================================
!=====  tensile-fault contribution  =====
!========================================
    if(disl3/=F0) then
      du( 1)=-alp(4)*(sd/R+qY*cd)        -alp(5)*(z*Y11-q2*Z32)
      du( 2)= alp(4)*F2*xy*sd+dtilde*X11 -alp(5)*C*(X11-q2*X32)
      du( 3)= alp(4)*(ytilde*X11+xy*cd)  +alp(5)*q*(C*et*X32+xi*Z32)
      du( 4)= alp(4)*xi/R3*sd+xi*q*Y32*cd+alp(5)*xi*(F3*C*et/R5-F2*Z32-Z0)
      du( 5)= alp(4)*F2*Y0*sd-dtilde/R3  +alp(5)*C/R3*(F1-F3*q2/R2)
      du( 6)=-alp(4)*YY0                 -alp(5)*(C*et*qR-q*Z0)
      du( 7)= alp(4)*(q/R3+Y0*sd*cd)     +alp(5)*(z/R3*cd+C*dtilde*qR-q*Z0*sd)
      du( 8)=-alp(4)*F2*xi*PPY*sd-ytilde*dtilde*X32 &
             +alp(5)*C*((ytilde+F2*q*sd)*X32-ytilde*q2*X53)
      du( 9)=-alp(4)*(xi*PPY*cd-X11+ytilde*ytilde*X32) &
             +alp(5)*(C*((dtilde+F2*q*cd)*X32-ytilde*et*q*X53)+xi*QQY)
      du(10)=  -et/R3+Y0*cd*cd           -alp(5)*(z/R3*sd-C*ytilde*qR-Y0*sd*sd+q*Z0*cd)
      du(11)= alp(4)*F2*xi*PPZ*sd-X11+dtilde*dtilde*X32 &
             -alp(5)*C*((dtilde-F2*q*cd)*X32-dtilde*q2*X53)
      du(12)= alp(4)*(xi*PPZ*cd+ytilde*dtilde*X32) &
             +alp(5)*(C*((ytilde-F2*q*sd)*X32+dtilde*et*q*X53)+xi*QQZ)
      do 444 i=1,12
        u(i)=u(i)+disl3/pi2*du(i)
 444  enddo
    endif
    return
end subroutine uC
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  dccon0(dip,sd,cd)
!*******************************************************************
!*****   calculate medium constants and fault-dip constants    *****
!*******************************************************************
!
!***** input
!*****   dip   : dip-angle (degree)
!***** output
!*****   sd,cd : sin, cos of dip-angle
    implicit none
    real*8, intent(in) :: dip
    real*8, intent(out) :: sd,cd

    real*8, parameter :: F0=0.0d0, F1=1.0d0
    real*8, parameter :: pi2=6.283185307179586d0
    real*8, parameter :: eps=1.0d-6
    real*8, parameter :: p18=pi2/360.d0
!-----
    sd=dsin(dip*p18)
    cd=dcos(dip*p18)
!### caution ### if cos(dip) is sufficiently small, it is set to zero
    if(dabs(cd)<eps) then
      cd=F0
      if(sd>F0) sd= F1
      if(sd<F0) sd=-F1
    endif
    return
end subroutine dccon0
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  dccon1(x,y,d,sd,cd,R,p,q,s,t)
!**********************************************************************
!*****   calculate station geometry constants for point source    *****
!**********************************************************************
!***** input
!*****   x,y,d : station coordinates in fault system
!*****   sd,cd : sin, cos of dip-angle
    implicit none
    real*8, intent(in) :: x,y,d
    real*8, intent(in) :: sd,cd
    real*8, intent(out) :: R
    real*8, intent(out) :: p,q,s,t
!-----
    p=y*cd+d*sd
    q=y*sd-d*cd
    s=p*sd+q*cd
    t=p*cd-q*sd
    R=dsqrt(x*x+y*y+d*d)
    return
end subroutine dccon1
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine  dccon2(xi,et,q,sd,cd,R,ytilde,dtilde)
!**********************************************************************
!*****   calculate station geometry constants for finite source   *****
!**********************************************************************
!***** input
!*****   xi,et,q : station coordinates in fault system
!*****   sd,cd   : sin, cos of dip-angle
    implicit none
    real*8, intent(in) :: xi,et,q
    real*8, intent(in) :: sd,cd
    real*8, intent(out) :: R,ytilde,dtilde
!-----
    R = dsqrt(xi*xi+et*et+q*q)
    ytilde = et*cd+q*sd
    dtilde = et*sd-q*cd
    return
end subroutine dccon2
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine math_singularity(xi,et,q,R,kxi,ket,THETA,ALX,ALE)
!***** input
!*****   xi,et,q : station coordinates in fault system
!*****   kxi,ket : kxi=1, ket=1 means R+xi<eps, R+et<eps, respectively
    implicit none
    real*8, intent(in) :: xi,et,q,R
    integer, intent(in) :: kxi,ket
    real*8, intent(out) :: THETA,ALX,ALE

    ! p.1034 (i)
    if(q==0.0d0) then
      THETA=0.0d0
    else
      THETA=datan(xi*et/(q*R))
    endif
    ! p.1034 (iii)
    if(kxi==1) then
      ALX=-dlog(R-xi)
    else
      ALX=dlog(R+xi)
    endif
    ! p.1034 (iv)
    if(ket==1) then
      ALE=-dlog(R-et)
    else
      ALE=dlog(R+et)
    endif

    return
end subroutine math_singularity
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine eq14_xy(xi,et,kxi,ket,R,X11,X32,Y11,Y32)
!***** input
!*****   xi,et   : station coordinates in fault system
!*****   kxi,ket : kxi=1, ket=1 means R+xi<eps, R+et<eps, respectively
    implicit none
    real*8, intent(in) :: xi,et
    integer, intent(in) :: kxi,ket
    real*8, intent(in) :: R
    real*8, intent(out) :: X11,X32,Y11,Y32

    X11 = 0.0d0
    X32 = 0.0d0
    Y11 = 0.0d0
    Y32 = 0.0d0
!-----
    if(kxi/=1) then
      X11=1.0d0/(R*(R+xi))        ! Eq.(14), p.1026
      X32=(2.0d0*R+xi)*X11*X11/R  ! Eq.(14), p.1026
    endif
!-----
    if(ket/=1) then
      Y11=1.0d0/(R*(R+et))        ! Eq.(14), p.1026
      Y32=(2.0d0*R+et)*Y11*Y11/R  ! Eq.(14), p.1026
    endif
!-----
    return
end subroutine eq14_xy
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine auxvar_derivative_yz(xi,q,sd,cd,R,ytilde,dtilde,X11,X32,Y32, &
                                EY,EZ,FY,FZ,GY,GZ,HY,HZ)
!***** input
!*****   xi,q  : station coordinates in fault system
!*****   sd,cd : sin, cos of dip-angle
    implicit none
    ! arguments
    real*8, intent(in) :: xi,q,sd,cd,R,ytilde,dtilde
    real*8, intent(in) :: X11,X32,Y32
    real*8, intent(out) :: EY,EZ,FY,FZ,GY,GZ,HY,HZ
    ! local variables
    real*8 :: R3
    R3 = R**3.0d0

    EY=sd/R-ytilde*q/R3           ! E  in Tab.8, p.1032
    EZ=cd/R+dtilde*q/R3           ! E' in Tab.9, p.1033
    FY=dtilde/R3+xi*xi*Y32*sd     ! F  in Tab.8, p.1032
    FZ=ytilde/R3+xi*xi*Y32*cd     ! F' in Tab.9, p.1033
    GY=2.0d0*X11*sd-ytilde*q*X32  ! G  in Tab.8, p.1032
    GZ=2.0d0*X11*cd+dtilde*q*X32  ! G' in Tab.9, p.1033
    HY=dtilde*q*X32+xi*q*Y32*sd   ! H  in Tab.8, p.1032
    HZ=ytilde*q*X32+xi*q*Y32*cd   ! H' in Tab.9, p.1033
    return
end subroutine auxvar_derivative_yz
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine assign_alpha(alpha,alp)
!***** input
!*****   alpha : medium constant  (lambda+myu)/(lambda+2*myu)
!***** output
!*****   alp   : coefficients related to alpha value
    implicit none
    real*8, intent(in) :: alpha
    real*8, intent(out) :: alp(5)
    alp(1)=(1.0d0-alpha)/2.0d0
    alp(2)= alpha/2.0d0
    alp(3)=(1.0d0-alpha)/alpha
    alp(4)=1.0d0-alpha
    alp(5)=alpha
    return
end subroutine assign_alpha
!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine varout(u,ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz)
    implicit none
    real*8, intent(in) :: u(12)
    real*8, intent(out) :: ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    ux=u(1)
    uy=u(2)
    uz=u(3)
    uxx=u(4)
    uyx=u(5)
    uzx=u(6)
    uxy=u(7)
    uyy=u(8)
    uzy=u(9)
    uxz=u(10)
    uyz=u(11)
    uzz=u(12)
    return
end subroutine varout
!----------------------------------------------------------------------

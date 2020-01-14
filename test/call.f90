! This is a conversion file from FORTRAN77 to Fortran90
! The original source is available online:
!   http://www.bosai.go.jp/study/application/dc3d/download/DC3Dmanual.pdf
program test
    implicit none

    real(8) :: X, Y, Z, DEPTH, DIP, AL1, AL2, AW1, AW2, DISL1, DISL2, DISL3, ALPHA
    real(8) :: UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
    integer :: IRET

    real(8) :: POT1, POT2, POT3, POT4

    !----- Check DC3D
    X=10.0d0
    Y=20.0d0
    Z=-30.0d0
    !Z=30.0d0
    DEPTH=50.0d0
    DIP=70.0d0
    AL1=-80.0d0
    AL2=120.0d0
    AW1=-30.0d0
    AW2=25.0d0
    DISL1=200.0d0
    DISL2=-150.0d0
    DISL3=100.0d0
    !-----
    write(*,*) '*** OUTPUT OF DC3D ***'
    write(*,'(a,2(f8.2,","))') 'DEPTH, DIP =',DEPTH,DIP
    write(*,'(a,4(f8.2,","))') 'AL1, AL2, AW1, AW2 =',AL1,AL2,AW1,AW2
    write(*,'(a,4(f8.2,","))') 'DISL1, DISL2, DISL3 =',DISL1,DISL2,DISL3
    write(*,'(3(a,f8.2,","))') 'X =',X,' Y =',Y,' Z =',Z
    !-----
    ALPHA=2.0d0/3.0d0
    call DC3D(ALPHA,X,Y,Z,DEPTH,DIP, &
              AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3, &
              UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
    write(*,*)
    write(*,'(a,I1)') 'IRET = ',IRET
    write(*,'(a,3(f11.6,","))') 'UX, UY, UZ =',UX,UY,UZ
    write(*,*) '   ANSWER = -37.8981,    63.1789,    14.9607'
    stop
end program test

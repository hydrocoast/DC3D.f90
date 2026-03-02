! This is a conversion file from FORTRAN77 to Fortran90
! The originAL source is available online:
!   http://www.bosai.go.jp/study/application/dc3d/download/DC3DmanuAL.pdf
program test
    implicit none

    reAL(8) :: x, y, z, depth, dip, AL1, AL2, AW1, AW2, disl1, disl2, disl3, alpha
    reAL(8) :: ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    integer :: iret

    reAL(8) :: pot1, pot2, pot3, pot4

    !----- Check DC3D
    x=10.0d0
    y=20.0d0
    z=-30.0d0
    !z=30.0d0
    depth=50.0d0
    dip=70.0d0
    AL1=-80.0d0
    AL2=120.0d0
    AW1=-30.0d0
    AW2=25.0d0
    disl1=200.0d0
    disl2=-150.0d0
    disl3=100.0d0
    !-----
    write(*,*) '*** output of DC3D ***'
    write(*,'(a,2(f8.2,","))') 'depth, dip =',depth,dip
    write(*,'(a,4(f8.2,","))') 'AL1, AL2, AW1, AW2 =',AL1,AL2,AW1,AW2
    write(*,'(a,4(f8.2,","))') 'disl1, disl2, disl3 =',disl1,disl2,disl3
    write(*,'(3(a,f8.2,","))') 'x =',x,' y =',y,' z =',z
    !-----
    alpha=2.0d0/3.0d0
    call DC3D(alpha,x,y,z,depth,dip, &
              AL1,AL2,AW1,AW2,disl1,disl2,disl3, &
              ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)
    write(*,*)
    write(*,'(a,I1)') 'iret = ',iret
    write(*,'(a,3(f11.6,","))') 'ux, uy, uz =',ux,uy,uz
    write(*,*) '   answer = -37.8981,    63.1789,    14.9607'
    stop
end program test

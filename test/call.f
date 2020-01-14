C Appendix-2 : Test routine to check displacement calculation by DC3D
C   available at http://www.bosai.go.jp/study/application/dc3d/download/DC3Dmanual.pdf
C----- Check DC3D
      X=10.
      Y=20.
      Z=-30.
      DEPTH=50.
      DIP=70.
      AL1=-80.
      AL2=120.
      AW1=-30.
      AW2=25.
      DISL1=200.
      DISL2=-150.
      DISL3=100.
C-----
      WRITE(6,*) '*** OUTPUT OF DC3D ***'
      WRITE(6,*)
      WRITE(6,*) 'DEPTH,DIP=',DEPTH,DIP
      WRITE(6,*) 'AL1,AL2,AW1,AW2=',AL1,AL2,AW1,AW2
      WRITE(6,*) 'DISL1,DISL2,DISL3=',DISL1,DISL2,DISL3
      WRITE(6,*) 'X=',X,' Y=',Y,' Z=',Z
C-----
      ALPHA=2./3.
      CALL DC3D(ALPHA,X,Y,Z,DEPTH,DIP,
     * AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,
     * UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
      WRITE(6,*)
      WRITE(6,*) 'IRET=',IRET
      WRITE(6,*) 'UX,UY,UZ=',UX,UY,UZ
      WRITE(6,*) ' ANSWER = -37.8981 63.1789 14.9607'
      STOP
      END

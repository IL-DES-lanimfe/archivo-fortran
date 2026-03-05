      PROGRAM TF2D
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (N=4000)
      PARAMETER (NK=600)

      DIMENSION R(N),GR(N),RK(NK),FK(NK)
      DIMENSION XKR(N),XJ0(N),XIN(N)
      DIMENSION X(N),Y(N)

      OPEN (15, STATUS = 'OLD', FILE = 'interpo2.dat')
      OPEN (16, STATUS = 'UNKNOWN', FILE = 'ski2.dat')

C     ***********************************************
C     DATOS PARA LOS COEFICIENTES DE J_(0)
C     ***********************************************

      DATA p1b,p2b,p3b,p4b,p5b/1.d0,-.1098628627d-2,
     *.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1b,q2b,q3b,q4b,q5b/
     *-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,
     *-.934945152d-7/
      DATA r1b,r2b,r3b,r4b,r5b,r6b/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,
     *-184.9052456d0/,s1b,s2b,
     *s3b,s4b,s5b,s6b/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/

C     ************************************************
C     EL NUMERO DE PUNTOS EN R ES N Y EL NUMERO DE 
C     PUNTOS EN K ES NK
C     ************************************************

C     ************************************************
C     PARAMETROS GENERALES
C     ************************************************

      PI=4.D0*ATAN(1.D0)
      DENS=0.8836

C     ************************************************
C     LECTURA DEL ARCHIVO DE ENTRADA (R vs.GR)
C     ************************************************

      DO 5 LL=1,N,1
      READ(15,*)R(LL),GR(LL)
C     TRANSFORMANDO G(R) EN G(R)-1
      GR(LL)=GR(LL)-1.D0
5     CONTINUE

      RMAX=R(N)
      DR=R(2)-R(1)

C     ***********************************************
C     LA MALLA EN K SE CONSTRUYE PROBANDO KMAX Y 
C     DELTA K
C     ***********************************************

      KMAX=60
      DK=KMAX/(1.D0*NK)
      DO 1 JJ=1,NK,1
      RK(JJ)=JJ*DK
1     CONTINUE

C     **********************************************
C     EVALUACION DE J_(0)(rk)
C     **********************************************

      DO 2 J=1,NK,1
      DO 3 K=1,N,1
      XKR(K)=RK(J)*R(K)
      XKR2=XKR(K)
       if(abs(XKR2).lt.8.)then
         yy=XKR2**2
         XJ0(K)=(r1b+yy*(r2b+yy*(r3b+yy*(r4b+yy*(r5b+yy*r6b)))))/(
     *s1b+yy*(s2b+
     *yy*(s3b+yy*(s4b+yy*(s5b+yy*s6b)))))
       else
         ax=abs(XKR2)
         z=8./ax
         yy=z**2
         xx=ax-.785398164
         XJ0(K)=sqrt(.636619772/ax)*(cos(xx)*(p1b+yy*(p2b+yy*(
     *p3b+yy*(p4b+
     *yy*p5b))))-z*sin(xx)*(q1b+yy*(q2b+yy*(q3b+yy*(q4b+yy*q5b)))))
       endif

C     ***************************************************
C     EVALUACION DEL INTEGRANDO DE FK
C     ***************************************************

      XIN(K)=R(K)*GR(K)*XJ0(K)
3     CONTINUE

C     ***************************************************
C     INTEGRACION POR SIMPSON
C     ***************************************************

      CALL SIMPSON(XIN,DR,N,RES1)

      FK(J)=1.D0+2.D0*PI*DENS*RES1

      WRITE(16,*)SNGL(RK(J)),SNGL(FK(J))
2     CONTINUE

      STOP
      END

C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON
C     ****************************************************

      SUBROUTINE SIMPSON(XIN,DR,N,RES1)
c      IMPLICIT REAL*4(A-H,O-Z)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	PARAMETER (M=4000)
      DIMENSION XIN(N)
	
	RES1=0.0
	S0=0.0
	S1=0.0
	S2=0.0

	DO 100 I=2,M-1,2
	S0=S0+XIN(I-1)
	S1=S1+XIN(I)
	S2=S2+XIN(I+1)
100	CONTINUE
	RES1=DR*(S0+4.0*S1+S2)/3.0
	IF(MOD(M,2).EQ.0) RES1=RES1+DR*(5.0*XIN(M)+8.0*XIN(M-1)-
     #	XIN(M-2))/12.0
	RETURN
	END

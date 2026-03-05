C     PROGRAMA PARA CALCULAR LOS MOMENTOS DE F(K,T) EN DOS DIMENSIONES
c	CASO DE DISCOS BLANDOS (CASO JESUS) 
C	PARA EL CASO JESUS EN DICIEMBRE DE 2005 (1-ENERO-2005, BITACORA AZUL)
c	con correccion del tercer momento (2-I-06)
c	Enviado a Pedro el 31 de oct. de 2006.

      PROGRAM COLJESUS1
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N=1000)
      PARAMETER (NL=1000)

      DIMENSION R(N),FR(N),RK(NL),FK(NL),XM39K(NL),F11(N)
      DIMENSION XM1K(NL),XM21K(NL),XM22K(NL),CTKA(NL),F12(N)
      DIMENSION XM23K(NL),XM2K(NL),BETAK(NL),BETAK2(NL),F13(N)
      DIMENSION XM31K(NL),XM32K(NL),XM33K(NL),ALFAK(NL),F14(N)
      DIMENSION XM34K(NL),XM35K(NL),XM36K(NL),BETAK1(NL)
      DIMENSION XM37K(NL),XM3K(NL),F12K(NL),BETAK3(NL),XM38K(NL)
      DIMENSION XU10(N),XU11(N),XU12(N),XU13(N),CTK(NL),XM2KK(NL)
      DIMENSION XU20(N),XU21(N),XU22(N),XU23(N),XU24(N),CTK1(NL)
      DIMENSION XU25(N),XU26(N),XU27(N),XKR(N),XJ0(N),XJ1(N)
      DIMENSION F1(N),F2(N),F3(N),F4(N),F5(N),F1KS(NL),F2KS(NL)
      DIMENSION F6(N),F7(N),F8(N),F9(N),F10(N),XD1KS(NL),XD2KS(NL)
      DIMENSION XD0K(NL),XAA(NL),XBB(NL),XD1K(NL),XD2K(NL)
      DIMENSION F1K(NL),F2K(NL),FTK(NL),XD00K(NL),XAA1(NL)
      DIMENSION XBBB(NL),XAAA(NL)


      OPEN (15, STATUS = 'OLD', FILE = 'grsspy.dat')
      OPEN (16, STATUS = 'unknown', FILE = 'abcjpy.dat')
      OPEN (18, STATUS = 'OLD', FILE = 'skddpy.dat')

      PI=4.0*ATAN(1.0)

c	WRITE(*,*)'DAME IT (tiempo en segs)'
c	READ(*,*)IT
c	TT=IT/50.75

      RHO=0.5
      XA=6.0
      IZ=50
      DO 5 LL=1,N,1
      READ(15,*)R(LL),FR(LL)
5     CONTINUE

      DO 6 MM=1,NL,1
      READ(18,*)RK(MM),FK(MM)
6     CONTINUE

      DR=R(2)-R(1)
      DO 7 JJ=1,N,1
      XU10(JJ)=XA/R(JJ)**IZ
      XU11(JJ)=-IZ*XU10(JJ)/R(JJ)
c      XU12(JJ)=IZ*(IZ-1)*XU10(JJ)/R(JJ)**2
c      XU13(JJ)=-IZ*(IZ-1)*(IZ-2)*XU10(JJ)/R(JJ)**3
c	con correccion lapsus indicada por marco
      XU12(JJ)=IZ*(IZ+1)*XU10(JJ)/R(JJ)**2
      XU13(JJ)=-IZ*(IZ+1)*(IZ+2)*XU10(JJ)/R(JJ)**3
      XU20(JJ)=R(JJ)*XU12(JJ)+XU11(JJ)
      XU21(JJ)=R(JJ)*XU12(JJ)
      XU22(JJ)=R(JJ)*XU12(JJ)-XU11(JJ)
      XU23(JJ)=R(JJ)*(XU12(JJ))**2
      XU24(JJ)=(XU11(JJ)**2)/R(JJ)
      XU25(JJ)=R(JJ)*XU13(JJ)
      XU26(JJ)=(R(JJ)*XU13(JJ))-(3.0*XU12(JJ))+(3.0*XU11(JJ)/R(JJ))
      XU27(JJ)=((XU11(JJ)**2)/R(JJ))-(R(JJ)*XU12(JJ)**2)
 7    CONTINUE

      DO 8 KK=1,N,1
      F1(KK)=FR(KK)*XU20(KK)
      F4(KK)=F1(KK)
      F5(KK)=FR(KK)*XU23(KK)
      F6(KK)=FR(KK)*XU24(KK)

      F11(KK)=FR(KK)**2*XU23(KK)
      F12(KK)=FR(KK)**2*XU24(KK)

 8    CONTINUE
      CALL SIMPSON(F1,DR,N,XINT1)
      CALL SIMPSON(F4,DR,N,XINT4)
      CALL SIMPSON(F5,DR,N,XINT5)
      CALL SIMPSON(F6,DR,N,XINT6)

      CALL SIMPSON(F11,DR,N,XINT11)
      CALL SIMPSON(F12,DR,N,XINT12)



      DO 10 II=1,NL,1

      XM1K(II)=-RK(II)**2

      XM21K(II)=RHO*PI*RK(II)**2*XINT1

      XM31K(II)=(3.0*PI*RHO*RK(II)**4)*XINT4
      XM32K(II)=(2.0*PI*RHO*RK(II)**2)*XINT5
      XM33K(II)=(2.0*PI*RHO*RK(II)**2)*XINT6 

      DO 9 I=1,N,1
      XKR(I)=RK(II)*R(I)
      XKR2=XKR(I)
      CALL XJOTA0(XKR2,XXJ1)
      XJ0(I)=XXJ1
      CALL XJOTA1(XKR2,XXJ2)
      XJ1(I)=XXJ2
      F2(I)=FR(I)*XU21(I)*XXJ1
      F3(I)=FR(I)*XU22(I)*(XXJ2/XKR2)
      F7(I)=FR(I)*XU25(I)*XXJ2
      F8(I)=FR(I)*XU26(I)*((XXJ1/XKR2)-(2.0*XXJ2/(XKR2**2)))
      F9(I)=FR(I)*XU23(I)*XXJ1
      F10(I)=FR(I)*XU27(I)*(XXJ2/XKR2)

      F13(I)=FR(I)**2*XU23(I)*(XXJ1-XXJ2/XKR2)
      F14(I)=FR(I)**2*XU24(I)*(XXJ2/XKR2)


 9    CONTINUE

      CALL SIMPSON(F2,DR,N,XINT2)
      CALL SIMPSON(F3,DR,N,XINT3)
      CALL SIMPSON(F7,DR,N,XINT7)
      CALL SIMPSON(F8,DR,N,XINT8)
      CALL SIMPSON(F9,DR,N,XINT9)
      CALL SIMPSON(F10,DR,N,XINT10)

      CALL SIMPSON(F13,DR,N,XINT13)
      CALL SIMPSON(F14,DR,N,XINT14)

      XINT15=XINT11+XINT12-2.0*XINT13-2.0*XINT14
      XM22K(II)=2.0*PI*RHO*RK(II)**2*XINT2
      XM23K(II)=2.0*PI*RHO*RK(II)**2*XINT3

      XM2K(II)=RK(II)**4+XM21K(II)-XM22K(II)
     #+XM23K(II)

      XM2KK(II)=RK(II)**8+2.0*(XM21K(II)-XM22K(II)
     #+XM23K(II))*RK(II)**4
      XM34K(II)=4.0*RHO*PI*(RK(II)**3)*XINT7
      XM35K(II)=4.0*RHO*PI*(RK(II)**3)*XINT8
      XM36K(II)=4.0*PI*RHO*(RK(II)**2)*XINT9
      XM37K(II)=4.0*PI*RHO*(RK(II)**2)*XINT10
      XM38K(II)=(XM2K(II)-RK(II)**4)**2/RK(II)**2
      XM39K(II)=2.0*PI*RHO**2*RK(II)**4*XINT15

      XM3K(II)=-RK(II)**6-XM31K(II)-XM32K(II)
     #-XM33K(II)-XM34K(II)-XM35K(II)+XM36K(II)+XM37K(II)
     #-XM38K(II)
c     #+XM39K(II)
c      CTK(II)=XM2KK(II)+XM3K(II)*RK(II)**2
      CTK1(II)=XM2K(II)**2+XM3K(II)*RK(II)**2
      XAA1(II)=(XM1K(II)**2)/FK(II)
 
      XAA(II)=XM2K(II)-XAA1(II)
      CTKA(II)=CTK1(II)/XAA(II)
      CTK(II)=4.0*CTKA(II)/FK(II)
      XBB(II)=(-XM3K(II)-(XM1K(II)**3)/(FK(II)**2)
     #+2.0*XM1K(II)*XM2K(II)/(FK(II)))/XAA(II)

C     COEFICIENTES DE LA LAURA

      XAAA(II)=XAA(II)/RK(II)**2
      XBBB(II)=XBB(II)-XAAA(II)


      WRITE(16,*)SNGL(RK(II)),SNGL(XAAA(II)),SNGL(XBBB(II))

      BETAK1(II)=4.0*(XAA(II)-XBB(II)*RK(II)**2)/FK(II)
      BETAK2(II)=XBB(II)+RK(II)**2/FK(II)
      BETAK3(II)=CTK(II)/BETAK2(II)**2
      BETAK(II)=DSQRT(1.0+BETAK3(II))
      XD00K(II)=(XBB(II)-(RK(II)**2)/FK(II))**2+
     #4.0*XAA(II)/FK(II)
      XD011=XD00K(II)
      XD0K(II)=DSQRT(XD011)
      XD1K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
     #+XD0K(II)/2.0
      XD2K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
     #-XD0K(II)/2.0

      ALFAK(II)=((RK(II)**2/FK(II))+XBB(II))/2.0

      F1K(II)=FK(II)*(XD1K(II)-XBB(II))/(XD1K(II)-XD2K(II))
      F2K(II)=FK(II)*(XD2K(II)-XBB(II))/(XD2K(II)-XD1K(II))
c      FTK(II)=F1K(II)*EXP(-XD1K(II)*TT)+F2K(II)
c     #*EXP(-XD2K(II)*TT)
c      WRITE(16,*)SNGL(RK(II)),SNGL(FTK(II))
10    CONTINUE

      STOP
      END
C     ***********************************************
C     SUBRUTINA PARA EL CALCULO DE LA J_(0). 
C     BASE: ABRAMOWITZ, p.369.
C     ***********************************************
      SUBROUTINE XJOTA0(XKR2,XXJ0)
      IMPLICIT REAL*8(A-H,O-Z)
       if(abs(XKR2).le.3.0)then
         yy=(XKR2/3.0)**2

         XXJ0=1.0-2.2499997*yy+1.2656208*yy**2-
     *0.3163866*yy**3+0.0444479*yy**4-0.0039444*
     *yy**5+0.0002100*yy**6    

       else
         yy=(3.0/XKR2)

         TETA0=XKR2-0.78539816-0.04166397*yy-0.00003954*
     *yy**2+0.00262573*yy**3-0.00054125*yy**4-
     *0.00029333*yy**5+0.00013558*yy**6

         XF0=0.79788456-0.00000077*yy-0.00552740*yy**2-
     *0.00009512*yy**3+0.00137237*yy**4-0.00072805*yy**5 +
     *0.00014476*yy**6

         XXJ0=XF0*COS(TETA0)/sqrt(XKR2)

       endif

       RETURN
       END

C     ***********************************************
C     SUBRUTINA PARA EL CALCULO DE LA J_(1). 
C     BASE: ABRAMOWITZ p.370.
C     ***********************************************
      SUBROUTINE XJOTA1(SPL,XXJ1)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(abs(SPL).le.3.0)then
        yy=(SPL/3.0)**2

        XXJ1=(0.5-0.56249985*yy+0.21093573*yy**2-
     *0.03954289*yy**3+0.00443319*yy**4-0.00031761*
     *yy**5+0.00001109*yy**6)*SPL

      else
        yy=(3.0/SPL)

         TETA1=SPL-2.35619449+0.12499612*yy+0.00005650*
     *yy**2-0.00637879*yy**3+0.00074348*yy**4+0.00079824*
     *yy**5-0.00029166*yy**6    

         XF1=0.79788456+0.00000156*yy+0.01659667*yy**2+
     *0.00017105*yy**3-0.00249511*yy**4+0.00113653*yy**5-
     *0.00020033*yy**6

        XXJ1=(XF1*COS(TETA1))/sqrt(SPL)

      endif
      RETURN
      END

C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON 
C     ****************************************************
      SUBROUTINE SIMPSON(XINI,DR,N,RES1)
      IMPLICIT REAL*8(A-H,O-Z)
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XINI(N)
	
	RES1=0.0
	S0=0.0
	S1=0.0
	S2=0.0

	DO 100 I=2,N-1,2
	S0=S0+XINI(I-1)
	S1=S1+XINI(I)
	S2=S2+XINI(I+1)
100	CONTINUE
	RES1=DR*(S0+4.0*S1+S2)/3.0
	IF(MOD(N,2).EQ.0) RES1=RES1+DR*(5.0*XINI(N)+8.0*XINI(N-1)-
     #	XINI(N-2))/12.0
	RETURN
	END












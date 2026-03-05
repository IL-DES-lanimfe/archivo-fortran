	Program dia
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Programa para calcular diametros de esfera dura equivalente			C
C																	C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 implicit none !toda variable se tiene que decalrar
	 integer i
	 double precision n,d,diametro
	 open(20,file='diametros.dat')
	 DO i=1,100
	  n=dble(i)
	  d=diametro(n)
	  write(20,*)n,d
	 ENDDO
	 close(20)
	 stop
	END

      Double precision function diametro(nu)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Este programa sirve para calcular el diametro de esfera dura equi-  C
C valente mediante la solución numérica de la ecuación 6.3.11 del  C
C libro "Theory of simple liquids", de Hansen y Macdonald.            C
C																	C
C Enero 2007															C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 implicit none !toda variable se tiene que declarar
CC
C Variables para hacer la integral
CC
	 Integer Nr
	 Parameter (Nr=2**11)
	 Double precision Iu(Nr),d
CC
C variables para el potencial en particular que estemos trabajando
CC
	 Double precision nu
	 Double precision R,pot
CC
C variables auxiliares
CC
	 Integer i
	 Double precision dr

	 dr=1.D-3
	 Do i=1,Nr
	  R=dble(i)*dr
	  pot=0.D0
CC
C Definicion del potencial
CC
	  If (R.le.1.D0)then
	   pot= (1.D0/R**nu)**2-2.D0*(1.D0/R**nu)+1.D0
	  endif
CC
C integrando para calcular d/sigma
CC
	  Iu(i)=1.D0-dexp(-pot)
!	 write(50,*)R,Iu(i)
	 enddo
! 	stop
CC
C integracion por simpson
CC
	 call SIMPSON(Iu,dr,Nr,d)
CC
!	write(*,*)'el diametro de esfera dura equivalente es:',d
	 diametro=d
	 Return
      END

	SUBROUTINE SIMPSON(XIN,DR,M,RES1)
C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON
C     ****************************************************

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XIN(M)
	
	RES1=0.D0
	S0=0.D0
	S1=0.D0
	S2=0.D0
   
	DO 100 I=2,M-1,2
	S0=S0+XIN(I-1)
	S1=S1+XIN(I)
	S2=S2+XIN(I+1)
100	CONTINUE
	RES1=DR*(S0+4.D0*S1+S2)/3.D0
        
	IF(MOD(M,2).EQ.0) RES1=RES1+DR*(5.D0*XIN(M)+8.D0*XIN(M-1)-
     #  XIN(M-2))/12.D0
	RETURN
	END
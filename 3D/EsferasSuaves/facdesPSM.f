      PROGRAM SDKPSM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Este programa calcula el factor de estructura de una sistema con un 	C
C potencial de interaccion tipo esferas penetrables (PSM por sus siglas C
C en ingles). Solucion analitica usando la aproximacion  de campo medio C
C (PRE 63, 031206 (2001) LIKOS ET AL).					C
C									C
C Julio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
	integer i
	double precision dy,y,Sk,fv,T
	double precision SdekPSM
	dy=0.1D-1
	fv=9.D0
	T=6.D0
	open(20,file='facdesPSMfv9T6.dat')
	do i=1,2**12
	 y=dble(i)*dy
	 Sk=SdekPSM(y,fv,T)
	 write(20,*)y,Sk
	enddo
	close(20)
	stop
      END
CC
C funcion de S(k) para PSM
CC
      DOUBLE PRECISION FUNCTION SdekPSM(k,eta,tem)
	implicit none
	Double precision k,eta,tem
	Double precision Senos
	Senos=(dsin(k)-k*dcos(k))/k**3
	SdekPSM=1.D0/(1.D0+24.D0*eta*Senos/tem)
	return
      END
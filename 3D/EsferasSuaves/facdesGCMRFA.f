      PROGRAM GaussianCoreRPA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para calcular el factor de estructura de un sistema con pot. C
C de interaccion tipo gaussiano. En este caso se usa la aproximacion	C
C de campo medio.							C
C									C
C Agosto 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
	Integer i
	Double precision dy,rho,eps,y,S
	Double precision SdkGCRPA
	dy=1.D-2
	rho=0.25D0
	eps=1.D0/0.008D0
	write(*,*)eps
	open(20,file='facdesGCRPArho0,25T0,008.dat')
	DO i=1,5000
	 y=dble(i)*dy
	 S=SdkGCRPA(y,rho,eps)
	 write(20,*)y,S
	ENDDO
	close(20)
	stop
      END
CC
C Funcion para calcular factor de estructura GCM
CC
      Double precision function SdkGCRPA(k,den,Tin)
	implicit none
	double precision k,den,Tin
	double precision alpha,pi
	pi=dacos(-1.D0)
	alpha=den*Tin*pi**(3.D0/2.D0)
	SdkGCRPA=1.D0/(1.D0+alpha*dexp(-k**2/4.D0))
	Return
      End
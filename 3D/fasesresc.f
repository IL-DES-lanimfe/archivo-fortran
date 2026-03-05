      PROGRAM fasesrescal
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Este programa es para cambiar el diagrama de fases escalando con el   C
C valor de esfera dura.							C
C									C
C Junio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
	Integer i
	Double precision fv,fvHS,nu,xx
	fvHS=0.552402675D0
	open(20,file='vidrioLJMRYDS.dat')
	open(30,file='vidrioHSLJMRYDS.dat')
	do i=1,8
	 read(20,*)fv,xx,nu
	 write(30,*)fv/fvHS,nu
	enddo
	close(20)
	close(30)
      END
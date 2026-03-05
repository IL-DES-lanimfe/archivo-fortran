      PROGRAM INVDIA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para invertir y elevar al cubo cociente d/sigma		C
C									C
C Junio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
	Integer i
	Double precision nu,d
	open(20,file='diametrosblip.dat')
	open(30,file='diametrosblipinv.dat')
	DO i=1,100
	 read(20,*)nu,d
	 write(30,*)1.D0/d,(1.D0/d)**3,0.5627D0*(1.D0/d)**3,nu
	ENDDO
	close(20)
	close(19)
	stop
      END
      PROGRAM DESVIACIONES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para claular el cociente n(nu)/n(HS) en la curva de vitrifi_ C
C cación para esferas suaves.						C
C									C
C Junio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none
	Integer i,j
	Double precision fv,nu,d,Dia(99,2)
	open(20,file='vidrioHSLJMRYDS.dat')
	open(30,file='diametroISO.dat')
	open(40,file='isoestructuraRYS.dat')
	Do i=1,99
	 read(30,*)Dia(i,1),Dia(i,2)
!	 write(60,*)Dia(i,1),Dia(i,2)
	Enddo
	close(30)
	Do i=1,8
	 read(20,*)fv,nu
	 do j=1,99
	  if(nu.eq.Dia(j,1))d=Dia(j,2)
	 enddo
	 write(40,*)fv*d**3,nu
	Enddo
	close(20)
	close(40)
	stop
      END
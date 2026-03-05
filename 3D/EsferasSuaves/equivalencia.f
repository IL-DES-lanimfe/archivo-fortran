	PROGRAM equiv
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Este progama sirve para convertir las fracciones de volumen de esf. C
C suave a la de esfera dura equivalente.								C
C																	C
C Febrero 2008														C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 implicit none !toda variable se tiene que declarar
	 integer i
	 double precision nu,nu2,xx,d,fv
	 open(20,file='diametrosD.dat')
	 open(30,file='vidrioLJMRY.dat')
	 open(40,file='vidrioLJMRYr.dat')
	 do i=1,100
	  read(20,*)nu,d
	  if((nu.ge.10).and.(mod(i,5).eq.0))then
	   write(*,*)nu
	   read(30,*)fv
	   write(40,*)fv*d**3,nu
	  endif
	 enddo
	 stop
	END


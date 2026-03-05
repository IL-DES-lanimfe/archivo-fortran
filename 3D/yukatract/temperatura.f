      PROGRAM convierte
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa es para cambiar de ejes de K a 1/K, la cual es propor-C
C cional a la temperatura.					     C
C								     C
C Enero 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se debe de declarar
	integer i
	double precision x,xx
	open(50,file='cristalyukatr.dat')
	open(51,file='cristalyukatrT.dat')
	do i=1,103
	 read(50,*)x,xx
	 write(51,*)x,1.D0/xx,xx
	enddo
	stop
      END
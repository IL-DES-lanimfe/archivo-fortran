      PROGRAM TENTO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa sirve para transformar el eje temporal de los datos  C
C experimentales reportados en PRE 49,4206, para la funcion de dis-  C
C persion intermedia de un sistema de esfera dura		     C
C								     C
C Marzo 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none !toda variable se tiene que declarar
CC
	Integer i,Np
	Parameter(Np=61)
	Double precision T,Fkt
	open(50,file='fktq4,10fv587.dat')
	open(60,file='fktent0q4,10fv587.dat')
	do i=1,Np
	 read(50,*)T,Fkt
	 write(60,*)10.D0**T*3.5D-6,Fkt
	enddo
	stop
      END
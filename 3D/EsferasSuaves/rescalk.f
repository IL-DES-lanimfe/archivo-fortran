      PROGRAM cambiodeunidades
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para cambiar las unides en la k de los factores no ergodicos C
C de esfera suave.							C
C									C
C Julio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none
	Integer i
	Double Precision d,y,f
	d=0.858444771728931D0
	open(20,file='facnoergoLJMRYnu5.dat')
	open(30,file='facnoergHSLJMRYnu5.dat')
	do i=1,2**12
	 read(20,*)y,f
	 write(30,*)y*d,f
	enddo
	close(20)
	close(30)
	stop
      END
     SUBROUTINE maxmin(Sk,Np,imax,imin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Esta subrutina calcula el valor maximo de una lista,esta diseñado  C
C para usarlo con el factor de estructura y por eso se usa un arregloC
C de dos columnas pero se puede modificar facilmente para alguna otraC
C lista								     C
C Tambien calcula el minimo que este despues del maximo		     C
C Noviembre 2006						     C
C parametros:							     C
C		*Sk:arreglo de 2 columnas			     C
C		*Np:numero de renglones del arreglo		     C
C		*imax: es el numero del renglon donde esta el maximo C
C		*imin: es el numero del renglon donde esta el minimo C
C								     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none !toda variable se tiene que declarar
CC
C variables para los parametros
CC
	Integer Nk,imax,imin
	Double precision Sk(Nk,2)
CC
C variables internas
CC
	integer i
	Double precision sminima,smaxima
CC
	smaxima=Sk(1,2)
	Do i=2,Nk
	 if(Sk(i,2).gt.smaxima) imax=i
	enddo
	sminima=S(imax,2)
	Do i=imax+1,Nk
	 if(Sk(i,2).lt.sminima) imin=i
	enddo
	Return
     END
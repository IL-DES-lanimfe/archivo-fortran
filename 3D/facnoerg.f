      Subroutine  facnoerg(gamma,Sk,gcero,kmin,nk,f)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Esta subrutiuna evalua los factores no ergodicos predichos por el  C
C criterio de atrapamiento.					     C
C Parametros: *gamma: valor del parametro gamma del criterio	     C
C	      *Sk:factor de estructura				     C
C	      *gcero:funcion relacionada con Vineyard		     C
C	      *kmin:posicion del primer minimo de S(k)		     C
C	      *nk:numero de puntos en k				     C
C	      *f:arreglo de salida				     C
C								     C
C Agosto de 2006						     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Implicit none !toda variable se tiene que declarar
CC
C Parametros
CC
	Integer nk
	Double precision kmin,gamma
	Double precision Sk(nk,2),gcero(nk),f(nk,2)
CC
C variables internas
CC
	Integer i
	Double precision l,den,num
	Double precision lambda
CC
C evaluando el factor no ergodico
CC
	DO i=1,nk
	 l=lambda(Sk(i,1),kmin)
	 den=Sk(i,2)*gcero(i)*l
	 num=-gamma*(Sk(i,1))**2
	 f(i,1)=Sk(i,1)
	 f(i,2)=1.D0/(1.D0+ (num/den))
	ENDDO
	Return
      End

     Double precision function lambda(x,xminima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								   C
C Funcion interpoladora lambda(k),se usa en la cerradura para      C 
C Cs(k,z) en la Teoria Autoconsistente. Ver detalles en la tesis   C
C de Laura Yeomans Reyna.                                          C
C Parametros: *x: valor de k en que queremos evaluar la funcion    C
C	      *xminima: valor de k en donde esta el primer minimo  C
C			de S(k)					   C
C								   C
C Agosto de 2006						   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       Implicit none

       Double precision x,xminima

       lambda= 1.D0/(1.D0 + (x/xminima)**2)
       Return

      End
      SUBROUTINE interpola(nk,Ski,gri,kk,Sk,gr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta subrutina es para interpolar el factor de estructura,se nece-C
C sita para hacer coincidir el arreglo de k que se lee de la funcionC
C hidrodinamica con el arreglo en k de S(k)			    C
C Parametros:							    C
C		*nk: numero de puntos en el arreglo		    C
C		*Ski:arreglo bidimensional original		    C
C		*gri:arreglo de g(r) original			    C
C		*kk:arreglo de k de la funcion hidrodinamica	    C
C		*Sk:arreglo bidimensional ajustado		    C
C		*gr:arreglo de g(r) ajustado			    C
C								    C
C Noviembre 2006					    	    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer nk,i
	Double precision Ski(nk,2),gri(nk,2),kk(nk)
	Double precision Sk(nk,2),gr(nk,2),dr,pi
	pi=dacos(-1.D0)
CC
C inicia la interpolacion de S(k)
CC
	Do i=1,nk
	 Sk(i,1)=kk(i)
	 if(Ski(i,1).EQ.kk(i)) Sk(i,2)=Ski(i,2)
	 if(Ski(i,1).GT.kk(i)) Sk(i,2)=(Ski(i-1,2)+Ski(i,2))/2.D0
	 if(Ski(i,1).LT.kk(i)) Sk(i,2)=(Ski(i+1,2)+Ski(i,2))/2.D0
	enddo
CC
C termina interpolacion de S(k)
C inicia intepolacion de g(r)
CC
	dr=pi/(dble(nk)*kk(nk))
	Do i=1,nk
	 gr(i,1)=dble(i)*dr
	 if(gri(i,1).EQ.gr(i,1)) gr(i,2)=gri(i,2)
	 if(gri(i,1).GT.gr(i,1)) gr(i,2)=(gri(i-1,2)+gri(i,2))/2.D0
	 if(gri(i,1).LT.gr(i,1)) gr(i,2)=(gri(i+1,2)+gri(i,2))/2.D0
	enddo
CC
C termina interpolacion de g(r)
CC
	return
      end
	  
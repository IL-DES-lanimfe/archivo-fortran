      PROGRAM facdesSW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa sirve para calcular el factor de estructura S(k) paraC
C un potencial tipo pozo cuadrado basandose en la solucion analitica C
C para PY usada en la siguiente referencia: PRE, 66, 021403 (2002).  C
C								     C
C Julio 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer i,Nk
	Parameter (Nk=2**12)
	Double Precision Sdek(Nk,2),k,dk,fv,eps,Tinv,SkSW
	k=0.01D0
	dk=0.01D0
	fv=0.5D0
	eps=0.03D0
	Tinv=0.001D0
	do i=1,Nk
	 Sdek(i,1)=k
	 Sdek(i,2)=SkSW(k,fv,eps,Tinv)
	 write(50,*)Sdek(i,1),Sdek(i,2)
	 k=k+dk
	enddo
      END
CC
      Double precision function SkSW(y,fv,rng,prf)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Esta funcion calcula el factor de estructura para un potencial de  C
C Pozo cuadrado. Basado en la solucion analitica que se usa en la sigC
C referencia: PRE, 66, 021403 (2002).				     C
C Potencial:	V(r)= infinito para r entre 0 y R'		     C
C		V(r)= -u para r entre R' y R			     C
C		V(r)= 0 para r mayor que R			     C
C Parametros:	* y: vector de onda				     C
C		* fv: fraccion de volumen			     C
C		* rng: (R-R')/R					     C
C		* prf: u/KT					     C
C								     C
C Julio 2007							     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
CC
C parametros
CC
	Double precision y,fv,rng,prf
CC
C variables internas
CC
	Double precision Delta,Gama,lambda,mu,alpha,beta,Saux
CC
C Funciones que se usan
CC
	Double precision f1,f2,f3,f5
CC
C Cantidades que dependen de los parametros del potencial y la fracc.
C de volumen
CC
	Delta=dexp(-prf)/(12.D0*rng)+fv/(1.D0-fv)
	Gama=fv*(1.D0+fv/2.D0)/(3.D0*(1.D0-fv)**2)
	lambda=6.D0*(Delta-dsqrt(Delta**2-Gama))/fv
	mu=lambda*fv*(1.D0-fv)
	alpha=(1.D0+2.D0*fv-mu)**2/(1.D0-fv)**4
	beta=-(3.D0*fv*(2.D0+fv)**2-2.D0*mu*(1.D0+7.D0*fv+fv**2)+mu**2
     #         *(2.D0+fv))/(2.D0*(1.D0-fv)**4)
CC
C la siguiente cantidad es la que viene en el articulo y es 1/S(Q)-1
CC
	Saux=24.D0*fv*(alpha*f2(y)+beta*f3(y)+fv*alpha*f5(y)/2.D0)
     #      +4.D0*(fv*lambda*rng)**2*(f2(rng*y)-f3(rng*y)/2.D0)
     #      +2.D0*(fv*lambda)**2*(f1(y)-rng**2*f1(rng*y))
     #      -2.D0*fv*lambda*(f1(y)-(1.D0-rng)**2*f1((1.D0-rng)*y))/rng
     #      -24.D0*fv*(f2(y)-(1.D0-rng)**3*f2((1.D0-rng)*y))
CC
C paso intermedio para despejar
CC
	Saux=Saux+1.D0
CC
C Valor final para S(k)
CC
	SkSW=1.D0/Saux
	return
      end

      Double precision function f1(x)
	implicit none
	double precision x
	f1=(1.D0-dcos(x))/x**2
	return
      end

      Double precision function f2(x)
	implicit none
	double precision x
	f2=(dsin(x)-x*dcos(x))/x**3
	return
      end

      Double precision function f3(x)
	implicit none
	double precision x
	f3=(2.D0*x*dsin(x)-(x**2-2.D0)*dcos(x)-2.D0)/x**4
	return
      end

      Double precision function f5(x)
	implicit none
	double precision x
	f5=((4.D0*x**3-24.D0*x)*dsin(x)-(x**4-12.D0*x**2+24.D0)*dcos(x)
     #      +24.D0)/x**6
	return
      end
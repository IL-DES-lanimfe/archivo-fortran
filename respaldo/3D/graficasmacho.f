      Program graficas
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C 
C Este programa evalua la integral de la funcional del parametro    C 
C gamma y la integral del criterio de atrapamiento.Esta diseñado    C
C para dos parametros de control (por ejemplo frac. de vol. y tempe-C
C ratura), pero se puede expandir facilmente para mas parametros.   C
C Para cambiar de sistema solo es necesario cambiar subrutina y     C
C variables para el factor de estrucutura.                          C
C								    C
C Agosto de 2006                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none
C indices para los ciclos en general
        integer iint,ig0 
C variables para la integral sobre vectores de onda
        double precision  y,dy,S  
C variables para evaluar la funcional de gamma
        double precision gamma,dgamma,gammamax,gammamin 
C numero de puntos a evaluar en gamma
	integer Ngam
C variables para las integrales de funcional (temporales)
        double precision num,den1,den2 
C variables las evaluaciones de la funcional de gamma
        double precision integral
C variables para parametros del criterio (labda,gcero,1er min S(k))
        double precision l,g,kmin
C Numero de puntos en vecs. de onda
	Integer nk
	parameter (nk=2**12)
C arreglos para S(k) y g(r)
        double precision  Sdek(nk,2),T(nk),g0(nk)
C parametros S(k) y g(r)
	Double precision parametro1,parametro2 
C variables para evaluar en los parametros de control
	double precision par2max,par2min,dpar2
C variables auxiliares
        Double precision pi,sga
        double precision criterio,gamreal
CC
CCCCCCCCCCCCCCCC Funciones que utiliza el programa CCCCCCCCCCCCCCCCCCC
CC
        double precision  lambda,kminima,sgamma
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
CC Esta parte siempre acompaña a la subrutiuna de Sdek diseñada por C
CCCCCCCCCCCCCCCCCCCC   ?????????????   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCC Variables,commons y decalraciones necesarias para S(k) CCCCCCC
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Archivos de salida CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C archivo de salida,gamma,funcional de gamma
        open(8,file='fungama.dat')
C archivo de salida funcional criterio 
        open(9,file='funcrit.dat')
CC
        pi= dacos(-1.D0)
CC
C saludo de beinvenida y entrada de datos
CC
	write(*,*)' '
	write(*,*)'Bienvenido a graficas, vamos a calcular'
	write(*,*)'las integrales del criterio y de la funcional de'
	write(*,*)'gamma para asi tener una idea de donde esta la'
	write(*,*)'transición vítrea. Necesito unos parametros'
	write(*,*)' '
	write(*,*)'Es el caso de (anote el sistema)'
	write(*,*)' '
CC
C El programa te hace las graficas de la funcional del criterio como 
C funcion del parametro 2 manteniendo el parametro 1 fijo, es decir 
C cada corrida del programa te genera una iso-parametro1 por decirlo 
C de alguna manera.
CC
        write(*,*)'¿cuanto vale (parametro1)?'
        read (*,*)parametro1
        write(*,*)'¿cuanto vale (parametro2) max.?'
        read (*,*)par2max
	write(*,*)'¿cuanto vale (parametro 2) min?'
        read (*,*)par2min
	write(*,*)'¿Longitud del incremento en (parametro 2)?'
        read (*,*)dpar2
	write(*,*)'¿cuantos puntos quieres en la fun. de gamma?'
        read (*,*)Ngam
CC
C Empieza el ciclo para barrer sobre parametro 2
CC
	DO parametro2=par2min,par2max,dpar2
          write(*,*)'parametros',parametro1,parametro2
CC
CCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculado S(k)'
CC Aqui se inserta la subrutina del factor de estrucutura
          call 
CC
C se concoce S(k) y g(r)
CC
	  dy=Sdek(2,1)-Sdek(1,1)
          write(*,*)'calculando la k del primer minimo de S(k)'
          kmin=kminima(SdeK,nk)
CC
C ahora se conoce la posicion del primer minimo de S(k)
CC
  	  write(*,*)'evaluando la funcion g-cero'
	  do ig0=1,nk
	   g0(ig0)=1.D0
CC
C la definicion de g-cero depende de la vineyard, para variar gcero se
C puede hacer una funcion
CC
C          y=SdeK(ig0)
C          g0(ig0)=gcero(y,parametro1,parametro2)
	  enddo
CC
CCCCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculando gamma del criterio'
          gamreal=sgamma(SdeK,nk,fv,g0,kmin)
          if (gamreal.eq.1.D20)then
           write(*,*)'no se encontro solucion para gamma',parametro2
           goto 200
          endif
	  write(*,*)'gamma=',gamreal
CC	  
C  ya se conoce gamma
CC
	   write(*,*)'calculando la funcional del criterio'
           criterio= 0.D0
CC
C evaluando el integrando del criterio
CC
           do iint=1, nk
	    y=Sdek(iint,1)
            S=Sdek(iint,2)
            g= g0(iint)
            l= lambda(y,kmin)
            sga=gamreal
          T(iint)=((S-1)*l*y)**2/(36.D0*pi*fv*(l*(S+1/g)-2*y**2*sga/g))
          enddo
CC
C evaluando la integral del criterio
CC
	   call SIMPSON(T,dy,Nk,criterio)
           write(9,*)parametro2,criterio
           write(*,*)'criterio listo',parametro2,criterio
              
	   write(*,*)'evaluando funcional de gamma'
CC
C Aqui se genera una grafica para la funcional de gamma como funcion de
C gamma (esta deberia de ser una fucion de 3 variables gamma,parametro1
C , parametro2). En cada corrida del programa se generan tantas curvas
C como puntos en el intervalo de evaluacion del parametro2
CC 
	    gammamin=-10.D0
	    gammamax=0.D0
	    dgamma=dabs(gammamax-gammamin)/dble(Ngam)
 200        continue
            Do gamma=gammamin,gammamax,dgamma
              integral= 0.D0
CC
C evaluacion del integrando de la funcional de gamma
CC
              do iint=1, Nk
		 y=Sdek(iint,1)
                 S=SdeK(iint,2)
                 g= g0(iint)
                 l= lambda(y,kmin)
                 num= (-1.D0)*((S-1.D0)*l)**2
                 den1= 36.D0*pi*fv*(l*(S+1.D0/g)-(2.D0*gamma/g)*y**2)
                 den2= 1.D0/(l-gamma*y**2) + 1.D0/(l*g*S - gamma*y**2)
                 T(iint)= num/(den1*den2)
              enddo
CC
C evaluacion de la integral
CC
              call SIMPSON(T,dy,Nk,integral)
	      write(8,*)gamma,integral
CC
CC Termina la evaluacion de la  integral de la funcional de gamma
CC
C aumentamos gamma para el siguiente punto
CC
          EndDo
          write(*,*)'funcional de gamma lista',parametro2
CC Termina el ciclo de evaluacion de la funcional
       ENDDO
CC Termina ciclo de parametros
       close(8)
       close(9)
CCC Fin del programa CCCC
       stop
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutinas del factor de estructura, proporcionadas por ???????  C
C 	 							   C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                    C
C Las siguentes funciones se utilizan en el criterio C
C                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Double precision function kminima(Sk,N)

       Implicit none
       
       Integer N
C       Double precision emp

       Integer h,m,j
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax,kmin,kmax

       Double precision facdes

C       y=3.D0
C       dy=0.001D0
       
C       Do i=1,N

C          S(i,1)=y
C          S(i,2)= facdes()
C          y= y+ dy

C       enddo

       smax=Sk(1,2)
               
       Do j=2, N
          if (Sk(j,2).GT.smax)then
            
             smax=Sk(j,2)
             m= j
          endif

       enddo
       
       smin= Sk(m,2)

       Do h=m , N-m

          if (Sk(h,2).LT.smin)then

             smin= Sk(h,2)
             kmin= Sk(h,1)

          endif

       enddo

       kminima=kmin
       Return

       End

       SUBROUTINE SIMPSON(XIN,DR,M,RES1)
C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON
C     ****************************************************

	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION XIN(M)
	
	RES1=0.D0
	S0=0.D0
	S1=0.D0
	S2=0.D0
   
	DO 100 I=2,M-1,2
	S0=S0+XIN(I-1)
	S1=S1+XIN(I)
	S2=S2+XIN(I+1)
100	CONTINUE
	RES1=DR*(S0+4.D0*S1+S2)/3.D0
        
	IF(MOD(M,2).EQ.0) RES1=RES1+DR*(5.D0*XIN(M)+8.D0*XIN(M-1)-
     #  XIN(M-2))/12.D0
	RETURN
	END

      Double precision function sgamma(Sk,maxlist,emp,gcero,kmin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Esta funcion sirve para resolver la ecuacion para el parametro     C
C gamma que aparece en el criterio de atrapamiento de colapso de la  C
C movilidad.							     C
C Parametros: *Sk: Factor de estrucutura			     C
C	      *maxlist: numero de puntos de S(k)		     C
C	      *emp:fraccion de volumen				     C
C	      *g0:funcion g cero				     C
C	      *kmin:posicion del primer minimo de S(k)	             C
C para los detalles ver tesis de Laura Yeomans Reyna.                C
C agosto 2006							     C
C								     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  
       Implicit none

       Integer maxlist
       Double precision Sk(maxlist,2),gcero(maxlist),emp

       Integer in,jn
       Double precision gami, gamf, gamma 
       Double precision T(maxlist),y,dy
       Double precision sol,num,den1,den2
       Double precision S,l,g,pi,kmin

       Double precision lambda

       pi= dacos(-1.D0)

       dy= Sk(2,1)-Sk(1,1)
                 
       gami= 0.D0
CC
C metodo de iteracion de punto fijo (100 iteraciones)
CC
       Do jn=1,100

         gamf=0.D0
CC
C evaluacion del integrando
CC
          do in=1, maxlist
	     y=Sk(in,1)
             S=Sk(in,2)!facdes(Sk,maxlist,y)
             g=gcero(in) !gcero(y,ancho,tem,gr,k)
             l= lambda(y,kmin)
             num= (-1.D0)*((S-1.D0)*l)**2
             den1= 36.D0*pi*emp*(l*(S+1.D0/g)-(2.D0*gami/g)*y**2)
             den2= 1.D0/(l-gami*y**2) + 1.D0/(l*g*S - gami*y**2)
C             T= (((S+1.D0)*y**2)*num)/(S*(den1*den2))
             T(in)= num/(den1*den2)
          enddo
CC
C subrutina para hacer la integral
CC
	  call SIMPSON(T,dy,maxlist,gamf)
CC
C criterio de convergencia
CC
          sol= dabs((gamf-gami)/gami)
          if (sol.LE.1.D-4)then
             sgamma= gamf
             Return
	  endif
CC
C por si no converge
CC
	  if (jn.eq.100)then
	     sgamma=1.D20
	     Return
          endif

          gami=gamf
       Enddo
      End

      Double precision function lambda(x,xminima)

       Implicit none

       Double precision x,xminima

       lambda= 1.D0/(1.D0 + (x/xminima)**2)
       Return

      End

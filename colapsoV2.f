       Program critcolapso

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Este programa aplica el criterio de atrapamiento al estilo         C
C "colapso de la difusion" desarrollado por Laura Yeomans y el Dr.   C
C Magdaleno Medina.Esta diseñado para dos parametros de control,     C
C la aplicacion del criterio a sistemas diferentes se da cambiando   C
C el factor de  estructura para el caso especifico que estemos       C
C tratando.                                                          C
C								     C
C Agosto 2006                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none !Toda variable se tiene que declarar
CC
C Declaraciones
CC
C Parametros de control
CC
	Integer Npar2
	Double precision parametro2
	Double precision par1max,par2max
	Double precision par1min,par2min
	Double precision dpar2
	Double precision fv
CC
C variables para bisecciones sucesivas
CC
	Double precision par1post,par1ant,par1med
	Double precision integralpost,integralant,integralmed
	Double precision solpost,solant,solmed
	Double precision solprod,prodant,prodpost
CC
C variables de estatica 
CC
	Integer Nk,ig0
	parameter(Nk=2**12)
	Double precision SdeK(Nk,2),g0(Nk)
	Double precision kmin
CC
C variables de integracion
CC
	Integer iint
	Double precision y,dy,S,g,l,gam,dzeta
	Double precision T(Nk)
CC
C variables axiliares
CC
	Double precision pi
CC
C funciones que se usan
CC
	Double precision sgamma,lambda,yminima!,gcero
CC
C Saludo inicial y lectura de datos
CC
	write(*,*)' '
	write(*,*)'Este programa es para generar un diagrama de fases'
	write(*,*)'de vitrificacion con dos parametros con CM.'
	write(*,*)'Espero que ya sepas cuales son y que significan,'
        write(*,*)'echa de una vez los valores'
	write(*,*)' '
 	write(*,*)'diga el intervalo de busqueda del parametro 1'
	write(*,*)'¿maximo?'
	read(*,*)par1max
	write(*,*)'¿minimo?'
	read(*,*)par1min
	write(*,*)'diga el intervalo de busqueda del parametro 2'
	write(*,*)'¿maximo?'
	read(*,*)par2max
	write(*,*)'¿minimo?'
	read(*,*)par2min
	write(*,*)'¿numero de puntos?'
	read(*,*)Npar2
CC
C ya entraron todos los datos obligatorios
CC	
C Si ninguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente
CC
C	 write(*,*)'¿fraccion de volumen?'
C	 read(*,*)fv
CC
	write(*,*)' '
	write(*,*)'ya rugistes leon, ahora a trabajar'
	write(*,*)' '
	pi= dacos(-1.D0)
	dpar2=dabs(par2max-par2min)/dble(Npar2)
        open(8,file='fasescm.dat') 
        open(9,file='facnoergocm.dat')
        open(10,file='S(k)maxcm.dat')
        DO parametro2=par2min, par2max,dpar2
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 par1post=par1max
         par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCC
	 write(*,*)'calculando S(k) cota inferior'
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
	 dy=Sdek(2,1)-Sdek(1,1)
	 write(*,*)'calculando k del primer minimo de S(k)'
	 kmin=yminima(Sdek,Nk)
	 write(*,*)'calculando la funcion g0(k)'
	 do ig0=1,Nk
	  g0(ig0)=1.D0
CC
C Aqui se puede insertar otra propuesta para g0(k) la cual dependera
C de la Vineyard
CC
	 enddo
CCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCC 
	 write(*,*)'evaluando criterio en cota inferior'
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
C	  fv=par1ant
C    	  fv=parametro2
         gam= sgamma(Sdek,Nk,fv,g0,kmin)
         dzeta=1.D0/gam
C        gam=0.D0
         if (gam.eq.1.D20)then
          write(*,*)'no se encontro solucion para gamma',par1ant
          goto 200
         endif
CC
CC inicia la evaluacion de la integral del criterio
CC
         do iint=1, Nk
	  y=Sdek(iint,1)
          S=Sdek(iint,2)
          g= g0(iint)
          l= lambda(y,kmin)
        T(iint)=((S-1)*l*y)**2/(36.D0*pi*fv*(l*(S+1/g)-2*y**2*gam/g))
         enddo
	 call SIMPSON(T,dy,Nk,integralant)
	 write(*,*)'la evaluacion de la integral en',par1ant
	 write(*,*)' es',integralant
	 solant=1.D0 - integralant
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C pequeño 
CC
C Criterio de convergencia
CC
         if (dabs(solant).lt.1.D-4)then
          write(8,*)par1ant,parametro2,integralant,gam
          write(*,*)'se encontro solucion para: ',par1ant
          write(*,*)'evaluando factor no ergodico'
          call facnoerg(gam,Sdek,g0,kmin,Nk,fnoer)
          do i=1,Nk
           write(9,*)fnoer(i,1),fnoer(i,2)
          enddo
          write(*,*)'calculando el maximo de S(k)'
          Smx=Smax(Sdek,Nk)
          write(10,*)par1ant,parametro2,Smx
          goto 200
         endif
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCC
         write(*,*)'calculando S(k) en la cota superior'
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
	 dy=Sdek(2,1)-Sdek(1,1)
	 write(*,*)'calculando k del primer minimo de S(k)'
	 kmin=yminima(Sdek,Nk)
	 write(*,*)'calculando la funcion g0(k)'
	 do ig0=1,Nk
	  g0(ig0)=1.D0
CC
C Aqui se puede insertar otra propuesta para g0(k) la cual dependera
C de la Vineyard
CC
	 enddo
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	 write(*,*)'evaluando criterio en la cota superior'
CC
C Se inicia el calculo del parametro gamma
CC
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
C	  fv=par1post
C    	  fv=parametro2
         gam= sgamma(Sdek,Nk,fv,g0,kmin)
         dzeta=1.D0/gam
C        gam=0.D0
         if (gam.eq.1.D20)then
          write(*,*)'no se encontro solucion para gamma',par1post
          goto 200
         endif
CC
C inicia la evaluacion de la integral del criterio
CC
         do iint=1, Nk
	  y=Sdek(iint,1)
          S=Sdek(iint,2)
          g= g0(iint)
          l= lambda(y,kmin)
        T(iint)=((S-1)*l*y)**2/(36.D0*pi*fv*(l*(S+1/g)-2*y**2*gam/g))
         enddo
	 call SIMPSON(T,dy,Nk,integralpost)
	 write(*,*)'evaluacion de la integral en',par1post
	 write(*,*)' es',integralpost
	 solpost=1.D0 - integralpost
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C grande 
CC
C Criterio de convergencia
CC
         if (dabs(solpost).lt.1.D-4)then
          write(8,*)par1post,parametro2,integralpost,gam
          write(*,*)'se encontro solucion para: ',par1post
          write(*,*)'evaluando factor no ergodico'
          call facnoerg(gam,Sdek,g0,kmin,Nk,fnoer)
          do i=1,Nk
           write(9,*)fnoer(i,1),fnoer(i,2)
          enddo
          write(*,*)'calculando el maximo de S(k)'
          Smx=Smax(Sdek,Nk)
          write(10,*)par1post,parametro2,Smx
          goto 200
         endif
CC
C Termina la evaluacion de la  integral para el valor del parametro1
C mas grande
CC
         solprod=solant*solpost
         IF(solprod.lt.0.D0)then
 100	  write(*,*)'bisecciones sucesivas'
          par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) en el valor intermedio'
C llame aqui a la subrutina de S(k) evaluandola en par1med y
C parametro2
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=yminima(Sdek,Nk)
	  write(*,*)'calculando la funcion g0(k)'
	  do ig0=1,Nk
	   g0(ig0)=1.D0
CC
C Aqui se puede insertar otra propuesta para g0(k) la cual dependera
C de la Vineyard
CC
	  enddo
CCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCC 
	  write(*,*)'evaluando criterio en valor intermedio'
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
C	  fv=par1med
C    	  fv=parametro2
          gam= sgamma(Sdek,Nk,fv,g0,kmin)
          dzeta=1.D0/gam
C          gam=0.D0
          if (gam.eq.1.D20)then
           write(*,*)'no se encontro solucion para gamma',par1med
           goto 200
          endif
CC
CC inicia la evaluacion de la integral del criterio
CC
          do iint=1, Nk
	   y=Sdek(iint,1)
           S=Sdek(iint,2)
           g= g0(iint)
           l= lambda(y,kmin)
        T(iint)=((S-1)*l*y)**2/(36.D0*pi*fv*(l*(S+1/g)-2*y**2*gam/g))
          enddo
	  call SIMPSON(T,dy,Nk,integralmed)
	  write(*,*)'evaluacion de la integral en',par1med
	  write(*,*)' es',integralmed
	  solmed=1.D0 - integralmed
CC
C Termina la evaluacion de la  integral para valor de parametro 1 
C intermedio 
CC
C Criterio de convergencia
CC
          if (dabs(solmed).lt.1.D-4)then
           write(8,*)par1med,parametro2,integralmed,gam
           write(*,*)'se encontro solucion para: ',par1med
           write(*,*)'evaluando factor no ergodico'
           call facnoerg(gam,Sdek,g0,kmin,Nk,fnoer)
           do i=1,Nk
            write(9,*)fnoer(i,1),fnoer(i,2)
           enddo
           write(*,*)'calculando el maximo de S(k)'
           Smx=Smax(Sdek,Nk)
           write(10,*)par1med,parametro2,Smx
           goto 200
          endif
CC
C Decidiendo el intervalo en donde esta la solucion
CC
          prodpost= solmed*solpost
          prodant=solmed*solant
          if(prodant.lt.0.D0)then
           par1post=par1med
	   solpost=solmed
           goto 100
          endif
          if(prodpost.lt.0.D0)then
           par1ant=par1med
	   solant=solmed
           goto 100
          endif
         ENDIF
         write(*,*)'la solucion no esta en ese intervalo'
 200     continue
        ENDDO 
        close(8)
	close(9)
	close(10)
	stop
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutina del factor de estructura.Para evitar broncas con la    C
C precision de la integracion que se necesita para resolver el     C
C criterio es necesario que el valor maximo de k para el cual      C
C tienes datos de S(k) sea al menos 100, aunque siempre puedes     C
C probar para ver cual es el valor optimo.			   C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Hasta aqui todo lo necesario para calcular S(k).                 C
C Las siguentes funciones se utilizan en el criterio               C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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
      Double precision function yminima(Sk,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C Esta funcion sirve para calcular el primer minimo de S(k), pero    C
C puede adapatarse para otra funcion solo hay que ponerla en un      C
C arreglo bidimensional                                              C
C Parametros: *Sk: arreglo en que esta S(k)                          C
C                 (primera columna k y segunda S(k))                 C
C	      *N: Tamaño del arreglo                                 C
C								     C
C Agosto del 2006
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Implicit none
       
       Integer N
C       Double precision emp

       Integer h,m,j
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax,kmin,kmax

c       Double precision facdes

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

       yminima=kmin
       Return

       End


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
             S=Sk(in,2)
             g=gcero(in)
             l= lambda(y,kmin)
             num= (-1.D0)*((S-1.D0)*l)**2
             den1= 36.D0*pi*emp*(l*(S+1.D0/g)-(2.D0*gami/g)*y**2)
             den2= 1.D0/(l-gami*y**2) + 1.D0/(l*g*S - gami*y**2)
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

       Double precision function Smax(Sk,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C Esta funcion sirve para calcular el maximo de S(k), pero           C
C puede adapatarse para otra funcion solo hay que ponerla en un      C
C arreglo bidimensional                                              C
C Parametros: *Sk: arreglo en que esta S(k)                          C
C                 (primera columna k y segunda S(k))                 C
C	      *N: Tamaño del arreglo                                C
C								     C
C Agosto del 2006                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        Implicit none
       
        Integer N
        Integer j
        Double precision Sk(N,2)
        Smax=Sk(1,2)
        Do j=2, N
         if (Sk(j,2).GT.smax)then
          Smax=Sk(j,2)
         endif
        enddo
        Return
       End


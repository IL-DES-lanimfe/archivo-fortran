      Program critresergo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Este programa aplica el criterio de atrapamiento al estilo         C
C "restauracion de la egodicidad" desarrollado por Marco Chavez Rojo,C
C Pedro E. Ramirez Gonzalez y el Dr. Magdaleno Medina, basandose en  C
C la forma en que MCT obtiene el suyo.Esta diseñado para dos para-   C
C metros de control,la aplicacion del criterio a sistemas diferentes C
C se da cambiando el factor de  estructura para el caso especifico   C
C que estemos tratando.                                              C
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
	 Integer Nk
	 parameter(Nk=2**12)
	 Double precision SdeK(Nk,2)
	 Double precision kmin
CC
C variables de integracion
CC
	 Integer iint
	 Double precision y,dy,S,g,l,gama,num,den,integral
	 Double precision T(Nk)
CC
C variables axiliares
CC
	 Integer i
	 Double precision pi,a,gamam
	 Double precision Smx
	 Double precision fnoer(Nk,2)
CC
C funciones que se usan
CC
	 Double precision sgamma,lambda,yminima,Smax 
CC
C Variables del factor de estructura
CC
	 Double precision 
CC
C Saludo inicial y lectura de datos
CC
	 write(*,*)' '
	 write(*,*)'Este programa es para generar un diagrama de fases'
	 write(*,*)'de vitrificacion con dos parametros con RE.'
	 write(*,*)'Espero que ya sepas cuales son y que significan,'
       write(*,*)'echa de una vez los valores'
	 write(*,*)'Estamos en el caso de '
	 write(*,*)'parametro 1:,parametro 2:'
 	 write(*,*)'diga el intervalo de busqueda del parametro 1'
	 write(*,*)'maximo?'
	 read(*,*)par1max
	 write(*,*)'minimo?'
	 read(*,*)par1min
 	 write(*,*)'diga el intervalo de busqueda del parametro 2'
	 write(*,*)'maximo?'
	 read(*,*)par2max
	 write(*,*)'minimo?'
	 read(*,*)par2min
	 write(*,*)'numero de puntos?'
	 read(*,*)Npar2
CC
C ya entraron todos los datos obligatorios
CC	
	 pi= dacos(-1.D0)
C Si ninguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente
CC
C	 write(*,*)'fraccion de volumen?'
C	 read(*,*)fv
CC
	 write(*,*)' '
	 write(*,*)'ya rugistes leon, ahora a trabajar'
	 write(*,*)' '
       open(8,file='fasesre.dat') 
       open(9,file='facnoergore.dat')
       open(10,file='Skmaxre.dat')
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 dpar2=(par2max-par2min)/dble(Npar2)
	 DO parametro2=par2min,par2max,dpar2
	  par1post=par1max
        par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) cota inferior'
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
CC
C inicia el calculo de S(k) y g(r)
CC
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=yminima(Sdek,Nk)
	  write(*,*)'calculando la funcion g0(k)'
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	  write(*,*)'evaluando criterio en cota inferior'
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
C	  fv=par1ant
C    	  fv=parametro2
CC
CC inicia la evaluacion de la integral del criterio
CC
	  integralant=0.D0
	  Do gama=0.D0,0.1D0,1.D-5
	   a=gama
         do iint=1, Nk
	    y=Sdek(iint,1)
          S=Sdek(iint,2)
          l= lambda(y,kmin)
	    num=a*((S-1)*l*y**2)**2
	    den=(a*y**2+S*l)*(a*y**2+l)
          T(iint)=num/(36.D0*pi*fv*den)
         enddo
	   call SIMPSON(T,dy,Nk,integral)
	   if(integral.gt.integralant)then
	    integralant=integral
	    gamam=a
	   endif
	  EndDo
	  solant=1.D0 - integralant
CC
C Criterio de convergencia
CC
	  if (dabs(solant).lt.1.D-4)then
         write(8,*)sngl(par1ant),sngl(parametro2),sngl(gamam)
         write(*,*)'se encontro solucion para: ',par1ant,parametro2
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	   write(9,*)'param',sngl(par1ant),sngl(parametro2)
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1ant),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
         goto 200
        endif
CC
C fin criterio de convergencia
CC
	  write(*,*)'la evaluacion de la integral en',par1ant
	  write(*,*)' es',integralant
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C pequeño 
CC
CCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(*,*)'calculando S(k) en la cota superior'
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCC 
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
CC
C inicia la evaluacion de la integral del criterio
CC
	  integralpost=0.D0
	  Do gama=0.D0,0.1D0,1.D-5
	   a=gama
         do iint=1, Nk
	    y=Sdek(iint,1)
          S=Sdek(iint,2)
          l= lambda(y,kmin)
	    num=a*((S-1)*l*y**2)**2
	    den=(a*y**2+S*l)*(a*y**2+l)
          T(iint)=num/(36.D0*pi*fv*den)
         enddo
	   call SIMPSON(T,dy,Nk,integral)
	   if(integral.gt.integralpost)then
	    integralpost=integral
	    gamam=a
	   endif
	  EndDo
	  solpost=1.D0 - integralpost
CC
C Criterio de convergencia
CC
	  if (dabs(solpost).lt.1.D-4)then
         write(8,*)sngl(par1post),sngl(parametro2),sngl(gamam)
         write(*,*)'se encontro solucion para: ',par1post,parametro2
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	   write(9,*)'param',sngl(par1post),sngl(parametro2)
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1post),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
         goto 200
        endif
CC
C Fin criterio de convergencia
CC
	  write(*,*)'la evaluacion de la integral en',par1post
	  write(*,*)' es',integralpost
CC
C Termina la evaluacion de la  integral para valor de parametro 1 mas
C grande 
CC
        solprod=solant*solpost
        IF(solprod.lt.0.D0)then
 100     continue
	   write(*,*)'bisecciones sucesivas'
         par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   write(*,*)'calculando S(k) en el valor intermedio'
C llame aqui a la subrutina de S(k) evaluandola en par1med y
C parametro2
	   dy=Sdek(2,1)-Sdek(1,1)
	   write(*,*)'calculando k del primer minimo de S(k)'
	   kmin=yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	   write(*,*)'evaluando criterio en valor intermedio'
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
C	   fv=par1med
C    	   fv=parametro2
CC
CC inicia la evaluacion de la integral del criterio
CC
	   integralmed=0.D0
	   Do gama=0.D0,0.1D0,1.D-5
	    a=gama
          do iint=1, Nk
	     y=Sdek(iint,1)
           S=Sdek(iint,2)
           l= lambda(y,kmin)
	     num=a*((S-1)*l*y**2)**2
	     den=(a*y**2+S*l)*(a*y**2+l)
           T(iint)=num/(36.D0*pi*fv*den)
          enddo
	    call SIMPSON(T,dy,Nk,integral)
	    if(integral.gt.integralmed)then
	     integralmed=integral
	     gamam=a
	    endif
	   EndDo
	   solmed=1.D0 - integralmed
CC
C Criterio de convergencia
CC
	   if (dabs(solmed).lt.1.D-4)then
          write(8,*)sngl(par1med),sngl(parametro2),sngl(gamam)
          write(*,*)'se encontro solucion para: ',par1med,parametro2
          write(*,*)'evaluando factor no ergodico'
          call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
	    write(9,*)'param',sngl(par1med),sngl(parametro2)
          do i=1,Nk
           write(9,*)fnoer(i,1),fnoer(i,2)
          enddo
          write(*,*)'calculando el maximo de S(k)'
          Smx=Smax(Sdek,Nk)
          write(10,*)sngl(par1med),sngl(parametro2),sngl(Smx),
     &	        sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
          goto 200
         endif
CC
C Fin criterio de convergencia
CC
	   write(*,*)'la evaluacion de la integral en',par1med
	   write(*,*)' es',integralmed

CC
C Termina la evaluacion de la  integral para valor de parametro 1 
C intermedio 
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
 200    continue
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
	
      Subroutine  facnoerg(gamma,Sk,kmin,nk,f)
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
	 Double precision Sk(nk,2),f(nk,2)
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
	  den=Sk(i,2)*l
	  num=gamma*(Sk(i,1))**2
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

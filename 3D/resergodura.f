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
	 Double precision y,dy,S,l,gama,num,den,integral
	 Double precision T(Nk)
CC
C variables axiliares
CC
	 Integer i
	 Double precision pi,a,gamam,mc
	 Double precision Smx
	 Double precision fnoer(Nk,2)
CC
C funciones que se usan
CC
	 Double precision lambda,yminima,Smax
CC
C Variables del factor de estructura
CC
!	Double precision Gder(Nk,2)
CC
C Saludo inicial y lectura de datos
CC
	 write(*,*)' '
	 write(*,*)'Este programa es para generar un diagrama de fases'
	 write(*,*)'de vitrificacion con dos parametros con RE.'
	 write(*,*)'Espero que ya sepas cuales son y que significan,'
       write(*,*)'echa de una vez los valores'
	 write(*,*)'Estamos en el caso de esferas duras'
	 write(*,*)'parametro 1:fracc. de volumen,parametro 2 no hay'
 	 write(*,*)'diga el intervalo de busqueda del parametro 1'
	 write(*,*)'maximo?'
	 read(*,*)par1max
	 write(*,*)'minimo?'
	 read(*,*)par1min
!  	write(*,*)'diga el intervalo de busqueda del parametro 2'
! 	write(*,*)'¿maximo?'
! 	read(*,*)par2max
! 	write(*,*)'¿minimo?'
! 	read(*,*)par2min
! 	write(*,*)'¿numero de puntos?'
! 	read(*,*)Npar2
CC
C ya entraron todos los datos obligatorios
CC	
	 pi= dacos(-1.D0)
C Si ninguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente
CC
C	 write(*,*)'¿fraccion de volumen?'
C	 read(*,*)fv
CC
	 write(*,*)' '
	 write(*,*)'ya rugistes leon, ahora a trabajar'
	 write(*,*)' '
       open(8,file='vidrioduraex.dat') 
       open(9,file='facnoergoduraex.dat')
       open(10,file='Skmaxduraex.dat')
	 mc=0.73D0
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
! 	dpar2=(par2max-par2min)/dble(Npar2)
! 	DO parametro2=par2min,par2max,dpar2
	 par1post=par1max
       par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 write(*,*)'calculando S(k) cota inferior'
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
CC
C inicia el calculo de S(k) y g(r)
CC
       call Sdeked(Sdek,par1ant,Nk,100.D0,10.D0)
	 dy=Sdek(2,1)-Sdek(1,1)
	 write(*,*)'calculando k del primer minimo de S(k)'
	 kmin=mc*yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	 write(*,*)'evaluando criterio en cota inferior'
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	 fv=par1ant
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
	   num=a*((S-1.D0)*l*y**2)**2
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
        write(8,*)sngl(par1ant),sngl(gamam)
        write(*,*)'se encontro solucion para: ',par1ant,gamam
        write(*,*)'evaluando factor no ergodico'
        call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
        do i=1,Nk
         write(9,*)fnoer(i,1),fnoer(i,2)
        enddo
        write(*,*)'calculando el maximo de S(k)'
        Smx=Smax(Sdek,Nk)
        write(10,*)sngl(par1ant),sngl(Smx),
     &            sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
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
       call Sdeked(Sdek,par1post,Nk,100.D0,10.D0)
	 dy=Sdek(2,1)-Sdek(1,1)
	 write(*,*)'calculando k del primer minimo de S(k)'
	 kmin=mc*yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCC 
	 write(*,*)'evaluando criterio en la cota superior'
CC
C Se inicia el calculo del parametro gamma
CC
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	 fv=par1post
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
        write(8,*)sngl(par1post),sngl(gamam)
        write(*,*)'se encontro solucion para: ',par1post,gamam
        write(*,*)'evaluando factor no ergodico'
        call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
        do i=1,Nk
         write(9,*)fnoer(i,1),fnoer(i,2)
        enddo
        write(*,*)'calculando el maximo de S(k)'
        Smx=Smax(Sdek,Nk)
        write(10,*)sngl(par1post),sngl(Smx),
     &	         sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
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
 100    continue
	  write(*,*)'bisecciones sucesivas'
        par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) en el valor intermedio'
Cllame aqui a la subrutina de S(k) evaluandola en par1med y parametro2
        call Sdeked(Sdek,par1med,Nk,100.D0,10.D0)
	  dy=Sdek(2,1)-Sdek(1,1)
	  write(*,*)'calculando k del primer minimo de S(k)'
	  kmin=mc*yminima(Sdek,Nk)
CCCCCCCCCCCCCCCCCCCC Fin modulo de estaticaCCCCCCCCCCCCCCCCCCCCCCCC 
	  write(*,*)'evaluando criterio en valor intermedio'
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	  fv=par1med
C    	  fv=parametro2
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
         write(8,*)sngl(par1med),sngl(gamam)
         write(*,*)'se encontro solucion para: ',par1med,gamam
         write(*,*)'evaluando factor no ergodico'
         call facnoerg(gamam,Sdek,kmin,Nk,fnoer)
         do i=1,Nk
          write(9,*)fnoer(i,1),fnoer(i,2)
         enddo
         write(*,*)'calculando el maximo de S(k)'
         Smx=Smax(Sdek,Nk)
         write(10,*)sngl(par1med),sngl(Smx),
     &              sngl(dsqrt(gamam)*(6.d0*fv/pi)**(1.d0/3.D0))
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
 200   continue
!	ENDDO
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

      Subroutine  Sdeked(Sdk,fracvol,Np,kmaxima,rmaxima)
!Sdeked(Sdk,Gr,fracvol,Np,kmaxima,rmaxima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Subrutina para evaluar el factor de estructura (S(k)), la función  C
C de distribución radial (g(r))  sistema con potencial de interac-   C
C ción tipo esfera dura utilizando la estrucutura de esfera dura     C
C obtenida con la solución exacta propuesta por M.S. Wertheim (PRL   C
C  vol.10,no.8, pp.321,1963) para la cerradura de Perckus-Yevick  y  C
C la expresion para S(k) en funcion de la transformada de Lapalace   C
C de g(r) (Physica 149A (1988),123-163).                             C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Parametros: *Sdk:arreglo para el factor de estructura (S(k))       C
C             *Gr: arreglo para la funcion de dist. radial (g(r))    C
C             *Yr: arreglo para la funcion Y(r)                      C
C             *fracvol:fraccion de volumen                           C
C             *Np:numero de puntos a evaluar en S(k) y g(r)          C
C             *kmaxima:valor de k maximo para evaluar S(k)           C
C             *rmaxima:valor de r maximo para evaluar g(r)           C
C             *n: exponente del potencial                            C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        integer Np,Npg
        double precision fracvol,fv,rmaxima,kmaxima
        double precision Sdk(Np,2)!,Gr(Np,2),Yr(Np,2)
!parametro para extende k para ganar buena g(r)
!        parameter(Npg=50)
        integer i
        double precision K,dK
!arreglo auxiliar para buena g(r)
!        double precision Sdkg(Npg*Np,2)
        double precision C1,C2
C funciones auxiliares para evaluar S(k)
        complex G
        double precision facdes
C Solo por comodidad
        fv=fracvol
C primero evaluamos una buena S(k) usado los resultados exactos 
C mediante la función facdes
        K=0.001D0
        dK=kmaxima/dble(Np)
        do i=1,Np
         Sdk(i,1)=K
         Sdk(i,2)=facdes(K,fv)                 
         K=K+dK
        enddo
CC
C definición del espaciado para una buena g(r)
CC
!         dK=1000.D0/Dble(Np)
!         K=0.001D0
CC
C Se usa el arreglo Sdkg para guardar información a k's
C grandes y poder hacer una buena evaluación de g(r)
CC                
!          do i=1,Npg*Np
!           Sdkg(i,1)=K
!           Sdkg(i,2)=facdes(K,fv)
!           K=K+dK
!          enddo
CC
C ya tenemos S(k) en hasta una k grande ahora le sacamos la
C transformada inversa de Fourier para obtener g(r)
CC
!         call  TIF3D(Npg*Np,Sdkg,Np,Gr,fv,rmaxima)
CC        
C ya se calcularon S(k) y g(r)
C Calcualmos Y(r)
CC
!         Do i=1,Np
!          Yr(i,1)=Gr(i,1)
!          if(Gr(i,1)/lam.le.1.D0) then
!                   
!            C1=(1.D0+2.D0*fve)**2-6.D0*fve*(1.D0+0.5D0*fve)**2*Gr(i,1)
!              C2=fve*(1.D0+2.D0*fve)**2*0.5D0*Gr(i,1)**3           
!              Yr(i,2)=(C1+C2)/(1.D0-fve)**4
!              Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)
!            else
!             
!              Yr(i,2)=Gr(i,2)
! 	     Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)           
!            endif
!            
!         Enddo
        Return
       End

       Complex function G(x,emp)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta funcion calcula la transformada de Lapalce de la funcion de  C
C distibuion radial. Como auxliar a calculo de factor de estructura C
C de esfera dura                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar

C declaracion de argumentos
        Complex x
        double precision emp

C variables locales
        Complex a,b,c,d,num,den1,den2,den3,den

C calculos auxiliares

         a= cmplx(1.D0+emp/2.D0,0.D0)
         b= cmplx(1.D0+2.D0*emp,0.D0)
         c= cmplx(1.D0-emp,0.D0)
         d= cmplx(18.D0,0.D0)

         num= x*(a*x+b)
         den1= cmplx(12.D0*emp,0.D0)*(a*x+b)
         den2= (c**2*x**3+cmplx(6.D0*emp,0.D0)*c*x**2)*exp(x)
         den3=(d*cmplx(emp,0.D0)**2*x-cmplx(12.D0*emp,0.D0)*b)*exp(x)
         den= den1+den2+den3

C valor que regresa la funcion

         G= num/den

         return

       end

      
       
       double precision function facdes(z,fv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Esta funcion calcula el factor de estructura.C
C de esfera dura.Necesario usar la funcion G   C
C auxiliar.                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
        Implicit none !aguas toda variable se tiene que declarar

C decalracion de argumentos

        double precision  z, fv,corte

C decalraciones locales 
        Complex cz, mcz

C declaracion funcion auxiliar      
        complex  G
        
C variables auxiliares

        double precision T0,T2,F1,F2 
        
        parameter(corte=0.5D0) !parametro para fijar k's pequeñas
        
C Evitando problemas a k's chicas

        if(z.lt.corte) then
        
            T0=(2.D0*(-1.D0+fv)**4*(0.5D0+fv))/(1.D0+2.D0*fv)**3
            F1=0.4D0*(-1.D0+fv)**4*fv*(0.5D0+fv)
            F2=(4.D0+(-2.75D0+fv)*fv)/((1.D0+2.D0*fv)**5) 
            T2=F1*F2
            
            facdes=T0+T2*z**2
           
            Return
          
        else

          cz= dcmplx(0.D0,z)
          mcz= dcmplx(0.D0,(-1.D0)*z)

C valor de la funcion
          facdes= 1.D0 + 12.D0*fv*(dble((G(mcz,fv)-G(cz,fv))/cz))
      
          Return
          
        endif
        
       End
        

      Subroutine TIF3D(NK,SDK,N,GDR,fv,rmx)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION SDK(NK,2),GDR(N,2)
      DIMENSION R(N),HR(N),RK(NK),HK(NK)
      DIMENSION XKR(NK),XIN(NK)

      
C     ************************************************
C     PARAMETROS GENERALES
C     ************************************************

       PI=4.D0*ATAN(1.D0)
       DENS=6.D0*fv/PI

C     ************************************************
C     LECTURA DEL ARCHIVO DE ENTRADA (R,GR)
C     ************************************************

      DO LL=1,NK,1

       RK(LL)=SDK(LL,1)
       HK(LL)=(SDK(LL,2)-1.D0)/DENS
      

      ENDDO

      KMAX=RK(NK)
      DRK=RK(2)-RK(1)

      RMAX=rmx
      DR=RMAX/(1.D0*N)
   
      DO  JJ=1,N,1
      R(JJ)=real(JJ)*DR
      ENDDO

      DO  J=1,N,1
         Do K=1,NK,1
            XKR(K)=R(J)*RK(K)
            
C     ***************************************************
C     EVALUACION DEL INTEGRANDO DE FK
C     ***************************************************

            XIN(K)=RK(K)*RK(K)*HK(K)*dSIN(XKR(K))/XKR(K)
          
         Enddo

C     ***************************************************
C     INTEGRACION POR SIMPSON
C     ***************************************************

         CALL SIMPSON(XIN,DRK,NK,RES1)

         HR(J)=4.D0*PI*RES1 /(2.D0*PI)**3  
     
         GDR(J,1)=R(J)
         GDR(J,2)=HR(J)+1.D0
         

      ENDDO

      RETURN
      END


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

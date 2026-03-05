       Program critcolapso

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Este programa aplica el criterio de atrapamiento al estilo         C
C "colapso de la difusion" desarrollado por Laura Yeomans y el Dr.   C
C Magdaleno Medina.Esta diseñado para dos parametros de control,    C 
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
	Double precision pi,Smx,fnoer(Nk,2)
	Double precision Gr(Nk,2),Yr(Nk,2)
	Integer i
CC
C funciones que se usan
CC
	Double precision sgamma,lambda,yminima,Smax!,gcero
CC
C Saludo inicial y lectura de datos
CC
	write(*,*)' '
	write(*,*)'Este programa es para generar un diagrama de fases'
	write(*,*)'de vitrificacion con dos parametros con CM.'
	write(*,*)'Espero que ya sepas cuales son y que significan,'
	write(*,*)'este es el caso de potencial (1/r)**n'
	write(*,*)'parametro1=fracc. de vol., parametro2=n'
        write(*,*)'echa de una vez los valores'
	write(*,*)' '
 	write(*,*)'establezca el intervalo de busqueda del parametro 1'
	write(*,*)'maximo?'
	read(*,*)par1max
	write(*,*)'minimo?'
	read(*,*)par1min
	write(*,*)'establezca el intervalo de busqueda del parametro 2'
	write(*,*)'maximo?'
	read(*,*)par2max
	write(*,*)'minimo?'
	read(*,*)par2min
	write(*,*)'numero de puntos?'
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
        open(8,file='fasescmsu.dat') 
        open(9,file='facnoergocmsu.dat')
        open(10,file='S(k)maxcmsu.dat')
        DO parametro2=par2min, par2max,dpar2
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 par1post=par1max
         par1ant= par1min
 100     continue 
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 write(*,*)'calculando S(k) cota inferior'
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
	 call Sdekes(Sdek,Gr,Yr,par1ant,Nk,100.D0,10.D0,parametro2)
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
CCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCC 
	 write(*,*)'evaluando criterio en cota inferior'
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
         fv=par1ant
C    	  fv=parametro2
         gam= sgamma(Sdek,Nk,fv,g0,kmin)
         dzeta=1.D0/gam
C        gam=0.D0
         if (gam.eq.1.D20)then
          write(*,*)'no se encontro solucion para gamma',par1ant
          par1ant=par1ant+1.D0
          goto 100
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
 300     continue
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         write(*,*)'calculando S(k) en la cota superior'
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
	 call Sdekes(Sdek,Gr,Yr,par1post,Nk,100.D0,10.D0,parametro2)
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
CCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCC 
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
         gam= sgamma(Sdek,Nk,fv,g0,kmin)
         dzeta=1.D0/gam
C        gam=0.D0
         if (gam.eq.1.D20)then
          write(*,*)'no se encontro solucion para gamma',par1post
          par1post=par1post-1.D0
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
	  write(*,*)'bisecciones sucesivas'
          par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) en el valor intermedio'
C llame aqui a la subrutina de S(k) evaluandola en par1med y parametro2
	  call Sdekes(Sdek,Gr,Yr,par1med,Nk,100.D0,10.D0,parametro2)
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
CCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCC 
	  write(*,*)'evaluando criterio en valor intermedio'
CC
C Si alguno de los parametros es la fraccion de volumen entonces 
C descomenta lo siguiente segun el caso
CC
	  fv=par1med
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
           goto 300
          endif
          if(prodpost.lt.0.D0)then
           par1ant=par1med
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

      Subroutine  Sdekes(Sdk,Gr,Yr,fracvol,Np,kmaxima,rmaxima,n)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Subrutina para evaluar el factor de estructura (S(k)), la función de C
C distribución radial (g(r)) y la función " " (Y(r)) para un sistema   C
C con potencial de interacción tipo (1/r)**n utilizando la estrucutura C
C de esfera dura obtenida con la solución exacta propuesta por M.S.    C
C Wertheim (PRL vol.10,no.8, pp.321,1963) para la cerradura de         C
C Perckus-Yevick  y la expresion para S(k) en funcion de la transfor-  C
C mada de Lapalace de g(r) (Physica 149A (1988),123-163).              C
C El proceso de reescalamiento para obtener la estructura de esfera    C
C sueve a partir de esfera dura viene descrito en Theory of simple     C
C liquids de J.P. Hansen y I.R. McDonald secc. 6.3, pp. 155            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C Parametros: *Sdk:arreglo para el factor de estructura (S(k))         C
C             *Gr: arreglo para la funcion de dist. radial (g(r))      C
C             *Yr: arreglo para la funcion Y(r)
C             *fracvol:fraccion de volumen                             C
C             *Np:numero de puntos a evaluar en S(k) y g(r)            C
C             *kmaxima:valor de k maximo para evaluar S(k)             C
C             *rmaxima:valor de r maximo para evaluar g(r)             C
C             *n: exponente del potencial                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        integer Np,Npg
        double precision fracvol,fv,rmaxima,kmaxima,n
        double precision Sdk(Np,2),Gr(Np,2),Yr(Np,2)
        parameter(Npg=50)!parametro para extende k para ganar buena g(r)
        integer i
        double precision K,dK
        double precision Sdkg(Npg*Np,2)!arreglo auxiliar para buena g(r)
        double precision C1,C2
	double precision lam,fve,Eul !variables para esfera dura equiv.
	parameter (Eul=0.5772D0)
   
C funciones auxiliares para evaluar S(k)

        complex G
        double precision facdes

C Solo por comodidad

        fv=fracvol
C primero evaluamos una buena S(k) usado los resultados exactos 
C mediante la función facdes
        
C razon entre diametro equiv. de esfera dura y diametro de esfera suave.
	lam=1.D0+Eul/n
	fve=fv*lam**3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
        K=0.001D0
        dK=kmaxima/dble(Np)
                        
        do i=1,Np
           
           Sdk(i,1)=lam*K
           Sdk(i,2)=facdes(lam*K,fve)                 
           K=K+dK

        enddo

C definición del espaciado para una buena g(r)

        dK=1000.D0/Dble(Np)
        K=0.001D0
        
C Se usa el arreglo Sdkg para guardar información a k's
C grandes y poder hacer una buena evaluación de g(r)
                
         do i=1,Npg*Np
           
           Sdkg(i,1)=lam*K
           Sdkg(i,2)=facdes(lam*K,fve)
           K=K+dK

         enddo
C ya tenemos S(k) en hasta una k grande ahora le sacamos la transformada
C inversa de Fourier para obtener g(r)

         call  TIF3D(Npg*Np,Sdkg,Np,Gr,fve,lam*rmaxima)
        
C ya se calcularon S(k) y g(r)
C Calcualmos Y(r)

        Do i=1,Np
            
            Yr(i,1)=Gr(i,1)
            
            if(Gr(i,1)/lam.le.1.D0) then
                  
           C1=(1.D0+2.D0*fve)**2-6.D0*fve*(1.D0+0.5D0*fve)**2*Gr(i,1)
             C2=fve*(1.D0+2.D0*fve)**2*0.5D0*Gr(i,1)**3           
             Yr(i,2)=(C1+C2)/(1.D0-fve)**4
             Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)
           else
            
             Yr(i,2)=Gr(i,2)
	     Gr(i,2)=Yr(i,2)*dexp(-(1.D0/(lam*Gr(i,1)))**n)           
           endif
           
        Enddo
        
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


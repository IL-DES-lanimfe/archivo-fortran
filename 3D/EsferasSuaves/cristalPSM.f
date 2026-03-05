       Program cristal

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    	C
C Este programa aplica el criterio de Verlet para determinar el punto	C
C de cristalizacion en un liquido.Esta diseñado para dos parametros  	C
C de control,la aplicacion del criterio a sistemas diferentes se da  	C
C cambiando el factor de  estructura para el caso especifico que     	C
C estemos tratando.                                                  	C
C								     	C
C Septiembre 2006                                                    	C
C								     	C
C Adaptacion a Esferas Penetrables (PSM)				C
C									C
C Julio 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
	Double precision solpost,solant,solmed
	Double precision solprod,prodant,prodpost
	Double precision Smx
CC
C variables de estatica 
CC
	Integer Nk
	parameter(Nk=2**12)
	Double precision SdeK(Nk,2)
CC
C variables axiliares
CC
	Integer i
	Double precision pi,dy
!	Double precision Rk(Nk),Sk(Nk)
CC
C funciones que se usan
CC
	Double precision Smax,SdekPSM
CC
C Saludo inicial y lectura de datos
CC
	write(*,*)' '
	write(*,*)'Este programa es para generar un diagrama de fases'
	write(*,*)'de cristalizacion con dos parametros'
	write(*,*)'Este es el caso del pot. de Yukawa atractivo'
        write(*,*)'parametro1:fraccion de volumen'
	write(*,*)'parametro2:T(temperatura)'
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
! 	write(*,*)'numero de puntos?'
! 	read(*,*)Npar2
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
!	z=100.D0
	dpar2=0.1D0 !dabs(par2max-par2min)/dble(Npar2)
        open(8,file='cristalPSM.dat') 
	dy=0.25D-2
      DO parametro2=par2min, par2max,dpar2
CC
C se usa el metodo de bisecciones sucesivas. Parametro 2 fijo y 
C bisecciones sucesivas en parametro 1 para localizar la transicion
C luego se cambia el parametro 2 y se repite la operacion
CC
	 par1post=par1max
       par1ant= par1min
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 write(*,*)'calculando S(k) cota inferior'
CC
C ponga aqui la subrutina de S(k) evaluandola en par1ant y parametro2
CC
!	 call msa(par1ant,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	 do i=1,Nk
	  Sdek(i,1)=dble(i)*dy
	  Sdek(i,2)=SdekPSM(Sdek(i,1),par1ant,parametro2)
	 enddo
CC
       write(*,*)'calculando el maximo de S(k)'
       Smx=Smax(Sdek,Nk)
	 solant=2.85D0-Smx
	 write(*,*)'maximo de S(k)',Smx,par1ant
CC
C Criterio de convergencia
CC
 	 if (dabs(solant).lt.1.D-4)then
          write(8,*)par1ant,parametro2,Smx
          write(*,*)'se encontro solucion para: ',par1ant,parametro2
          goto 200
         endif	
CC
C fin criterio de convergencia
CC
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         write(*,*)'calculando S(k) en la cota superior'
CC
C ponga aqui la subrutina de S(k) evaluandola en par1post y parametro2
CC
!	 call msa(par1post,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	 do i=1,Nk
	  Sdek(i,1)=dble(i)*dy
	  Sdek(i,2)=SdekPSM(Sdek(i,1),par1post,parametro2)
	 enddo
CC
	write(*,*)'calculando el maximo de S(k)'
        Smx=Smax(Sdek,Nk)
	solpost=2.85D0-Smx
	write(*,*)'maximo de S(k)',Smx,par1post
CC
C Criterio de convergencia
CC
        if (dabs(solpost).lt.1.D-4)then
          write(8,*)par1post,parametro2,Smx
          write(*,*)'se encontro solucion para: ',par1post,parametro2
          goto 200
         endif 
CC
C Termina la evaluacion de la  integral para el valor del parametro1
C mas grande
CC
         solprod=solant*solpost
         IF(solprod.lt.0.D0)then
	  write(*,*)'bisecciones sucesivas'
 100	  continue
          par1med=(par1post+par1ant)/2.D0
CCCCCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  write(*,*)'calculando S(k) en el valor intermedio'
CC
C llame aqui a la subrutina de S(k) evaluandola en par1med y parametro2
CC
!	  call msa(par1med,parametro2,z,0.01D0,Nk,Rk,Sk,xi)
	  do i=1,Nk
	  Sdek(i,1)=dble(i)*dy
	  Sdek(i,2)=SdekPSM(Sdek(i,1),par1med,parametro2)
	  enddo
CC
	  write(*,*)'calculando el maximo de S(k)'
          Smx=Smax(Sdek,Nk)
	  solmed=2.85D0 - Smx
	  write(*,*)'maximo de S(k)',Smx,par1med
CC
C Criterio de convergencia
CC
          if (dabs(solmed).lt.1.D-4)then
           write(8,*)par1med,parametro2,Smx
           write(*,*)'se encontro solucion para: ',par1med,parametro2
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
	stop
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutina del factor de estructura 			           C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION SdekPSM(k,eta,tem)
	implicit none
	Double precision k,eta,tem
	Double precision Senos
	Senos=(dsin(k)-k*dcos(k))/k**3
	SdekPSM=1.D0/(1.D0+24.D0*eta*Senos/tem)
	return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Hasta aqui todo lo necesario para calcular S(k).                 C
C Las siguentes funciones se utilizan en el criterio               C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


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
CC
	Smax=Sk(251,2)
        Do j=252, N
         if (Sk(j,2).GE.Smax)then
          Smax=Sk(j,2)
         endif
        enddo
        Return
       End


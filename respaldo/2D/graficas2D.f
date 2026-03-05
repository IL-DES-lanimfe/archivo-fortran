      Program graficasdz

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C	 
C Este programa evalua la integral de la funcional del criterio de	C	 
C atrapamiento para el caso de "restauracion de la ergdicidad".Esta	C
C hecho para dos parametros de control (por ejemplo frac. de vol.		C 
C y temperatura), pero se puede expandir facilmente para mas para-	C
C metros.Para cambiar de sistema solo es necesario cambiar subrutina	C
C y variables para el factor de estrucutura.							C	  
C																	C
C Agosto de 2006														C
C																	C
C Modificado a dos dimensiones.										C
C																	C
C Enero 2008															C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Implicit none
C indices para los ciclos en general
       integer iint,i
C variables para la integral sobre vectores de onda
       double precision  y,dy,S  
C variables para evaluar la funcional del criterio
       double precision gama,dgama,gamamax,gamamin 
C numero de puntos a evaluar en gama
	 integer Ngama
C variables para las integrales de funcional (temporales)
       double precision num,den
C variables las evaluaciones de la funcional del criterio
       double precision integral
C variables para parametros del criterio (labda,1er min S(k))
       double precision l,kmin
C Numero de puntos en vecs. de onda
	 Integer nk
	 parameter (nk=1000)
C arreglos para S(k) y g(r)
       double precision  Sdek(nk,2),T(nk)
C variables para evaluar en los parametros de control
	 double precision par2max,par2min,dpar2
C variables auxiliares
       Double precision pi,fa 
C parametros S(k) y g(r)
	 Double precision parametro1,parametro2
CC
CCCCCCCCCCCCCCCC Funciones que utiliza el programa CCCCCCCCCCCCCCCCCCC
CC
       double precision  lambda,kminima
CCCCCCCCCCCCCCCCCCC Archivos de salida CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C archivo de salida funcional criterio 
       open(9,file='funcrit.dat')
	 open(10,file='solsgam.dat')
CC
       pi= dacos(-1.D0)
CC
C saludo de beinvenida y entrada de datos
CC
	 write(*,*)' '
	 write(*,*)'Bienvenido a graficas, vamos a calcular'
	 write(*,*)'las integrales del criterio y de la funcional de'
	 write(*,*)'gamma para asi tener una idea de donde esta la'
	 write(*,*)'transicion vitrea'
	 write(*,*)' '
CC
C El programa te hace las graficas de la funcional del criterio como 
C funcion del parametro 2 manteniendo el parametro 1 fijo, es decir 
C cada corrida del programa te genera una iso-parametro1 por decirlo 
C de alguna manera.
CC
       write(*,*)'cuanto vale parametro1?'
       read (*,*)parametro1
       write(*,*)'cuanto vale parametro 2?'
       write(*,*)'maximo?'
       read (*,*)par2max
	 write(*,*)'minimo?'
       read (*,*)par2min
	 write(*,*)'Longitud del incremento en parametro2?'
       read (*,*)dpar2

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
CC
C se concoce S(k) y g(r)
CC
	  dy=Sdek(2,1)-Sdek(1,1)
        write(*,*)'calculando la k del primer minimo de S(k)'
        kmin=kminima(SdeK,nk)
CC
C ahora se conoce la posicion del primer minimo de S(k)
CC
CCCCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCC
CCe
	  write(*,*)'calculando la funcional del criterio'
	  gamamin=0.D0
	  gamamax=0.1D0
	  dgama=1.D-5
CC
C Segun sea el caso
CC
C	  fa=parametro1
C	  fa=parametro2

        Do gama=gamamin,gamamax,dgama
         integral= 0.D0
CC
C evaluacion del integrando de la funcional del criterio
CC
         do iint=1, Nk
		y=Sdek(iint,1)
          S=SdeK(iint,2)
          l= lambda(y,kmin)
          num=gama*y**3*((S-1.D0)*l)**2
          den=16.D0*fa*(gama*y**2+S*l)*(gama*y**2+l)
          T(iint)= num/den
         enddo
CC
C evaluacion de la integral
CC
         call SIMPSON(T,dy,Nk,integral)
	   write(9,*)gama,integral
CC
C Buscamos soluciones en el lado vitreo
CC
	   if((gama.le.0.01).and.(dabs(integral-1.D0).le.1.D-4))then
	    write(10,*)parametro2,gama,1.D0/gama
	   endif
CC
CC Termina la evaluacion de la  integral de la funcional de gamma
CC
C aumentamos gamma para el siguiente punto
CC
        EndDo
        write(*,*)'funcional del criterio lista',parametro2
CC Termina el ciclo de evaluacion de la funcional
       ENDDO
CC Termina ciclo de parametros
       close(9)
	 close(10)
CCC Fin del programa CCCC
       stop
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                  C
C Subrutinas del factor de estructura							     C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Hasta aqui factor de estrucutura									C
C																	C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                    C
C Las siguentes funciones se utilizan en el criterio C
C                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Double precision function kminima(Sk,N)

       Implicit none
       
       Integer N
       Integer h,m,j
       Double precision Sk(N,2)
       Double precision smin,smax,kmin,kmax

       Double precision facdes

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

      Double precision function lambda(x,xminima)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Funcion interpoladora lambda(k),se usa en la cerradura para			C 
C Cs(k,z) en la Teoria Autoconsistente. Ver detalles en la tesis		C
C de Laura Yeomans Reyna.												C
C Parametros: *x: valor de k en que queremos evaluar la funcion		C
C	      *xminima: valor de k en donde esta el primer minimo		C
C			de S(k)													C
C																	C
C Agosto de 2006														C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       Implicit none

       Double precision x,xminima

       lambda= 1.D0/(1.D0 + (x/xminima)**2)
       Return

      End

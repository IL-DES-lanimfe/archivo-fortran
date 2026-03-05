      Program graficas

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C 
C Este programa evalua la integral de la funcional del parametro    C 
C gamma y la integral del criterio de atrapamiento.Esta diseñado    C
C para dos parametros de control (por ejemplo frac. de vol. y tempe-C
C ratura), pero se puede expandir facilmente para mas parametros.   C
C Para cambiar de sistema solo es necesario cambiar subrutina y     C
C variables para el factor de estrucutura.                          C  C								    C
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
        double precision  Sdek(nk,2),T(nk),g0(nk),Gr(nk,2),yr(nk,2)
C parametros S(k) y g(r)
	Double precision parametro1,fv 
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
	write(*,*)'Es el caso de pot(u=1/r**n)'
	write(*,*)' '
CC
C El programa te hace las graficas de la funcional del criterio como 
C funcion del parametro 2 manteniendo el parametro 1 fijo, es decir 
C cada corrida del programa te genera una iso-parametro1 por decirlo 
C de alguna manera.
CC
        write(*,*)'¿cuanto vale (n)?'
        read (*,*)parametro1
        write(*,*)'¿cuanto vale (fv) max.?'
        read (*,*)par2max
	write(*,*)'¿cuanto vale (fv) min?'
        read (*,*)par2min
	write(*,*)'¿Longitud del incremento en (fv)?'
        read (*,*)dpar2
	write(*,*)'¿cuantos puntos quieres en la fun. de gamma?'
        read (*,*)Ngam
CC
C Empieza el ciclo para barrer sobre parametro 2
CC
	DO fv=par2min,par2max,dpar2
          write(*,*)'parametros',parametro1,fv
CC
CCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculado S(k)'
CC Aqui se inserta la subrutina del factor de estrucutura
          call Sdekes(Sdek,Gr,Yr,fv,Nk,100.d0,10.d0,parametro1)
cccccccccccccccc   OJO: PARAMETRO1 = FRACCION DE VOLUMEN fv
cccccccccccccccc        PARAMETRO2 = EXPONENTE n DEL POTENCIAL
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
           write(*,*)'no se encontro solucion para gamma',fv
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
           write(9,*)fv,criterio
           write(*,*)'criterio listo',fv,criterio
              
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
          write(*,*)'funcional de gamma lista',fv
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

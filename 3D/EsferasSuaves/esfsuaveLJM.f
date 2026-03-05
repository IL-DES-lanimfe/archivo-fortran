      Program estructuraesfsuave

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C Este programa sirve para probar una subrutina para calcular el C
C el factor de estructura,la funcion de distribucion radial y la C
C funcion y(r) para un potencial de esfera suave (1/r)**n. Usa   C
C los resultados analiticos de M.S. Wertheim (PRL vol.10,no.8,   C
C pp.321,1963) para la solucion de la ec. de Ornstein-Zernik con C
C una cerradura de Percus-Yevick y un proceso de equivalencia    C
C entre esfera dura y esfera suave descrito en Hansen y McDonald.C
C 17/Julio/2006                                                  C
C															   C
C Adaptado a LJM												   C
C															   C
C Enero 2007                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Implicit none !aguas toda se tiene que declarar
CC
       integer Nk,i
       parameter (Nk=2**10)
       double precision Gder(Nk,2),Yder(Nk,2),Sdek(Nk,2)
       double precision n,fv,kmax,rmax,nu
       double precision siginv,diametro,pi
CC
       open(10,file='facdesLJMVWfv5146nu18.dat')
       open(11,file='fdisradLJMVWfv5146nu18.dat')
!        open(12,file='YderLYMVMnu50.dat')
CC 
       write(*,*)'Fracc de volumen?'
       read(*,*)fv
       write(*,*)'valor max. de k para el fac. de estructura?'
       read(*,*)kmax
       write(*,*)'valor max. de r para la fun. de dist. radial?'
       read(*,*)rmax
       write(*,*)'exponente del potencial?'
       read(*,*)nu 
CC
       pi=4.d0*datan(1.d0)
       siginv=diametro(nu)
!        fv=pi*n/(6.D0*siginv**3)
       call Sdekes(Sdek,Gder,Yder,fv,Nk,kmax,rmax,siginv,nu)
CC
       do i=1,Nk
        write(10,*)Sdek(i,1),Sdek(i,2)
        write(11,*)Gder(i,1),Gder(i,2)
!         write(12,*)Yder(i,1),Yder(i,2)
       enddo
       close(10)
       close(11)
!        close(12)
       stop
      END

      Subroutine  Sdekes(Sdk,Gr,Yr,fracvol,Np,kmaxima,rmaxima,lam,n)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C Subrutina para evaluar el factor de estructura (S(k)), la función   C
C de distribución radial (g(r)) y la función " " (Y(r)) para un sis-  C
C tema con potencial de interacción tipo (1/r)**n utilizando la       C
C estrucutura de esfera dura obtenida con la solución exacta propuestaC
C por M.S. Wertheim (PRL vol.10,no.8, pp.321,1963) para la cerradura  C
C de Perckus-Yevick  y la expresion para S(k) en funcion de la	      C
C transformada de Lapalace de g(r) (Physica 149A (1988),123-163).     C
C El proceso de reescalamiento para obtener la estructura de esfera   C
C suave a partir de esfera dura viene descrito en Theory of simple    C
C liquids de J.P. Hansen y I.R. McDonald secc. 6.3, pp. 155           C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        Implicit none !aguas toda variable se tiene que declarar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C Parametros: *Sdk:arreglo para el factor de estructura (S(k))       C
C             *Gr: arreglo para la funcion de dist. radial (g(r))    C
C             *Yr: arreglo para la funcion Y(r)
C             *fracvol:fraccion de volumen                           C
C             *Np:numero de puntos a evaluar en S(k) y g(r)          C
C             *kmaxima:valor de k maximo para evaluar S(k)           C
C             *rmaxima:valor de r maximo para evaluar g(r)           C
C	      *lam:diam. duro/diam. suave
C             *n: exponente del potencial                            C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        integer Np,Npg
        double precision fracvol,fv,rmaxima,kmaxima,n
        double precision Sdk(Np,2),Gr(Np,2),Yr(Np,2)
!parametro para extende k para ganar buena g(r)
        parameter(Npg=10)
        integer i
        double precision K,dK,r
!arreglo auxiliar para buena g(r)
        double precision Sdkg(Npg*Np,2)
        double precision C1,C2
!variables para esfera dura equiv.
	double precision lam,fve,pot 
	double precision diametro
C funciones auxiliares para evaluar S(k)
        complex G
        double precision facdes
C variables necesarias para correccion de Verlet-Weis
	double precision Anum,Aden,A,mu,fvw,gc
C Solo por comodidad

        fv=fracvol
C primero evaluamos una buena S(k) usado los resultados exactos 
C mediante la función facdes
        
C razon entre diametro equiv. de esfera dura y diametro de esfera suave.
!	lam=diametro(n)!1.D0+Eul/n
	fve=fv*lam**3
CC
C correccion de Verlet-Weis
CC
	fvw=fve-fve**2/16.D0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
        K=0.001D0
        dK=kmaxima/dble(Np)
C
        do i=1,Np
CC
C La K que vamos a meter esta en unidades de dw, portanto la k correcta
C de esfera suave es como sigue
CC
           Sdk(i,1)=K!*((fve/fvw)**(1.D0/3.D0)/lam)
CC
C correccion de Verlet-Weis
CC
           Sdk(i,2)=facdes(K,fvw)
           K=K+dK
C
        enddo

C definición del espaciado para una buena g(r)
!	goto 100
        dK=500.D0/Dble(Np)
        K=0.001D0
        
C Se usa el arreglo Sdkg para guardar información a k's
C grandes y poder hacer una buena evaluación de g(r)
         do i=1,Npg*Np
           Sdkg(i,1)=K
           Sdkg(i,2)=facdes(K,fvw)
           K=K+dK
         enddo
C ya tenemos S(k) en hasta una k grande ahora le sacamos la
C transformada inversa de Fourier para obtener g(r)
         call  TIF3D(Npg*Np,Sdkg,Np,Gr,fvw,rmaxima,gc)
C ya se calcularon S(k) y g(r)
C Calcualmos Y(r)
CC
C Correccion de Verlet-Weis
CC
	Anum=3.D0*(fvw**2*(1.D0-0.7117D0*fvw-0.114D0*fvw**2))
	Aden=4.D0*(1.D0-fvw)**4
	A=Anum/Aden
	mu=24.D0*A/(fvw*gc)
CC
        Do i=1,Np
CC
C la g(r) que salio de la transformada esta en unidades de dw
C debemos transformarla a las unidades de d
CC
	    Yr(i,1)=Gr(i,1)*(fvw/fve)**(1.D0/3.D0)
	    r=Yr(i,1)
            if(Yr(i,1).le.1.D0) then
             C1=(1.D0+2.D0*fvw)**2-6.D0*fvw*(1.D0+0.5D0*fvw)**2*r
             C2=fve*(1.D0+2.D0*fvw)**2*0.5D0*r**3           
             Yr(i,2)=(C1+C2)/(1.D0-fvw)**4
           else
             Yr(i,2)=Gr(i,2)
     #               +(A/r)*dexp(-mu*(r-1.D0))*dcos(mu*(r-1.D0))
           endif
CC
C Salida final de la g(r) de esfera suave
CC
	   Gr(i,1)=Yr(i,1)*lam
	   pot=0.D0
	   if(Gr(i,1).le.1.D0)then
	    pot=(1.D0/Gr(i,1)**n)**2-2.D0*(1.D0/Gr(i,1)**n)+1.D0
	   endif
	   Gr(i,2)=dexp(-pot)*Yr(i,2)
!	   Gr(i,1)=Yr(i,1)/(fvw/fve)**(1.D0/3.D0)
	
        Enddo
! 100	continue
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
        

      Subroutine TIF3D(NK,SDK,N,GDR,fv,rmx,gc)
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
	 if((GDR(J,1).GE.1.D0).and.(GDR(J,1).LT.1.D0+DR))gc=GDR(J,2)
         

      ENDDO

      RETURN
      END

C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON
C     ****************************************************

      SUBROUTINE SIMPSON(XIN,DR,M,RES1)

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

      Double precision function diametro(nu)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C Este programa sirve para calcular el diametro de esfera dura equi-C
C valente mediante la solución numérica de la ecuación 6.3.11 del   C
C libro "Theory of simple liquids", de Hansen y Macdonald.          C
C								    C
C Enero 2007							    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
CC
C Variables para hacer la integral
CC
	Integer Nr
	Parameter (Nr=2**11)
	Double precision Iu(Nr),d
CC
C variables para el potencial en particular que estemos trabajando
CC
	Double precision nu
	Double precision R,pot
CC
C variables auxiliares
CC
	Integer i
	Double precision dr

	dr=1.D-3
	Do i=1,Nr
	 R=dble(i)*dr
	 pot=0.D0
CC
C Definicion del potencial
CC
	 If (R.le.1.D0)then
	  pot= (1.D0/R**nu)**2-2.D0*(1.D0/R**nu)+1.D0
	 endif
CC
C integrando para calcular d/sigma
CC
	 Iu(i)=1.D0-dexp(-pot)
!	 write(50,*)R,Iu(i)
	enddo
! 	stop
CC
C integracion por simpson
CC
	call SIMPSON(Iu,dr,Nr,d)
CC
!	write(*,*)'el diametro de esfera dura equivalente es:',d
	diametro=d
	Return
      END
CC
 
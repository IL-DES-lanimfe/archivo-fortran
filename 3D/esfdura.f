      Program estructuraesfdura

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                C
C Este programa sirve para probar una subrutina para calcular el C
C el factor de estructura,la funcion de distribucion radial para C
C un potencial de esfera dura. Usalos resultados analiticos de   C
C M.S. Wertheim (PRL vol.10,no.8,pp.321,1963) para la solucion deC
C la ec. de Ornstein-Zernik con una cerradura de Percus-Yevick.  C
C								 C
C Octubre 2006                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       Implicit none !aguas toda se tiene que declarar
       integer Nk,i
       parameter (Nk=2**12)
       double precision Sdek(Nk,2),Gder(Nk,2)
       double precision fv !,kmax,rmax
       character*1 dec
       write(*,*)'fraccion de volumen?'
       read(*,*)fv
!        write(*,*)'¿valor max. de k para el fac. de estructura?'
!        read(*,*)kmax
!        write(*,*)'¿valor max. de r para la fun. de dist. radial?'
!        read(*,*)rmax
!       write(*,*)'quieres calcular g(r) tambien?(s o n)'
!       read(*,*)dec
	 dec='n'
CC
C Correccion de Verlet -Weis
CC
	! fv=(0.491/0.494)*fv
       call Sdeked(Sdek,Gder,fv,Nk,100.D0,10.D0,dec)
       open(10,file="facdesdurafv501.dat")
       do i=1,Nk
        write(10,*)Sdek(i,1),Sdek(i,2)
       enddo
       close(10)
       if(dec.eq.'s')then
        open(11,file='fdisraddurafv74.dat')
        do i=1,Nk
	 write(11,*)Gder(i,1),Gder(i,2)
        enddo
        close(11)
       endif
       stop
      END

      Subroutine  Sdeked(Sdk,Gr,fracvol,Np,kmaxima,rmaxima,dec)
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
C             *fracvol:fraccion de volumen                           C
C             *Np:numero de puntos a evaluar en S(k) y g(r)          C
C             *kmaxima:valor de k maximo para evaluar S(k)           C
C             *rmaxima:valor de r maximo para evaluar g(r)           C
C             *dec: si vale 's' se calcula g(r),'n' para no calcular C
C		    g(r)		                             C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        integer Np,Npg
        double precision fracvol,fv,rmaxima,kmaxima
        double precision Sdk(Np,2),Gr(Np,2)
	character*1 dec
!parametro para extende k para ganar buena g(r)
        parameter(Npg=50)
        integer i
        double precision K,dK
!arreglo auxiliar para buena g(r)
        double precision Sdkg(Npg*Np,2)
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
C por si no se necesita la funcion de distribucion radial
CC
	if (dec.eq.'n') Return
CC
C definición del espaciado para una buena g(r)
CC
        dK=1000.D0/Dble(Np)
        K=0.001D0
CC
C Se usa el arreglo Sdkg para guardar información a k's
C grandes y poder hacer una buena evaluación de g(r)
CC                
         do i=1,Npg*Np
          Sdkg(i,1)=K
          Sdkg(i,2)=facdes(K,fv)
          K=K+dK
         enddo
CC
C ya tenemos S(k) en hasta una k grande ahora le sacamos la
C transformada inversa de Fourier para obtener g(r)
CC
         call  TIF3D(Npg*Np,Sdkg,Np,Gr,fv,rmaxima)
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

 

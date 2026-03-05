      PROGRAM spinoPSM
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C									C
C Programa para calcular la espinodal para el modelo de esferas 	C
C penetrables (PSM por sus siglas en ingles). 				C
C									C
C Agosto de 2008							C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer i,j
	Double Precision T,fv
	Double Precision dy,dfv
	Double Precision Smx,S,Smi,y,SdekPSM

CC
C Archivo de Salida
CC
	 open(20,file='spinodalPSM.dat')
CC
	dy=0.5D-2
	DO T=3.D0,14.D0,1.D0
CC
C parametros iniciales
CC
	 if(T.eq.3.D0)fv=4.286D0
	 if(T.eq.4.D0)fv=5.745D0
	 if(T.eq.5.D0)fv=7.205D0
	 if(T.eq.6.D0)fv=8.657D0
	 if(T.eq.7.D0)fv=10.112D0
	 if(T.eq.8.D0)fv=11.57D0
	 if(T.eq.9.D0)fv=13.019D0
	 if(T.eq.10.D0)fv=14.472D0
	 if(T.eq.11.D0)fv=15.925D0
	 if(T.eq.12.D0)fv=17.378D0
	 if(T.eq.13.D0)fv=18.83D0
	 if(T.eq.14.D0)fv=20.282D0
	 dfv=1.D-4
CC
C ciclo para acercarse a la espinodal
CC
	 DO i=1,10000
CC
C calculando el maximo de S(k)
CC
	  Smx=0.D0
	  do j=1,5000
	   y=dble(j)*dy
	   S=SdekPSM(y,fv,T)
	   if(Smx.LT.S)Smx=S
	   if(S.lt.0.D0)then
	    write(*,*)'ya se paso'
	    fv=fv-dfv
	    dfv=dfv/10.D0
	    goto 1000
	   endif
	  enddo
CC
	  write(*,*)'Calculado el maximo de S(k)',T,fv,Smx
CC
C criterio de convergencia
CC
	  Smi=1.D0/Smx
	  if(Smi.lt.1.D-5)then
	   write(20,*)fv,T
	   write(*,*)'se encontro solucion',fv,T
	   goto 2000
	  endif
CC
C por si nos pasamos
CC
 1000	  continue
CC
C Aumentamos la fraccion de volumen para seguir buscando
CC
	  fv=fv+dfv
	 ENDDO
CC
C Para cuando ya encontro la solucion
CC
 2000	 continue
	ENDDO
	close(20)
	stop
      END
CC
C Factor de Estructura
CC
      DOUBLE PRECISION FUNCTION SdekPSM(k,eta,tem)
	implicit none
	Double precision k,eta,tem
	Double precision Senos
	Senos=(dsin(k)-k*dcos(k))/k**3
	SdekPSM=1.D0/(1.D0+24.D0*eta*Senos/tem)
	return
      END
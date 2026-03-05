      PROGRAM desplazamientocuadratricomedio
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa es para calcular el dezplazamiento cuadratico medio  C
C como funcion del tiempo					     C
C								     C
C Diciembre 2006						     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none !toda variable se tiene que declarar
	Integer i,nt,it,niter,iiter,ic
	Integer met1,met2
	Parameter(nt=3000)
	Parameter(niter=100)
	Double precision T(nt),DELTZ(nt),CO(nt),f(nt),W(nt)
	Double precision dt,conv,fvieja,solu,Aes,Bes,fa,xxx
	Double precision WW(nt),WWI(nt),conv2,solu2,WWvieja
	Double precision D0tau,De,Ddt(nt),tet0,tb
	Double precision convr1,convr2,conv2r1,conv2r2
	parameter(tb=63.39D0)
CC
C Leemos el archivo de entrada de la deltaz(t) es necesario que el 
C pimer punto este cercano a cero,ademas damos la propuesta inicial
C para f(t)
CC
	open(21,file='coefdespcuadb.dat')
	read(21,*)Aes,Bes
	close(21)
	open(20,file='deltazHppppb.dat')
	open (32,file='Ddetsexppppb.dat')
	Do i=1,nt
	 read(20,*)T(i),xxx,DELTZ(i)
	 tet0=T(i)
	 D0tau=Aes/(2.D0*Bes+Aes**2)
	 De=1.D0-Aes*D0tau
	 f(i)=(Aes*D0tau)*dexp(-T(i))
	 CO(i)=0.D0
	 W(i)=(Aes*D0tau)*dexp(-T(i))
!	 DELTZ(i)=Aes*dexp(-De*T(i)/D0tau)
         Ddt(i)=De-
     #   (De-1.D0)*D0tau*(1.D0-dexp(-tet0/D0tau))/tet0
	 write(32,*)T(i),Ddt(i)*T(i),Ddt(i)
	Enddo
	close(20)
	close(32)
	dt=T(2)-T(1)
CC
C Iniciamos el ciclo para correr sobre los tiempos
CC
	open(30,file='Ddetppppb.dat')
 	open(31,file='Ddet2ppppb.dat')
	DO it=1,nt
	 write(*,*)'tiempo:',T(it)
CC
C interruptores para los dos metodos
CC
	 met1=0
	 met2=0
CC
C Iniciamos el ciclo de iteraciones
CC
	 Do iiter=1,niter
	  write(*,*)'iteracion:',iiter
CC
C empezamos con la convolucion
CC
	  do ic=1,it
	   CO(ic)=DELTZ(it-ic)*f(ic)
 	   WWI(ic)=DELTZ(it+1-ic)*WW(ic)
	  enddo
	   call SIMPSON(CO,dt,it,conv)
	   call SIMPSON(WWI,dt,it,conv2)
! 	  if((it.gt.500).and.(it.le.750))then
! 	   call SIMPSON(CO,dt,500,convr1)
! 	   call SIMPSON(WWI,dt,500,conv2r1)
! 	   call SIMPSON(CO,2.D0*dt,it-500,convr2)
! 	   call SIMPSON(WWI,2.D0*dt,it-500,conv2r2)
! 	   conv=convr1+convr2
! 	   conv2=conv2r1+conv2r2
! 	  endif
CC
C fin de la convolucion
CC
	  f(it)=DELTZ(it)-conv
	  WW(it)=T(it)-conv2
CC
C criterios de convergencia
CC
	  if(met1.eq.1)goto 100
	  solu=dabs((f(it)-fvieja)/fvieja)
	  if(solu.LE.1.D-6)then
	   do ic=1,it
	    CO(ic)=(T(it)-T(ic))*f(ic)
	   enddo
	   call SIMPSON(CO,dt,it,conv)
	   W(it)=T(it)-conv
	   write(30,*)real(T(it)),real(tb*T(it)),real(W(it)),
     #         real(W(it)/T(it)),real(f(it))
	   write(*,*)'metodo 1 listo'
	   met1=1
	   goto 100	
	  endif
	  fvieja=f(it)
 100	  continue
CC
	  if(met2.eq.1)goto 200
	  solu2=dabs((WW(it)-WWvieja)/WWvieja)
	  if(solu2.LE.1.D-6)then
	   write(31,*)real(T(it)),real(tb*T(it)),real(WW(it)),
     #                real(WW(it)/T(it))
	   write(*,*)'metodo 2 listo'
	   met2=1
	   goto 200	
	  endif
	  WWvieja=WW(it)
 200	  continue
CC
	  if((met1.eq.1).and.(met2.eq.1))goto 300
	 Enddo
CC
C Fin de iteraciones
CC
 300	 continue
	ENDDO
	close(30)
	close(31)
      END
C     ****************************************************
C     SUBRUTINA DE INTEGRACION POR SIMPSON 
C     ****************************************************
      SUBROUTINE SIMPSON(XINI,DR,N,RES1)
      IMPLICIT REAL*8(A-H,O-Z)
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XINI(N)
	
	RES1=0.0
	S0=0.0
	S1=0.0
	S2=0.0

	DO 100 I=2,N-1,2
	S0=S0+XINI(I-1)
	S1=S1+XINI(I)
	S2=S2+XINI(I+1)
100	CONTINUE
	RES1=DR*(S0+4.0*S1+S2)/3.0
	IF(MOD(N,2).EQ.0) RES1=RES1+DR*(5.0*XINI(N)+8.0*XINI(N-1)-
     #	XINI(N-2))/12.0
	RETURN
	END

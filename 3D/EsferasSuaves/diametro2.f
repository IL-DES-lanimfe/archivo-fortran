	PROGRAM diametros2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Este programa calcula los diametros efectivos de un sistema de esfe-C
C ras suaves fijandose en que coincida el maximo de S(k).             C
C																	C
C Octubre 2007														C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 implicit none !toda variable se tiene que declarar
	 Integer Nk
	 parameter (Nk=2**12)
	 Double precision nu,phi,dphi,Sdek(Nk,2),Gder(Nk,2)
	 Double precision SmxR,Smx,d
	 Double precision diametro,Smax
CC
C primero calculamos la S(k) de referencia
CC
	 nu=50.D0
	 phi=0.52D0
c... Reads inpud data
       call input(0.52D0,50.D0,0.28D0)
c... perform calculations
       call oz2(Sdek,Gder)
CC
C calculamos el maximo
CC
	 SmxR=Smax(Sdek,Nk)
CC
	 open(20,file='diametrosS.dat')
	 Do nu=6.D0,15.D0,1.D0
	  dphi=0.1D0
CC
C valor de nu para el que se quiere calcular
CC
!	 nu=10.D0
CC
C Primera aproximacion al diametro efectivo con la funcion blip
CC
	  d=diametro(nu)
	  phi=0.5D0/d**3
CC
C Ahora calculamos S(k) para esa nu y vemos si la funcion blip la levanta
C
 123	  continue
	  write(*,*)phi,nu
c... Reads inpud data
       call input(phi,nu,0.5D0)
c... perform calculations
       call oz2(Sdek,Gder)
CC
C Calculamos maximo
CC
	  Smx=Smax(Sdek,Nk)
CC
C ahora vemos si son equivalentes
CC
	  if(dabs(SmxR-Smx).le.1.D-3)then
	   d=(0.5D0/phi)**(1.D0/3.D0)
	   write(20,*)nu,d,phi
	   goto 234
	  endif
	  if(Smx.gt.SmxR)then
	   phi=phi-dphi
	   dphi=dphi/10.D0
	  endif
	    
CC
C Si no ha levantado
CC
	  phi=phi+dphi
	  goto 123
 234    continue
       enddo
	 close(20)
	 stop
	END

      Double precision function diametro(nu)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C																	C
C Este programa sirve para calcular el diametro de esfera dura equi-  C
C valente mediante la solución numérica de la ecuación 6.3.11 del  C
C libro "Theory of simple liquids", de Hansen y Macdonald.            C
C																	C
C Enero 2007															C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine input(fv,xnu,alpha)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      implicit integer*4(i-n), real*8(a-h,o-z)
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
c
      pi=4.d0*datan(1.d0)
c
c--- rho=densidad numerica
c--- x(i)=fraccion molar de la especie i
c--- closure=tipo de cerradura (2=HNC,3=RY)     
c
      fi = fv
      rho = (6.0/pi)*fi
!      rho = 0.5
c      
      x(1) = 1.0
      rmax = 65.0  ! rmax=800 para yukawas
		   !  rmax = 10-50 para corto alcance
c
      alfa = alpha
      closure = 3.0
      nr   = nm
      nq   = nr
      nrho = 100
      ez   = 1.d-4
c
      sigma1 = 1.0
      sigma2 = 1.0
      sigma2 = sigma2 / sigma1
      sigma1 = sigma1 / sigma1
      sig(1) = sigma1
      sig(2) = (sigma1 + sigma2 )/2.0
      sig(3) = sigma2
      do i = 1, 3
      sigg(i) = sig(i)/ 2.0
      enddo
c
      SIG3B=X(1)*SIGMA1**3+X(2)*SIGMA2**3
      x(2)=1.d0-x(1)
c
      Dq = PI/RMAX
      DR = RMAX/(1.D0*NR)
      DO I = 1, Nr
      R(I) = (I-1)*DR
      q(I) = (I-1)*Dq
      enddo
      qmax = q(nr)
c
      CALL POT(nr,xnu)
c

      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE POT(n,ex)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER(nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
c
c
      ipote = 3
c
      dmed = rho**(-1.0/3.0)
      write(*,*)rho,dmed
c
c! potencial $u(r)=\epsilon(1/r)^{\nu}$
c
      if (ipote.eq.1) then
      rlamb = 10
      E(1) = 1.0
      E(3) = 1.0
      E(2) = SQRT(e(1)*e(3))
      z(1) = rlamb
      z(2) = rlamb
      z(3) = rlamb
      DO K=1,3
      DO I=1,N
      if(r(i).lt.sigg(k))then
      U(i,k)  = 0.0
      UP(I,K) = 0.0
      else
      U(I,K) = E(K)*(sig(k)/r(i))**(rlamb)
      UP(I,K) = U(i,k)*(rlamb)          ! = -f(r)*r
c      UP(I,K) = U(i,k)*(rlamb/r(i))    ! = -f(r)
      endif
      ENDDO
      ENDDO
      endif
c
c! potencial $u(r)=4\epsilon_{ij}[(\sigma_{ij}/r)^{2\nu}$-(\sigma_{ij}/r)^{\nu}+1/4$
c! $\sigma_{ij}=(\sigma_i+\sigma_j)/2$, $\epsilon_{ij}=(\epsilon_i\epsilon_j)^{1/2}$
c
      if (ipote.eq.2) then
      RLAMB = 6.0
      E(1)=1.0
      E(3)=1.0
      E(2)=SQRT(e(1)*e(3))
      Z(1)=RLAMB
      Z(2)=RLAMB
      Z(3)=RLAMB
c
      DO K=1,3
      DO I=1,N
      arg4 = sig(k)*2.0**(1.0/rlamb)
      if((r(i).lt.sigg(k)).or.(r(i).gt.arg4))then
      U(i,k)  = 0.0
      UP(i,k) = 0.0
      else
      arg1 = (sig(k)/r(i))**(rlamb)
      arg2 = arg1 * arg1
      arg3 = 0.25
c      arg3 = 0.0
      U(I,K)  = (4.0*E(K)) * (arg2 - arg1 + arg3)
      UP(I,K) = 4.0 * (2.0*arg2 - arg1) * (rlamb)         ! = -f(r)*r
c      UP(I,K) = (2*arg2-arg1)*(rlamb/r(i))   ! = -f(r)
      endif
      ENDDO
      ENDDO
      endif
c
c! potencial $u(r)=4\epsilon_{ij}[(\sigma_{ij}/r)^{2\nu}$-(\sigma_{ij}/r)^{\nu}+1/4$
c! $\sigma_{ij}=(\sigma_i+\sigma_j)/2$, $\epsilon_{ij}=(\epsilon_i\epsilon_j)^{1/2}$
c
      if (ipote.eq.3) then
      RLAMB = ex
      E(1)=1.0
      E(3)=1.0
      E(2)=SQRT(e(1)*e(3))
      Z(1)=RLAMB
      Z(2)=RLAMB
      Z(3)=RLAMB
c
      DO K=1,3
      DO I=1,N
      if((r(i).lt.sigg(k)).or.(r(i).gt.sig(k)))then
      U(i,k)  = 0.0
      UP(i,k) = 0.0
      else
      arg1 = (sig(k)/r(i))**(rlamb)
      arg2 = arg1 * arg1
      arg3 = 1.0
c      arg3 = 0.0
      U(I,K)  = E(K) * (arg2 - 2*arg1 + arg3)
      UP(I,K) = (2.0*arg2 - 2*arg1) * (rlamb)         ! = -f(r)*r
c      UP(I,K) = (2*arg2-arg1)*(rlamb/r(i))   ! = -f(r)
      endif
      ENDDO
      ENDDO
      endif
c
C! potencial $u(r)=\epsilon\exp[-\kappa(r-\sigma)]/r$
c
      if (ipote.eq.4) then
      rlamb = 0.15
      E(1) = 500.0
      E(3) = 800.0
      E(2) = SQRT(e(1)*e(3))
      z(1) = rlamb
      z(2) = rlamb
      z(3) = rlamb
      DO K=1,3
      DO I=1,N
      if(r(i).lt.sig(k))then
      U(i,k)  = 0.0
      UP(i,k) = 0.0
      else
      U(I,K)  = E(K)*EXP(-Z(K)*(R(I)-1.0))/R(I)
      UP(I,K) = U(i,k)*(1.0+rlamb*r(i))         ! = -f(r)*r
c      UP(I,K) = U(i,k)*(1.0+rlamb*r(i))/r(i)   ! = -f(r)
      endif
      ENDDO
      ENDDO
      endif
c
c... DATOS A CAMBIAR ---- estar=gamma* --- zstar=kappa* 
c
      if (ipote.eq.5) then
      Zstar = 8.0
      Astar = 3.84
c
      Zyuk = zstar/dmed
      Ayuk = Astar*exp(-zyuk)*dmed*exp(zstar)
c
      WRITE(*,*)'dmed=',dmed
      WRITE(*,*)'Zyuk,Z*',zyuk,zstar
      WRITE(*,*)'Ayuk,Astar',Ayuk,Astar
c      PAUSE
c
      RLAMB = Zyuk
      E(1)=(SQRT(Ayuk))
      E(3)=0.0
      E(2)=e(1)*e(3)
      e(1)=e(1)**2
      e(3)=e(3)**2
      Z(1)=RLAMB
      Z(2)=RLAMB
      Z(3)=RLAMB
c
      DO K=1,3
      DO I=1,N
      if(r(i).lt.sig(k))then
c      if(r(i).lt.sig1)then
      U(i,k)  = 0.0
      UP(i,k) = 0.0
      write(*,*)'sigk',sig(k)
      else
      U(I,K)  = E(K)*EXP(-Z(K)*(R(I)-1.0))/R(I)
      UP(I,K) = U(i,k)*(1.0+rlamb*r(i))          ! = -f(r)*r
c      UP(I,K) = U(i,k)*(1.0+rlamb*r(i))/r(i)    ! = -f(r)
      endif
      ENDDO
      ENDDO
      endif
c
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine oz2(Sk,Gr)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION g(nm,3),c(nm,3),GAM(NM,3),GAM0(NM,3)
      DIMENSION GAM1(NM,3),R1(NM),DIF(2),rk(nm),c1(nm,3)
      DIMENSION GAMAUX(NM,3)
      dimension gh(nm,3),gh2(nm,3),ck(nm,3),s(nm,3)
	dimension Sk(nm,2),Gr(nm,2)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
c
      tfin = 0.0
      pi=4.d0*datan(1.d0)
c
      do k=1,3
      do i=1,nm
      GAM0(I,K)=0.D0
      enddo
      enddo
c -----------------------------------Anfang.
      rhoa=rho
      DT=1.0/(1.D0*NRHO)
      DRHO=RHO/(1.d0*nrho)
c
      KJ=1
      RHO=KJ*DRHO
      T=DT*KJ
c
12    CALL Ng(NR,GAM0,GAM,C)
c
      IF (KJ.GT.1) GOTO 20
      DO 15 K=1,3
         DO 15 I=1,Nr
15        GAM0(I,K)=GAM(I,K)
      KJ=KJ+1
      RHO=KJ*DRHO
      T=DT*KJ
      GOTO 12
20    CONTINUE
      CALL EXTRAP(GAM0,GAM,Nr,RHO,DRHO)
22    DO 25 K=1,3
         DO 25 I=1,Nr
25       GAM1(I,K)=GAM(I,K)
      KJ=KJ+1
      RHO=KJ*DRHO
      T=DT*KJ
      CALL Ng(NR,GAM0,GAM,C)
      IF (KJ.EQ.(NRHO)) GOTO 30
      CALL EXTRAP(GAM1,GAM,Nr,RHO,DRHO)
      DO 27 K=1,3
         DO 27 I=1,NR
27          GAM0(I,K)=GAM1(I,K)
      GOTO 22
 30   continue
c ----------------------------------- PY o HNC
      T=1.D0
      RHO=RHOA
      CALL Ng(NR,GAM0,GAM,C)
      write(*,*)'Salida HNC:'
c ----------------------------------- RY
      if (closure.lt.2.5) goto 1000
c ----------------------------------- RY
      ddrho=rhoa/100.d0
      dalfa=-alfa/50.d0
500   T=1.D0
      RHO=RHOA
      CALL Ng(NR,GAM0,GAM,C)
      do 100 k=1,3
      do 100 i=1,nr
      GAMAUX(I,K)=GAM(I,K)
100   GAM0(I,K)=GAM(I,K)      
      CALL TERMO(NR,GAM,C,PV0,CHIC0,ENER0)
      CHIC = CHIC0
c------
      RHO=RHO-DDRHO
      CALL Ng(NR,GAM0,GAM,C)
      CALL TERMO(NR,GAM,C,PV1,CHIC1,ENER1)
c------
      RHO=RHO+2.D0*DDRHO
      CALL Ng(NR,GAM0,GAM,C)
      CALL TERMO(NR,GAM,C,PV2,CHIC2,ENER2)
      CALL RY(PV1,PV2,CHIC,DDRHO,ALFA,DALFA,IRY)
      IF(IRY.EQ.1) GOTO 500    
c ------------------------
      write(*,*)'Salida RY:'
1000  T=1.d0
      TFIN = 1.0
      write(*,*)'calculo final'
      RHO=RHOA
      CALL Ng(NR,GAM0,GAM,C)
c --------------------------
      CALL ESCRIBE(NR,R,GAM,C,Sk,Gr)
      CALL TERMO(NR,GAM,C,PV,CHIC,ENER)
      pexv=pv/rho-1.d0
c
c      write(*,*)'  '
c      WRITE(*,*)'sqmax=',sqmax
c      write(*,*)'  '
      WRITE(*,*)'RHO=',SNGL(rho)
      WRITE(*,*)'CHIC=',SNGL(CHIC)
      WRITE(*,*)'PEXV=',SNGL(PexV)
      WRITE(*,*)'EEXE=',SNGL(ENER)
      return
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE Termo(NR,GAM,C,PV1,CHIC,ENER)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION GAM(NM,3),C(NM,3),R1(NM),G(NM,3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq

c
      pi=4.d0*datan(1.d0)
c      write(*,*)'dr',dr,rho,rlamb,ipote
c      dr = RMAX/(1.D0*NR)
c
      do i = 1, nr
         r1(i)=x(1)**2*c(i,1)+2.d0*x(1)*x(2)*c(i,2)+x(2)**2*c(i,3)
         R1(I)=R1(I)*R(I)**2
      enddo
c
      chic=0.d0
      do i = 2, nr - 1
      chic = chic + r1(i)
      enddo
      chic = dr*(chic+(r1(1)+r1(nr))/2.d0)
      chic = 1.d0-4.d0*pi*rho*chic
c
      do k = 1, 3
      do i = 1, nr
      g(i,k)=gam(i,k)+c(i,k)+1.d0
      enddo
      enddo
c
c      do i = 1, nr
c      r1(i)=x(1)**2*g(i,1)*u(i,1)+2.0*x(1)*x(2)*g(i,2)*u(i,2)
c      r1(i)=(r1(i)+x(2)**2*g(i,3)*u(i,3))*(rlamb)*r(i)**2
c      enddo
      do i = 1, nr
      ru1 = x(1)**2*g(i,1)*up(i,1)
      ru2 = 2.0*x(1)*x(2)*g(i,2)*up(i,2)
      ru3 = x(2)**2*g(i,3)*up(i,3)
      r1(i)= (ru1+ru2+ru3)*r(i)**2
c      ru1 = x(1)**2*g(i,1)*u(i,1)
c      ru2 = 2.0*x(1)*x(2)*g(i,2)*u(i,2)
c      ru3 = x(2)**2*g(i,3)*u(i,3)
c      r1(i)= (ru1+ru2+ru3)*(1.0+rlamb*r(i))*r(i)**2
      enddo
c
      pv1=0.d0
      do  i = 2, nr - 1
      pv1 = pv1 + r1(i)
c      write(*,*)'pv1',pv1
      enddo
      pv1 = dr*(pv1+(r1(1)+r1(nr))/2.0)
c      write(*,*)'pv1',pv1,u(1,1),u(1,2),u(1,3)
      pv1 = rho*(1.0+2.0*pi*rho*pv1/3.0)
c      write(*,*)'pv1',pv1
c
      do i=1,nr
      r1(i)=x(1)**2*g(i,1)*u(i,1)+2.0*x(1)*x(2)*g(i,2)*u(i,2)
      r1(i)=(r1(i)+x(2)**2*g(i,3)*u(i,3))*r(i)**2
      enddo
      ENER=0.D0
      DO I=2,NR-1
      ENER=ener+r1(i)
      enddo
      ener=dr*(ener+(r1(1)+r1(nr))/2.0)
      ENER=2.D0*PI*RHO*ENER
      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE RY(PV1,PV2,CHIC,DDRHO,ALFA,DALFA,IRY)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      dimension DIF(2)
      data ix /1/
      save ix,dif
c
      IRY = 0
c      ddrho=rhoa/100.d0
      CHIV=(PV2-PV1)/(2.D0*DDRHO)
      DIF(IX) = CHIC - CHIV
      WRITE(*,*)'CHIC=',SNGL(CHIC),'  CHIV=',SNGL(CHIV)
      WRITE(*,*)'ALFA=',sngl(ALFA),'  DIF=',sngl(DIF(IX))
      WRITE(*,*)' '
      IF (IX.LT.2) THEN
      IX=IX+1
      ALFA=ALFA+DALFA
      IRY = 1
      RETURN
      ENDIF
      PROD=DIF(1)*DIF(2)
      IF (PROD.LT.0.D0) THEN
      ALFA=ALFA-DALFA
      A=(DIF(2)-DIF(1))/DALFA
      B=DIF(1)-A*ALFA
      ALFA=-B/A
      WRITE(*,*)'ALFA=',sngl(ALFA)
      IRY = 0
      RETURN
      ELSE
      DIF(1)=DIF(2)
      ALFA=ALFA+DALFA
      IRY = 1
      RETURN
      ENDIF
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ESCRIBE(NR,RR,GAM,C,Sk,Gr)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION rr(nm),g(nm,3),c(nm,3),GAM(NM,3)
      DIMENSION GAM1(NM,3),R1(NM),rk(nm),c1(nm,3)
      dimension gh(nm,3),gh2(nm,3),ck(nm,3),s(nm,3)
	dimension Sk(nm,2),Gr(nm,2)
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
c
      do k=1,3
      G(1,K)=0.D0
      DO I=2,NR
      g(i,k)=gam(i,k)+c(i,k)+1.d0
      enddo
      enddo
c --- salen las funciones de correlacion gij(r)
      open(10,file='Gr.out',status='unknown')
      open(20,file='Cr.out',status='unknown')
         do i=1,nr
         do k=1,3
         gh(i,k)=gam(i,k)+c(i,k)
         enddo
         write(10,402)rr(i),gh(i,1)+1.0,gh(i,2)+1.0,gh(i,3)+1.0
         write(20,402)rr(i),c(i,1),c(i,2),c(i,3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  Gr(i,1)=rr(i)
	  Gr(i,2)=gh(i,1)+1.0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         enddo
       close(10)
       close(20)
c --- --- --- --- --- --- --- ---

c === calculo de factores de estructura ===
      nk=nr
      rkmax=qmax/2.0
      DK=RKMAX/(1.D0*NK)
      DO I=1,NK
      RK(I)=I*DK
      ENDDO
      CALL FT(GH,GH2,RR,RK,Nr,Nk,DR)
      CALL FT(C,C1,RR,RK,Nr,Nk,DR)

c === salen las cij(k) y la ck(k) para el cilindro ==========
      open(10,file='hkij.out',status='unknown')
      open(11,file='ckij.out',status='unknown')
      open(12,file='ck.out',status='unknown')
       do i=1,nk
       write(10,402)rk(i),gh2(i,1),gh2(i,2),gh2(i,3)
       write(11,402)rk(i),c1(i,1),c1(i,2),c1(i,3)
       write(12,*)rk(i),c1(i,1)
       enddo
 402   format(4(f15.5))
      close(10)
      close(11)
      close(12)
c === === === === ===
       DO I=1,Nr
       DO K=1,3
       CK(I,K)=C1(i,k)
       ENDDO
       ENDDO
       DO I=1,Nr
       DELTA=(1.D0-RHO*x(1)*Ck(I,1))*(1.D0-RHO*x(2)*Ck(I,3))
       delta=delta-rho**2*x(1)*x(2)*ck(i,2)**2
       s(I,1)=(1.d0-rho*x(2)*ck(i,3))*ck(i,1)
       s(I,1)=(s(i,1)+rho*x(2)*ck(i,2)**2)/DELTA
       s(I,2)=Ck(I,2)/DELTA
       s(I,3)=(1.d0-RHO*x(1)*Ck(I,1))*ck(i,3)
       s(I,3)=(s(i,3)+RHO*x(1)*Ck(I,2)**2)/DELTA
       s(i,1)=x(1)+rho*x(1)**2*s(i,1)
       s(i,2)=rho*x(1)*x(2)*s(i,2)
       s(i,3)=x(2)+rho*x(2)**2*s(i,3)
       enddo

      sqmax = 0.0
      write(*,*)'salida para Sij(k): Sqij.out'
      write(*,*)'orden: S11, S12, S22'
      open(20,file='Sq.out',status='unknown')
      do i=1,nr
      write(20,66)rk(i),s(i,1),s(i,2),s(i,3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	Sk(i,1)=rk(i)
	Sk(i,2)=s(i,1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(s(i,1).gt.sqmax) sqmax = s(i,1)
      enddo
      close(20)
      write(*,*)'  '
      WRITE(*,*)'sqmax=',sqmax
      write(*,*)'  '
 66   format(4(f16.12,5x)) 
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE Ng(N,GAM0,GAM,C)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION GAM0(NM,3),GAM(NM,3),C(NM,3)
      dimension g1(nm,3),d1(nm,3),g2(nm,3),d2(nm,3),g3(nm,3)
      dimension d3(nm,3),d01(nm,3),d02(nm,3),d01d01(3),d01d02(3)
      dimension d3d01(3),d02d02(3),d3d02(3),const1(3),const2(3)
      dimension f(nm,3)
      dimension fng(nm,3,3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
c
      pi=4.d0*datan(1.d0)
c
      do k=1,3
      do i=1,n
      F(I,K)=GAM0(I,K)
      enddo
      enddo
c
      call ONg(n,f,g1,c)
c
      do k=1,3
      do i=1,n
      d1(i,k)=g1(i,k)-f(i,k)
      enddo
      enddo
c
      call ONg(n,g1,g2,c)
c
      do k=1,3
      do i=1,n
      d2(i,k)=g2(i,k)-g1(i,k)
      enddo
      enddo
c
300   continue
c
      call ONg(n,g2,g3,c)
c
      IF (KJ.LT.2) GOTO 200
c
      do k=1,3
      do i=1,n
        d3(i,k)=g3(i,k)-g2(i,k)
      enddo
      enddo
c
      do k=1,3
      do i=1,n
         D01(I,K)=D3(I,K)-D2(I,K)
         d02(i,k)=d3(i,k)-d1(i,k)
      enddo
      enddo
c
      call pp(n,dr,d01,d01,d01d01)
      call pp(n,dr,d01,d02,d01d02)
      call pp(n,dr,d3,d01,d3d01)
      call pp(n,dr,d02,d02,d02d02)
      call pp(n,dr,d3,d02,d3d02)
c
      V=1.D-50
      DO 50 K=1,3
         IF (D01D01(K).LE.V.OR.D01D02(K).LE.V) GOTO 100
         IF ((D02D02(K)-D01D02(K)**2/D01D01(K)).LE.V) THEN
            GOTO 100
         ELSE
            CONST2(K)=D3D02(K)-D01D02(K)*D3D01(K)/D01D01(K)
            CONST2(K)=CONST2(K)/(D02D02(K)-D01D02(K)**2/D01D01(K))
            CONST1(K)=(D3D02(K)-D02D02(K)*CONST2(K))/D01D02(K)
            DO 55 I=1,N
               F(I,K)=(1.D0-CONST1(K)-CONST2(K))*G3(I,K)
55             F(I,K)=F(I,K)+CONST1(K)*G2(I,K)+CONST2(K)*G1(I,K)
         ENDIF
50    CONTINUE
c
      call ONg(n,f,g3,c)
c
      do 60 k=1,3
         do 60 i=1,n
60          d3(i,k)=g3(i,k)-f(i,k)
c
100   CALL PRES(D3,DR,N,ETA)
c      WRITE(*,*)sngl(RHO),sngl(EZ),sngl(ETA),N
      IF (ETA.LE.Ez) goto 200
      do 70 k=1,3
         do 70 i=1,n
            G1(I,K)=G2(I,K)
            D1(I,K)=D2(I,K)
            G2(I,K)=G3(I,K)
70          d2(i,k)=d3(i,k)
      GOTO 300
200   do 80 k=1,3
         do 80 i=1,n
80          GAM(I,K)=G3(I,K)
      WRITE(*,*)'ITERACION',kj,sngl(ETA),N
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine pp(n,dr,f,g,prod)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      implicit integer*4(i-n), real*8(a-h,o-z)
      dimension f(nm,3),g(nm,3),prod(3)
      do 10 k=1,3
         prod(k)=0.d0
         do 20 i=2,n-1
20       prod(k)=prod(k)+f(i,k)*g(i,k)
10       prod(k)=(prod(k)+(f(1,k)*g(1,k)+f(n,k)*g(n,k))/2.d0)*dr
      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine ONg(n,gamma0,gamma1,c1)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION GAMMA0(NM,3),GAMMA1(NM,3),C1(NM,3)
      DIMENSION S(NM,3),CK(NM,3)
      DIMENSION A0(3),a1(3),a2(3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO

      call closrel(n,gamma0,c1)
      CALL FFTM(C1,RMAX,N,1)
         DO I=1,N
         DO K=1,3
         CK(I,K)=C1(i,k)
         ENDDO
         ENDDO
      do 10 i=1,n
         DELTA=(1.D0-RHO*X(1)*C1(I,1))*(1.D0-RHO*X(2)*C1(I,3))
         DELTA=DELTA-RHO**2*X(1)*X(2)*C1(I,2)**2
         GAMMA1(I,1)=(1.D0-RHO*X(2)*C1(I,3))*C1(I,1)
         GAMMA1(I,1)=(GAMMA1(I,1)+RHO*X(2)*C1(I,2)**2)/DELTA
         GAMMA1(I,1)=GAMMA1(I,1)-C1(I,1)
         GAMMA1(I,2)=C1(I,2)/DELTA-C1(I,2)
         GAMMA1(I,3)=(1.D0-RHO*X(1)*C1(I,1))*C1(I,3)
         GAMMA1(I,3)=(GAMMA1(I,3)+RHO*X(1)*C1(I,2)**2)/DELTA
10       GAMMA1(I,3)=GAMMA1(I,3)-C1(I,3)
      call fftm(gamma1,rmax,n,2)
      call closrel(n,gamma1,c1)
c
c
c
       IF(TFIN.LT.1.0) RETURN
       DO I=1,N
       DELTA=(1.D0-RHO*x(1)*Ck(I,1))*(1.D0-RHO*x(2)*Ck(I,3))
       delta=delta-rho**2*x(1)*x(2)*ck(i,2)**2
       s(I,1)=(1.d0-rho*x(2)*ck(i,3))*ck(i,1)
       s(I,1)=(s(i,1)+rho*x(2)*ck(i,2)**2)/DELTA
       s(I,2)=Ck(I,2)/DELTA
       s(I,3)=(1.d0-RHO*x(1)*Ck(I,1))*ck(i,3)
       s(I,3)=(s(i,3)+RHO*x(1)*Ck(I,2)**2)/DELTA
       s(i,1)=x(1)+rho*x(1)**2*s(i,1)
       s(i,2)=rho*x(1)*x(2)*s(i,2)
       s(i,3)=x(2)+rho*x(2)**2*s(i,3)
       enddo
      sqmax = 0.0
      ich = 20
      open(ich,file='Sq2.out',status='unknown')
      do i=1,n
      write(ich,66)q(i),s(i,1),s(i,2),s(i,3)
      if(s(i,1).gt.sqmax) sqmax = s(i,1)
      enddo
      write(*,*)'sqmax==',sqmax
      close(ich)
 66   format(4(f16.12,5x))
      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine closrel(n,gamma,c)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION GAMMA(NM,3),C(NM,3),sigaux(3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
C ----------------------------- PY 
      do k = 1, 3
      if(ipote.eq.1.or.ipote.eq.2.or.ipote.eq.3) then
      sigaux(k) = sigg(k)
      else
      sigaux(k) = sig(k)
      endif
      enddo
      do k = 1, 3
      do i = 1, n
      C(I,k) = GAMMA(I,k) + 1.D0
      enddo
      enddo
c...PY
      if (closure.lt.1.5) then
      do k=1,3
      DO I=1,N
         if (r(i).lt.sigaux(k)) then
         c(i,k) = -c(i,k)
         else
           arg=u(i,k)*t
           if(arg.gt.70.0)then
           c(i,k) = -c(i,k)
           else
           c(i,k)=-(dexp(-arg)-1.0)*c(i,k)
           endif
         endif
       enddo
       enddo
       return
       endif
c...HNC
      if(closure.gt.1.5 .and. closure.lt.2.5)then        
      do k=1,3
      DO I=1,N
         if (r(i).lt.sigaux(k)) then
         c(i,k) = -c(i,k)
         else
           arg = u(i,k)*t - gamma(i,k)
           if(arg.gt.70.0)then
           c(i,k) = -c(i,k)
           else
           c(i,k) = dexp(-arg) - c(i,k)
           endif
         endif
       enddo
       enddo
       return
       endif
C -------------------------------- RY
      do k=1,3
      do i=1,n
      if (r(i).lt.sigaux(k)) then
      c(i,k) = -c(i,k)
      else
      f=1.d0-dexp(-alfa*r(i))
      arg1 = (dexp(gamma(i,k)*f)-1.d0)/f
      c(i,k)=dexp(-u(i,k)*T)*(1.d0+arg1)-c(i,k)
      endif
      enddo
      enddo
      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE EXTRAP(GAM0,GAM,M,RHO,DRHO)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION GAM0(NM,3),GAM(NM,3)
      DO 5 K=1,3
         DO 5 I=1,M
            A=(GAM(I,K)-GAM0(I,K))/DRHO
            B=GAM(I,K)-A*RHO
5           GAM0(I,K)=A*(RHO+DRHO)+B
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE INTERP(M,N,R1,GAM,GAM0,R)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION R1(NM),GAM(NM,3),GAM0(NM,3)
      DIMENSION R(NM),A(NM),B(NM)
      DO 90 K=1,3
         DO 10 I=1,M-1
            A(I)=(GAM(I+1,K)-GAM(I,K))/(R1(I+1)-R1(I))
10          B(I)=GAM(I+1,K)-A(I)*R1(I+1)
         DO 20 J=1,M-1
            DO 20 I=(J-1)*N/M+1,J*N/M
20             GAM0(I,K)=A(J)*R(I)+B(J)
         DO 40 I=N-N/M+1,N
40          GAM0(I,K)=GAM0(N-N/M,K)
90    CONTINUE
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PRES(F,DR,N,ETA)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      REAL*8 F,DR,ETA,AUX
      DIMENSION F(NM,3)
      ETA=0.D0
      DO 20 K=1,3
         AUX=0.D0
         DO 10 I=2,N-1
10          AUX=AUX+F(I,K)**2
         AUX=(AUX+(F(1,K)**2+F(N,K)**2)/2.D0)*DR
20       ETA =ETA+AUX
      ETA=DSQRT(ETA)
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE FFTM(FA,RMAX,N,II)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N), REAL*8(A-H,O-Z)
      DIMENSION FA(NM,3),F(NM)
      DO 10 I=1,3
         DO 20 K=1,N
20          F(K)=FA(K,I)
         CALL FFT(F,RMAX,N,II)
         DO 30 K=1,N
30          FA(K,I)=F(K)
10    CONTINUE
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE CALINT(F,DR,N,RINT)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION F(NM)
      RINT=0.D0
      DO 10 I=2,N-1
10       RINT=RINT+F(I)
      RINT=DR*(RINT+(F(1)+F(N))/2.D0)
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE INTT(H,DR,N,SFT)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N), REAL*8(A-H,O-Z)
      DIMENSION H(NM,3),SFT(3),F(NM)
      DO 10 K=1,3
         DO 20 I=1,N
20          F(I)=H(I,K)
         CALL CALINT(F,DR,N,RINT)
10       SFT(K)=RINT
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE FT(C,C1,R,RK,N,NK,DR)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
      IMPLICIT INTEGER*4(I-N), REAL*8(A-H,O-Z)
      DIMENSION R(NM),RK(NM),C(NM,3),C1(NM,3),SFT(3),CA(NM,3)
      PI=4.D0*DATAN(1.D0)
      DO 10 I=1,NK
      DO 20 K=1,3
      DO 20 J=1,N
20    CA(J,K)=R(J)*C(J,K)*DSIN(RK(I)*R(J))

      CALL INTT(CA,DR,N,SFT)
      DO 30 K=1,3
30    C1(I,K)=4.D0*PI*SFT(K)/RK(I)
10    CONTINUE
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE FFT(F,RMAX,N,II)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (nm=2**12)
C SI II=1 TF, SI II=2 TFI.
      REAL*8 F,A,RMAX
      DIMENSION F(NM)
      DO 10 I=1,N
         F(I)=(I-1)*F(I)
10    CONTINUE
      CALL SINFT(F,N,NM)
      IF (II.EQ.1) GOTO 20
      A=N*(1.D0/2.D0/RMAX**3)
      GOTO 30
20    A=4*RMAX**3/(1.D0*N**2)
      F(1)=0.D0
30    DO 40 I=2,N
      F(I)=A*F(I)/(1.D0*(I-1))
40    CONTINUE
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE SINFT(Y,N,NMAX)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,Y
      DIMENSION Y(NMAX)
      PI=4*DATAN(1.D0)
      THETA=PI/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.D0
      M=N/2
      DO 11 J=1,M
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
      Y1=WI*(Y(J+1)+Y(N-J+1))
      Y2=0.5*(Y(J+1)-Y(N-J+1))
      Y(J+1)=Y1+Y2
      Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL REALFT(Y,M,+1,NMAX)
      SUM=0.0
      Y(1)=0.5*Y(1)
      Y(2)=0.D0
      DO 12 J=1,N-1,2
      SUM=SUM+Y(J)
      Y(J)=Y(J+1)
      Y(J+1)=SUM
12    CONTINUE
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE REALFT(DATA,N,ISIGN,NMAX)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,DATA
      DIMENSION DATA(2*NMAX)
      PI=4*DATAN(1.D0)
      THETA=PI/DBLE(N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
      C2=-0.5
      CALL FOUR1(DATA,N,+1,NMAX)
      ELSE
      C2=0.5
      THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3
c     DO 11 I=2,N/2+1
      do 11 i=2,n/2
      I1=2*I-1
      I2=I1+1
      I3=N2P3-I2
      I4=I3+1
      WRS=SNGL(WR)
      WIS=SNGL(WI)
      H1R=C1*(DATA(I1)+DATA(I3))
      H1I=C1*(DATA(I2)-DATA(I4))
      H2R=-C2*(DATA(I2)+DATA(I4))
      H2I=C2*(DATA(I1)-DATA(I3))
      DATA(I1)=H1R+WRS*H2R-WIS*H2I
      DATA(I2)=H1I+WRS*H2I+WIS*H2R
      DATA(I3)=H1R-WRS*H2R+WIS*H2I
      DATA(I4)=-H1I+WRS*H2I+WIS*H2R
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
      H1R=DATA(1)
      DATA(1)=H1R+DATA(2)
      DATA(2)=H1R-DATA(2)
      ELSE
      H1R=DATA(1)
      DATA(1)=C1*(H1R+DATA(2))
      DATA(2)=C1*(H1R-DATA(2))
      CALL FOUR1(DATA,N,-1,NMAX)
      ENDIF
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE FOUR1(DATA,NN,ISIGN,NMAX)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,PI,DATA
      DIMENSION DATA(2*NMAX)
      PI=4*DATAN(1.D0)
      N=2*NN
      J=1
      DO 11 I=1,N,2
      IF (J.GT.I) THEN
      TEMPR=DATA(J)
      TEMPI=DATA(J+1)
      DATA(J)=DATA(I)
      DATA(J+1)=DATA(I+1)
      DATA(I)=TEMPR
      DATA(I+1)=TEMPI
      ENDIF
      M=N/2
1     IF ((M.GE.2).AND.(J.GT.M)) THEN
      J=J-M
      M=M/2
      GOTO 1
      ENDIF
      J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
      ISTEP=2*MMAX
      THETA=2*PI/(ISIGN*MMAX)
      WPR=-2.D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.D0
      WI=0.D0
      DO 13 M=1,MMAX,2
      DO 12 I=M,N,ISTEP
      J=I+MMAX
      TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
      TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
      DATA(J)=DATA(I)-TEMPR
      DATA(J+1)=DATA(I+1)-TEMPI
      DATA(I)=DATA(I)+TEMPR
      DATA(I+1)=DATA(I+1)+TEMPI
12    CONTINUE
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
13    CONTINUE
      MMAX=ISTEP
      GOTO 2
      ENDIF
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
      Program graficasdz

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   	C 
C Este programa evalua la integral de la funcional del criterio de  	C 
C atrapamiento para el caso de "restauracion de la ergdicidad".Esta 	C
C diseñado para dos parametros de control (por ejemplo frac. de vol.	C
C y temperatura), pero se puede expandir facilmente para mas para-  	C
C metros.Para cambiar de sistema solo es necesario cambiar subrutina	C
C y variables para el factor de estructura. 			    	C
C								    	C
C Enero de 2007                                                     	C
C									C
C Adaptado a Gaussian core model					C
C									C
C agosto 2008								C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
C variables para parametros del criterio (labda,gcero,1er min S(k))
        double precision l,kmin
C Numero de puntos en vecs. de onda
	Integer nk
	parameter (nk=2**12)
C arreglos para S(k) y g(r)
        double precision  Sdek(nk,2),T(nk)
C parametros S(k) y g(r)
	Double precision parametro1,parametro2 
C variables para evaluar en los parametros de control
	double precision par2max,par2min,dpar2
C variables auxiliares
	Integer apro
        Double precision pi,fv
	Double precision Gder(Nk,2),eps
CC
CCCCCCCCCCCCCCCC Funciones que utiliza el programa CCCCCCCCCCCCCCCCCCC
CC
        double precision  lambda,kminima

CCCCCCCCCCCCCCCCCCC Archivos de salida CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C archivo de salida funcional criterio 
        open(9,file='funcritGC.dat')
	open(10,file='solucionesGC.dat')
	open(20,file='SsdekGC.dat')
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
	write(*,*)'Es el caso de pot. de esferas penetrables'
	write(*,*)'parametro 1: Densidad reducida'
	write(*,*)'parametro 2: T (temperatura)'
	write(*,*)' '
CC
C El programa te hace las graficas de la funcional del criterio como 
C funcion del parametro 2 manteniendo el parametro 1 fijo, es decir 
C cada corrida del programa te genera una iso-parametro1 por decirlo 
C de alguna manera.
CC
        write(*,*)'¿cuanto vale el parametro 1?'
        read (*,*)parametro1
        write(*,*)'¿entre que intervalo evaluamos el parametro 2?'
        write(*,*)'¿maximo?'
        read (*,*)par2max
	write(*,*)'¿minimo?'
        read (*,*)par2min
	write(*,*)'¿Longitud del incremento en el parametro 2?'
        read (*,*)dpar2
! 	write(*,*)'RPA(1), HNC(2)'
! 	read(*,*)apro
! 	write(*,*)'¿cuantos puntos quieres en la fun. del criterio?'
!         read (*,*)Ngama
CC
C Empieza el ciclo para barrer sobre parametro 2
CC
	fv=pi*parametro1/6.D0
	DO parametro2=par2min,par2max,dpar2
          write(*,*)'parametros',parametro1,parametro2
CC
CCCCCCCCCCCCCCCCCCC Modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculado S(k)'
CC Aqui se inserta la subrutina del factor de estrucutura
CC
	 eps=1.D0/parametro2
! 	 if(apro.eq.1)then
! 	  dy=1.5D-2
! 	  do i=1,Nk
! 	   Sdek(i,1)=dble(i)*dy
! 	   Sdek(i,2)=SdkGCRPA(Sdek(i,1),parametro1,eps)
! 	  enddo
! 	 endif
! 	 if(apro.eq.2)then
c... Reads inpud data
	  call input(parametro1,eps,1.D0)
c... perform calculations
	  call oz2(Sdek,Gder)
! 	 endif

C se concoce S(k) y g(r)
CC

	  dy=Sdek(2,1)-Sdek(1,1)
          write(*,*)'calculando la k del primer minimo de S(k)'
          kmin=kminima(SdeK,nk)
	  write(*,*)'kmin',kmin
CC
C ahora se conoce la posicion del primer minimo de S(k)
CC
CCCCCCCCCCCCCCCCCCCCCC Fin modulo de estatica CCCCCCCCCCCCCCCCCCCCCCCCC
CC
	  write(*,*)'calculando la funcional del criterio'
CC
C dependiendo del caso
CC
! 	    fv=parametro2
	    gamamin=0.0D0
	    gamamax=0.3D0
	    dgama=5.D-5
            Do gama=gamamin,gamamax,dgama
              integral= 0.D0
CC
C evaluacion del integrando de la funcional del criterio
CC
              do iint=1, nk
		 y=Sdek(iint,1)
                 S=SdeK(iint,2)
                 l= lambda(y,kmin)
                 num= gama*((S-1.D0)*l*y**2)**2
                 den= 36.D0*pi*fv*(gama*y**2+S*l)*(gama*y**2+l)
                 T(iint)= num/den
              enddo

CC
C evaluacion de la integral
CC
              call SIMPSON(T,dy,nk,integral)
	      write(9,*)gama,1.D0/gama,integral
CC
C buscamos soluciones dentro del lado vitreo
CC
	if((gama.le.0.01).and.(dabs(1.D0-integral).le.1.D-4))then
	 write(10,*)parametro2,gama,1.D0/gama
	endif
	if((gama.gt.0.01).and.(dabs(1.D0-integral).le.1.D-5))then
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
C Subrutinas del factor de estructura,				   C
C 	 							   C
C                                                                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC 
C Subrutinas para S(k) GC con HNC
CC
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine input(fv,epsi,ex)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
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
      rho = fi !(6.0/pi)*fi
c      rho = 0.8
c      
      x(1) = 1.0
      rmax = 160.0  ! rmax=800 para yukaws; rmax = 10-50 para corto alcance
c
      alfa = 0.5
      closure = 2.0
      nr   = 2**12
      nq   = nr
      nrho = 100
      ez   = 1.d-3
c
      sigma1 = 1.0
      sigma2 = 1.0
      sigma2 = sigma2 / sigma1
      sigma1 = sigma1 / sigma1
      sig(1) = sigma1
      sig(2) = (sigma1 + sigma2 )/2.0
      sig(3) = sigma2
      do i = 1, 3
      sigg(i) = sig(i) * 0.7
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
      CALL POT(nr,epsi,ex)
c

      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE POT(n,Tinv,expo)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER(NM=2**14)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
c
c
      ipote = 6
c
      dmed = rho**(-1.0/3.0)
      write(*,*)rho,dmed
c
c
c...Potencial de Lennard-Jones
c
      if (ipote.eq.6) then
      RLAMB = expo
      E(1)=Tinv
      E(3)=Tinv
      E(2)=SQRT(e(1)*e(3))
      Z(1)=RLAMB
      Z(2)=RLAMB
      Z(3)=RLAMB
c
      DO K=1,3
      DO I=1,N
      arg1 = (r(i)/sig(k))**2
      U(I,K)  = E(K)*dexp(-arg1)
      UP(I,K) = 2.D0*U(I,K)*(arg1)**2        ! = -f(r)*r
      ENDDO
      ENDDO
      endif
c
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine oz2(Sk,Gr)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION g(nm,3),c(nm,3),GAM(NM,3),GAM0(NM,3)
      DIMENSION GAM1(NM,3),R1(NM),DIF(2),rk(nm),c1(nm,3)
      DIMENSION GAMAUX(NM,3)
      dimension gh(nm,3),gh2(nm,3),ck(nm,3),s(nm,3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
      Dimension Sk(nr,2),Gr(nr,2)
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
      PARAMETER (NM=2**14)
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
c      PARAMETER (NM=2**14)
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
      SUBROUTINE ESCRIBE(NR,RR,GAM,C,Sdk,Gdr)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION rr(nm),g(nm,3),c(nm,3),GAM(NM,3)
      DIMENSION GAM1(NM,3),R1(NM),rk(nm),c1(nm,3)
      dimension gh(nm,3),gh2(nm,3),ck(nm,3),s(nm,3)
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
      Dimension Sdk(nr,2),Gdr(nr,2)
c
      do k=1,3
      G(1,K)=0.D0
      DO I=2,NR
      g(i,k)=gam(i,k)+c(i,k)+1.d0
      enddo
      enddo
c --- salen las funciones de correlacion gij(r)
!       open(10,file='Gr.out',status='unknown')
!       open(20,file='Cr.out',status='unknown')
         do i=1,nr
         do k=1,3
         gh(i,k)=gam(i,k)+c(i,k)
         enddo
!          write(10,402)rr(i),gh(i,1)+1.0,gh(i,2)+1.0,gh(i,3)+1.0
!          write(20,402)rr(i),c(i,1),c(i,2),c(i,3)
CC
	 Gdr(i,1)= rr(i)
	 Gdr(i,2)= gh(i,1)+1.0
CC
         enddo
!        close(10)
!        close(20)
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
!       open(10,file='hkij.out',status='unknown')
!       open(11,file='ckij.out',status='unknown')
!       open(12,file='ck.out',status='unknown')
!        do i=1,nk
!        write(10,402)rk(i),gh2(i,1),gh2(i,2),gh2(i,3)
!        write(11,402)rk(i),c1(i,1),c1(i,2),c1(i,3)
!        write(12,*)rk(i),c1(i,1)
!        enddo
!  402   format(4(f15.5))
!       close(10)
!       close(11)
!       close(12)
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
!       write(*,*)'salida para Sij(k): Sqij.out'
!       write(*,*)'orden: S11, S12, S22'
!       open(20,file='Sq.out',status='unknown')
      do i=1,nr
!       write(20,66)rk(i),s(i,1),s(i,2),s(i,3)
CC
	Sdk(i,1)= rk(i)
	Sdk(i,2)= s(i,1)
CC
      if(s(i,1).gt.sqmax) sqmax = s(i,1)
      enddo
!       close(20)
      write(*,*)'  '
      WRITE(*,*)'sqmax=',sqmax
      write(*,*)'  '
!  66   format(4(f16.12,5x)) 
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE Ng(N,GAM0,GAM,C)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
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
      !WRITE(*,*)'ITERACION',kj,sngl(ETA),N
      RETURN
      END
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine pp(n,dr,f,g,prod)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
!       ich = 20
   !   open(ich,file='Sq2.out',status='unknown')
      do i=1,n
  !    write(ich,66)q(i),s(i,1),s(i,2),s(i,3)
      if(s(i,1).gt.sqmax) sqmax = s(i,1)
      enddo
   !   write(*,*)'sqmax==',sqmax
   !   close(ich)
 !66   format(4(f16.12,5x))
      return
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine closrel(n,gamma,c)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=2**14)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION GAMMA(NM,3),C(NM,3),sigaux(3)
      COMMON /OZ2C/T,TFIN,KJ
      common /syst1/rho,sig(3),sigg(3),x(2),eps
      COMMON /PARAMETERN/CLOSURE,ALFA,NR,NQ,EZ,NRHO
      common /pote/E(3),z(3),u(nm,3),up(nm,3),rlamb,ipote
      common /syst2/r(nm),q(nm),rmax,qmax,dr,dq
C ----------------------------- PY 
      do k = 1, 3
      if(ipote.eq.1.or.ipote.eq.2.or.ipote.eq.3
     &    .or.ipote.eq.6) then
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
      PARAMETER (NM=2**14)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C Hasta aqui factor de estrucutura				    C
C							            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

       Do h=m , N

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


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c            PROGRAMA: P-D.F
c Se resuelve la ecuacion de O-Z para un sistema de discos
c que interactuan a traves de un potencial de Yukawa. Los
c datos que definen al sistema se dan en la sunrutina INPUT.
c Se usan las cerraduras PY, HNC y RY. S(q) se calcula 
c mediante un algoritmos externos e internos.
c PRUEBA DE ADAPTACION PARA DISCOS SUAVES...LYR OCT/01
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      PARAMETER (NM=1500)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION R(NM),G(NM),C(NM)
      common /struktur/r,g,c
      COMMON /syst/RHO,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      common /pote/zk,z,u(nm)
      call input
c ------------------Se calcula g(r) y c(r).
      CALL OZ1(G,C)
c -----------------------------------S(q).
      stop
      end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine input
      implicit integer*4(i-n), real*8(a-h,o-z)
      PARAMETER (NM=1500)
      COMMON /syst/RHO,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      common /pote/zk,z,u(nm)
      pi=4.d0*datan(1.d0)

c     Atencion!!!: escala=diametro.'

      write(*,*)'elige la cerradura: py(1), hnc(2), ry(3)'
      READ(*,*)CLOSURE
c       closure=2.d0
        IF (CLOSURE.GT.2.5) THEN
        write(*,*)'Parametro de mezcla (escalado)?'
        READ(*,*)ALFA
c        alfa=0.3
        endif
      diam=1.d0
      write(*,*)'densidad reducida'
      read(*,*)rho
c      rho=8.4d-4
      RMAX=50.D0
      NR=1000
      write(*,*)'dame nrho'
      read(*,*)nrho
c      nrho=10
c --- zk= amplitud de la Yukawa, z= screening del potendial ---
c --- El potencial es de la forma Beta*U(r)=ZK*Dexp(-Z*r)/r ---
      write(*,*)'Amplitud del potencial'
      read(*,*) zk
      write(*,*)'Apantallamiento del potencial'
      read(*,*)z
c      zk=zk*dexp(z)
c      z=0.15d0
c      zk=500.d0*dexp(z)
      deps=1.D-3
      return
      end

c -------------------------------------------
c --- se calcula el factor de estructura ---
c --- con el algoritmo interno ---
      subroutine estructura(cc)
      parameter (nm=1500)
      implicit double precision (a-h,o-z)
      dimension cc(nm),c(nm),s(nm)
      COMMON /syst/RHO,diam
      COMMON /OZ1C/R(NM),Q(NM),TAM(NM),KJ
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      pi=4.d0*datan(1.d0)
c ---------------------------s(k) se calcula.

      do i=1,nr
      c(i)=cc(i)
      enddo

c -------------------------------
      CALL FFT(C,RMAX,QMAX,Nr,1)                                     
c ----------------------------------
      do j=1,nr
      s(j)=1.d0-rho*c(j)
      s(j)=1.d0/s(j)
      enddo
      write(*,*)'archivo de salida'
      open(20,file='skdd.dat',status='unknown')
      do i=1,nr
      write(20,*)q(i),s(i)
      enddo
      close(20)
      write(*,*)'salida para s(q)(este es el bueno)= ysin.out'
      return
      end

c --------------------------------------
c --- se calcula el factor de estructura ---
c --- con algoritmo externo ---
      subroutine estructura2(cc,r,nr)
      parameter (nm=1500)
      implicit double precision (a-h,o-z)
      external bessj0
      dimension cc(nm),c(nm),c1(nm),c2(nm),s(nm),r(nm),q(nm)
      COMMON /syst/RHO,diam
      pi=4.d0*datan(1.d0)
c ---------------------------s(k) se calcula.
      do i=1,nr
      c(i)=cc(i)
      enddo

c   -----------
      nq=300
      dq=0.1
c   ----------

      do j=1,nq
      q(j)=dq*j
      do i=1,nr
      x=q(j)*r(i)
      c1(i)=r(i)*bessj0(x)
      enddo
      call pp(nr,r,c1,c,prod)
      s(j)=2.d0*pi*prod
      enddo

      do j=1,nq
      s(j)=1.d0-rho*s(j)
      s(j)=1.d0/s(j)
      enddo
      pen=(s(1)-s(2))/dq
      b=s(1)-pen*q(1)
c      s(0)=b
      s0=b
c      q(0)=0.d0
      q0=0.d0
      write(*,*)'salida suplementaria para s(q)=ysex.out'
c      open(20,file='out2.dat',status='unknown')
      write(20,*)q0,s0
      do i=1,nq
c      write(20,*)q(i),s(i)
      enddo
      close(20)
      return
      end
c -------------------------
      SUBROUTINE OZ1(G,C)
      PARAMETER (NM=1500)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION G(NM),c(nm),GAM(NM),GAM0(NM)
      DIMENSION GAM1(NM),R1(NM),DIF(2),u1(nm)
      COMMON /OZ1C/R(NM),Q(NM),TAM(NM),KJ
      COMMON /syst/RHO,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      common /pote/zk,z,u(nm)
      PI=4.D0*DATAN(1.D0)

      DO I=1,NM
      GAM0(I)=0.d0
      enddo

      rhoa=rho
      DRHO=RHO/(1.d0*nrho)
      CALL CEROS(nr)
      qmax=tam(nr)/rmax

      DO  I=1,nr
         R(I)=tam(I)/qmax
         q(I)=tam(I)/rmax
      enddo

      call pot(nr,r)

      do i=1,nr
      if(r(i).lt.diam)then
      u1(i)=0.d0
      else
      u1(i)=u(i)*(z*r(i))
c      u1(i)=u(i)*(z*r(i)**2+r(i))
      endif
      enddo

      KJ=1
      RHO=KJ*DRHO
12    CALL Ng(nr,GAM0,GAM,C)
      IF (KJ.Gt.1) GOTO 20
      DO 15 I=1,nr
15       GAM0(I)=GAM(I)
      KJ=KJ+1
      RHO=KJ*DRHO
      GOTO 12
20    CALL EXTRAP(GAM0,GAM,nr,RHO,DRHO)
22    DO  I=1,nr
         GAM1(I)=GAM(I)
      enddo

      KJ=KJ+1
      RHO=KJ*DRHO

      CALL Ng(nr,GAM0,GAM,C)
       IF (KJ.EQ.(NRHO-1)) goto 5
 10    continue
       IF (KJ.EQ.(NRHO-1)) goto 14
      CALL EXTRAP(GAM1,GAM,nr,RHO,DRHO)
      DO 27 I=1,nr
27       GAM0(I)=GAM1(I)
      GOTO 22

 5    DO  I=1,nr
      GAM0(I)=GAM(I)
      enddo

C ------------------------------ PY- o HNC.
      IF (CLOSURE.LT.2.5) GOTO 1000
C ------------------------------ RY.
 33   DDRHO=RHOA/100.D0
      DALFA=-ALFA/20.D0
      IX=1
500   RHO=RHOA
      CALL Ng(NR,GAM0,GAM,C)
      DO 100 I=1,NR
      gam1(i)=gam0(i)
100      GAM0(I)=GAM(I)
      call pp(nr,r,c,r,chic)
      CHIC=1.D0-2.D0*PI*RHO*CHIC
      RHO=RHO-DDRHO
      CALL Ng(NR,GAM0,GAM,C)
c ---------------------------g(sigma+) se calcula.
      do i=1,nr
      g(i)=gam(i)+c(i)+1.d0
      enddo
      call pp(nr,r,u1,g,pexv1)
      pexv1=pexv1+dexp(-z*rmax)*(rmax+2.d0/z)*zk
      PEXV1=pi*rho*pexv1/2.d0
      pv1=rho*(1.d0+pexv1)
      write(16,*)'pexv1',pexv1
      RHO=RHO+2.D0*DDRHO
      CALL Ng(NR,GAM0,GAM,C)
c ---------------------------g(sigma+) se calcula.
      do i=1,nr
      g(i)=gam(i)+c(i)+1.d0
      enddo
      call pp(nr,r,u1,g,pexv2)
      pexv2=pexv2+dexp(-z*rmax)*(rmax+2.d0/z)*zk
      PEXV2=pi*rho*pexv2/2.d0
      pv2=rho*(1.d0+pexv2)
      write(16,*)'pexv2',pexv2
c ----------------------------
      CHIV=(PV2-PV1)/(2.D0*DDRHO)
      DIF(IX)=CHIC-CHIV
      WRITE(*,*)'CHIC=',sngl(chic),'CHIV=',sngl(chiv)
      WRITE(*,*)'ALFA=',SNGL(ALFA),'DIF=',SNGL(DIF(IX))
      IF (IX.LT.2) THEN
      IX=IX+1
      ALFA=ALFA+DALFA
      GOTO 500
      ENDIF
      PROD=DIF(1)*DIF(2)
      IF (PROD.LT.0.D0) THEN
      ALFA=ALFA-DALFA
      A=(DIF(2)-DIF(1))/DALFA
      B=DIF(1)-A*ALFA
      ALFA=-B/A
      WRITE(*,*)'ALFA=',SNGL(ALFA)
      GOTO 1000
      ELSE
      DIF(1)=DIF(2)
      ALFA=ALFA+DALFA
      GOTO 500
      ENDIF
1000  RHO=RHOA
      CALL Ng(NR,GAM0,GAM,C)
c --------------------------------------g(r) se calcula.
      write(*,*)'salida para g(r)=ygr.out'
      open(10,file='grss.dat',status='unknown')
      DO I=1,NR
         G(I)=GAM(I)+C(I)+1.D0
         write(10,11) r(i),g(i)
      enddo
      close(10)
c ---------------------------
      call pp(nr,r,c,r,chic)
      CHIC=1.D0-2.D0*PI*RHO*CHIC
      chic=1.d0/chic
      call pp(nr,r,u1,g,pexv)
      pexv=pexv+dexp(-z*rmax)*(rmax+2.d0/z)*zk
      PEXV=pi*rho*pexv/2.d0
      write(*,*)'salida para termodinamica=yres.out'
      open(10,file='yres.out',access='append',status='unknown')
      WRITE(10,13)RHO,PEXV,1.d0+pexv,CHIC,alfa
c      WRITE(*,13)RHO,PEXV,CHIC,alfa
      close(10)
      call estructura(c)
      call estructura2(c,r,nr)
      write(*,*)'  '
      goto 10
 11   format(2(2x,f10.4))
 13   FORMAT(5(F12.6,3x))
 14   return
      END
  
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE Ng(N,GAM0,GAM,C)
      PARAMETER (NM=1500)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION GAM0(NM),GAM(NM),C(NM),f(nm)
      DIMENSION G1(NM),G2(NM),G3(NM),g4(nm),g5(nm)
      DIMENSION d1(nm),d2(nm),D3(NM),d4(nm),d5(nm)
      dimension D01(NM),D02(NM),d03(nm),d04(nm)
      COMMON /OZ1C/R(NM),Q(NM),TAM(NM),KJ
      COMMON /syst/RHO,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      do 5 i=1,n
5        F(I)=GAM0(I)
      call ONg(n,f,g1,c)
      do 10 i=1,n
10       D1(I)=G1(I)-F(I)
      call ONg(n,g1,g2,c)
      do 20 i=1,n
20       D2(I)=G2(I)-G1(I)
      call ONg(n,g2,g3,c)
      do 21 i=1,n
21      d3(i)=g3(i)-g2(i)
 300    call ONg(n,g3,g4,c)
      if(kj.lt.2) goto 200
      do 23 i=1,n
23       D4(I)=G4(I)-G3(I)
      do 40 i=1,n
         D01(I)=D4(I)-D3(I)
         d02(i)=d4(i)-d2(i)
 40      d03(i)=d4(i)-d1(i)
      call pp(n,r,d01,d01,d11)
      call pp(n,r,d01,d02,d12)
      call pp(n,r,d01,d03,d13)
      call pp(n,r,d02,d02,d22)
      call pp(n,r,d02,d03,d23)
      call pp(n,r,d03,d03,d33)
      d21=d12
      d31=d13
      d32=d23
      call pp(n,r,d4,d01,bb1)
      call pp(n,r,d4,d02,bb2)
      call pp(n,r,d4,d03,bb3)
c --- Se analiza el sistema de ecuaciones ----
      ddet=d11*(d22*d33-d23*d32)-d12*(d21*d33-d31*d23)
      ddet=ddet+d13*(d21*d32-d31*d22)     
      deter=ddet
      auxx=dabs(deter)
      auxx2=1.d-50
      if(auxx.lt.auxx2) goto 100
      ddet=bb1*(d22*d33-d23*d32)-d12*(bb2*d33-bb3*d23)
      ddet1=ddet+d13*(bb2*d32-bb3*d22)     
      ddet=d11*(bb2*d33-d23*bb3)-bb1*(d21*d33-d31*d23)
      ddet2=ddet+d13*(d21*bb3-d31*bb2)     
      ddet=d11*(d22*bb3-bb2*d32)-d12*(d21*bb3-d31*bb2)
      ddet3=ddet+bb1*(d21*d32-d31*d22)     
      coef1=ddet1/deter
      coef2=ddet2/deter
      coef3=ddet3/deter
         DO 55 I=1,N
            F(I)=(1.D0-coef1-coef2-coef3)*G4(I)
55          F(I)=F(I)+coef1*g3(i)+coef2*g2(i)+coef3*g1(i)
      call ONg(n,f,g4,c)
      do 60 i=1,n
60       D4(I)=G4(I)-F(I)
100   CALL PP(n,r,D4,D4,D44)
      ETA=DSQRT(D44)
      write(6,*)'eta=',eta
      IF (ETA.LE.deps) goto 200
      do 70 i=1,n
         G1(I)=G2(I)
         g2(i)=g3(i)
         g3(i)=g4(i)
         D1(I)=D2(I)
         d2(i)=d3(i)
 70      d3(i)=d4(i)
      GOTO 300
200    WRITE(6,80)n,rho,ETA
 80    FORMAT(2x,'N=',I5,1x,'DENSIDAD=',F10.6,1x,'PRECISION=',E15.7)
       DO  I=1,N
       GAM(I)=G4(I)
      enddo
      RETURN
      END
c ******************************
      subroutine pp(n,r,f,g,prod)
      PARAMETER (NM=1500)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION F(NM),G(NM),r(nm),aux(nm)
      do i=1,n
         aux(i)=f(i)*g(i)
      enddo
      PROD=0.D0
      DO  I=1,N-1
         PROD=PROD+(aux(I+1)+aux(I))*(R(I+1)-R(I))
      ENDDO
      PROD=PROD/2.D0
      return
      end
c *************************************
      subroutine ONg(n,gamma0,gamma1,c1)
      PARAMETER (NM=1500)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION GAMMA0(NM),GAMMA1(NM),C1(NM)
      COMMON /OZ1C/R(NM),Q(NM),TAM(NM),KJ
      COMMON /syst/RHO,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      call closrel(n,r,gamma0,c1)
      CALL FFT(C1,RMAX,QMAX,N,1)
      DO 10 I=1,N
         DELTA=(1.D0-RHO*C1(I))
10       GAMMA1(I)=RHO*C1(I)**2/DELTA
      CALL FFT(GAMMA1,RMAX,QMAX,N,2)
      call closrel(n,r,gamma1,c1)
      return
      end
C ******************************
      SUBROUTINE EXTRAP(GAM0,GAM,M,RHO,DRHO)
      PARAMETER (NM=1500)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION GAM0(NM),GAM(NM)
         DO 5 I=1,M
            A=(GAM(I)-GAM0(I))/DRHO
            B=GAM(I)-A*RHO
5           GAM0(I)=A*(RHO+DRHO)+B
      RETURN
      END

c --------------------------------------
      SUBROUTINE CLOSREL(n,r,gamma,c)
      PARAMETER (NM=1500)
      implicit integer*4(i-n), real*8(a-h,o-z)
      DIMENSION R(NM),GAMMA(NM),C(NM)
       common /syst/rho,diam
      common /parametern/closure,alfa,rmax,qmax,deps,nr,nrho
      common /pote/zk,z,u(nm)
C ----------------------------- PY.
      if (closure.lt.1.5) then
      DO I=1,N
      c(i)=gamma(i)*(dexp(-u(i))-1.d0)+dexp(-u(i))-1.d0
      enddo
      return
      endif
c ----------------------------- HNC.
      if (closure.gt.1.5.and.closure.lt.2.5) then
      do i=1,n
      x=gamma(i)-u(i)
      c(i)=dexp(x)-gamma(i)-1.d0
      enddo
      return
      endif
C ----------------------------- RY.
      do i=1,n
      f=1.d0-dexp(-alfa*r(i))
      c(i)=(dexp(gamma(i)*f)-1.d0)/f
      c(i)=dexp(-u(i))*(1.d0+c(i))-gamma(i)-1.d0
      enddo
      return
      end

C ------------------------------------------------
c ---  Esta subrutina tabula al potencial tipo Yukawa 
c ---  adimensionado con la energia termica (KB*T)
      subroutine pot(n,r)
      implicit integer*4(i-n), real*8(a-h,o-z)
      parameter(nm=1500) 
      dimension r(nm)
      common /pote/zk,z,u(nm)
      do i=1,n
      if(r(i).lt.1.d0) then
      u(i)=1.d30
      else
c      x=z*r(i)
c      u(i)=zk*dexp(-x)/r(i)
       u(i)=zk/r(i)**z
      endif
      enddo
      close(10)
      write(17,*)'u(rmax)=',u(n)
      return
      end
C *************************************************************
C PROGRAMA PARA LA TRANSFORMADA DE FOURIER EN DOS DIMENSIONES *
C *************************************************************
      SUBROUTINE FFT(FR,RMAX,QMAX,N,IT)
      PARAMETER(NM=1500)
      implicit double precision (a-h,o-z)
      DIMENSION FR(nm),FK(nm),DELK(NM),DELR(NM)
      EXTERNAL BESSJ0,BESSJ1
      COMMON/OZ1C/R(NM),Q(NM),TAM(NM),KJ
      pi=4.0d0*datan(1.0d0)
      IF(IT.EQ.1) THEN
        DO I=1,n-1
          DELR(I)=(QMAX*BESSJ1(tam(I)))**2
          DELR(I)=2.0D0/DELR(I)
        ENDDO
        DO J=1,N-1
         FK(J)=0.0d0
         DO I=1,N-1
         X=Q(J)*R(I) 
         FK(J)= FK(J)+DELR(I)*FR(I)*BESSJ0(X)
        ENDDO
         FK(J)=2.0d0*PI*FK(J)
        ENDDO
         DO I=1,N-1
         FR(I)=FK(I)
        ENDDO
      ELSE
        DO I=1,n-1
          DELK(I)=(RMAX*BESSJ1(tam(I)))**2
          DELK(I)=2.0D0/DELK(I)
        ENDDO
        DO I=1,N-1
         FK(I)=0.0D0
        DO J=1,N-1
         X=Q(J)*R(I)
         FK(I)=FK(I)+DELK(J)*FR(J)*BESSJ0(X)
        ENDDO
         FK(I)=FK(I)/(2.0D0*PI)
        ENDDO
         DO I=1,N-1
         FR(I)=FK(I)
        ENDDO
      ENDIF
      RETURN
      END
c ***********************************************************************
C     *** PROGRAMA PARA CALCULAR LAS RAICES DE LA ***
C     ***     FUNCION DE BESSEL DE ORDEN CERO     ***
      SUBROUTINE CEROS(N)
      implicit double precision (a-h,o-z)
      PARAMETER(NM=1500)
      DIMENSION BE(NM) 
      COMMON/OZ1C/R(NM),Q(NM),TAM(NM),KJ
      TAM(1)=2.4048255577d0
      TAM(2)=5.5200781103d0
      TAM(3)=8.6537279129d0
      TAM(4)=11.7915344391d0
      TAM(5)=14.9309177086d0
      TAM(6)=18.0710639679d0
      TAM(7)=21.2116366299d0
      TAM(8)=24.3524715308d0
      TAM(9)=27.4934791320d0
      PI=4.0d0*datan(1.0d0)
      DO  I=10,N
       BE(I)=(FLOAT(I) -0.25d0)*PI
       Z=8.0d0*BE(I)
       A=1.0d0/Z
       C=124.0d0/(3.0d0*(Z**3))
       D=32.0d0*3779.0d0/(15.0d0*(Z**5))
       E=64.0d0*6277237.0d0/(105.0d0*(Z**7))
       TAM(I)=BE(I)+A-C+D-E
      enddo
      return
      END
c ***********************************************************************
      FUNCTION BESSJ0(X)
      implicit double precision (a-h,o-z)
c      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
c     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D
     *0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(dABS(X).LT.8.0d0)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=dABS(X)
        Z=8.0d0/AX
        Y=Z**2
        XX=AX-.785398164d0
        BESSJ0=SQRT(.6366197720d0/AX)*(dCOS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*dSIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
c ***********************************************************************
      FUNCTION BESSJ1(X)
      implicit double precision (a-h,o-z)
c      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
c     *    S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(dABS(X).LT.8.0d0)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=dABS(X)
        Z=8.0d0/AX
        Y=Z**2
        XX=AX-2.356194491d0
        BESSJ1=SQRT(.636619772d0/AX)*(dCOS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*dSIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *SIGN(1.0d0,X)
      ENDIF
      RETURN
      end
c ***********************************************************************

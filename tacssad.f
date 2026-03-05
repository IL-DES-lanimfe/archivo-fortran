      PROGRAM TAC3SS
C     PROGRAMA BASE: TAC3.f
C     ESTE PROGRAMA CALCULA LAS FUNCIONES DE DISPERSION INTERMEDIA F(k,t) y 
C     SU PARTE SELF Fs(k,t), DE LA TEORIA AUTOCONSISTENTE PROPUESTA EN LA 
C     TESIS LYR. ADICIONALMENTE CALCULA LAS FUNCIONES MEMORIA C(t), C(k,t) y 
C     Cs(k,t).
C     REQUIERE COMO INSUMOS: EL FACTOR DE ESTRUCTURA S(k) Y LOS COEFICIENTES 
C     DE LA MEMORIA EXPONENCIAL SIMPLE COLECTIVA Y SELF.
C     ESTE PROGRAMA CORRESPONDE A CASO DE SISTEMAS COLOIDALES TRIDIMENSIONALES.
C     PRUEBAS CON SS-12: 15-ABRIL-2001 (SANTA01)
C     PROBANDO SOFT/SPHERES, VISITA M.M.N. (26/I/02)
C
C****************************************************************************
C	CON APROXIMACION VINEYARD ADITIVA INCLUIDA:
C	C(K,Z)=Cs(K,Z)+Csexp(K,Z)-CSsexp(K,Z)
C	CS(K,Z)=CSsexp(K,Z)+LAMDA(K)*(DELTA(Z)-CSsexp(K,Z)
C	donde: C(K,Z)=> memoria colectiva y CS(K,Z)=> memoria self
C	agosto 22 de 2006
C****************************************************************************
C
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      PARAMETER (N=1024)
      PARAMETER (M=3000)

      DIMENSION Y(N),T(M),S(N),C(M),FS(N,M),XC(N,M)
      DIMENSION XINC(N),DIFC(M),CFS(M),FSI(M),COLD(M)
      DIMENSION CXC(M),XCI(M)
      DIMENSION CSC(M),CS(N,M),CI(M),XX(N)
      DIMENSION AS(N),BS(N),AT(N),BT(N),XLAM(N)
      DIMENSION CCI(M),CC(N,M),CCC(M),CTC(M)
C
      dimension Sdek(N,2)

      OPEN (10,FILE='Self.dat',status='unknown')
      OPEN (16,FILE='colectiva.dat',status='unknown')
      OPEN (11,FILE='facdes.dat',status='unknown')
      OPEN (12,FILE='ctss26.dat')
      OPEN (15,FILE='pss26.dat')
      OPEN (13,FILE='ckss26.dat')

      PI=4.*ATAN(1.0)

C	SELECCION DE TIPO DE APROXIMACION VINEYARD:
C	IVIN=1 PARA MULTIPLICATIVA
C	IVIN=2 PARA ADITIVA
C	PEDRO: AQUI PODRIMOS EJECUTAR EL PROGRAMA PARA LA VINEYARD DE TU INTERES O
C	BIEN QUE SIEMPRE CALCULE AMBAS. SEGUN TU ELECCION Y GUSTO. YO TE 
C	INCLUYO EN ESTE PROGRAMA LA OPCION DE EJECUTAR CADA OPCION. 
C	LO DIFERENTE EN EL PROGRAMA ESTA ENTRE LAS LINEAS 135 A 143

C	EJEMPLO: VINEYARD MULTIPLICATIVA
	IVIN=1

C     DATOS DE ENTRADA:
      DO  I=1,N
      READ(11,*)Sdek(I,1),Sdek(I,2)
      enddo   
      YMIN=yminima(Sdek,N)
C      YMAX=6.8555
c      RHO=0.5
c      PHI=PI*RHO/6.0
      PHI=0.55
      RHO=6.0*PHI/PI

C     AQUI VOY AQUI VOY.....

C     ESCALA DE DT SEGUN LAS COMPARACIONES CON DB
c      DT=0.0005
      DT=0.0001

      TOL=0.0005

      DO 1 I=1,N
       Y(I)=Sdek(I,1)
       S(I)=Sdek(I,2)
C     FUNCION LAMDA
      XLAM(I)=1.0/(1.0+(Y(I)/YMIN)**2)
1     CONTINUE   

      DY=Y(2)-Y(1)

C     COEFICIENTES DE LA MEMORIA EXPONENCIAL SELF (AS(K) Y BS(K)):

      DO 66 JJ=1,N
         READ(10,*)XX(JJ),AS(JJ),BS(JJ)
66    CONTINUE

C     COEFICIENTES DE LA MEMORIA EXPONENCIAL COLECTIVA AT(K) Y BT(K):

      DO 77 LL=1,N
         READ(16,*)XX(LL),AT(LL),BT(LL)
77    CONTINUE


C     MALLA DEL EJE TEMPORAL:

      T(1)=0.0
      DO 4 L=2,M
      T(L)=T(L-1)+DT
4     CONTINUE 

C     INSUMOS PARA FS, XC Y CC:

      DO 2 J=1,N
      DO 3 K=1,M
         FS(J,K)=EXP(-Y(J)**2*T(K))
         XC(J,K)=EXP(-Y(J)**2*T(K)/S(J))
         CC(J,K)=0.0
3     CONTINUE
2     CONTINUE

C     BANDERIN ITERATIVO:

      ITER=1

C     CALCULO DE C(T):

14    DO 7 KK=1,M
      DO 8 LL=1,N
      XINC(LL)=Y(LL)**4*(S(LL)-1.0)**2*XC(LL,KK)*FS(LL,KK)/S(LL)
8     CONTINUE
      CALL TRAPE(XINC,DY,N,RES)
      C(KK)=RES/(36.0*PHI*PI)
7     CONTINUE

C     CALCULO DE Cs(Y,T) Y Cc(Y,T):

      DO 21 I=1,M
      DO 30 J=1,N
      DO 40 K=1,I
      CCC(K)=CC(J,K)
      CTC(K)=C(K)    
40    CONTINUE

      CALL TRAPE(CCC,DT,I,RES00)
      CALL TRAPE(CTC,DT,I,RES01)

C     FUNCION MEMORIA SELF:

      CS(J,I)=C(I)*XLAM(J)+(1.0-XLAM(J))*AS(J)*EXP(-BS(J)*T(I))

C     FUNCION MEMORIA COLECTIVA SEGUN LA VINEYARD DESEADA.
C	VINEYARD MULTIPLICATIVA IVIN=1,VINEYARD ADITIVA IVIN=2
	IF(IVIN.EQ.1)THEN
      CC(J,I)=AT(J)*(1-XLAM(J))+XLAM(J)*(AT(J)/AS(J))*C(I)+
     #(AT(J)/AS(J))*BS(J)*XLAM(J)*RES01-BT(J)*RES00
	ENDIF
	IF(IVIN.EQ.2)THEN
      CC(J,I)=CS(J,I)+AT(J)*EXP(-BT(J)*T(I))-AS(J)*EXP(-BS(J)*T(I))
	ENDIF

30    CONTINUE
21    CONTINUE

C     COMPARANDO:

      IF(ITER.NE.1)THEN
      DO 15 III=1,M
      DIFC(III)=ABS(C(III)-COLD(III))
15    CONTINUE
      CALL TRAPE(DIFC,DT,M,RESD)
      CALL TRAPE(C,DT,M,RESNEW)
      ERROR=RESD/RESNEW
      IF(ERROR.LT.TOL)THEN
      GO TO 100
      ENDIF

C     VERIFICANDO CONVERGENCIA
      WRITE(*,*)ITER,ERROR
      ENDIF

C     CALCULO ITERATIVO DE FS Y XC:

      DO 9 IJ=1,M
      IF(IJ.EQ.1)THEN
      DO 10 IK=1,N
      FS(IK,IJ)=1.0
      XC(IK,IJ)=1.0
10    CONTINUE
      GO TO 9
      ENDIF

      DO 11 IL=1,N
      DO 12 IM=1,IJ

C     CONVOLUCIONES INVOLUCRADAS

      CFS(IM)=CS(IL,IJ+1-IM)*FS(IL,IM)
      CXC(IM)=CC(IL,IJ+1-IM)*XC(IL,IM)

C     FUNCIONES INVOLUCRADAS

      CI(IM)=CS(IL,IM)
      FSI(IM)=FS(IL,IM)
      CCI(IM)=CC(IL,IM)
      XCI(IM)=XC(IL,IM)

12    CONTINUE

C     INTEGRALES DE Fs(k,t)

      CALL TRAPE(CI,DT,IJ,RES1)
      CALL TRAPE(FSI,DT,IJ,RES2)
      CALL TRAPE(CFS,DT,IJ,RES3)

C     INTEGRALES DE Xc(k,t)=F(k,t)/F(k,0)=F(k,t)/S(k)

      CALL TRAPE(CCI,DT,IJ,RES4)
      CALL TRAPE(XCI,DT,IJ,RES5)
      CALL TRAPE(CXC,DT,IJ,RES6)

C     EVALUANDO 

      FS(IL,IJ)=1.0+RES1-Y(IL)**2*RES2-RES3
      XC(IL,IJ)=1.0+RES4-(Y(IL)**2/S(IL))*RES5-RES6

11    CONTINUE
9     CONTINUE

C     GUARDANDO LA C(T) VIEJA:

      DO 13 IN=1,M
      COLD(IN)=C(IN)
13    CONTINUE

C     ACUMULANDO EL BANDERIN Y CERRANDO EL PROCESO ITERATIVO

      ITER=ITER+1
      GO TO 14

C     ESCRIBIR LA C(T) QUE CONVERGIO:

100   CONTINUE
      DO 16 JJJ=1,M
      WRITE(12,*)SNGL(T(JJJ)),SNGL(C(JJJ))
16    CONTINUE

C     ESCRIBIR LAS FUNCIONES DE DISPERSION (F(k,t) Y Fs(k,t)) Y LAS FUNCIONES 
C     MEMORIA (C(k,t) Y Cs(k,t)) CON QUE CONVERGIO C(T), PARA LOS TIEMPOS TAU 
C     DE INTERES (DB POR EJEMPLO).

      DO 26 I=100,M,100
         TAU=I*DT
         WRITE(15,*)I,'TAU=',TAU
         WRITE(13,*)I,'TAU=',TAU
      DO 27 J=1,N
      XC(J,I)=XC(J,I)*S(J)
      WRITE(15,*)SNGL(Y(J)),SNGL(FS(J,I)),SNGL(XC(J,I))
      WRITE(13,*)SNGL(Y(J)),SNGL(CS(J,I)),SNGL(CC(J,I))
27    CONTINUE
26    CONTINUE

C     SI EL TIEMPO ES SUFICIENTEMENTE LARGO, CALCULAR EL COEFICIENTE DE DIFUSION
C     A TIEMPOS LARGOS

C      CALL TRAPE(C,DT,M,RESDL)

C      DLDO=1.0/(1.0+RESDL)
C      WRITE(*,*)PHI,DLDO

222   CONTINUE
      close (10)
      close (16)
      close (11)
      close (12)
      close (15)
      close (13)

      STOP
      END


C     SUBRUTINA DE INTEGRACION CON TRAPECIO:

      SUBROUTINE TRAPE(F,HACHE,NN,R)
      IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
      DIMENSION F(NN)

      IFIN=NN-1
      IF(IFIN.EQ.1)THEN
      R=HACHE*(F(1)+F(2))/2.0
      GO TO 3
      ENDIF
      R=0.0
      DO 1 I=2,IFIN
      R=R+F(I)
1     CONTINUE
      R=HACHE*((F(1)+F(NN))/2.0+R)
3     RETURN
      END


      Double precision function yminima(Sk,N)

       Implicit none
       
       Integer N
C       Double precision emp

       Integer h,m,j
       Double precision Sk(N,2)
C       Double precision y,dy
       Double precision smin,smax,kmin,kmax

c       Double precision facdes

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

       yminima=kmin
       Return

      End








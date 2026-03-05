      PROGRAM DECIMATIONH
C     PROGRAMA BASE: TAC3.f
C     ESTE PROGRAMA CALCULA LAS FUNCIONES DE DISPERSION INTERMEDIA 
C     F(k,t) y SU PARTE SELF Fs(k,t), DE LA TEORIA AUTOCONSISTENTE 
C     PROPUESTA EN LA TESIS LYR. ADICIONALMENTE CALCULA LAS 
C     FUNCIONES MEMORIA C(t), C(k,t) y Cs(k,t).
C     REQUIERE COMO INSUMOS: EL FACTOR DE ESTRUCTURA S(k) Y LOS 
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL SIMPLE COLECTIVA Y 
C     SELF.
C     ESTE PROGRAMA CORRESPONDE A CASO DE SISTEMAS COLOIDALES 
C     TRIDIMENSIONALES.
C     workshop 2006, tratando de implementar algoritmo de Fuchs para 
C     resolver a tiempos largotes.

       IMPLICIT INTEGER*4(I-N),REAL*8(A-H,O-Z)
!tamaño de la malla en K     
       PARAMETER (NK=1000)
       
!NUMERO DE PUNTOS EN TIEMPOS CORTOS           
       PARAMETER (NT2=10)
             
!tamaño de la malla temporal      
       PARAMETER (MT=300)   
      
!PARAMETRO QUE ME MULTIPLICA EL kmin         
       PARAMETER (RN=1.) 
      
!NUMERO DE REESCALAMIENTOS      
C      PARAMETER (NRES=2)  
      
!NUMERO DE ITERACIONES           
       PARAMETER (NITER=50000)
      
!MALLA K, MALLA T, FACTOR DE ESTRUCTURA            
       DIMENSION Y(NK),T(0:MT+1),S(NK) 
      
!DELTAZ EN CADA REESCALAMIENTO      
       DIMENSION DELTAZ(MT),DELTAZOLD(MT)
      
!FUNCIONES EN CADA REESCALAMIENTO      
       DIMENSION FS(NK,MT),XC(NK,MT)
      
!MEMORIA SELF PARA CADA REESCALAMIENTO      
       DIMENSION CS(NK,MT)
      
!COEFICIENTES SEXP E INTERPOLADORA      
       DIMENSION AS(NK),BS(NK),AT(NK),BT(NK),XLAM(NK)
      
!MEMORIA COLECTIVA PARA CADA REESCALAMIENTO      
       DIMENSION CC(NK,MT)
      
!SEXP EN CADA REESCALAMIENTO      
       DIMENSION CSSEXP(NK,MT),CCSEXP(NK,MT)
      
!INTEGRALES DE COLA EN CADA REESCALAMIENTO      
       DIMENSION DFS(NK,MT),DXC(NK,MT)
       DIMENSION DCS(NK,MT),DCC(NK,MT)
      
!SUMAS EN CADA REESCALAMIENTO      
       DIMENSION A(MT),ASS(MT)
      
!G VINEYARD, GG CANTIDADES ESTATICAS      
       DIMENSION G(NK,MT),GG(NK)
      
!SALIDAS      
       DIMENSION XCOUT(NK,MT),FSOUT(NK,MT),DELTAZOUT(MT)
      
C auxiliares
       CHARACTER*2 DEC
       DIMENSION Sdek(NK,2),Gder(NK,2),GG0(NK),FH(NK)
       PARAMETER (XNU=50.D0)
       OPEN (12,FILE='grsspy.dat')
       OPEN (13,FILE='skddpy.dat')
       OPEN (14,FILE='absjpy.dat')	
       OPEN (15,FILE='abcjpy.dat')
!       OPEN (50,FILE='FUNCION HIDRODINAMICA.DAT')
       OPEN (25,FILE='fktmaxpy.dat')
       OPEN (26,FILE='deltazpy.dat')
       OPEN (27,FILE='fktpy.dat')
       OPEN (28,FILE='sueltospy.dat')
CC Mensaje inicial
       write(*,*)' '
       write(*,*)'Este programa use el metodo de Fuchs para resolver'
       write(*,*)'el esquema autoconsistente para el caso de un'
       write(*,*)'potencial (1/r)**n con n igual a: ',XNU
       write(*,*)' '
       PI=4.*ATAN(1.0)
       IRES=0
C     DATOS DE ENTRADA:
       WRITE(*,*)'¿QUE LAMBDA QUIERES?(1:=1,2:LORENZIANA)'
       READ(*,*)LAMBEL
       WRITE(*,*)'¿QUE VINEYARD QUIERES?(1:MULTIPLICATIVA,2:ADITIVA)'
       READ(*,*)LVINE
       WRITE(*,*)'¿FRACCION DE AREA?'
       READ(*,*)PHI
       WRITE(*,*)'¿CUANTOS REESCALAMIENTOS?'
       READ(*,*)NRES
       WRITE(*,*)'¿QUIERES CHECAR SOLO LA ESTATICA?(s o n)'
       READ(*,*)DEC
c      PHI=0.55D0
       RHO=6.0*PHI/PI
C     ESCALA DE DT
       DT=9.9999999D-5
       TOL=0.00001
C      dk=100.D0/dble(NK)
C      Y(1)=0.001d0
       WRITE(*,*)'Leyendo insumos estaticos'
       DO IK=1,NK
!	READ(50,*)Y(IK)!,FH(IK)
	FH(IK)=1.D0
 	READ(12,*)Gder(IK,1),Gder(IK,2)
 	READ(13,*)Sdek(IK,1),Sdek(IK,2)
	READ(14,*)XX,AS(IK),BS(IK)
	READ(15,*)XX,AT(IK),BT(IK)
	
! 	write(8,*)Sdek(IK,1),Sdek(IK,2)	
       ENDDO
       CLOSE(12)
       CLOSE(13)
       CLOSE(14)
       CLOSE(15)
       write(*,*)'ya esta S(k), g(r) y H(k)'
       WRITE(*,*)'Calculando insumos estaticos.Tenga paciencia'
CCCCCCCCCCC FACTOR DE ESTRUCTURA Y CANTIDAES ESTATICASCCCCCCCCCCCCCCC
       write(*,*)'ya esta S(k), g(r) y Y(r)',PHI
CC minimo de S(k) para lambda y maximo para graficar
       YMIN=yminima(Sdek,NK)
       write(*,*)'el minimo de S(k)',YMIN
       IKM=Lmaximo(Sdek,NK)
       write(*,*)'el maximo de S(k)',IKM,Sdek(IKM,1),Sdek(IKM,2)
       write(28,*)'el maximo de S(k)',IKM,Sdek(IKM,1),Sdek(IKM,2)
! CC
! C inicia el calculo de los coeficientes As y Bs
! CC
!        call SELFMOM(PHI,XNU,6.D0,NK-1,Gder,NK,Sdek,AS,BS)
!        write(*,*)'ya esta As(k) y Bs(k)',PHI
! CC
! C inicia el calculo de los coeficientes Ac y Bc
! CC
!        call COLMOM(PHI,XNU,1.D0,Sdek,NK,Gder,NK-1,FH,AT,BT) 
!        write(*,*)'ya esta Ac(k) y Bc(k)',PHI
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO 1 IK=1,NK
C      READ(12,*)Y(IK),S(IK)
	Y(IK)=Sdek(IK,1)
        S(IK)=Sdek(IK,2)
CC
C Valor inicial de F(k,t) y Correlador
CC
	WRITE(27,*)Y(IK),S(IK),1.D0
	
C	FH(IK)=S(IK)+1.d0
C     FUNCION LAMDA
        IF(LAMBEL.EQ.1)THEN
         XLAM(IK)=1.D0
        ELSE
         XLAM(IK)=1.0/(1.0+(Y(IK)/(RN*YMIN))**2)
        ENDIF
C     FUNCION G(k)
C      GG(IK)=(Y(IK)**4*(S(IK)-1.0)**2)/S(IK)
        GG(IK)=(Y(IK)**3*(S(IK)-1.0)**2)/S(IK)
        GG0(IK)=GG(IK)/(32.D0*PI*PHI)
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL SELF (AS(K) Y BS(K)):
C      WRITE(*,*)AS(IK),BS(IK)
C     COEFICIENTES DE LA MEMORIA EXPONENCIAL COLECTIVA AT(K) Y BT(K):
C      READ(11,*)XX,AT(IK),BT(IK)
1      CONTINUE 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     DIFERENCIAL DE K
        DY=Y(2)-Y(1)
C     calculo de deltaz(t=0)
        CALL SIMPSON(GG0,DY,NK,DEZ0)
        WRITE(*,*)'DELTAZ EN t IGUAL A CERO ES:',DEZ0
        WRITE(28,*)'DELTAZ EN t IGUAL A CERO ES:',DEZ0
CC
C ESCRIBIMOS EL VALOR INICIAL DEL CORRELADOR EN EL MAXIMO
CC
	WRITE(25,*)0.D0,1.D0
C     MALLA DEL EJE TEMPORAL:
        T(0)=0.0
        DO IT=1,MT
         T(IT)=T(IT-1)+DT
 !     READ(13,*)T(IT)!,DELTAZ(IT)
         GRANDO = 0.
         DO IK=1,NK
          FS(IK,IT)=EXP(-Y(IK)**2*T(IT))
C         XC(IK,IT)=EXP(-Y(IK)**2*T(IT)/S(IK))
	  XC(IK,IT)=EXP(-Y(IK)**2*FH(IK)*T(IT)/S(IK))
          CCSEXP(IK,IT)=AT(IK)*EXP(-BT(IK)*T(IT))                    
          CSSEXP(IK,IT)=AS(IK)*EXP(-BS(IK)*T(IT)) 
	  CS(IK,IT)=CSSEXP(IK,IT)
	  CC(IK,IT)=CCSEXP(IK,IT)
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
          GRANDO = GRANDO + GG(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO
          DELTAZ(IT)=GRANDO*DY/(32.D0*PI*PHI)
	  if(IT.le.NT2/2)then
	   write(26,*)T(IT),DELTAZ(IT)
	   XCMAX=XC(IKM,IT)
	   WRITE(25,*)T(IT),XCMAX
	  endif
        ENDDO
	if(DEC.eq.'s')goto 1000
C  AQUI EMPIEZA LO SABROSO>
	write(*,*)' '
	write(*,*)'Por fin termina estatica, empieza lo SABROSO'
	H=DT
        DO IT=NT2/2+1,MT,1  !COMIENZA INTERVALO TEMPORAL
         !write(*,*)"tiempo",it
         DELTAZ(IT)=DELTAZ(IT-1) !IMPUT INICIAL DE DELTAZ
         DO IA=1,NITER         !CICLO DE ITERACIONES        
          DELTAZOLD(IT)=DELTAZ(IT)
          !write(*,*)"iteracion",IA
          GRANDO=0.
          DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
           CS(IK,IT)=CSSEXP(IK,IT)
     ^               +(DELTAZ(IT)-CSSEXP(IK,IT))*XLAM(IK)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   IF(LVINE.EQ.1)THEN
            CONVOL=0.
            DO ITT=1,IT
             CONVOL=CONVOL+DELTAZ(ITT)*EXP(BT(IK)*(T(ITT)-T(IT)))
            ENDDO
            CC(IK,IT)=CS(IK,IT)+
     ^      (1.-XLAM(IK))*(AS(IK)*CCSEXP(IK,IT)/AT(IK)-CSSEXP(IK,IT))+
     ^      (BS(IK)-BT(IK))*XLAM(IK)*CONVOL*DT
            CC(IK,IT)=AT(ik)*CC(IK,IT)/AS(IK)
	   ELSE
	     CC(IK,IT)=CS(IK,IT)+CCSEXP(IK,IT)-CSSEXP(IK,IT)
	   ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           SUMAS=0.
           SUMAC=0.
           IF(MOD(IT,2).EQ.0.) THEN
            DO IS=2,IT/2    !CICLO DE SUMATORIA
             SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
             SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^            +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
            ENDDO           !TERMINA CICLO DE SUMATORIA
            ASS(IT)=-CS(IK,IT/2)*FS(IK,IT/2)+DFS(IK,1)*CS(IK,IT-1)
     ^              +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
            DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	    FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
            A(IT)=-CC(IK,IT/2)*XC(IK,IT/2)+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	    DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	    XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
           ELSE
           DO IS=2,IT/2    !CICLO DE SUMATORIA
            SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
            SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^            +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
           ENDDO           !TERMINA CICLO DE SUMATORIA
           SUMAS=SUMAS+DFS(IK,IT-IT/2)*(CS(IK,IT-IT/2)-CS(IK,IT/2))/2.
     ^           +DCS(IK,IT-IT/2)*(FS(IK,IT-IT/2)-FS(IK,IT/2))/2.
           SUMAC=SUMAC+DXC(IK,IT-IT/2)*(CC(IK,IT-IT/2)-CC(IK,IT/2))/2.
     ^          +DCC(IK,IT-IT/2)*(XC(IK,IT-IT/2)-XC(IK,IT/2))/2.
      PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)+FS(IK,IT/2))/4.
           ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^              +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
           DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	   FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
      PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)+XC(IK,IT/2))/4.
           A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	   DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	   XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ENDIF
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	  GRANDO=GRANDO+GG(IK)*FS(IK,IT)*XC(IK,IT)
         ENDDO             !CICLO SOBRE VECTORES DE ONDA
	 DELTAZ(IT)=GRANDO*DY/(36.*PHI*PI)
         GRANDO=0.
!CRITERIO DE CONVERGENCIA
	 IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	  CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   XCMAX=XC(IKM,IT)
	   WRITE(25,*)T(IT),XCMAX
           GOTO 1001
	  ENDIF
         ENDIF
	 IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	  CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   XCMAX=XC(IKM,IT)
	   WRITE(25,*)T(IT),XCMAX
           GOTO 1001
	  ENDIF
         ENDIF
!FIN CRITERIO DE CONVERGENCIA
         IF(IA.EQ.NITER) then
          WRITE(*,*) "NO CONVIRGIO",T(IT)
	  stop
         endif
         ENDDO              !CICLO DE ITERACIONES         
 1001    continue		
        ENDDO                !TERMINA INTERVALO TEMPORAL

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
COMIENZAN LOS BRINCOS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 201    IRES=IRES+1
        DT=DT+DT
        H=DT
C   REDEFINO FUNCIONES Y DIFERENCIALES ("decimation")
        DO IT=1,MT/2
         DELTAZ(IT)=DELTAZ(2*IT)
         T(IT)=T(2*IT)
         DO IK=1,NK
          FS(IK,IT)=FS(IK,2*IT)
          XC(IK,IT)=XC(IK,2*IT)
          CS(IK,IT)=CS(IK,2*IT)
          CC(IK,IT)=CC(IK,2*IT)
          DFS(IK,IT)=.5*(DFS(IK,2*IT)+DFS(IK,2*IT))
          DXC(IK,IT)=.5*(DXC(IK,2*IT)+DXC(IK,2*IT))
          DCS(IK,IT)=.5*(DCS(IK,2*IT)+DCS(IK,2*IT))
          DCC(IK,IT)=.5*(DCC(IK,2*IT)+DCC(IK,2*IT))
         ENDDO
        ENDDO
        DO IT=MT/2+1,MT,1  !COMIENZA INTERVALO TEMPORAL
         T(IT)=T(IT-1)+DT
CINPUT INICIAL DE DELTAZ
         DELTAZ(IT)=2.*DELTAZ(IT-1)-DELTAZ(IT-2) 
         DO IA=1,NITER         !CICLO DE ITERACIONES        
          DELTAZOLD(IT)=DELTAZ(IT)
          GRANDO=0.
          DO IK=1,NK         !CICLO SOBRE VECTORES DE ONDA
           CSSEXP(IK,IT)=AS(IK)*EXP(-BS(IK)*T(IT))
           CCSEXP(IK,IT)=AT(IK)*EXP(-BT(IK)*T(IT))
           CS(IK,IT)=CSSEXP(IK,IT)
     ^                  +(DELTAZ(IT)-CSSEXP(IK,IT))*XLAM(IK)
c              if(IT.EQ.184)WRITE(*,*)'sexp',CSSEXP(IK,IT),CCSEXP(IK,IT)
c              if(IT.EQ.184)WRITE(*,*)'self',CS(IK,IT),DELTAZ(IT)
c              convol=grando2(ik)+CS(IK,IT)exp(-BT(ik)*T(it))*H     !ojo con los reescal
CCCCCCC   vineyard CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	   IF(LVINE.EQ.1)THEN
             CONVOL=0.
            DO ITT=1,IT
             CONVOL=CONVOL+DELTAZ(ITT)*EXP(BT(IK)*(T(ITT)-T(IT)))
            ENDDO
            CC(IK,IT)=CS(IK,IT)+
     ^      (1.-XLAM(IK))*(AS(IK)*CCSEXP(IK,IT)/AT(IK)-CSSEXP(IK,IT))+
     ^      (BS(IK)-BT(IK))*XLAM(IK)*CONVOL*DT
	    CC(IK,IT)=AT(ik)*CC(IK,IT)/AS(IK)
	   ELSE
            CC(IK,IT)=CS(IK,IT)+CCSEXP(IK,IT)-CSSEXP(IK,IT)
	   ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           SUMAS=0.
           SUMAC=0.
           IF(MOD(IT,2).EQ.0.) THEN
C OJO! FALTA UTILIZAR LOS VALORES YA CALCULADOS EN LUGAR DE LOS PROMEDIOS
           DO IS=2,IT/2    !CICLO DE SUMATORIA
            SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
            SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^            +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
           ENDDO           !TERMINA CICLO DE SUMATORIA
           ASS(IT)=-CS(IK,IT/2)*FS(IK,IT/2)+DFS(IK,1)*CS(IK,IT-1)
     ^              +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
           DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	   FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
           A(IT)=-CC(IK,IT/2)*XC(IK,IT/2)+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C              DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	   DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	   XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ELSE
          DO IS=2,IT/2    !CICLO DE SUMATORIA
           SUMAS=SUMAS+DFS(IK,IS)*(CS(IK,IT-IS)-CS(IK,IT-IS+1))
     ^            +DCS(IK,IS)*(FS(IK,IT-IS)-FS(IK,IT-IS+1))
           SUMAC=SUMAC+DXC(IK,IS)*(CC(IK,IT-IS)-CC(IK,IT-IS+1))
     ^           +DCC(IK,IS)*(XC(IK,IT-IS)-XC(IK,IT-IS+1))
          ENDDO           !TERMINA CICLO DE SUMATORIA
          SUMAS=SUMAS+DFS(IK,IT-IT/2)*(CS(IK,IT-IT/2)-CS(IK,IT/2))/2.
     ^             +DCS(IK,IT-IT/2)*(FS(IK,IT-IT/2)-FS(IK,IT/2))/2.
          SUMAC=SUMAC+DXC(IK,IT-IT/2)*(CC(IK,IT-IT/2)-CC(IK,IT/2))/2.
     ^             +DCC(IK,IT-IT/2)*(XC(IK,IT-IT/2)-XC(IK,IT/2))/2.
      PROM=(CS(IK,IT-IT/2)+CS(IK,IT/2))*(FS(IK,IT-IT/2)+FS(IK,IT/2))/4.

          ASS(IT)=-PROM+DFS(IK,1)*CS(IK,IT-1)
     ^             +(1./H + DCS(IK,1))*FS(IK,IT-1)+SUMAS
          DEN1=(1./H+Y(IK)**2+DCS(IK,1))
      	  FS(IK,IT)=(CS(IK,IT)*(1.-DFS(IK,1))+ASS(IT))/DEN1
      PROM=(CC(IK,IT-IT/2)+CC(IK,IT/2))*(XC(IK,IT-IT/2)+XC(IK,IT/2))/4.
          A(IT)=-PROM+DXC(IK,1)*CC(IK,IT-1)
     ^            +(1./H + DCC(IK,1))*XC(IK,IT-1)+SUMAC
C         DEN2=(1./H+(Y(IK)**2)/S(IK)+DCC(IK,1))
	  DEN2=(1./H+(FH(IK)*Y(IK)**2)/S(IK)+DCC(IK,1))
	  XC(IK,IT)=(CC(IK,IT)*(1.-DXC(IK,1))+A(IT))/DEN2
          ENDIF
c              if(IT.EQ.184)WRITE(*,*)'sum',SUMAS,SUMAC,ASS(IT),A(IT)
          DFS(IK,IT)=FS(IK,IT)
          DXC(IK,IT)=XC(IK,IT)
          DCS(IK,IT)=CS(IK,IT)
          DCC(IK,IT)=CC(IK,IT)
	  GRANDO=GRANDO+GG(IK)*FS(IK,IT)*XC(IK,IT)
c              if(IT.EQ.184)WRITE(*,*)'grando',GRANDO
         ENDDO             !CICLO SOBRE VECTORES DE ONDA
	 DELTAZ(IT)=GRANDO*DY/(36.*PHI*PI)
         GRANDO=0.
!CRITERIO DE CONVERGENCIA
	 IF((DELTAZ(IT)-DELTAZOLD(IT)).GE.0.D0)THEN
	  CRITERIO=(DELTAZ(IT)-DELTAZOLD(IT))/DELTAZ(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   XCMAX=XC(IKM,IT)
	   WRITE(25,*)T(IT),XCMAX
           GOTO 2006
	  ENDIF
         ENDIF
	 IF((DELTAZOLD(IT)-DELTAZ(IT)).GT.0.D0)THEN
	  CRITERIO=(DELTAZOLD(IT)-DELTAZ(IT))/DELTAZOLD(IT)
	  IF(CRITERIO.LT.TOL)THEN
	   if(mod(it,100).eq.0.)then
	    WRITE(*,*)"ad",t(it)
	    do IK=2,NK
	     WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
	    enddo
	   endif
	   write(26,*)T(IT),DELTAZ(IT)
	   XCMAX=XC(IKM,IT)
	   WRITE(25,*)T(IT),XCMAX
           GOTO 2006
	  ENDIF
         ENDIF
!FIN CRITERIO DE CONVERGENCIA
         IF(IA.EQ.NITER) then
          WRITE(*,*) "NO CONVIRGIO",T(IT)
	  Stop
         endif
        ENDDO              !CICLO DE ITERACIONES         
 2006  continue		
       ENDDO                !TERMINA INTERVALO TEMPORAL

       write(*,*)ires,t(mt)
       write(28,*)'valor final de deltaz:',T(MT),DELTAZ(MT)
       IF(IRES.LE.NRES)GOTO 201
!        DO IT=1,MT,10
!         DO IK=1,NK
!          WRITE(27,*)Y(IK),sngl(XC(IK,IT)*S(IK)),sngl(XC(IK,IT))
!          !WRITE(14,*)sngl(Y(IK)),sngl(XCOUT(IK,IT+1)*S(IK))
!         ENDDO
!        ENDDO

      
 1000  close(25)
       close(26)
       close(27)
       STOP
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C SUBRUTINAS PARA EL FACTOR DE ESTRUCUTURA                          C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								    C
C FIN DE SUBRUTINAS PARA EL FACTOR DE ESTRUCTURA. INICIAN SUBRUTI- C
C NAS PARA CALCULAR AS Y BS					    C
C								    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!       SUBROUTINE SELFMOM(PHI,XZK,A,N,GR,NL,SK,AS,BS)
! C     PROGRAMA PARA CALCULAR LA Fs(K,T) TEORICA EN TRES DIMENSIONES 
! C	CON LA APROXIMACION SEXP, PARA EL CASO DE SOFT-SPHERES Y LOS 
! C	COEFICIENTES As(K) Y Bs(K) DE LA TEORIA AUTOCONSISTENTE.
! C     SALIDA: PARA CADA T ARCHIVO K VS Fs(K,T) O As(K) Y Bs(K)
! C     ADAPTADO PARA SELF DEL PROGRAMA COL3SS.F (11-IV-01)
! C     RUMBO A CANCUN-01
! C     MODIFICADO PARA SOFT-SPHERES EN VISITA DE M.M.N. (26-I-02)
! C     AUTOR: LAURA YEOMANS REINA
! C     TRANSFORMADO A SUBRUTINA POR PEDRO EZEQUIEL RAMIREZ 19/JULIO/2006
! 
!       IMPLICIT REAL*8(A-H,O-Z)
! C      PARAMETER (N=1023)
! C      PARAMETER (NL=1024)
! 
!       DIMENSION R(N),FR(N),RK(NL),FK(NL)
!       DIMENSION XM1K(NL),XM21K(NL),XM22K(NL),XM2KA(NL)
!       DIMENSION XM23K(NL),XM2K(NL),XAAA(NL),XBBB(NL)
!       DIMENSION XM31K(NL),XM32K(NL),XM33K(NL)
!       DIMENSION XM34K(NL),XM35K(NL),XM36K(NL)
!       DIMENSION XM37K(NL),XM3K(NL),F12K(NL),XM38K(NL)
!       DIMENSION XU10(N),XU11(N),XU12(N),XU13(N),XM2KK(NL)
!       DIMENSION XU20(N),XU21(N),XU22(N),XU23(N),XU24(N)
!       DIMENSION XU25(N),XKR(N),XJ0(N),XJ1(N)
!       DIMENSION F1(N),F2(N),F3(N),F4(N),F5(N),F2PP(N)
!       DIMENSION F2P(N),F3P(N),F8P(N),F9P(N),XM39K(NL)
!       DIMENSION F6(N),F7(N),F8(N),F9(N),XD1KS(NL),XD2KS(NL)
!       DIMENSION XD0K(NL),XAA(NL),XBB(NL),XD1K(NL),XD2K(NL)
!       DIMENSION F1K(NL),F2K(NL),FTK(NL),XD00K(NL)
!       DIMENSION XM32PK(NL),XM33PK(NL),XM36PK(NL),XM37PK(NL)
! 
! C ARREGLOS NECESARIOS PARA MODIFICAR A SUBRUTINA
! 
! 	DIMENSION GR(NL,2),SK(NL,2),AS(NL),BS(NL)
!       OPEN (15,FILE = 'momentosselfnu10fv805.dat')
!       OPEN (16,FILE = 'coefselfnu10fv805.dat')
! C      OPEN (18, STATUS = 'OLD', FILE = 'skss17p.dat')
! 
!       PI=4.0*ATAN(1.0)
! c	WRITE(*,*)'DAME S0,RK0'
! c	READ(*,*)S0,RK0
! c      WRITE(*,*)'DAME DT'
! c      READ(*,*)DT
! c	WRITE(*,*)'DAME RHO'
! c	READ(*,*)RHO
! C	WRITE(*,*)'DAME PHI'
! C	READ(*,*)PHI
! c      TT=DT*0.00448903
! c	TT=DT*0.0005
! C      PHI=0.00044
! C	RHO=0.003
!        RHO=6.0*PHI/PI
! C      PHI=PI*RHO/6.0
! C      RHO=0.0008403
! C       WRITE(*,*)'DAME Z'
! C       READ(*,*)XZK
! C      XZK=0.15
! C	WRITE(*,*)'DAME A'
! C      READ(*,*)A
! C      XK=XK*EXP(XZK)
!       DO 5 LL=1,N,1
! C      READ(15,*)R(LL),FR(LL)
! 	R(LL)=GR(LL,1)
! 	FR(LL)=GR(LL,2)
! 	if(FR(LL).LE.1.D-14) ICRP=LL
! 5     CONTINUE
!       DO 6 MM=1,NL,1
! C      READ(18,*)RK(MM),FK(MM)
! 	RK(MM)=SK(MM,1)
! 	FK(MM)=SK(MM,2)
! c     EN REALIDAD DE LOS FACTORES DE ESTRUCTURA SOLO LA MALLA (CASO SELF):
!       FK(MM)=1.0
! 6     CONTINUE
! 
!       DR=R(2)-R(1)
!       DO 7 JJ=ICRP,N,1
! C     PARA EL POTENCIAL SOFT-SPHERES-XZK
!       XU10(JJ)=A/R(JJ)**XZK
!       XU11(JJ)=-XZK*XU10(JJ)/R(JJ)
!       XU12(JJ)=(XZK)*(XZK+1)*XU10(JJ)/R(JJ)**2
!       XU13(JJ)=-(XZK)*(XZK+1)*(XZK+2)*XU10(JJ)/R(JJ)**3
! 
!       XU20(JJ)=XU12(JJ)+2.0*XU11(JJ)/R(JJ)
!       XU21(JJ)=XU12(JJ)-XU11(JJ)/R(JJ)
! 
!       XU22(JJ)=XU11(JJ)/R(JJ)
!       XU23(JJ)=XU21(JJ)/R(JJ)
!       XU24(JJ)=XU12(JJ)**2-(XU11(JJ)/R(JJ))**2
!       XU25(JJ)=(XU11(JJ)/R(JJ))**2
!       IF(XU10(JJ).LE.1.D-14)THEN
!        ICRG=JJ
!        NI=ICRG-ICRP
!        GOTO 5000
!       ENDIF
!  7    CONTINUE
! 
!  5000 CONTINUE
!       DO 8 KK=ICRP,ICRG,1
!       F1(KK)=FR(KK)*XU20(KK)*R(KK)**2
!       F2(KK)=FR(KK)*XU24(KK)*R(KK)**2
!       F3(KK)=FR(KK)*XU25(KK)*R(KK)**2
!       F2P(KK)=FR(KK)*F2(KK)
!       F3P(KK)=FR(KK)*F3(KK)
!       F2PP(KK)=FR(KK)*XU12(KK)*R(KK)**4
!  8    CONTINUE
! 
!       CALL SIMPSON(F1,DR,NI,XINT1)
!       CALL SIMPSON(F2,DR,NI,XINT2)
!       CALL SIMPSON(F3,DR,NI,XINT3)
!       CALL SIMPSON(F2P,DR,NI,XINT2P)
!       CALL SIMPSON(F3P,DR,NI,XINT3P)
! 
!       DO 10 II=1,NL,1
! 
! C   PRIMER MOMENTO:
! 
!       XM1K(II)=-RK(II)**2
! 
! C   CHECANDO EL COMPORTAMIENTO DEL PRIMER MOMENTO
! c      WRITE(16,*)RK(II),XM1K(II)
! 
!       XM21K(II)=4.0*RHO*PI*RK(II)**2*XINT1/3.0
! 
!       XM31K(II)=-(4.0*PI*RHO*RK(II)**4)*XINT1
!       XM32K(II)=-(8.0*PI*RHO*RK(II)**2)*XINT2/3.0
!       XM33K(II)=-(8.0*PI*RHO*RK(II)**2)*XINT3 
!       XM32PK(II)=-(8.0*PI*RHO*RK(II)**2)*XINT2P/3.0
!       XM33PK(II)=-(8.0*PI*RHO*RK(II)**2)*XINT3P 
! 
!       DO 9 I=ICRP,ICRG,1
!       XKR(I)=RK(II)*R(I)
!       XKR2=XKR(I)
! 
!       F4(I)=R(I)**2*(FR(I)*XU12(I)*SIN(XKR2)/XKR2)
!       F5(I)=R(I)**2*(FR(I)*XU21(I)*(SIN(XKR2)/XKR2**3-
!      %COS(XKR2)/XKR2**2))
! C      F6(I)=R(I)**2*(FR(I)*XU13(I)*((6.0/XKR2**3-
! C     %1.0/XKR2)*COS(XKR2)+
! C     $3.0*(1.0/XKR2**2-2.0/XKR2**4)*SIN(XKR2)))
! C      F7(I)=R(I)**2*(FR(I)*(XU21(I)/R(I))*(SIN(XKR2)/XKR2**2+
! C     &3.0*(COS(XKR2)/XKR2**3-SIN(XKR2)/XKR2**4)))
! C      F8(I)=R(I)**2*(FR(I)*XU24(I)*((1.0/XKR2-
! C     %2.0/XKR2**3)*SIN(XKR2)+
! C     &2.0*COS(XKR2)/XKR2**2))
! C      F9(I)=R(I)**2*(FR(I)*XU25(I)*SIN(XKR2)/XKR2)
!       F8P(I)=F8(I)*FR(I)
!       F9P(I)=F9(I)*FR(I)
!  9    CONTINUE
! 
!       CALL SIMPSON(F4,DR,NI,XINT4)
!       CALL SIMPSON(F5,DR,NI,XINT5)
!       CALL SIMPSON(F6,DR,NI,XINT6)
!       CALL SIMPSON(F7,DR,NI,XINT7)
!       CALL SIMPSON(F8,DR,NI,XINT8)
!       CALL SIMPSON(F9,DR,NI,XINT9)
!       CALL SIMPSON(F8P,DR,NI,XINT8P)
!       CALL SIMPSON(F9P,DR,NI,XINT9P)
!       CALL SIMPSON(F2PP,DR,NI,XINT2PP)
! 
!       XM22K(II)=-4.0*PI*RHO*RK(II)**2*XINT4
!       XM23K(II)=8.0*PI*RHO*RK(II)**2*XINT5
! 
! C   SEGUNDO MOMENTO:
!       XM2K(II)=XM1K(II)**2+XM21K(II)
! 
! C   SEGUNDO MOMENTO ASINTOTICA HASTA ORDEN RK**4
! c      XM2KA(II)=(1.0+2.0*PI*RHO*XINT2PP/3.0)*XM1K(II)**2
! 
! C   CHECANDO EL COMPORTAMIENTO SEGUNDO MOMENTO
! c      WRITE(16,*)RK(II),XM2K(II)
!       
! C      XM34K(II)=-8.0*PI*RHO*RK(II)**3*XINT6
! C      XM35K(II)=48.0*RHO*PI*RK(II)**3*XINT7
! C     OJO SEGUN MI REVISION EL SIGNO EN LAS SIGUIENTES CUATRO LINEAS DEBE
! C     SER +
!       XM36K(II)=8.0*PI*RHO*RK(II)**2*XINT8
!       XM37K(II)=8.0*PI*RHO*RK(II)**2*XINT9
! C      XM36PK(II)=8.0*PI*RHO*RK(II)**2*XINT8P
! C      XM37PK(II)=8.0*PI*RHO*RK(II)**2*XINT9P
! 
! 
! C   INCLUYENDO TERMINOS APROXIMADOS
! 
! C      S0=(FK(1)-1.0)
! 
!       XM38K(II)=-(XM2K(II)-RK(II)**4)**2/RK(II)**2
! C      XM39K(II)=XM32PK(II)+XM33PK(II)+XM36PK(II)+XM37PK(II)
! 
! C   TERCER MOMENTO:(CON 2 TERMINOS APROXIMADOS)
! 
!       XM3K(II)=XM1K(II)**3+XM31K(II)+XM32K(II)+XM33K(II)
!      #+XM38K(II)
! c     %+XM39K(II)
! 
! C   CHECANDO EL COMPORTAMIENTO  TERCER MOMENTO
!       WRITE(15,*)Real(RK(II)),Real(XM1K(II)),Real(XM2K(II))
!      &           ,Real(XM3K(II))
! 
! 
!  
! C   COEFICIENTE A(K):
! 
!       XAA(II)=XM2K(II)-(XM1K(II)**2)/FK(II)
! 
! C   CHECANDO EL COMPORTAMIENTO  
! c      WRITE(16,*)RK(II),XAA(II)
! 
! C   COEFICIENTE B(K):
! 
!       XBB(II)=(-XM3K(II)-(XM1K(II)**3)/FK(II)**2
!      $+2.0*XM1K(II)*XM2K(II)/FK(II))/XAA(II)
! 
! C   CHECANDO EL COMPORTAMIENTO  
! c      WRITE(16,*)RK(II),XBB(II)
! 
! 
! C   COEFICIENTES DE LA FUNCION DE FRICCION QUE USA LA LAURA
! 
!       XAAA(II)=XAA(II)/RK(II)**2
!       XBBB(II)=XBB(II)-XAAA(II)
! 
!       WRITE(16,*)RK(II),XAAA(II),XBBB(II)
! CCCCCCCCCCCC SALIDA DE LA SUBRUTINA CCCCCCCCCCCCCCCCCCCCC
! 
! 	AS(II)=XAAA(II)
! 	BS(II)=XBBB(II)
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C   RADICAL
! 
!       XD00K(II)=(XBB(II)-(RK(II)**2)/FK(II))**2+
!      #4.0*XAA(II)/FK(II)
!       XD011=XD00K(II)
!       XD0K(II)=SQRT(XD011)
! 
! C   D1(K):
! 
!       XD1K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
!      #+XD0K(II)/2.0
! 
! C   D2(K):
! 
!       XD2K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
!      #-XD0K(II)/2.0
! 
! 
! C   F1(K) EN SEXP:
! 
!       F1K(II)=FK(II)*(XD1K(II)-XBB(II))/(XD1K(II)-XD2K(II))
! 
! C   F2(K) EN SEXP:
! 
!       F2K(II)=FK(II)*(XD2K(II)-XBB(II))/(XD2K(II)-XD1K(II))
! 
! C      WRITE(16,*)SNGL(RK(II)),SNGL(XAAA(II)),SNGL(XBBB(II))
! 
!        F12K(II)=F1K(II)+F2K(II)
! 
!        FTK(II)=F1K(II)*EXP(-XD1K(II)*TT)+F2K(II)*EXP(-XD2K(II)*TT)
! 
! c       WRITE(16,*)SNGL(RK(II)),SNGL(FTK(II))
! 
! c      WRITE(16,*)SNGL(RK(II)),SNGL(XD1K(II)),SNGL(XD2K(II))
! c      WRITE(16,*)SNGL(RK(II)),SNGL(F1K(II)),SNGL(F2K(II))
! 
! 10    CONTINUE
!       CLOSE(15)
!       CLOSE(16)
!       RETURN
!       END
! 
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! C								     C
! C TERMINAN SUBRUTINAS PARA CALCULAR AS Y BS. INICIAN SUBRUTINAS PARA C
! C CALCUALR AC Y BC.                                                  C
! C								     C
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! 
!       SUBROUTINE COLMOM(PHI,XZK,A,SK,NL,GR,N,FH,AC,BC)
! C     PROGRAMA PARA CALCULAR LA F(K,T) TEORICA EN TRES DIMENSIONES 
! C     ENTRADAS: ARCHIVO R VS G(R) SIMULADAS CON DB Y
! C     K VS S(K) TAMBIEN CON DB,PARA EL CASO DE SOFT-SPHERES
! C     SALIDA: PARA CADA T ARCHIVO K VS F(K,T)
! C     ADAPTADO PARA 3D DEL PROGRAMA COL3SEX.F DE ACU� (27-II-01)
! C     RUMBO A CANCUN-01
! C     ADECUACION PARA SOFT-SPHERES, VISITA M.M.N. (26-I-02)
! C     MODIFICADO A SUBRUTINA POR PEDRO EZEQUIEL RAMIREZ GONZALEZ
! C     18/JULIO/2006
!       
!       IMPLICIT REAL*8(A-H,O-Z)
! 
!     
!       DIMENSION R(N),FR(N),RK(NL),FK(NL)
!       DIMENSION XM1K(NL),XM21K(NL),XM22K(NL),XM2KA(NL)
!       DIMENSION XM23K(NL),XM2K(NL),XAAA(NL),XBBB(NL)
!       DIMENSION XM31K(NL),XM32K(NL),XM33K(NL)
!       DIMENSION XM34K(NL),XM35K(NL),XM36K(NL)
!       DIMENSION XM37K(NL),XM3K(NL),F12K(NL),XM38K(NL)
!       DIMENSION XU10(N),XU11(N),XU12(N),XU13(N),XM2KK(NL)
!       DIMENSION XU20(N),XU21(N),XU22(N),XU23(N),XU24(N)
!       DIMENSION XU25(N),XKR(N),XJ0(N),XJ1(N)
!       DIMENSION F1(N),F2(N),F3(N),F4(N),F5(N),F2PP(N)
!       DIMENSION F2P(N),F3P(N),F8P(N),F9P(N),XM39K(NL)
!       DIMENSION F6(N),F7(N),F8(N),F9(N),XD1KS(NL),XD2KS(NL)
!       DIMENSION XD0K(NL),XAA(NL),XBB(NL),XD1K(NL),XD2K(NL)
!       DIMENSION F1K(NL),F2K(NL),FTK(NL),XD00K(NL)
!       DIMENSION XM32PK(NL),XM33PK(NL),XM36PK(NL),XM37PK(NL)
!       
! C ARREGLOS QUE SE NECESITAN PARA MODIFICAR A SUBRUTINA
!       DIMENSION SK(NL,2),GR(NL,2),AC(NL),BC(NL),FH(NL)
! CC
! C Arreglos para moementos con reescalamiento hidrodinamico
! CC
!       DIMENSION XM1KH(NL),XM2KH(NL),XM3KH(NL)
!       OPEN (15,FILE = 'momentoscolectivonu10fv805.dat')
!       OPEN (16,FILE = 'coefcolectivonu10fv805.dat')
!       OPEN (17,FILE = 'momentosH.dat')
! 
!       PI=4.0*ATAN(1.0)
! c	WRITE(*,*)'DAME S0,RK0'
! c	READ(*,*)S0,RK0
! C      WRITE(*,*)'DAME DT'
! C      READ(*,*)DT
! c      TT=DT*0.00448903
! C	TT=DT*0.0005
! C      PHI=0.465
! C	RHO=0.003
!       RHO=6.0*PHI/PI
! C      PHI=PI*RHO/6.0
! C      RHO=0.0008403
! C      XZK=50.0
! C      A=1.0
! C      XK=XK*EXP(XZK)
!       
!       DO 5 LL=1,N,1
! C      READ(15,*)R(LL),FR(LL)
! 	R(LL)=GR(LL,1)
! 	FR(LL)=GR(LL,2)
! 	if(FR(LL).LE.1.D-14) ICRP=LL
! 5     CONTINUE
!       DO 6 MM=1,NL,1
! C      READ(18,*)RK(MM),FK(MM)
! 	RK(MM)=SK(MM,1)
! 	FK(MM)=SK(MM,2)
! c      FK(MM)=1.0
! 6     CONTINUE
! 
!       DR=R(2)-R(1)
!       DO 7 JJ=1,N,1
! C     PARA EL POTENCIAL SOFT-SPHERES
!       XU10(JJ)=A/R(JJ)**XZK
!       XU11(JJ)=-XZK*XU10(JJ)/R(JJ)
!       XU12(JJ)=(XZK)*(XZK+1)*XU10(JJ)/R(JJ)**2
!       XU13(JJ)=-(XZK)*(XZK+1)*(XZK+2)*XU10(JJ)/R(JJ)**3
! 
!       XU20(JJ)=XU12(JJ)+2.0*XU11(JJ)/R(JJ)
!       XU21(JJ)=XU12(JJ)-XU11(JJ)/R(JJ)
! 
!       XU22(JJ)=XU11(JJ)/R(JJ)
!       XU23(JJ)=XU21(JJ)/R(JJ)
!       XU24(JJ)=XU12(JJ)**2-(XU11(JJ)/R(JJ))**2
!       XU25(JJ)=(XU11(JJ)/R(JJ))**2
!       IF(XU10(JJ).LE.1.D-14)THEN
!        ICRG=JJ
!        NI=ICRG-ICRP
!        GOTO 5000
!       ENDIF
!  7    CONTINUE
! 
!  5000 CONTINUE
!       DO 8 KK=ICRP,ICRG,1
!       F1(KK)=FR(KK)*XU20(KK)*R(KK)**2
!       F2(KK)=FR(KK)*XU24(KK)*R(KK)**2
!       F3(KK)=FR(KK)*XU25(KK)*R(KK)**2
!       F2P(KK)=FR(KK)*F2(KK)
!       F3P(KK)=FR(KK)*F3(KK)
!       F2PP(KK)=FR(KK)*XU12(KK)*R(KK)**4
!  8    CONTINUE
! 
!       CALL SIMPSON(F1,DR,N,XINT1)
!       CALL SIMPSON(F2,DR,N,XINT2)
!       CALL SIMPSON(F3,DR,N,XINT3)
!       CALL SIMPSON(F2P,DR,N,XINT2P)
!       CALL SIMPSON(F3P,DR,N,XINT3P)
! 
! 
!       DO 10 II=1,NL,1
! 
! C   PRIMER MOMENTO:
! 
!       XM1K(II)=-RK(II)**2
! CC
! C PRIMER MOMENTO CON REESCALAMIENTO HIDORDINAMICO
! CC
!       XM1KH(II)=FH(II)*XM1K(II)
! 
! C   CHECANDO EL COMPORTAMIENTO DEL PRIMER MOMENTO
! c      WRITE(16,*)RK(II),XM1K(II)
! 
!       XM21K(II)=4.0*RHO*PI*RK(II)**2*XINT1/3.0
! 
!       XM31K(II)=-(4.0*PI*RHO*RK(II)**4)*XINT1
!       XM32K(II)=-(8.0*PI*RHO*RK(II)**2)*XINT2/3.0
!       XM33K(II)=-(8.0*PI*RHO*RK(II)**2)*XINT3 
!       XM32PK(II)=-(8.0*PI*RHO*RK(II)**2)*XINT2P/3.0
!       XM33PK(II)=-(8.0*PI*RHO*RK(II)**2)*XINT3P 
! 
!       DO 9 I=ICRP,ICRG,1
!       XKR(I)=RK(II)*R(I)
!       XKR2=XKR(I)
! 
!       F4(I)=R(I)**2*(FR(I)*XU12(I)*SIN(XKR2)/XKR2)
!       F5(I)=R(I)**2*(FR(I)*XU21(I)*(SIN(XKR2)/XKR2**3-
!      %COS(XKR2)/XKR2**2))
!       F6(I)=R(I)**2*(FR(I)*XU13(I)*((6.0/XKR2**3-
!      %1.0/XKR2)*COS(XKR2)+
!      $3.0*(1.0/XKR2**2-2.0/XKR2**4)*SIN(XKR2)))
!       F7(I)=R(I)**2*(FR(I)*(XU21(I)/R(I))*(SIN(XKR2)/XKR2**2+
!      &3.0*(COS(XKR2)/XKR2**3-SIN(XKR2)/XKR2**4)))
!       F8(I)=R(I)**2*(FR(I)*XU24(I)*((1.0/XKR2-
!      %2.0/XKR2**3)*SIN(XKR2)+
!      &2.0*COS(XKR2)/XKR2**2))
!       F9(I)=R(I)**2*(FR(I)*XU25(I)*SIN(XKR2)/XKR2)
!       F8P(I)=F8(I)*FR(I)
!       F9P(I)=F9(I)*FR(I)
!  9    CONTINUE
! 
!       CALL SIMPSON(F4,DR,NI,XINT4)
!       CALL SIMPSON(F5,DR,NI,XINT5)
!       CALL SIMPSON(F6,DR,NI,XINT6)
!       CALL SIMPSON(F7,DR,NI,XINT7)
!       CALL SIMPSON(F8,DR,NI,XINT8)
!       CALL SIMPSON(F9,DR,NI,XINT9)
!       CALL SIMPSON(F8P,DR,NI,XINT8P)
!       CALL SIMPSON(F9P,DR,NI,XINT9P)
!       CALL SIMPSON(F2PP,DR,NI,XINT2PP)
! 
!       XM22K(II)=-4.0*PI*RHO*RK(II)**2*XINT4
!       XM23K(II)=8.0*PI*RHO*RK(II)**2*XINT5
! 
! C   SEGUNDO MOMENTO:
!       XM2K(II)=XM1K(II)**2+XM21K(II)+XM22K(II)+XM23K(II)
! CC
! C SEGUNDO MOMENTO CON REESCALAMIENTO HIDRODINAMICO
! CC
!       XM2KH(II)=XM2K(II)*FH(II)**2
! 
! C   SEGUNDO MOMENTO ASINTOTICA HASTA ORDEN RK**4
! c      XM2KA(II)=(1.0+2.0*PI*RHO*XINT2PP/3.0)*XM1K(II)**2
! 
! C   CHECANDO EL COMPORTAMIENTO SEGUNDO MOMENTO
! c      WRITE(16,*)RK(II),XM2K(II)
!       
!       XM34K(II)=-8.0*PI*RHO*RK(II)**3*XINT6
!       XM35K(II)=48.0*RHO*PI*RK(II)**3*XINT7
! C     OJO SEGUN MI REVISION EL SIGNO EN LAS SIGUIENTES CUATRO LINEAS DEBE
! C     SER +
!       XM36K(II)=8.0*PI*RHO*RK(II)**2*XINT8
!       XM37K(II)=8.0*PI*RHO*RK(II)**2*XINT9
!       XM36PK(II)=8.0*PI*RHO*RK(II)**2*XINT8P
!       XM37PK(II)=8.0*PI*RHO*RK(II)**2*XINT9P
! 
! 
! C   INCLUYENDO TERMINOS APROXIMADOS
! 
! C      S0=(FK(1)-1.0)
! 
!       XM38K(II)=-(XM2K(II)-RK(II)**4)**2/RK(II)**2
!       XM39K(II)=XM32PK(II)+XM33PK(II)+XM36PK(II)+XM37PK(II)
! 
! C   TERCER MOMENTO:(CON 2 TERMINOS APROXIMADOS)
! 
!       XM3K(II)=XM1K(II)**3+XM31K(II)+XM32K(II)+
!      #XM33K(II)+XM34K(II)+XM35K(II)+XM36K(II)+XM37K(II)
!      #+XM38K(II)
! c     %+XM39K(II)
! CC
! C TERCER MOMENTO CON REESCALAMIENTO TERMODINAMICO
! CC
!       XM3KH(II)=XM3K(II)*FH(II)**3
! 
! C   CHECANDO EL COMPORTAMIENTO MOMENTOS
!       WRITE(15,*)Real(RK(II)),Real(XM1K(II)),Real(XM2K(II))
!      &           ,Real(XM3K(II))
! CC
! C CHECANDO EL COMPORTAMIENTO MOMENTOS CON REESCALAMIENTO
! CC
!       WRITE(17,*)Real(RK(II)),Real(XM1KH(II)),Real(XM2KH(II))
!      &           ,Real(XM3KH(II))
! 
! C   COEFICIENTE A(K):
! 
!       XAA(II)=XM2K(II)
!      $-(XM1K(II)**2)/FK(II)
! 
! C   CHECANDO EL COMPORTAMIENTO  
! c      WRITE(16,*)RK(II),XAA(II)
! 
! C   COEFICIENTE B(K):
! 
!       XBB(II)=(-XM3K(II)-(XM1K(II)**3)/FK(II)**2
!      $+2.0*XM1K(II)*XM2K(II)/FK(II))/XAA(II)
! 
! C   CHECANDO EL COMPORTAMIENTO  
! c      WRITE(16,*)RK(II),XBB(II)
! 
! 
! C   COEFICIENTES DE LA FUNCION DE FRICCION QUE USA LA LAURA
!       XAAA(II)=XAA(II)/RK(II)**2
!       XBBB(II)=XBB(II)-XAAA(II)
!       WRITE(16,*)RK(II),XAAA(II),XBBB(II)
! C SALIDA DE LA SUBRUTINA
! 
! 	AC(II)=XAAA(II)
! 	BC(II)=XBBB(II)
! 	
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
! C   RADICAL
! 
!       XD00K(II)=(XBB(II)-(RK(II)**2)/FK(II))**2+
!      #4.0*XAA(II)/FK(II)
!       XD011=XD00K(II)
!       XD0K(II)=SQRT(XD011)
! 
! C   D1(K):
! 
!       XD1K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
!      #+XD0K(II)/2.0
! 
! C   D2(K):
! 
!       XD2K(II)=(RK(II)**2/FK(II)+XBB(II))/2.0
!      #-XD0K(II)/2.0
! 
! 
! C   F1(K) EN SEXP:
! 
!       F1K(II)=FK(II)*(XD1K(II)-XBB(II))/(XD1K(II)-XD2K(II))
! 
! C   F2(K) EN SEXP:
! 
!       F2K(II)=FK(II)*(XD2K(II)-XBB(II))/(XD2K(II)-XD1K(II))
! 
! C      WRITE(16,*)SNGL(RK(II)),SNGL(XAAA(II)),SNGL(XBBB(II))
! 
!        F12K(II)=F1K(II)+F2K(II)
! 
!        FTK(II)=F1K(II)*EXP(-XD1K(II)*TT)+F2K(II)*EXP(-XD2K(II)*TT)
! 
! c       WRITE(16,*)SNGL(RK(II)),SNGL(FTK(II))
! 
! c      WRITE(16,*)SNGL(RK(II)),SNGL(XD1K(II)),SNGL(XD2K(II))
! c      WRITE(16,*)SNGL(RK(II)),SNGL(F1K(II)),SNGL(F2K(II))
! 
! 10    CONTINUE
!       CLOSE(15)
!       CLOSE(16)
!       CLOSE(17)
!       RETURN
!       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C SUBRUTINAS AUXILIARES, SIMPSON SE USA EN LAS TRES ANTERIORES Y LAS C
C OTRAS DOS EN SELFMOM Y COLMOM					     C
C							             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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


      SUBROUTINE XJOTA0(XKR2,XXJ0)
C     ***********************************************
C     SUBRUTINA PARA EL CALCULO DE LA J_(0). 
C     BASE: ABRAMOWITZ, p.369.
C     ***********************************************
      IMPLICIT REAL*4(A-H,O-Z)
       if(abs(XKR2).le.3.0)then
         yy=(XKR2/3.0)**2

         XXJ0=1.0-2.2499997*yy+1.2656208*yy**2-
     *0.3163866*yy**3+0.0444479*yy**4-0.0039444*
     *yy**5+0.0002100*yy**6    

       else
         yy=(3.0/XKR2)

         TETA0=XKR2-0.78539816-0.04166397*yy-0.00003954*
     *yy**2+0.00262573*yy**3-0.00054125*yy**4-
     *0.00029333*yy**5+0.00013558*yy**6

         XF0=0.79788456-0.00000077*yy-0.00552740*yy**2-
     *0.00009512*yy**3+0.00137237*yy**4-0.00072805*yy**5 +
     *0.00014476*yy**6

         XXJ0=XF0*COS(TETA0)/sqrt(XKR2)

       endif

       RETURN
       END


      SUBROUTINE XJOTA1(SPL,XXJ1)
C     ***********************************************
C     SUBRUTINA PARA EL CALCULO DE LA J_(1). 
C     BASE: ABRAMOWITZ p.370.
C     ***********************************************
      IMPLICIT REAL*4 (A-H,O-Z)
      if(abs(SPL).le.3.0)then
        yy=(SPL/3.0)**2

        XXJ1=(0.5-0.56249985*yy+0.21093573*yy**2-
     *0.03954289*yy**3+0.00443319*yy**4-0.00031761*
     *yy**5+0.00001109*yy**6)*SPL

      else
        yy=(3.0/SPL)

         TETA1=SPL-2.35619449+0.12499612*yy+0.00005650*
     *yy**2-0.00637879*yy**3+0.00074348*yy**4+0.00079824*
     *yy**5-0.00029166*yy**6    

         XF1=0.79788456+0.00000156*yy+0.01659667*yy**2+
     *0.00017105*yy**3-0.00249511*yy**4+0.00113653*yy**5-
     *0.00020033*yy**6

        XXJ1=(XF1*COS(TETA1))/sqrt(SPL)

      endif
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C							             C
C La siguiente funcion sirve para calcular el primer minimo de S(k)  C
C para la funcion interpoladora.                                     C
C								     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								      C
C Para observar el plateu localizamos el primer pico del correlador.  C
C la siguiente funcion hace esa chamba.				      C
C								      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Integer function Lmaximo(Sk,NK)
       Implicit none
       Integer NK
       Integer j,h
       Double precision Sk(NK,2)
       Double precision smax
     	smax=Sk(1,2)
       Do h=2,NK
	if (Sk(h,2).GT.smax)then
	  smax=Sk(h,2)
	  j=h
	endif
       enddo
       Lmaximo=j
       Return

       End


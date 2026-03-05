      PROGRAM cristalizacion
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C								     C
C Este programa aplica el criterio de Hansen-Verlet para ubicar la   C
C cristalizacion, evalua el factor de estructura y calcula el maximo C
C usa bisecciones sucesivas para determinar cuando el maximo vale 2.8C
C ,esta hecho para dos parametros de control, la aplicacion a un sis-C
C tema particular se da mediante una subrutina para S(k)	     C
C								     C
C Noviembre de 2006					 	     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	implicit none ! Toda variable se tiene que declarar
CC
C subrutinas y funciones que se incluyen desde archivos externos
CC
	include 'maxmin.f'

CC
C variables para los parametros de control
CC
	Integer Npar2
	Double precision par1max,par1min,par2max,par2min
	Double precision parametro2
	Double precision dpar2
CC
C variables para las bisecciones sucesivas
CC
	Double precision prodmaxmin,prodant
	Double precision par1ant,par1post,par1med
	Double precision solant,solpost,solmed
CC
C variables del factor de estructura
CC
	Integer Nk
	parameter(Nk=2**12)
	Double precision Sdek(Nk,2)
CC
C variables auxiliares
CC
	Integer imaximo,iminimo
CC
C funciones externas
CC
	external maxmin
CC
C Recoleccion de datos
CC
	write(*,*)'Programa cristalizacion'
	write(*,*)' '
	write(*,*)'diga el intervalo de busqueda en el parametro 1'
	write(*,*)'¿maximo?'
	read(*,*)par1max
	write(*,*)'¿minimo?'
	read(*,*)par1min
	write(*,*)'diga el intervalo de busqueda en el parametro 2'
	write(*,*)'¿maximo?'
	read(*,*)par2max
	write(*,*)'¿minimo?'
	read(*,*)par2min
	write(*,*)'¿numero de puntos para dividir el intervalo?'
	read(*,*)Npar2
CC
C Archivo de salida
CC
	open(10,file='cristalizacion.dat')
CC
C Fijamos el parametro 2 y hacemos bisecciones sucesivas en el 
C parametro 1
CC
	dpar2=(par2max-par2min)/dble(Npar2)
	DO parametro2=par2min,par2max,dpar2
	 write(*,*)'parametro 2=',parametro2
	 write(*,*)'calculando S(k) en la cota inferior'
CC
C aqui va la subrutina de S(k) evaluada en par1min y parametro2
CC
	 call
CC
C calculamos el maximo de S(k)
CC
	 call maxmin(Sdek,Nk,imaximo,iminimo)
	 smx=Sdek(imaximo,2)
	 write(*,*)'el maximo de S(k)',smx,par1min
CC
C guardamos el valor de la solucion y elegimos si ya le atinamos
CC
	 solant=2.8-smx
	 if(dabs(solant).EQ.1.D-4)then
	  write(10,*)real(par1min),real(parametro2),real(smx)
	  write(*,*)'Se halló solución',par1min,parametro2
	  goto 200
	 endif
CC
	 write(*,*)'calculando S(k) cota superior'
CC
C aqui va la subrutina de S(k) evaluada en par1max y parametro2
CC
	 call
CC
C calculamos el maximo de S(k)
CC
	 call maxmin(Sdek,Nk,imaximo,iminimo)
	 smx=Sdek(imaximo,2)
	 write(*,*)'el maximo de S(k)',smx,par1max
CC
C guardamos el valor de la solucion y elegimos si ya le atinamos
CC
	 solpost=2.8-smx
	 if(dabs(solpost).EQ.1.D-4)then
	  write(10,*)real(par1max),real(parametro2),real(smx)
	  write(*,*)'Se halló solución',par1max,parametro2
	  goto 200
	 endif
CC
C bisecciones sucesivas
CC
	 prodmaxmin=solant*solpost
	 IF (prodmaxmin.LT.0.D0)then
	  par1ant=par1min
	  par1post=par1max
 100	  continue
	  par1med=(par1post+par1ant)/2.D0
	  write(*,*)'calculando S(k) en valor intermedio',par1med
CC
C aqui va la subrutina de S(k) evaluada en par1med y parametro2
CC
	  call
CC
C calculamos el maximo de S(k)
CC
	  call maxmin(Sdek,Nk,imaximo,iminimo)
	  smx=Sdek(imaximo,2)
	  write(*,*)'el maximo de S(k)',smx,par1med
CC
C guardamos el valor de la solucion y elegimos si ya le atinamos
CC
	  solmed=2.8-smx
	  if(dabs(solmed).EQ.1.D-4)then
	   write(10,*)real(par1med),real(parametro2),real(smx)
	   write(*,*)'Se halló solución',par1med,parametro2
	   goto 200
	  endif
CC
C elegimos el intervalo donde esta la solucion
CC
	  prodant=solmed*solant
	  if(prodant.LT.0.D0)then
	   par1post=par1med
	   solpost=solmed
	   goto 100
	  endif
	  par1ant=par1med
	  solant=solmed
	  goto 100
	 ENDIF
CC
C si no se cumple la implicacion anterior
CC
	 write(*,*)'la solucion no esta en este intervalo'
CC
 200	 continue
	ENDDO
	close(10)
	stop
      END
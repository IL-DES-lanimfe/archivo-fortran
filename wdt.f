	program descuadmed
	implicit double precision (a-h,o-z),integer*4(i-n)
	parameter(np=7000,nt=10**6)
	dimension w(0:nt),dz(0:np),t(0:np),tv(0:nt),wold(0:nt)
c
! 	lectura de la delta zeta estrella
	open(1,file='deltazdurafv485.dat')
	do i=1,7000
	read(1,*)t(i),dz(i)
	enddo
	close(1)
! mallas iniciales para el tiempo y w de t
	dt=1.d-3
	do i=0,nt
	tv(i)=dble(i)*dt
	w(i)=tv(i)
	enddo
! 	stop
! 
	open(2,file='w.dat')
	do i=1,nt !ciclo temporal
	w(i)=w(i-1)
	do m=1,300 !ciclo de iteraciones
	wold(i)=w(i)
	suma=0.d0
	do j=1,i ! ciclo de la integral
	suma=suma+delta(tv(i)-tv(j),t,dz,np)*w(j)
	enddo
	w(i)=tv(i)-suma*dt
	if(dabs(w(i)-wold(i))/w(i).lt.1.d-4)goto 10
	enddo
10	continue
	print*,m,tv(i),w(i)
	write(2,*)tv(i),w(i)
	enddo
	close(2)

	stop
	end

	double precision function delta(tv,t,dz,np)
	implicit double precision (a-h,o-z),integer*4(i-n)
	dimension w(np),dz(np),t(np)
	do i=1,np
	if(tv.eq.t(i))then
	delta=dz(i)
	return
	endif
	if(tv.lt.t(1))then
	delta=dz(1)
	return
	endif
	if(tv.lt.t(i))then
	delta=(dz(i)+dz(i-1))/2.d0
	return
	endif
	enddo
	return
	end
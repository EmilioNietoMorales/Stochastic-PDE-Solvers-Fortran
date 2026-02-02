C-----------------------Empezamos el programa principal----------------------------

      program subproyecto1
      implicit none
      real*8 Rcil, Lcil, MCstep, RND, x, y, z
      integer*8 NMC, acc
      integer i
      character*10 nada
      dimension MCstep(1:5), acc(1:5)
      
C comprobamos que los datos están bien introducidos (imprimiendolos en la terminal): 
      
      open(1,FILE='datos_in.txt')
      
      read(1,*) Rcil, nada, Lcil,nada, NMC,nada
      
      do i=1,5
      	read(1,*) MCstep(i), nada
      end do
      
      close(1)
      
      print*, 'Radio cilindro = ',Rcil, 'Largo cilindro = ',Lcil,
     &'NMC = ',NMC,'MCstep = {', MCstep(1),',',MCstep(2),',',
     &MCstep(3),',',MCstep(4),',',MCstep(5),'}'
     
C Llamamos a la subrutina MC para obtener la simulacion. Tomamos de punto inicial x=y=0, z=RND

      call init_random_seed()
      call random_number(RND)

      x = 0
      y = 0
      z = RND*Lcil
      acc=0
      
     
      
      call MC(Rcil,Lcil,NMC,MCstep(1),x,y,z,acc(1))    
            

      
      open(45,FILE='MC_acc.txt')
      
      do i=1,5
      	write(45,*) MCstep(i), acc(i)/real(NMC)
      end do
      
      close(45)
      
      end program subproyecto1
      
      
      
C-------------------------------subrutina MC-----------------------------------------------------


      subroutine MC(Rcil,Lcil,NMC,pasoMC,x,y,z,acc)
      implicit none
      real*8 Rcil, Lcil, pasoMC, x, y, z,RND, R, xf, yf, zf, r2
      integer*8 NMC, acc
      integer i
      

      R = 0.5
      r2 = (Rcil - R)**2
C hacemos un bucle de 1 a NMC (intentos de movimientos) donde cada iteración sea un movimiento
C abrimos un fichero donde almacenaremos la trayectoria de la particula

      open(2,FILE='trayectoria1.txt')
      

      do i=1, NMC
      
C dentro del bucle programamos 1 movimiento de la particula      
      
      	call init_random_seed()
      	
      	call random_number(RND)
      	xf = x + (RND - R)*pasoMC
      	
      	call random_number(RND)
      	yf = y +(RND - R)*pasoMC
      	
      	call random_number(RND)
      	zf = z +(RND - R)*pasoMC
      	
C corregimos zf:
	zf = zf - Lcil * floor(zf/Lcil)     	
      	
      	if(yf**2+xf**2 .LE. r2) then
      		
      		x = xf
      		y = yf
      		z = zf
      		
      		write(2,*) x , y , z
      		
      		acc = acc + 1
      	end if 
      end do
      
      print*, 'La tasa de aceptación es: ', real(acc)/NMC * 100, '%'
      
      close(2)
      
     
      
      
      
      
      return 
      end
      
      
      
c-------------------------subrutina para inicializar el generador de numeros aleatorios      
      subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms

      call random_seed(size = n)
      allocate(seed(n))
! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream",
     +  form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
        read(un) seed
        close(un)
      else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
        call system_clock(count)
        if (count /= 0) then
          t = transfer(count, t)
        else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 
     -         + dt(2) * 31_8 * 24 * 60 * 60 * 1000 
     -         + dt(3) * 24 * 60 * 60 * 60 * 1000 
     -         + dt(5) * 60 * 60 * 1000 
     -         + dt(6) * 60 * 1000 + dt(7) * 1000 
     -         + dt(8)
          t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n.ge.3) then
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          end if
        else
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed

      

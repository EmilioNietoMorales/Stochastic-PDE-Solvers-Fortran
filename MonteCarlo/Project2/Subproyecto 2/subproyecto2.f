C-----------------------Empezamos el programa principal----------------------------

      program subproyecto2
      implicit none
      real*8 Rcil, Lcil, MCstep, RND, x, y, z,volbinz,volbinr,pi,r
      real*8 normalR, normalZ
      integer*8 NMC, accM,accP, Nmed, histoR, histoZ, Nbinr, Nbinz
      integer i,j
      character*10 nada
      parameter (pi=4.d0*atan(1.d0))
      
C dimensionamos los vectores histoR y histoZ con un numero suficientemente grande y los 
C inicializamos a 0
      
      dimension histoR(1:300), histoZ(1:300)
      
      do i=1,300
      	histoR(i) = 0
      	histoZ(i) = 0
      end do
      
C comprobamos que los datos están bien introducidos (imprimiendolos en la terminal): 
      
      open(1,FILE='datos_in2.txt')
      
      read(1,*) Rcil, nada, Lcil,nada, NMC,nada, Nmed, nada
     &,MCstep, nada, Nbinr, nada, Nbinz, nada
      
      
      close(1)
      
      print*, 'Radio cilindro = ',Rcil, 'Largo cilindro = ',Lcil,
     &'NMC = ',NMC,'Nmed = ',Nmed,'MCstep = ',MCstep,'Nbinr = ',Nbinr
     &,'Nbinz = ',Nbinz
     
C Llamamos a la subrutina precalentamiento para obtener la simulacion. Tomamos de punto inicial
C x=y=0, z=RND

      call init_random_seed()
      call random_number(RND)

      x = 0
      y = 0
      z = RND*Lcil
      accP=0
      accM = 0
      

      call precalentamiento(Rcil,Lcil,NMC,MCstep,x,y,z,accP)    
      
C una vez finalizado el precalentamiento, sabemos que los valores de x,y,z son los de la iteración NMC/2 de precalentamiento. Hacemos el bucle de medida

      do i=1, Nmed
      	do j=1, NMC/Nmed
      		call movimiento(Rcil,Lcil,MCstep,x,y,z,accM)
      	end do

      	call medidaR(x,y,z,Nbinr,histoR,Rcil)
      	call medidaZ(x,y,z,Nbinz,histoZ,Lcil)
      
      end do
      
      print*, 'La tasa de aceptación en el proceso de medida es: ',
     & real(accM)/NMC * 100, '%'
     
C Normalizamos la distribución y la escribimos en un archivo
     
      r = Rcil-0.5d0
      volbinr = 1.d0
      volbinz = 1.d0
      
      open(5,FILE='datoshistR2_3.txt')
      open(6,FILE='datoshistZ2_3.txt')
      
      do i =1,Nbinz
      	volbinz = pi*r**2*Lcil/Nbinz
  
       	normalZ = histoZ(i)/(volbinz*Nmed)
     
     	write(5,*) Lcil*i/Nbinz-Lcil/(2.d0*Nbinz),normalZ
      end do
      
      do i=1,Nbinr
      	volbinr = Lcil*pi*((r*real(i)/NbinR)**2-
     &(r*real(i-1)/Nbinr)**2)
     
     	normalR = histoR(i)/(volbinr*Nmed)
     	write(6,*) r*i/Nbinr-r/(2.d0*Nbinr),normalR
      end do
       
      close(5)
      close(6)
      
      
      end program subproyecto2
      
      
      
      
      
      
C--------------------subrutina movimiento (realiza un movimiento de la particula)----------------


      subroutine movimiento(Rcil,Lcil,MCstep,x,y,z,acc)
      implicit none
      real*8 Rcil, Lcil, MCstep, x, y, z,RND, R, xf, yf, zf, r2
      integer*8 acc
     
      R = 0.5
      r2 = (Rcil - R)**2
     
C programamos un movimiento de la particula (suponiendola una esfera dura)
      
      	
      call random_number(RND)
      xf = x + (RND - R)*MCstep
      	
      call random_number(RND)
      yf = y +(RND - R)*MCstep
      	
      call random_number(RND)
      zf = z +(RND - R)*MCstep
      	
C corregimos zf:

      zf = zf - Lcil * floor(zf/Lcil)   
      
C Un if para ver si la nueva posición de la partícula está dentro del cilindro (teniendo en cuenta el radio de esta)        	
      	
      if(yf**2+xf**2 .LE. r2) then
      		
      	x = xf
      	y = yf
      	z = zf   		
      	acc = acc + 1
      	
      end if     
      
      return 
      end
      
      
      
C-------------------------------subrutina precalentamiento------------------------------------


      subroutine precalentamiento(Rcil,Lcil,NMC,pasoMC,x,y,z,acc)
      implicit none
      real*8 Rcil, Lcil, pasoMC, x, y, z,RND, R, xf, yf, zf, r2
      integer*8 NMC, acc
      integer i
      

      R = 0.5
      r2 = (Rcil - R)**2
C hacemos un bucle de 1 a NMC (intentos de movimientos) donde cada iteración sea un movimiento

      do i=1, NMC/2
      
	call movimiento (Rcil,Lcil,pasoMC,x,y,z,acc)
	
      end do
      
      print*, 'La tasa de aceptación en el precalentamiento es: ',
     &2*real(acc)/NMC * 100, '%'
          
      
      return 
      end
      
      
      
      
C-----------------------Subrutina de medidaR(x,y,z,Nbinr,histoR,Rcil)------------------------


      subroutine medidaR(x,y,z,Nbinr,histoR, Rcil)
      implicit none
      real*8 x,y,z,h,Rcil
      integer*8 Nbinr, histoR, i
      dimension histoR(1:300)
      
      h = (Rcil-0.5d0)/Nbinr

C Calculamos en que bin esta, le sumamos uno porque dimensionamos el histoR en 1 
      
      i = int(sqrt(x**2+y**2)/h) + 1
      
      histoR(i) = histoR(i) + 1

      return 
      end
      
      
C-------------------------Subrutina de medidaZ(x,y,z,Nbinz,histoZ,Lcil)-----------------------

      subroutine medidaZ(x,y,z,Nbinz,histoZ,Lcil)
      implicit none
      real*8 x,y,z,h,Lcil
      integer*8 Nbinz, histoZ, i
      dimension histoZ(1:300)
      
      h = Lcil/Nbinz
      
C Calculamos en que bin esta, le sumamos uno porque dimensionamos el histoZ en 1

      i = int(z/h) + 1

      histoZ(i) = histoZ(i) + 1
      
      
      
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

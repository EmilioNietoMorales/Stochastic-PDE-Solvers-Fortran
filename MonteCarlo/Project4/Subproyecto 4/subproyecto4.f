C-----------------------Empezamos el programa principal----------------------------

      program subproyecto4
      implicit none
      real*8 Rcil, Lcil, MCstep, RND, x, y, z,volbinz,volbinr,pi,r
      real*8 normalR, normalZ, dx, dy ,dz, ix, iy, iz 
      real* 8 phiRND
      integer*8 NMC, accM,accP, Nmed, histoR, histoZ, Nbinr, Nbinz
      integer i,j,Npart,k,p,sol
      character*10 nada
      parameter (pi=4.d0*atan(1.d0))
      dimension x(1:100), y(1:100), z(1:100)
      dimension histoR(1:50), histoZ(1:50)

      do i=1,50
       histoR(i) = 0
       histoZ(i) = 0
      end do
C comprobamos que los datos están bien introducidos (imprimiendolos en la terminal): 
      
      open(1,FILE='datos_in4.txt')
      
      read(1,*) Npart,nada,Rcil, nada, Lcil,nada, NMC,nada, Nmed
     &,nada,MCstep, nada, Nbinr, nada, Nbinz, nada

      close(1)
      
      print*, 'Radio cilindro = ',Rcil, 'Largo cilindro = ',Lcil,
     &'NMC = ',NMC,'Nmed = ',Nmed,'MCstep = ',MCstep,'Nbinr = ',Nbinr
     &,'Nbinz = ',Nbinz
     
C Llamamos a la subrutina precalentamiento para obtener la simulacion. Tomamos de punto inicial
C x=y=0, z=RND

      call init_random_seed()
      
      do i=1,100
      	x(i) = 0
      	y(i) = 0
      	z(i) = 0
      end do

C Aqui tenemos que colocar las particulas en posiciones iniciales donde no se solapen unas 
C con otras, ya que tenemos el potencial de esferas duras

      sol = 1
      r = Rcil - 0.5d0
      
C a cada particula le damos una posicion aleatroria y comprobamos que se verifiquen las condiciones de solapamiento   

      do i = 1, Npart
      	do while(sol == 1)
      		
      		sol = 0
      		
      		call random_number(phiRND)
      		ix = r*cos(2*pi*phiRND)
      	
      		call random_number(phiRND)
     		iy = r*sin(2*pi*phiRND)
      	
     		call random_number(RND)
     		iz = Lcil*RND
     		
C No haria falta que comprobemos que esta dentro del cilindro, ya que tomando coordenadas cilindricas, nos aseguramos que las posiciones estan dentro de este. Aún asi lo haremos (por costumbre ya :)  )
     	
     		if(ix**2 + iy**2 .GT. r**2) then
     			sol = 1
     		end if
      	
     		do j = 1, i-1
    			dx = ix - x(j)
     			dy = iy - y(j)
     			dz = iz - z(j)
     			dz = dz - Lcil*nint(dz/Lcil)
     			if (dx**2+dy**2+dz**2.LT.(2*0.5d0)**2) then
     				sol  = 1
     			end if
      		
     		end do
      	
     	end do
      	
     	if(sol == 0) then
     		x(i) = ix
     		y(i) = iy
     		z(i) = iz
     		sol = 1
C     		print*, 'una particula colocada'
     	end if      	
      
      end do   
      
      
      print*, 'Terminada la colocacion de particuas iniciales'
      
      	
      accP=0
      accM = 0
      
C      print*, 'Llegamos aqui'

      call precalentamiento(Npart,Rcil,Lcil,NMC,MCstep,x,y,z,accP)    
      

      
C una vez finalizado el precalentamiento, sabemos que los valores de x,y,z son los de la iteración NMC/2 de precalentamiento. Hacemos el bucle de medida

      do i=1, Nmed
      	do j=1, NMC/Nmed
      		do k=1,Npart
      		   call movimiento(Npart,Rcil,Lcil,MCstep,x,y,z,accM)
      		end do
      	end do

      	call medidaR(Npart,x,y,z,Nbinr,histoR,Rcil)
      	call medidaZ(Npart,x,y,z,Nbinz,histoZ,Lcil)
 
      end do
      
      print*, 'La tasa de aceptación en el proceso de medida es: ',
     &real(accM)/NMC/Npart * 100,'%'
     
     
     
C Normalizamos la distribución y la escribimos en un archivo
      

      volbinr = 1.d0
      volbinz = 1.d0
   
      
C Escribimos los datos de los histogramas en un fichero     
      
      open(5,FILE='datoshistZ4_2.txt')
      open(6,FILE='datoshistR4_2.txt')
      
      do i=1,Nbinz
      
      	volbinz = pi*r**2*Lcil/Nbinz
      	normalZ = histoZ(i)/(volbinz*Nmed)
        
      	write(5,*) Lcil*i/Nbinz-Lcil/(2.d0*Nbinz),normalZ

      end do
      
      do i=1,Nbinr
      
      	volbinr = Lcil*pi*((r*real(i)/NbinR)**2-
     &(r*real(i-1)/Nbinr)**2)
     
     	normalR= histoR(i)/(volbinr*Nmed)
     	
     	write(6,*) r*i/Nbinr-r/(2.d0*Nbinr), normalR
      end do
      
      close(6)
      close(5)
      
      
      
      
      
      end program subproyecto4
      
      
      
      
      
      
C--------------------subrutina movimiento (realiza un movimiento de particula al azar)----------------

C Tenemos que hacer que el movimiento sea de una particula aleatoria

      subroutine movimiento(Npart,Rcil,Lcil,MCstep,x,y,z,acc)
      implicit none
      real*8 Rcil, Lcil, MCstep, x, y, z,RND, R, xf, yf, zf, r2
      real*8 dx,dy,dz,D
      integer*8 acc
      integer Npart,sol,j,p
      dimension x(1:100), y(1:100), z(1:100)
     
      R = 0.5
      r2 = (Rcil - R)**2
      D = 2*R
     
C programamos un movimiento de la particula (suponiendola una esfera dura)
C Elegimos una particula al azar y calculamos su posicion final:

      call random_number(RND)
      p = int(RND*Npart) + 1
      
      sol = 0
      	
      call random_number(RND)
      xf = x(p) + (RND - R)*MCstep
      	
      call random_number(RND)
      yf = y(p) +(RND - R)*MCstep
      	
      call random_number(RND)
      zf = z(p) +(RND - R)*MCstep
      	
C corregimos zf:

      zf = zf - Lcil * floor(zf/Lcil)   
      
C Un if para ver si la nueva posición de la partícula está dentro del cilindro (teniendo en cuenta el radio de esta)

      if(xf**2 + yf**2 .GT. r2) then
      	sol = 1
      end if   
      
C Comprobamos ahora que la nueva posicion de la particula aleatoria cumple que la distancia a cada una de las otras particulas es mayor que el diametro de una particula  y corregimos las distancias solo en z (ya que en x e y no tenemos condiciones de contorno periodicas)
      
      do j=1,Npart
      	if(j .NE. p) then
      		dx = xf - x(j)     		
      		dy = yf - y(j)
      		dz = zf - z(j)
      		dz = dz - Lcil*nint(dz/Lcil)
      		if(dx**2 +dy**2 + dz**2 .LT. D**2) then
      			sol = 1
      		end if 
      	end if     
      end do  
      
C Si verifica las condiciones anteriores (que no se salga del cilindro y no solape con ninguna otra particula) aceptamos el movimiento       	
      	
      if(sol == 0) then
      		
      	x(p) = xf
      	y(p) = yf
      	z(p) = zf   		
      	acc = acc + 1
      	
      end if     
      
      return 
      end
      
      
      
C-------------------------------subrutina precalentamiento------------------------------------


      subroutine precalentamiento(Npart,Rcil,Lcil,NMC,pasoMC,x,y,z,acc)
      implicit none
      real*8 Rcil, Lcil, pasoMC, x, y, z,RND
      integer*8 NMC, acc
      integer i,Npart,j,k
      dimension x(1:100), y(1:100), z(1:100)
      
C hacemos un bucle de 1 a NMC (intentos de movimientos) donde cada iteración sea un intento MC

      do i=1, NMC/2
      	do j=1,Npart
		call movimiento (Npart,Rcil,Lcil,pasoMC,x,y,z,acc)
	end do
      end do
      
      print*, 'La tasa de aceptación en el precalentamiento es: ',
     &2*real(acc)/NMC/Npart * 100, '%'
          
      
      return 
      end
      
      
      
      
C-----------------------Subrutina de medidaR(x,y,z,Nbinr,histoR,Rcil)------------------------


      subroutine medidaR(Npart,x,y,z,Nbinr,histoR, Rcil)
      implicit none
      real*8 x,y,z,h,Rcil
      integer*8 Nbinr, histoR
      integer i,j,Npart
      dimension histoR(1:50)
      dimension x(1:100), y(1:100), z(1:100)

      
      h = (Rcil-0.5d0)/Nbinr

C Calculamos en que bin esta cada particula, le sumamos uno porque dimensionamos el histoR en 1 
      
      do j=1,Npart
      
      	i = int(sqrt(x(j)**2+y(j)**2)/h) + 1
      
      	histoR(i) = histoR(i) + 1

      end do
      
      return 
      end
      
      
C-------------------------Subrutina de medidaZ(x,y,z,Nbinz,histoZ,Lcil)-----------------------

      subroutine medidaZ(Npart,x,y,z,Nbinz,histoZ,Lcil)
      implicit none
      real*8 x,y,z,h,Lcil
      integer*8 Nbinz, histoZ
      integer i,j,Npart
      dimension histoZ(1:50)
      dimension x(1:100), y(1:100), z(1:100)
      
      h = Lcil/Nbinz
      
C Calculamos en que bin esta, le sumamos uno porque dimensionamos el histoZ en 1

      do j=1,Npart
      
      	i = int(z(j)/h) + 1

      	histoZ(i) = histoZ(i) + 1
      
      end do
      
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


      PROGRAM main
     
      implicit none
      double precision, parameter :: pi=4.0*datan(1.0d0)
      double precision dt , ttotal, dx, D, X,alpha
      integer, parameter :: N=100  ! Número de celdas
      integer i,k,nitmax,nskip,j, info, ipiv(0:N) ! contadores espacial y temporal, número máximo de iteraciones, nskip
      character*10 nada
      
      double precision T(0:N),Told(0:N)  ! Told(i)=u_i^n, T(i)=u_i^(n+1)
      complex*16 A(0:N,0:N),b(0:N,1)
      
C Creamos tres vectores Tmax (0:N) Tmin(0:N) y Tmed(0:N) para almacenar los valores max min y med de cada profundidad:

      double precision Tmax(0:N), Tmin(0:N), Tmed(0:N)

C     Leemos parámetros de un fichero: intervalo, constante de difusión, paso temporal, tiempo, nskip
      open(12,FILE='param.dat')
      read(12,*) nada,dt,nada,ttotal,nada,nitmax,nada,nskip,nada,
     &X,nada,D,nada,dx
      close(12)
      
      alpha = D*dt/(dx**2)

C     Definimos algunos parámetros por comodidad: longitud del intervalo, paso espacial, alpha, número de pasos en el tiempo, ...
          
C     T(x,t=0): especificamos el estado inicial

      do i=0,N
        Told(i) = 10
      end do
      
C inicializamos los valores de Tmax, Tmin, para compararlos despues con cada paso temporal (si 
C sube o baja) y cambiarlos.
      
      do i =0,N
      	Tmax(i) = Told(i)
      	Tmin(i) = Told(i)    
      end do


c     Inicializamos la matriz de coeficientes A_ij = 0
      do j=0,N
        do i=0,N
          A(i,j)=0.0d0
        end do
      end do
      
c     Inicializamos el vector de coeficientes b_i = 0
      do i=0,N
          b(i,1)=0.0d0
      end do
      
           
C     Escribimos el estado inicial
      open (5,FILE='T1.dat')
      open (7, FILE='T2.dat')
      open (9, FILE = 'T3.dat')
      
      do i=0,N
        write(5,*) .0, i*dx, Told(i) ! 3 columnas: t, x, u
      end do

C     Bucle principal
      do k=1,nitmax ! bucle en el tiempo
      
c     Especificamos los coeficientes no nulos de A y B  
      	do i=1,N-1
      	  A(i,i)= 1 + alpha 
       	  A(i,i-1)= -alpha/2 
       	  A(i,i+1)= -alpha/2 
          b(i,1)=alpha/2*(Told(i-1)+Told(i+1))+(1-alpha)*Told(i)
      	end do

c     Coeficientes para imponer las condiciones de contorno de neumann
      	A(0,0)= 1.0d0
      	A(0,1)= 0.0d0
      	b(0,1)= 10 + 25 * sin(2*pi*k*dt) 
      	A(N,N)= 1.0d0 
      	A(N,N-1)= -1.0d0
      	b(N,1)= 0.0d0 
      	
      	call ZGESV(N+1,1,A,N+1,ipiv,b,N+1,info)
      	
        do i=0,N ! bucle en el espacio
          T(i)=b(i,1) ! avanzamos la solución en el tiempo
        end do
        
C Obtenemos la suma de las temperaturas para cada profundidad y en cada iteracion temporal:
	
	do i = 0, N
		Tmed(i)  = Tmed(i) + T(i)
	end do

        do i=0,N
          Told(i)= T(i)  ! actualizamos Told para la siguiente iteración
        end do
        if (mod(k,nskip).eq.0) then ! esto hay que entenderlo bien
          do i=0,N
            write(5,*) k * dt,i*dx,T(i) ! 3 columnas: t, x, T
          end do
          
C Escribimos las evoluciones temporales de T a las temperaturas x=0 (i = 0), x=0.6 (i=15)
C x= 1 (i=25) y a x = 2 (i = 50)
	
	  write (7,*) k * dt, T(0) , T(15), T(25), T(50)	
	
C comparamos Tmax y Tmin con el valor actual de T para ver si ha aumentado o disminuido
	
	  do i = 0, N
		if(T(i).LT.Tmin(i)) then
			Tmin(i) = T(i)
		end if
		
		if (T(i).GT.Tmax(i)) then
			Tmax(i) = T(i)
		end if
		
	  end do
	          
        end if      
      end do
      
      do i = 0, N
  
C Vamos a obtener la profundidad minima a la que colocar las tuberias para que el agua no se 
C congele, para ello obtenemos el valor de i para el cual Tmin este entre 0 y 0,5 ya que es un
C rango aceptable teniendo en cuenta que solo tenemos solo 100 puntos y la profundidad minima sera
C entonces, (i) * dx, ya que para profundidades mayores la temperatura será mayor que 0
	if (Tmin(i).GT.0 .AND. Tmin(i).LT.0.5) then
	  print*, 'La profundidad minima para que no congele es: ', 
     &(i) * dx
  
	end if 
      	Tmed(i) = Tmed(i)/nitmax
	write(9,*) i * dx, Tmax(i), Tmin(i), Tmed(i)     
      end do

      close (5)   
      close(7)  
      close (9)
      END



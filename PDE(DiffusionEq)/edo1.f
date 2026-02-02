c     Programa para obtener la solución de una ODE con condiciones de contorno de Dirichlet, Neumann o mixtas
c     linf, lsup : límites inferior y superior del intervalo
c     Nodos en los extremos de las celdas
c     (i=0,N) -> (x_0=linf,x_N=lsup) : puntos de frontera
c     A(i,j) -> Matriz de los coeficientes

      program MAIN
      integer, parameter :: N=200				
      double precision, parameter :: pi=4*atan(1.0d0)
      integer :: i,info
      double precision :: linf,lsup,x,h,yexact,error
      double precision :: A(0:N,0:N),b(0:N,1),y(0:N),ipiv(N+1)

      linf=0.0d0
      lsup=2.0d0
      h=(lsup-linf)/N ! tamaño de celda

c     Inicializamos la matriz de coeficientes A_ij = 0
      do j=0,N
        do i=0,N
          A(i,j)=0.0d0
        end do
      end do

c     Especificamos los coeficientes de no nulos de A y el vector b  
      do i=1,N-1
        x=linf+i*h ! relación entre  i (variable discreta) y x (variable continua)
        A(i,i)=... ! completar
        A(i,i-1)=... ! completar
        A(i,i+1)=... ! completar
        b(i,1)=... ! completar
      end do

c     Coeficientes para imponer las condiciones de contorno
      A(0,0)=... ! completar
      b(0,1)=... ! completar
      A(N,N)=... ! completar
      b(N,1)=... ! completar

      call DGESV(N+1,1,A,N+1,ipiv,b,N+1,info)

c     Escribimos la solución      
      open (1,file='sol1.dat')

      do i=0,N
        x=linf+i*h
        write(1,*) x,b(i,1)
      end do

      close(1)

      end program


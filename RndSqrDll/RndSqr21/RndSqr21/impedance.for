      subroutine impedance(n,dC,Um,em,dZ0)
! Результирующая матрица импедансов (Ом)
      implicit complex*16(c,w,z), real*8(a-b,d-h,o-v,x-y)
      dimension dZ0(n,n), dC(n*n), Um(n*n), em(n)
      dimension dV(n,n), dCn(n,n), dUn(n,n), dEmn(n)
      dimension dCe(n,n), dUe(n,n)
      k=1
      do i=1,n
          do j=1,n
              dCn(i,j) = dC(k)
              k=k+1
          end do
      end do
      k=1
      do i=1,n
          do j=1,n
              dUn(i,j) = Um(k)
              k=k+1
          end do
      end do
      do i=1,n
          dEmn(i) = em(i)
      end do
      dUe=0
      dV=0
      c=2.99792458
      do i=1,n
          dV(i,i) = sqrt(dEmn(i))/c
      end do
      dZ0=dCn
      call matr(dCn, dCe, n)
      dCn=dZ0
      dZ0=dUn
      call matr(dUn, dUe, n)
      dUn=dZ0
      dZ0=matmul(matmul(matmul(dUn, dV), dUe), dCe)*10e3
      return
      end
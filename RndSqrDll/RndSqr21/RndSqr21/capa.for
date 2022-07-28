      subroutine capa(bC,dC,er,n)
! Результирующая матрица емкостей (пФ/м)
      real*8 bC(1), dC(1), er
	do 1 i=1,n
	do 1 j=1,n
1	dC((j-1)*n+i) = 8.854 * ( bC((j-1)*n+i) + (er-1.) * dC((j-1)*n+i) )
      return
	end
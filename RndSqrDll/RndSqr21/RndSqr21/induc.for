      subroutine induc(bC, dL, n)
! –езультирующа€ матрица индуктивностей (мк√н/м)
      real*8 bC(1), dL(1), d
	pi = dacos(-1.d0)
	do 1 i=1,n
	do 1 j=1,n
1	dL((j-1)*n+i) = bC((j-1)*n+i)
	call dminv(dL,n,d)
	do 2 i=1,n
	do 2 j=1,n
2	dL((j-1)*n+i) = 0.4 * pi * dL((j-1)*n+i)
      return
	end
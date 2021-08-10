subroutine cal_N_3d(kesi,yita,zeta,N)
!This function calculates N, dNdkesi, J and the determinant of Jacobian matrix.
implicit none

real(kind=8),intent(in)::kesi,yita,zeta
real(kind=8),intent(out):: N(3,24)
real(kind=8) N1,N2,N3,N4,N5,N6,N7,N8
real(kind=8) o        !0
real(kind=8) ONE      !1
real(kind=8) EIG      !8

o   = 0.0D0
ONE = 1.0D0
EIG = 8.0D0

!Calculate N.
N1 = (ONE-kesi)*(ONE-yita)*(ONE-zeta)/EIG    !kesi=-1,yita=-1,zeta=-1
N2 = (ONE+kesi)*(ONE-yita)*(ONE-zeta)/EIG    !kesi= 1,yita=-1,zeta=-1
N3 = (ONE+kesi)*(ONE+yita)*(ONE-zeta)/EIG    !kesi= 1,yita= 1,zeta=-1
N4 = (ONE-kesi)*(ONE+yita)*(ONE-zeta)/EIG    !kesi=-1,yita= 1,zeta=-1
N5 = (ONE-kesi)*(ONE-yita)*(ONE+zeta)/EIG    !kesi=-1,yita=-1,zeta= 1
N6 = (ONE+kesi)*(ONE-yita)*(ONE+zeta)/EIG    !kesi= 1,yita=-1,zeta= 1
N7 = (ONE+kesi)*(ONE+yita)*(ONE+zeta)/EIG    !kesi= 1,yita= 1,zeta= 1
N8 = (ONE-kesi)*(ONE+yita)*(ONE+zeta)/EIG    !kesi=-1,yita= 1,zeta= 1

N(1,1:24) =[N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o,o]
N(2,1:24) =[o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8,o]
N(3,1:24) =[o,o,N1,o,o,N2,o,o,N3,o,o,N4,o,o,N5,o,o,N6,o,o,N7,o,o,N8]

return 
end SUBROUTINE cal_N_3d               

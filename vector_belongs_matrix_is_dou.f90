subroutine vector_belongs_matrix_is_dou(m,n,Matrix,Vector,Location,Yes)   
!If a Vector belongs to a Matrix(m,n).
!Author:  Fang Shi
!Email:   shifang@hyit.edu.cn       
!Website: http://phipsi.top, copyright(c) 2016-2021   
!Date:    2017-11-04.

integer i,m,n              
real(kind=8) Matrix(m,n)
real(kind=8) Vector(n),Vector_tem(n)
integer,intent(out)::Location
logical,intent(out)::Yes

Yes = .False.

do i=1,m
  Vector_tem = Matrix(i,1:n)
  if ((abs(MaxVal(Vector_tem-Vector)) < 1.0D-13).and. &        
      (abs(MinVal(Vector_tem-Vector)) < 1.0D-13)) then
      Yes = .True.
      Location = i
      exit
  end if
end do

return 
end subroutine vector_belongs_matrix_is_dou
    



      subroutine Tool_Read_File(Full_Name,m,n,Temp_DATA,Flag_Blank)
C     Read file.
      implicit none

      character*200 Full_Name
      integer i,j,m,n,istat
      logical Flag_Blank
      real(kind=8) Temp_DATA(m,n)

      open(112,file=Full_Name,status='old')
      Read(112,*,IOSTAT=istat)
      close(112)
      if ( istat /= 0 ) then   !blank
          Flag_Blank = .True.
      else                     !non-blank
          Flag_Blank = .False.
          open(122,file=Full_Name,status='old')
          read(122,*)((Temp_DATA(i,j),j=1,n),i=1,m)
          close(122)
      end if
      
      RETURN
      END SUBROUTINE Tool_Read_File
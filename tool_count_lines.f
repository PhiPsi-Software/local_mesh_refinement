      INTEGER FUNCTION tool_count_lines(Full_Name)
      !count the line number of a file.
      implicit none
      real line
      integer i
      character*200 Full_Name
      i=0  
      open(10,file=Full_Name,status='old')
      do while (.true.)
          read(10,*,end=100) line
          i=i+1
      enddo
  100 continue
      close(10)
      tool_count_lines = i
      
      RETURN
      END 
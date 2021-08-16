program local_grid_refinement
!Refine XFEM-enriched 3d hexahedral elements.
!author:  Fang Shi
!email:   shifang@hyit.edu.cn       
!website: http://phipsi.top, copyright(c) 2016-2021   
!date:    2021-08-09.

implicit none
integer i_e,i_n
integer n1,n2,n3,n4,n5,n6,n7,n8
integer nn(8),nn_l(9)
real(kind=8) x_nodes_l(9),y_nodes_l(9),z_nodes_l(9)
real(kind=8) x_nodes(8),y_nodes(8),z_nodes(8)
integer ele_enr_nd(1000,8)           
integer ele_enr_nodes_num(1000)         
integer c_added_node                        
integer c_added_ele                        
real(kind=8) c_node_coor(3)
real(kind=8) c_kesi,c_yita,c_eta
integer c_location
logical yes_exist
integer temp_added_nodes(100)
integer i_added_n
integer temp_elem_node(100000,8)              
integer temp_elem_node_old(100000,8)        
integer temp_elem_mat(100000)             
real(kind=8) temp_coor(100000,3)          
integer i
character (len=200) c_file_name_1,c_file_name_2,c_file_name_3
integer refine_template,edge_case,face_case
real(kind=8) added_node_coor(100,3),added_node_coor_14(20,3)
real(kind=8) added_node_coor_4(7,3)
real(kind=8) added_node_coor_5678(40,3)
integer new_orig_node(8)
integer j_ele
logical refined_elements(100000)
integer wanted_elem(100000,1:8),wanted_elem_mat(100000)
integer num_wanted_elem,i_check
integer vertex_case
real(kind=8),allocatable ::coor(:,:)
integer,allocatable ::elem_node(:,:),elem_mat(:),enriched_node_flag(:)
real(kind=8),allocatable ::temp_data(:,:)
real(kind=8),allocatable ::g_x_nodes(:,:),g_y_nodes(:,:),g_z_nodes(:,:)
integer num_node,num_elem
real(kind=8) one,two,thr,zr
integer tool_count_lines
logical flag_blank

!**************************
!   set the input files.
!**************************
c_file_name_1 = 'input.node'    !node x,y,z
c_file_name_2 = 'input.elem'    !element node number.
c_file_name_3 = 'input.ennd'    !enriched nodes.


!**************************
!    read *.node file.
!**************************
print *, "    Read node file..." 
num_node = tool_count_lines(c_file_name_1) !get number of nodes.
print *, "    Number of nodes:",num_node
allocate(coor(num_node,3))
allocate(temp_data(num_node,3))
call tool_read_file(c_file_name_1,num_node,3,temp_data,flag_blank)
coor = temp_data
deallocate(temp_data)

!**************************
!    read *.elem file.
!**************************
print *, "    Read element file..." 
num_elem = tool_count_lines(c_file_name_2) !get number of elements.
print *, "    Number of elements:",num_elem
allocate(elem_node(num_elem,8))
allocate(elem_mat(num_elem)) 
allocate(temp_data(num_elem,9))
call tool_read_file(c_file_name_2,num_elem,9,temp_data,flag_blank)
elem_node = int(temp_data(:,1:8))
elem_mat  = int(temp_data(:,9))
deallocate(temp_data)
allocate(g_x_nodes(num_elem,8))
allocate(g_y_nodes(num_elem,8))
allocate(g_z_nodes(num_elem,8))
do i=1,num_elem
  n1  = elem_node(i,1)                                                
  n2  = elem_node(i,2)                                              
  n3  = elem_node(i,3)                                             
  n4  = elem_node(i,4)            
  n5  = elem_node(i,5)   
  n6  = elem_node(i,6)   
  n7  = elem_node(i,7)   
  n8  = elem_node(i,8)   
  nn  = [n1,n2,n3,n4,n5,n6,n7,n8]                                                
  g_x_nodes(i,1:8) = coor(nn,1)
  g_y_nodes(i,1:8) = coor(nn,2) 
  g_z_nodes(i,1:8) = coor(nn,3) 
end do

!**************************
!    Read *.enrn file.
!**************************
print *, "    Read enrn file..." 
allocate(enriched_node_flag(num_node))
enriched_node_flag(1:num_node) = 0
allocate( temp_data(num_node,1))
call tool_read_file(c_file_name_3,num_node,1,temp_data,flag_blank)
enriched_node_flag  = int(temp_data(:,1))
deallocate(temp_data)
print *, "    Number of enriched nodes:",sum(enriched_node_flag)

!**************************
!   initialization 
!**************************
one = 1.0d0
two = 2.0d0
thr = 3.0d0
zr  = 0.0d0
ele_enr_nd(1:num_elem,1:8)          = 0   
ele_enr_nodes_num(1:num_elem)       = 0    
c_added_node                        = 0    
temp_coor(1:num_node*200,1:3)       = zr  
c_added_ele                         = 0    
temp_elem_node(1:num_elem*200, 1:8) = 0   
temp_elem_mat(1:num_elem*200)       = 0
refined_elements(1:num_elem*200)    = .false.    !used to mark the refined elements
c_added_node  = c_added_node + num_node 
c_added_ele   = c_added_ele  + num_elem

temp_elem_node(1:num_elem,1:8) = elem_node(1:num_elem,1:8)
temp_elem_mat(1:num_elem)      = elem_mat(1:num_elem)
temp_coor(1:num_node,1:3)      = coor(1:num_node,1:3)


!**************************
!  loop over elements.
!**************************
do i_e = 1,num_elem
  n1  = elem_node(i_e,1)
  n2  = elem_node(i_e,2)
  n3  = elem_node(i_e,3)
  n4  = elem_node(i_e,4)
  n5  = elem_node(i_e,5)
  n6  = elem_node(i_e,6)
  n7  = elem_node(i_e,7)
  n8  = elem_node(i_e,8)
  nn  = [n1,n2,n3,n4,n5,n6,n7,n8]
  !coordinates of nodes of the current element
  x_nodes = g_x_nodes(i_e,1:8)
  y_nodes = g_y_nodes(i_e,1:8)
  z_nodes = g_z_nodes(i_e,1:8)
  nn_l= [n1,n2,n3,n4,n5,n6,n7,n8,n1]
  x_nodes_l = coor(nn_l,1)
  y_nodes_l = coor(nn_l,2)
  z_nodes_l = coor(nn_l,3)
  !---------------------------------------------------
  ! determine the enriched nodes of each element.
  !---------------------------------------------------
  do i_n =1,8
      if(enriched_node_flag(nn(i_n))==1) then 
        ele_enr_nd(i_e,i_n) = 1
      endif
  end do
  
  ele_enr_nodes_num(i_e) = sum(ele_enr_nd(i_e,1:8))
  
  !mark the refined element.
  if (ele_enr_nodes_num(i_e)>=1) then
    refined_elements(i_e) = .true.
  end if
  !-----------------------------------------------------------
  ! determine refine template: refine_template=1,2,4,8.
  !-----------------------------------------------------------
  refine_template  = 0
  edge_case        = 0
  face_case        = 0
  !////////////////////
  ! 1 enriched nodes.
  !////////////////////
  if (ele_enr_nodes_num(i_e)==1)then
    refine_template = 1
    if(ele_enr_nd(i_e,4)==1) vertex_case =1
    if(ele_enr_nd(i_e,2)==1) vertex_case =2
    if(ele_enr_nd(i_e,1)==1) vertex_case =3
    if(ele_enr_nd(i_e,3)==1) vertex_case =4
    if(ele_enr_nd(i_e,6)==1) vertex_case =5
    if(ele_enr_nd(i_e,8)==1) vertex_case =6
    if(ele_enr_nd(i_e,5)==1) vertex_case =7
    if(ele_enr_nd(i_e,7)==1) vertex_case =8
  end if
  !////////////////////
  ! 2 enriched nodes.
  !////////////////////
  if (ele_enr_nodes_num(i_e)==2)then
    !if 2 nodes are on the same edge.
    if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,4)==1) then
      refine_template = 2;edge_case =1
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1)then
      refine_template = 2;edge_case =2
    else if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,8)==1)then
      refine_template = 2;edge_case =3
    else if(ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,7)==1)then
      refine_template = 2;edge_case =4
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1)then
      refine_template = 2;edge_case =5
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1)then
      refine_template = 2;edge_case =6
    else if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,6)==1)then
      refine_template = 2;edge_case =7
    else if(ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,8)==1)then
      refine_template = 2;edge_case =8
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,5)==1)then
      refine_template = 2;edge_case =9
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,6)==1)then
      refine_template = 2;edge_case =10
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,7)==1)then
      refine_template = 2;edge_case =11
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,8)==1)then
      refine_template = 2;edge_case =12
    !if 2 nodes are diagonal on a face.
    else if(ele_enr_nd(i_e,8)==1 .and. ele_enr_nd(i_e,6)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,7)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,2)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,3)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,7)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,6)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,5)==1)then
      refine_template = 4
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,8)==1)then
      refine_template = 4
    !other cases.
    else
      refine_template = 8
    end if
  end if
  !////////////////////
  ! 3 enriched nodes.
  !////////////////////
  if (ele_enr_nodes_num(i_e)==3)then
    if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,7)==1) then
      refine_template = 4; face_case =1
    else if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =1
    else if(ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =1
    else if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =1
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1) then
      refine_template = 4; face_case =2
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1) then
      refine_template = 4; face_case =2
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,4)==1) then
      refine_template = 4; face_case =2
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1) then
      refine_template = 4; face_case =2
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,7)==1) then
      refine_template = 4; face_case =3
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,6)==1) then
      refine_template = 4; face_case =3
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,6)==1) then
      refine_template = 4; face_case =3
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,6)==1) then
      refine_template = 4; face_case =3
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,6)==1) then
      refine_template = 4; face_case =4
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,5)==1) then
      refine_template = 4; face_case =4
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,5)==1) then
      refine_template = 4; face_case =4
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,5)==1) then
      refine_template = 4; face_case =4
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,5)==1) then
      refine_template = 4; face_case =5
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =5
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =5
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =5
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4; face_case =6
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,7)==1) then
      refine_template = 4; face_case =6
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,8)==1 .and. ele_enr_nd(i_e,7)==1) then
      refine_template = 4; face_case =6
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,8)==1 .and. ele_enr_nd(i_e,7)==1) then
      refine_template = 4; face_case =6
    !other cases.
    else
      refine_template = 8
    end if
  end if
  !////////////////////
  !  4 enriched nodes.
  !////////////////////
  if (ele_enr_nodes_num(i_e)==4)then
    !if 4 nodes are on the same face.
    if(ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,8)==1) then
      refine_template = 4;face_case =1
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1)then
      refine_template = 4;face_case =2
    else if(ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,7)==1 .and. ele_enr_nd(i_e,6)==1)then
      refine_template = 4;face_case =3
    else if(ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,2)==1 .and. ele_enr_nd(i_e,6)==1 .and. ele_enr_nd(i_e,5)==1)then
      refine_template = 4;face_case =4
    else if(ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,1)==1 .and. ele_enr_nd(i_e,5)==1 .and. ele_enr_nd(i_e,8)==1)then
      refine_template = 4;face_case =5
    else if(ele_enr_nd(i_e,3)==1 .and. ele_enr_nd(i_e,4)==1 .and. ele_enr_nd(i_e,8)==1 .and. ele_enr_nd(i_e,7)==1)then
      refine_template = 4;face_case =6
    !other cases.
    else
      refine_template = 8     
    end if
  end if
  !//////////////////////////////////////////
  !other cases:more than 4 enriched nodes.
  !//////////////////////////////////////////
  if (ele_enr_nodes_num(i_e)>=5)then
    refine_template = 8
  end if
  
  !print to screen.
  if (refine_template>0)then
     print *,'    > Element ',i_e,' refine template is',refine_template
  endif
  
  !---------------------------------------
  ! refine according to refine_template.
  !---------------------------------------
  added_node_coor(1:100,1:3) = zr
  select case (refine_template)
  
  case (1)
  !///////////////////////////////////////////////////
  !                 case  1: refine_template=1
  !!///////////////////////////////////////////////////
  added_node_coor_4(1,1:3)  = [-one,      one/thr,       -one];
  added_node_coor_4(2,1:3)  = [-one/thr,  one/thr,       -one];
  added_node_coor_4(3,1:3)  = [-one/thr,  one,           -one];
  added_node_coor_4(4,1:3)  = [-one,      one/thr,   -one/thr];
  added_node_coor_4(5,1:3)  = [-one/thr,  one/thr,   -one/thr];
  added_node_coor_4(6,1:3)  = [-one/thr,  one,       -one/thr];
  added_node_coor_4(7,1:3)  = [-one,      one,       -one/thr];
  if (vertex_case ==1) then
    added_node_coor(1:7,1) =  added_node_coor_4(1:7,1)
    added_node_coor(1:7,2) =  added_node_coor_4(1:7,2)
    added_node_coor(1:7,3) =  added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [1,2,3,4,5,6,7,8]
  else if(vertex_case ==2)then
    added_node_coor(1:7,1) =  added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =  added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =  added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [1,4,3,2,5,8,7,6] 
  else if(vertex_case ==3)then
    added_node_coor(1:7,1) =  -added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =  added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =  added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [2,3,4,1,6,7,8,5]
  else if(vertex_case ==4)then
    added_node_coor(1:7,1) =  added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =  -added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =  added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [4,1,2,3,8,5,6,7]
  else if(vertex_case ==5)then
    added_node_coor(1:7,1) =   added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =   added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =   -added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [5,8,7,6,1,4,3,2]
  else if(vertex_case ==6)then
    added_node_coor(1:7,1) =   added_node_coor_4(1:7,1)
    added_node_coor(1:7,2) =   added_node_coor_4(1:7,2)
    added_node_coor(1:7,3) =  -added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [5,6,7,8,1,2,3,4]
  else if(vertex_case ==7)then
    added_node_coor(1:7,1) =  -added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =   added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =  -added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [6,7,8,5,2,3,4,1]
  else if(vertex_case ==8)then
    added_node_coor(1:7,1) =   added_node_coor_4(1:7,2)
    added_node_coor(1:7,2) =  -added_node_coor_4(1:7,1)
    added_node_coor(1:7,3) =  -added_node_coor_4(1:7,3)
    new_orig_node(1:8) = [8,5,6,7,4,1,2,3]
  end if  
  !add nodes.
  do i_added_n = 1,7
    c_kesi = added_node_coor(i_added_n,1)
    c_yita = added_node_coor(i_added_n,2)
    c_eta  = added_node_coor(i_added_n,3)
    !calculate global coordinate system according to the local coordinate system.
    call cal_coor_by_kesiyita_3d(c_kesi,c_yita,c_eta,x_nodes,y_nodes,z_nodes, c_node_coor(1),c_node_coor(2),c_node_coor(3))
    !check if the newly added node already exist.
    call vector_belongs_matrix_is_dou(c_added_node,3,temp_coor(1:c_added_node,1:3), c_node_coor(1:3),c_location,yes_exist)
    !if not exist before then add it. 
    if(yes_exist .eqv. .false.)then
      c_added_node = c_added_node + 1
      temp_coor(c_added_node,1:3)=c_node_coor(1:3)
      temp_added_nodes(i_added_n) = c_added_node
    !if exists before.
    else
      temp_added_nodes(i_added_n) = c_location
    end if
  end do  
  !add elements.
  !element 1
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(1))
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(2))
  temp_elem_node(c_added_ele,3) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(5))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,7) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(4)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 2
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(2))
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(3))
  temp_elem_node(c_added_ele,3) = temp_added_nodes(3)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,7) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(5)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 3
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(7)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(5))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,7) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,8) = nn(new_orig_node(8))
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 4
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(3)
  temp_elem_node(c_added_ele,4) = nn(new_orig_node(4))
  temp_elem_node(c_added_ele,5) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(7)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !Make sure that the determinant detJ is positive 
  if (vertex_case == 2 .or. vertex_case == 6 .or. vertex_case == 7 .or. vertex_case == 8) then
    temp_elem_node_old = temp_elem_node
    temp_elem_node(c_added_ele-3:c_added_ele,1) = temp_elem_node_old(c_added_ele-3:c_added_ele,4)
    temp_elem_node(c_added_ele-3:c_added_ele,2) = temp_elem_node_old(c_added_ele-3:c_added_ele,3)
    temp_elem_node(c_added_ele-3:c_added_ele,3) = temp_elem_node_old(c_added_ele-3:c_added_ele,2)
    temp_elem_node(c_added_ele-3:c_added_ele,4) = temp_elem_node_old(c_added_ele-3:c_added_ele,1)
    temp_elem_node(c_added_ele-3:c_added_ele,5) = temp_elem_node_old(c_added_ele-3:c_added_ele,8)
    temp_elem_node(c_added_ele-3:c_added_ele,6) = temp_elem_node_old(c_added_ele-3:c_added_ele,7)
    temp_elem_node(c_added_ele-3:c_added_ele,7) = temp_elem_node_old(c_added_ele-3:c_added_ele,6)
    temp_elem_node(c_added_ele-3:c_added_ele,8) = temp_elem_node_old(c_added_ele-3:c_added_ele,5)
  end if  
  case (2)
  !///////////////////////////////////////////////////
  !                 case  2: refine_template=2
  !///////////////////////////////////////////////////
  added_node_coor_14(1,1:3)  = [-one,     -one/thr,     -one];
  added_node_coor_14(2,1:3)  = [-one,      one/thr,     -one];
  added_node_coor_14(3,1:3)  = [-one/thr, -one,         -one];
  added_node_coor_14(4,1:3)  = [-one/thr, -one/thr,     -one];
  added_node_coor_14(5,1:3)  = [-one/thr,  one/thr,     -one];
  added_node_coor_14(6,1:3)  = [-one/thr,  one,         -one];
  added_node_coor_14(7,1:3)  = [ one/two, -one/thr,     -one];
  added_node_coor_14(8,1:3)  = [ one/two,  one/thr,     -one];
  added_node_coor_14(9,1:3)  = [-one,     -one,     -one/thr];
  added_node_coor_14(10,1:3) = [-one,     -one/thr, -one/thr];
  added_node_coor_14(11,1:3) = [-one,      one/thr, -one/thr];
  added_node_coor_14(12,1:3) = [-one,      one,     -one/thr];
  added_node_coor_14(13,1:3) = [-one,     -one/thr,  one/two];
  added_node_coor_14(14,1:3) = [-one,      one/thr,  one/two];
  added_node_coor_14(15,1:3) = [-one/thr, -one,     -one/thr];
  added_node_coor_14(16,1:3) = [-one/thr, -one/thr, -one/thr];
  added_node_coor_14(17,1:3) = [-one/thr,  one/thr, -one/thr];
  added_node_coor_14(18,1:3) = [-one/thr,  one,     -one/thr];
  added_node_coor_14(19,1:3) = [ one/two, -one/thr,  one/two];
  added_node_coor_14(20,1:3) = [ one/two,  one/thr,  one/two];
  
  if (edge_case ==1) then
    added_node_coor(1:20,1) =  added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) =  added_node_coor_14(1:20,2)
    added_node_coor(1:20,3) =  added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [1,2,3,4,5,6,7,8]
  else if(edge_case ==2)then
    added_node_coor(1:20,1) = -added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) = -added_node_coor_14(1:20,2)
    added_node_coor(1:20,3) =  added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [3,4,1,2,7,8,5,6]
  else if(edge_case ==3)then
    added_node_coor(1:20,1) =  added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) = -added_node_coor_14(1:20,2)
    added_node_coor(1:20,3) = -added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [8,7,6,5,4,3,2,1]
  else if(edge_case ==4)then
    added_node_coor(1:20,1) = -added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) = -added_node_coor_14(1:20,2)
    added_node_coor(1:20,3) = -added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [7,8,5,6,3,4,1,2]
  else if(edge_case ==5)then
    added_node_coor(1:20,1) =  -added_node_coor_14(1:20,2)
    added_node_coor(1:20,2) =   added_node_coor_14(1:20,1)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [2,3,4,1,6,7,8,5]
  else if(edge_case ==6)then
    added_node_coor(1:20,1) =   added_node_coor_14(1:20,2)
    added_node_coor(1:20,2) =  -added_node_coor_14(1:20,1)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [4,1,2,3,8,5,6,7]
  else if(edge_case ==7)then
    added_node_coor(1:20,1) =  -added_node_coor_14(1:20,2)
    added_node_coor(1:20,2) =   added_node_coor_14(1:20,1)
    added_node_coor(1:20,3) =  -added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [6,7,8,5,2,3,4,1]
  else if(edge_case ==8)then
    added_node_coor(1:20,1) =  -added_node_coor_14(1:20,2)
    added_node_coor(1:20,2) =  -added_node_coor_14(1:20,1)
    added_node_coor(1:20,3) =  -added_node_coor_14(1:20,3)
    new_orig_node(1:8) = [7,6,5,8,3,2,1,4]
  else if(edge_case ==9)then
    added_node_coor(1:20,1) =   added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) =   added_node_coor_14(1:20,3)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,2)
    new_orig_node(1:8) = [1,2,6,5,4,3,7,8]
  else if(edge_case ==10)then
    added_node_coor(1:20,1) =   -added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) =   added_node_coor_14(1:20,3)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,2)
    new_orig_node(1:8) = [2,1,5,6,3,4,8,7]
  else if(edge_case ==11)then
    added_node_coor(1:20,1) =   -added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) =   -added_node_coor_14(1:20,3)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,2)
    new_orig_node(1:8) = [3,4,8,7,2,1,5,6]
  else if(edge_case ==12)then
    added_node_coor(1:20,1) =   added_node_coor_14(1:20,1)
    added_node_coor(1:20,2) =   -added_node_coor_14(1:20,3)
    added_node_coor(1:20,3) =   added_node_coor_14(1:20,2)
    new_orig_node(1:8) = [4,3,7,8,1,2,6,5]
  end if
  !add nodes.
  do i_added_n = 1,20
    c_kesi = added_node_coor(i_added_n,1)
    c_yita = added_node_coor(i_added_n,2)
    c_eta  = added_node_coor(i_added_n,3)
    !calculate global coordinate system according to the local coordinate system.
    call cal_coor_by_kesiyita_3d(c_kesi,c_yita,c_eta,x_nodes,y_nodes,z_nodes, c_node_coor(1),c_node_coor(2),c_node_coor(3))
    !check if the newly added node already exist.
    call vector_belongs_matrix_is_dou(c_added_node,3,temp_coor(1:c_added_node,1:3), c_node_coor(1:3),c_location,yes_exist)
    !if not exist before then add it. 
    if(yes_exist .eqv. .false.)then
      c_added_node = c_added_node + 1
      temp_coor(c_added_node,1:3)=c_node_coor(1:3)
      temp_added_nodes(i_added_n) = c_added_node
    !if exist before. 
    else
      temp_added_nodes(i_added_n) = c_location
    end if
  end do
  !add elements.
  !element 1
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(1))
  temp_elem_node(c_added_ele,2) = temp_added_nodes(3)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(9)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(15)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(10)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 2
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(10)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(11)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 3
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,4) = nn(new_orig_node(4))
  temp_elem_node(c_added_ele,5) = temp_added_nodes(11)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(18)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(12)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 4
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(3)
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(2))
  temp_elem_node(c_added_ele,3) = temp_added_nodes(7)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(15)
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,7) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(16)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 5
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(7)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(8)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(17)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 6
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(8)
  temp_elem_node(c_added_ele,3) = nn(new_orig_node(3))
  temp_elem_node(c_added_ele,4) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,7) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,8) = temp_added_nodes(18)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 7
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(9)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(15)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(10)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(5))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,7) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(13)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 8
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(5))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,7) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,8) = nn(new_orig_node(8))
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 9
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(11)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(18)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(12)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,7) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,8) = nn(new_orig_node(8))
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 10
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(10)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(11)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(14)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 11
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(2))
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(3))
  temp_elem_node(c_added_ele,3) = temp_added_nodes(8)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(7)
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(6))
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(7))
  temp_elem_node(c_added_ele,7) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(19)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !Make sure that the determinant detJ is positive 
  if (edge_case == 4 .or. edge_case == 7 .or.  &
        edge_case == 9 .or. edge_case == 11) then
    temp_elem_node_old = temp_elem_node
    temp_elem_node(c_added_ele-10:c_added_ele,1) =  temp_elem_node_old(c_added_ele-10:c_added_ele,4)
    temp_elem_node(c_added_ele-10:c_added_ele,2) =  temp_elem_node_old(c_added_ele-10:c_added_ele,3)
    temp_elem_node(c_added_ele-10:c_added_ele,3) =  temp_elem_node_old(c_added_ele-10:c_added_ele,2)
    temp_elem_node(c_added_ele-10:c_added_ele,4) =  temp_elem_node_old(c_added_ele-10:c_added_ele,1)
    temp_elem_node(c_added_ele-10:c_added_ele,5) =  temp_elem_node_old(c_added_ele-10:c_added_ele,8)
    temp_elem_node(c_added_ele-10:c_added_ele,6) =  temp_elem_node_old(c_added_ele-10:c_added_ele,7)
    temp_elem_node(c_added_ele-10:c_added_ele,7) =  temp_elem_node_old(c_added_ele-10:c_added_ele,6)
    temp_elem_node(c_added_ele-10:c_added_ele,8) =  temp_elem_node_old(c_added_ele-10:c_added_ele,5)
  end if
  
  case (4)
  !///////////////////////////////////////////////////
  !                 case  3: refine_template=4
  !///////////////////////////////////////////////////
  added_node_coor_5678(1,1:3)  = [-one/thr,  -one,one]
  added_node_coor_5678(2,1:3)  = [one/thr,   -one,one]
  added_node_coor_5678(3,1:3)  = [-one,  -one/thr,one]
  added_node_coor_5678(4,1:3)  = [-one/thr,-one/thr,one]
  added_node_coor_5678(5,1:3)  = [ one/thr,-one/thr,one]
  added_node_coor_5678(6,1:3)  = [-one,  one/thr, one]
  added_node_coor_5678(7,1:3)  = [-one/thr,one/thr, one]
  added_node_coor_5678(8,1:3)  = [one/thr, one/thr, one]
  added_node_coor_5678(9,1:3)  = [one,   one/thr, one]
  added_node_coor_5678(10,1:3) = [-one/thr,  one, one]
  added_node_coor_5678(11,1:3) = [one/thr,   one, one]
  added_node_coor_5678(12,1:3) = [-one,    -one,one/thr]
  added_node_coor_5678(13,1:3) = [-one/thr,  -one,one/thr]
  added_node_coor_5678(14,1:3) = [ one/thr,  -one,one/thr]
  added_node_coor_5678(15,1:3) = [ one,    -one,one/thr]
  added_node_coor_5678(16,1:3) = [-one,  -one/thr,one/thr]
  added_node_coor_5678(17,1:3) = [-one/thr,-one/thr,one/thr]
  added_node_coor_5678(18,1:3) = [ one/thr,-one/thr,one/thr]
  added_node_coor_5678(19,1:3) = [ one,  -one/thr,one/thr]
  added_node_coor_5678(20,1:3) = [-one,  one/thr,one/thr]
  added_node_coor_5678(21,1:3) = [-one/thr, one/thr,one/thr]
  added_node_coor_5678(22,1:3) = [one/thr,  one/thr,one/thr]
  added_node_coor_5678(23,1:3) = [one,    one/thr,one/thr]
  added_node_coor_5678(24,1:3) = [-one,     one,one/thr]
  added_node_coor_5678(25,1:3) = [-one/thr,   one,one/thr]
  added_node_coor_5678(26,1:3) = [one/thr,    one,one/thr]
  added_node_coor_5678(27,1:3) = [one,      one,one/thr]
  added_node_coor_5678(28,1:3) = [-one/thr,  -one,-one/two]
  added_node_coor_5678(29,1:3) = [ one/thr,  -one,-one/two]
  added_node_coor_5678(30,1:3) = [-one,   -one/thr,-one/two]
  added_node_coor_5678(31,1:3) = [ one,   -one/thr,-one/two]
  added_node_coor_5678(32,1:3) = [-one,    one/thr,-one/two]
  added_node_coor_5678(33,1:3) = [ one,     one/thr,-one/two]
  added_node_coor_5678(34,1:3) = [-one/thr,   one,-one/two]
  added_node_coor_5678(35,1:3) = [ one/thr,   one,-one/two]
  added_node_coor_5678(36,1:3) = [-one/thr,-one/thr,zr]
  added_node_coor_5678(37,1:3) = [ one/thr,-one/thr,zr]
  added_node_coor_5678(38,1:3) = [-one/thr, one/thr,zr]
  added_node_coor_5678(39,1:3) = [ one/thr, one/thr,zr]
  added_node_coor_5678(40,1:3) = [ one,  -one/thr,one]
  
  if (face_case == 1)then
    added_node_coor(1:40,1) =  added_node_coor_5678(1:40,1);
    added_node_coor(1:40,2) =  added_node_coor_5678(1:40,2);
    added_node_coor(1:40,3) =  added_node_coor_5678(1:40,3);
    new_orig_node(1:8) = [1,2,3,4,5,6,7,8];
  else if (face_case == 2)then
    added_node_coor(1:40,1) =  added_node_coor_5678(1:40,1);
    added_node_coor(1:40,2) =  added_node_coor_5678(1:40,2);
    added_node_coor(1:40,3) =  -added_node_coor_5678(1:40,3);
    new_orig_node(1:8) = [5,6,7,8,1,2,3,4];
  else if (face_case == 3)then
    added_node_coor(1:40,1) =  added_node_coor_5678(1:40,3);
    added_node_coor(1:40,2) =  added_node_coor_5678(1:40,2);
    added_node_coor(1:40,3) =  added_node_coor_5678(1:40,1);
    new_orig_node(1:8) = [1,5,8,4,2,6,7,3];
  else if (face_case == 4)then
    added_node_coor(1:40,1) =  -added_node_coor_5678(1:40,1);
    added_node_coor(1:40,2) =  -added_node_coor_5678(1:40,3);
    added_node_coor(1:40,3) =  -added_node_coor_5678(1:40,2);
    new_orig_node(1:8) = [7,8,4,3,6,5,1,2];
  else if (face_case == 5)then
    added_node_coor(1:40,1) =  -added_node_coor_5678(1:40,3);
    added_node_coor(1:40,2) =  added_node_coor_5678(1:40,1);
    added_node_coor(1:40,3) =  added_node_coor_5678(1:40,2);
    new_orig_node(1:8) = [2,3,7,6,1,4,8,5];
  else if (face_case == 6)then
    added_node_coor(1:40,1) =  -added_node_coor_5678(1:40,1);
    added_node_coor(1:40,2) =   added_node_coor_5678(1:40,3);
    added_node_coor(1:40,3) =  -added_node_coor_5678(1:40,2);
    new_orig_node(1:8) = [6,5,1,2,7,8,4,3];
  end if
  
  !add nodes.
  do i_added_n = 1,40
    c_kesi = added_node_coor(i_added_n,1)
    c_yita = added_node_coor(i_added_n,2)
    c_eta  = added_node_coor(i_added_n,3)
    !calculate global coordinate system according to the local coordinate system.
    call cal_coor_by_kesiyita_3d(c_kesi,c_yita,c_eta,x_nodes,y_nodes,z_nodes,c_node_coor(1),c_node_coor(2),c_node_coor(3))
    !check if the newly added node already exist.
    call vector_belongs_matrix_is_dou(c_added_node,3,temp_coor(1:c_added_node,1:3), c_node_coor(1:3),c_location,yes_exist)
    !if not exist before then add it. 
    if(yes_exist .eqv. .false.)then
      c_added_node = c_added_node + 1
      temp_coor(c_added_node,1:3)=c_node_coor(1:3)
      temp_added_nodes(i_added_n) = c_added_node
    !if exists before. 
    else
      temp_added_nodes(i_added_n) = c_location
    end if
  end do
  !add elements.
  !element 1
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(1));
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(2));
  temp_elem_node(c_added_ele,3) = nn(new_orig_node(3));
  temp_elem_node(c_added_ele,4) = nn(new_orig_node(4));
  temp_elem_node(c_added_ele,5) = temp_added_nodes(28);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(29);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(35);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(34);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 2
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(1));
  temp_elem_node(c_added_ele,2) = temp_added_nodes(28);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(34);
  temp_elem_node(c_added_ele,4) = nn(new_orig_node(4));
  temp_elem_node(c_added_ele,5) = temp_added_nodes(30);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(38);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(32);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 3
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(28);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(29);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(35);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(34);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(37);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(38);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 4
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(2));
  temp_elem_node(c_added_ele,2) = nn(new_orig_node(3));
  temp_elem_node(c_added_ele,3) = temp_added_nodes(35);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(29);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(31);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(33);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(37);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 5
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(1));
  temp_elem_node(c_added_ele,2) = temp_added_nodes(28);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(30);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(12);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(13);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(16);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 6
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(30);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(38);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(32);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(16);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(20);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 7
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(32);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(38);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(34);
  temp_elem_node(c_added_ele,4) = nn(new_orig_node(4));
  temp_elem_node(c_added_ele,5) = temp_added_nodes(20);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(25);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(24);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 8
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(28);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(29);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(37);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(13);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(14);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(17);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 9
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(36);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(37);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(38);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(21);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 10
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(38);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(35);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(34);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(26);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(25);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 11
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(2));
  temp_elem_node(c_added_ele,2) = temp_added_nodes(31);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(37);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(29);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(15);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(19);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(14);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 12
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(31);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(33);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(37);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(19);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(23);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(18);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 13
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(new_orig_node(3));
  temp_elem_node(c_added_ele,2) = temp_added_nodes(35);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(39);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(33);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(27);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(26);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(23);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 14
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(12);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(13);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(16);
  temp_elem_node(c_added_ele,5) = nn(new_orig_node(5));
  temp_elem_node(c_added_ele,6) = temp_added_nodes(1);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(4);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(3);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 15
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(16);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(20);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(3);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(4);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(7);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(6);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 16
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(20);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(25);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(24);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(6);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(7);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(10);
  temp_elem_node(c_added_ele,8) = nn(new_orig_node(8));
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 17
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(13);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(14);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(1);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(2);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(5);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(4);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 18
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(17);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(4);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(5);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(8);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(7);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 19
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(21);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(26);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(25);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(7);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(8);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(11);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(10);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 20
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(14);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(15);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(19);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(2);
  temp_elem_node(c_added_ele,6) = nn(new_orig_node(6));
  temp_elem_node(c_added_ele,7) = temp_added_nodes(40);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(5);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 21
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(18);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(19);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(23);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(5);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(40);
  temp_elem_node(c_added_ele,7) = temp_added_nodes(9);
  temp_elem_node(c_added_ele,8) = temp_added_nodes(8);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 22
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(22);
  temp_elem_node(c_added_ele,2) = temp_added_nodes(23);
  temp_elem_node(c_added_ele,3) = temp_added_nodes(27);
  temp_elem_node(c_added_ele,4) = temp_added_nodes(26);
  temp_elem_node(c_added_ele,5) = temp_added_nodes(8);
  temp_elem_node(c_added_ele,6) = temp_added_nodes(9);
  temp_elem_node(c_added_ele,7) = nn(new_orig_node(7));
  temp_elem_node(c_added_ele,8) = temp_added_nodes(11);
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  
  !Make sure that the determinant detJ is positive 
  if (face_case == 2 .or. face_case == 3 .or.  &
        face_case == 5 .or. face_case == 6) then
    temp_elem_node_old = temp_elem_node
    temp_elem_node(c_added_ele-21:c_added_ele,1) =  temp_elem_node_old(c_added_ele-21:c_added_ele,4)
    temp_elem_node(c_added_ele-21:c_added_ele,2) =  temp_elem_node_old(c_added_ele-21:c_added_ele,3)
    temp_elem_node(c_added_ele-21:c_added_ele,3) =  temp_elem_node_old(c_added_ele-21:c_added_ele,2)
    temp_elem_node(c_added_ele-21:c_added_ele,4) =  temp_elem_node_old(c_added_ele-21:c_added_ele,1)
    temp_elem_node(c_added_ele-21:c_added_ele,5) =  temp_elem_node_old(c_added_ele-21:c_added_ele,8)
    temp_elem_node(c_added_ele-21:c_added_ele,6) =  temp_elem_node_old(c_added_ele-21:c_added_ele,7)
    temp_elem_node(c_added_ele-21:c_added_ele,7) =  temp_elem_node_old(c_added_ele-21:c_added_ele,6)
    temp_elem_node(c_added_ele-21:c_added_ele,8) =  temp_elem_node_old(c_added_ele-21:c_added_ele,5)
  end if
  
  case (8)
  !///////////////////////////////////////////////////
  !                 case  4: refine_template=8
  !///////////////////////////////////////////////////
  added_node_coor(1,1:3)   = [-one/thr,   -one,   -one]
  added_node_coor(2,1:3)   = [ one/thr,   -one,   -one]
  added_node_coor(3,1:3)   = [-one,     -one, -one/thr]
  added_node_coor(4,1:3)   = [-one/thr,   -one, -one/thr]
  added_node_coor(5,1:3)   = [one/thr,    -one, -one/thr]
  added_node_coor(6,1:3)   = [ one,     -one, -one/thr]
  added_node_coor(7,1:3)   = [-one,     -one,  one/thr]
  added_node_coor(8,1:3)   = [-one/thr,   -one,  one/thr]
  added_node_coor(9,1:3)   = [ one/thr,   -one,  one/thr]
  added_node_coor(10,1:3)  = [ one,     -one,  one/thr]
  added_node_coor(11,1:3)  = [-one/thr,   -one,    one]
  added_node_coor(12,1:3)  = [ one/thr,   -one,    one]
  added_node_coor(13,1:3)   = [-one/thr,   -one/thr,   -one]
  added_node_coor(14,1:3)   = [ one/thr,   -one/thr,   -one]
  added_node_coor(15,1:3)   = [-one,     -one/thr, -one/thr]
  added_node_coor(16,1:3)   = [-one/thr,   -one/thr, -one/thr]
  added_node_coor(17,1:3)   = [one/thr,    -one/thr, -one/thr]
  added_node_coor(18,1:3)   = [ one,     -one/thr, -one/thr]
  added_node_coor(19,1:3)   = [-one,     -one/thr,  one/thr]
  added_node_coor(20,1:3)   = [-one/thr,   -one/thr,  one/thr]
  added_node_coor(21,1:3)   = [ one/thr,   -one/thr,  one/thr]
  added_node_coor(22,1:3)   = [ one,     -one/thr,  one/thr]
  added_node_coor(23,1:3)   = [-one/thr,   -one/thr,    one]
  added_node_coor(24,1:3)   = [ one/thr,   -one/thr,    one]
  added_node_coor(25,1:3)   = [-one/thr,   one/thr,   -one]
  added_node_coor(26,1:3)   = [ one/thr,   one/thr,   -one]
  added_node_coor(27,1:3)   = [-one,     one/thr, -one/thr]
  added_node_coor(28,1:3)   = [-one/thr,   one/thr, -one/thr]
  added_node_coor(29,1:3)   = [one/thr,    one/thr, -one/thr]
  added_node_coor(30,1:3)   = [ one,     one/thr, -one/thr]
  added_node_coor(31,1:3)   = [-one,     one/thr,  one/thr]
  added_node_coor(32,1:3)   = [-one/thr,   one/thr,  one/thr]
  added_node_coor(33,1:3)   = [ one/thr,   one/thr,  one/thr]
  added_node_coor(34,1:3)   = [ one,     one/thr,  one/thr]
  added_node_coor(35,1:3)   = [-one/thr,   one/thr,    one]
  added_node_coor(36,1:3)   = [ one/thr,   one/thr,    one]
  added_node_coor(37,1:3)   = [-one/thr,   one,   -one]
  added_node_coor(38,1:3)   = [ one/thr,   one,   -one]
  added_node_coor(39,1:3)   = [-one,     one, -one/thr]
  added_node_coor(40,1:3)   = [-one/thr,   one, -one/thr]
  added_node_coor(41,1:3)   = [one/thr,    one, -one/thr]
  added_node_coor(42,1:3)   = [ one,     one, -one/thr]
  added_node_coor(43,1:3)   = [-one,     one,  one/thr]
  added_node_coor(44,1:3)   = [-one/thr,   one,  one/thr]
  added_node_coor(45,1:3)   = [ one/thr,   one,  one/thr]
  added_node_coor(46,1:3)   = [ one,     one,  one/thr]
  added_node_coor(47,1:3)   = [-one/thr,   one,    one]
  added_node_coor(48,1:3)   = [ one/thr,   one,    one]
  added_node_coor(49,1:3)   = [ -one,   -one/thr,    -one]
  added_node_coor(50,1:3)   = [  one,   -one/thr,    -one]
  added_node_coor(51,1:3)   = [ -one,    one/thr,    -one]
  added_node_coor(52,1:3)   = [  one,    one/thr,    -one]
  added_node_coor(53,1:3)   = [ -one,   -one/thr,     one]
  added_node_coor(54,1:3)   = [  one,   -one/thr,     one]
  added_node_coor(55,1:3)   = [ -one,    one/thr,     one]
  added_node_coor(56,1:3)   = [  one,    one/thr,     one]
  !add nodes.
  do i_added_n = 1,56
    c_kesi = added_node_coor(i_added_n,1)
    c_yita = added_node_coor(i_added_n,2)
    c_eta  = added_node_coor(i_added_n,3)
    !calculate global coordinate system according to the local coordinate system.
    call cal_coor_by_kesiyita_3d(c_kesi,c_yita,c_eta,  &
        x_nodes,y_nodes,z_nodes, c_node_coor(1),c_node_coor(2),c_node_coor(3))
    !check if the newly added node already exist.
    call vector_belongs_matrix_is_dou(c_added_node,3,  &
        temp_coor(1:c_added_node,1:3), c_node_coor(1:3),c_location,yes_exist)
    !if not exist before then add it. 
    if(yes_exist .eqv. .false.)then
      c_added_node = c_added_node + 1
      temp_coor(c_added_node,1:3)=c_node_coor(1:3)
      temp_added_nodes(i_added_n) = c_added_node
    !if exists before. 
    else
      temp_added_nodes(i_added_n) = c_location
    end if
  end do
  !add elements.
  !element 1
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = nn(1)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(49)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(3)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(15)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 2
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(49)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(25)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(51)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(15)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(28)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(27)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 3
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(51)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(25)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(37)
  temp_elem_node(c_added_ele,4) = nn(4)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(27)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(28)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(40)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(39)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 4
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(1)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(4)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(16)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 5
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(13)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(26)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(25)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(16)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(29)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(28)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 6
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(25)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(26)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(38)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(37)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(28)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(29)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(41)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(40)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 7
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(2)
  temp_elem_node(c_added_ele,2) = nn(2)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(50)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(5)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(6)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(18)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(17)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 8
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(14)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(50)
  temp_elem_node(c_added_ele,3) = temp_added_nodes(52)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(26)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(17)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(18)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(30)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(29)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 9
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,1) = temp_added_nodes(26)
  temp_elem_node(c_added_ele,2) = temp_added_nodes(52)
  temp_elem_node(c_added_ele,3) = nn(3)
  temp_elem_node(c_added_ele,4) = temp_added_nodes(38)
  temp_elem_node(c_added_ele,5) = temp_added_nodes(29)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(30)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(42)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(41)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  do j_ele = 1,9
    temp_elem_node(c_added_ele+j_ele,1:4) =  temp_elem_node(c_added_ele+j_ele-9,5:8)
  end do
  !element 10
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(7)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(8)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(19)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 11
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(19)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(32)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(31)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 12
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(31)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(32)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(44)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(43)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 13
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(8)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(9)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(21)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(20)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 14
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(20)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(21)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(33)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(32)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 15
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(32)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(33)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(45)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(44)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 16
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(9)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(10)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(22)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(21)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 17
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(21)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(22)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(34)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(33)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 18
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(33)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(34)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(46)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(45)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  
  do j_ele = 1,9
    temp_elem_node(c_added_ele+j_ele,1:4) =  temp_elem_node(c_added_ele+j_ele-9,5:8)
  end do
  
  !element 19
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = nn(5)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(11)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(23)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(53)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 20
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(53)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(23)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(35)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(55)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 21
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(55)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(35)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(47)
  temp_elem_node(c_added_ele,8) = nn(8)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 22
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(11)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(12)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(24)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(23)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 23
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(23)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(24)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(36)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(35)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 24
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(35)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(36)
  temp_elem_node(c_added_ele,7) = temp_added_nodes(48)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(47)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 25
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) =temp_added_nodes(12)
  temp_elem_node(c_added_ele,6) =nn(6)
  temp_elem_node(c_added_ele,7) =temp_added_nodes(54)
  temp_elem_node(c_added_ele,8) =temp_added_nodes(24)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 26
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) =temp_added_nodes(24)
  temp_elem_node(c_added_ele,6) =temp_added_nodes(54)
  temp_elem_node(c_added_ele,7) =temp_added_nodes(56)
  temp_elem_node(c_added_ele,8) =temp_added_nodes(36)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
  !element 27
  c_added_ele = c_added_ele +1
  temp_elem_node(c_added_ele,5) = temp_added_nodes(36)
  temp_elem_node(c_added_ele,6) = temp_added_nodes(56)
  temp_elem_node(c_added_ele,7) = nn(7)
  temp_elem_node(c_added_ele,8) = temp_added_nodes(48)
  temp_elem_mat(c_added_ele)    = elem_mat(i_e)
end select
end do

!******************************************
!  Delete the refined original elemets.
!******************************************
num_wanted_elem = 0
do i_check = 1,c_added_ele
  !do not include refined original elemets.
  if(refined_elements(i_check) .eqv. .false.) then
    num_wanted_elem = num_wanted_elem +1
    wanted_elem(num_wanted_elem,1:8)=temp_elem_node(i_check,1:8)
    wanted_elem_mat(num_wanted_elem)=temp_elem_mat(i_check)
  else
  end if
end do
!
!**********************************
!  Update num_node and num_elem.
!**********************************
print *,' '
print *,'    Number of added nodes:', c_added_node - num_node
print *,'    Number of added elements:',num_wanted_elem-num_elem
num_node = c_added_node
num_elem = num_wanted_elem

!*********************************************
!   Reallocate and update global variables.
!*********************************************
print *,'    Reallocate and update variables...'
deallocate(coor);             allocate(coor(num_node,3))
deallocate(elem_node);        allocate(elem_node(num_elem,8))
deallocate(elem_mat);         allocate(elem_mat(num_elem))
coor(1:num_node,1:3)      = temp_coor(1:c_added_node,1:3)
elem_node(1:num_elem,1:8) = wanted_elem(1:num_wanted_elem,1:8)
elem_mat(1:num_elem)      = wanted_elem_mat(1:num_wanted_elem)

!***********************************************
!    save files *.node and *.elem.
!***********************************************
c_file_name_1  ='output.node'
c_file_name_2  ='output.elem'

print *, "    Save *.node file..."
open(201,file=c_file_name_1,status='unknown')
do i=1,num_node
  write(201, '(3e20.12)') coor(i,1:3)
end do
close(201)

print *, "    Save *.elem file..."
open(201,file=c_file_name_2,status='unknown')
do i=1,num_elem
  write(201, '(9i10)') elem_node(i,1:8),elem_mat(i)
end do
close(201)

print *,' '
print *,'    All done.'

return
end program local_grid_refinement

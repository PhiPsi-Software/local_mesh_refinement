% Written By: Fang Shi, 2021-08-10
% Website: phipsi.top
% Email: shifang@hyit.edu.cn

clear all; close all; clc; format compact;  format long;

c_figure = figure('units','normalized','position',[0.2,0.2,0.6,0.6],'Visible','on','Renderer','OpenGL');  
hold on;
axis off; 
axis equal;

view ([100,100,100])

title('Original Finite Element Mesh','FontSize',10)
axis off; axis equal;

%--------------------------
%       Read files.
%--------------------------
Node_Coor   = load(['input.node']);
Elem_Info   = load(['input.elem']);


Elem_Node = Elem_Info(:,1:8);
Elem_Material = Elem_Info(:,9);       % Material number of elements
num_of_Material = max(Elem_Material); % Number of materials
% Get the numbers of nodes, elements, boundaries and forces.
Num_Node = size(Node_Coor,1);
Num_Elem = size(Elem_Node,1);
% Get G_X_NODES,G_Y_NODES,G_Z_NODES
G_X_NODES = zeros(8,Num_Elem);
G_Y_NODES = zeros(8,Num_Elem);
G_Z_NODES = zeros(8,Num_Elem);
% Get the max and min value of node coordinates.
Max_X_Coor = max(max(Node_Coor(:,1)));
Min_X_Coor = min(min(Node_Coor(:,1)));
Max_Y_Coor = max(max(Node_Coor(:,2)));
Min_Y_Coor = min(min(Node_Coor(:,2)));
Max_Z_Coor = max(max(Node_Coor(:,3)));
Min_Z_Coor = min(min(Node_Coor(:,3)));	

%--------------------------
%   Draw the axes.
%--------------------------
axis([Min_X_Coor Max_X_Coor Min_Y_Coor Max_Y_Coor Min_Z_Coor Max_Z_Coor]);
c_X_Length = Max_X_Coor-Min_X_Coor;
c_Y_Length = Max_Y_Coor-Min_Y_Coor;
c_Z_Length = Max_Z_Coor-Min_Z_Coor;
Axes_length = (c_X_Length+c_Y_Length+c_Z_Length)/15;
plot3([0 Axes_length],[0 0],[0 0],'LineWidth',2,'Color','red')	
ts = text((c_X_Length+c_Y_Length+c_Z_Length)/14, 0, 0,"x",'Color','red','FontSize',15,'FontName','Consolas','FontAngle','italic');
plot3([0 0],[0 Axes_length],[0 0],'LineWidth',2,'Color','green')	
ts = text(0,(c_X_Length+c_Y_Length+c_Z_Length)/14,0,"y",'Color','green','FontSize',15,'FontName','Consolas','FontAngle','italic');
plot3([0 0],[0 0],[0 Axes_length],'LineWidth',2,'Color','blue')	
ts = text(0,0,(c_X_Length+c_Y_Length+c_Z_Length)/14,"z",'Color','blue','FontSize',15,'FontName','Consolas','FontAngle','italic');
	
%--------------------------
%   Plot nodes and mesh.
%--------------------------
%Plot mesh
c_plot_count = 0;
Ploted_Ele_lines =[0 0];
to_be_plot_count = 0;
to_be_plot_x = [];to_be_plot_y = [];to_be_plot_z = [];  
for iElem = 1:Num_Elem
	NN = [Elem_Node(iElem,1) Elem_Node(iElem,2) Elem_Node(iElem,3) Elem_Node(iElem,4) ...
		  Elem_Node(iElem,5) Elem_Node(iElem,6) Elem_Node(iElem,7) Elem_Node(iElem,8)];     
	for i=1:3
		if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	
			c_plot_count = c_plot_count + 1;
			Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
			to_be_plot_count = to_be_plot_count +1;
			to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+1),1)];
			to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+1),2)];
			to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+1),3)];		
		end			
	end
	for i=5:7
		if not(ismember(sort([NN(i),NN(i+1)]),Ploted_Ele_lines,'rows')) 	  
			c_plot_count = c_plot_count + 1;
			Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+1)]);
			to_be_plot_count = to_be_plot_count +1;
			to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+1),1)];
			to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+1),2)];
			to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+1),3)];		
			
		end					  
	end
	for i=1:4
		if not(ismember(sort([NN(i),NN(i+4)]),Ploted_Ele_lines,'rows')) 		
			c_plot_count = c_plot_count + 1;
			Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(i),NN(i+4)]);
			to_be_plot_count = to_be_plot_count +1;
			to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(i),1) Node_Coor(NN(i+4),1)];
			to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(i),2) Node_Coor(NN(i+4),2)];
			to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(i),3) Node_Coor(NN(i+4),3)];						
		end				  
	end	
	if not(ismember(sort([NN(1),NN(4)]),Ploted_Ele_lines,'rows')) 		
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(1),NN(4)]);		
		to_be_plot_count = to_be_plot_count +1;
		to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(1),1) Node_Coor(NN(4),1)];
		to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(1),2) Node_Coor(NN(4),2)];
		to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(1),3) Node_Coor(NN(4),3)];				
	end
	if not(ismember(sort([NN(5),NN(8)]),Ploted_Ele_lines,'rows')) 	
		c_plot_count = c_plot_count + 1;
		Ploted_Ele_lines(c_plot_count,1:2) = sort([NN(5),NN(8)]);	
		to_be_plot_count = to_be_plot_count +1;
		to_be_plot_x(to_be_plot_count,1:2) = [Node_Coor(NN(5),1) Node_Coor(NN(8),1)];
		to_be_plot_y(to_be_plot_count,1:2) = [Node_Coor(NN(5),2) Node_Coor(NN(8),2)];
		to_be_plot_z(to_be_plot_count,1:2) = [Node_Coor(NN(5),3) Node_Coor(NN(8),3)];					
	end
end 	
plot3(to_be_plot_x',to_be_plot_y',to_be_plot_z','LineWidth',0.5,'Color','blue')	
% Plot nodes.
x_node =[];y_node =[];z_node =[];
count = 0;
for j =1:Num_Node
	count = count +1;
	x_node(count) = Node_Coor(j,1);                                          
	y_node(count) = Node_Coor(j,2);  
	z_node(count) = Node_Coor(j,3);  		
end
plot3(x_node,y_node,z_node,'k.','MarkerSize',12.0)   
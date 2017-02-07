%Initialize Mesh
m=4; n=4;                   % of cells in each direction
TotCells=m*n;               % total number of cells
TotNodes=(m+1)*(n+1);       % total number of nodes
CellXWid=1.0; CellYWid=1.0; % side lengths for an individual cell
Nodes = zeros(TotNodes,2);  % matrix for Nodes coords
Cells = zeros(TotCells,4);  % index into Nodes for each cell node
                            % (lower left to upper left, counter-clockwise)
 
hold on;
 
% Set up the Nodes
for i=1:m+1
    for j=1:n+1
        k=(i-1)*(n+1) + j;
        Nodes(k,1)=j*CellXWid; Nodes(k,2)=i*CellYWid;
    end
end
 
% Plot the nodes
for k=1:TotNodes
    plot(Nodes(k,1),Nodes(k,2),'bO');
end
 
 
% Set up Cells, i.e. specify the 4 global node indices for each cell vertex
for i=1:m
    for j=1:n
        k = (i-1)*M+j;
        Cells(k,1)=k+i-1; Cells(k,2)=k+1+i-1; Cells(k,3)=Cells(k,2)+n+1; Cells(k,4)=Cells(k,2)+n;
    end
end
 
% Plot the cells
for k=1:TotCells
    for r=1:4
        x1=Nodes(Cells(k,r),1); y1=Nodes(Cells(k,r),2);
        if r<4
            x2=Nodes(Cells(k,r+1),1); y2=Nodes(Cells(k,r+1),2);
        else
            x2=Nodes(Cells(k,1),1); y2=Nodes(Cells(k,1),2);
        end
        x=[x1 x2]; y=[y1 y2];
        plot(x,y);
    end
end
        

%Calculate Lengths
L=zeros(n*m,4);
for i=1:n*m
    for j=1:3
        L(i,j)=norm(Nodes(Cells(i,j),:)-Nodes(Cells(i,j+1),:));
    end
    L(i,4)=norm(Nodes(Cells(i,4),:)-Nodes(Cells(i,1),:));
end

%Calculate Areas with shoelace formula
A=zeros(n*m,1);
for i=1:n*m
    A(i)=.5*(det(Nodes([Cells(i,1),Cells(i,2)],:))+det(Nodes([Cells(i,2),Cells(i,3)],:))+det(Nodes([Cells(i,3),Cells(i,4)],:))+det(Nodes([Cells(i,4),Cells(i,1)],:)));
end
        
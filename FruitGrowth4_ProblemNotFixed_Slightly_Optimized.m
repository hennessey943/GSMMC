% Couldn't fix the problem, but I optimized a couple for loops

%Fruit Growth Algorithm
%For each increment of time, compute dc(j)-> dw(k)->dA(k)
%Update network of cells to adjust cell area/volume

%Constants:
R=8.31;
T=293;
alpha=1;
D=1;
phi=1;
Z=1;
kappa=1;
Ciso=1;
Wiso=1;
Cout=.5*Ciso;
Pout=0;
tfinal=10;
dt=.001;

%Initialize Mesh
m=6; n=6;                   % of cells in each direction
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
        k = (i-1)*m+j;
        Cells(k,1)=k+i-1; Cells(k,2)=k+1+i-1; Cells(k,3)=Cells(k,2)+n+1; Cells(k,4)=Cells(k,2)+n;
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
Linit=L;
%Calculate Areas with shoelace formula
B=zeros(n*m,1);
for i=1:n*m
    B(i)=.5*(det(Nodes([Cells(i,1),Cells(i,2)],:))+det(Nodes([Cells(i,2),Cells(i,3)],:))+det(Nodes([Cells(i,3),Cells(i,4)],:))+det(Nodes([Cells(i,4),Cells(i,1)],:)));
end
A=zeros(n,m);
for i=1:n
    for j=1:m
        A(i,j)=B((j-1)*4+i);
    end
end

%Initial dx and dy
cx=zeros(n*m,1);
cy=zeros(n,m,1);
for i=1:n*m
    cx(i)=(L(i,2)+L(i,4))/2;
    cy(i)=(L(i,3)+L(i,1))/2;
end
dx=zeros(n,m);
dy=zeros(n,m);
for i=1:n
    for j=1:m
        dx(i,j)=cx((j-1)*4+i);
        dy(i,j)=cy((j-1)*4+i);
    end
end

%Initial concentrations
C=Ciso*ones(n,m);

%Initial water
W=Wiso*ones(n,m);

%Initial Pressure
P=zeros(n,m);

%Flux Functions Initializations

Fx=zeros(n,m);
Fy=zeros(n,m);
JL=zeros(n,m);
JR=zeros(n,m);
JT=zeros(n,m);
JB=zeros(n,m);
RHSC=zeros(n,m);
RHSW=zeros(n,m);
UB=zeros(n,m);
UT=zeros(n,m);
UL=zeros(n,m);
UR=zeros(n,m);


%% Update time step
Lnew=zeros(n,m);
%parpool
for t=1:tfinal
    %Move Vertices
    %Assume squares (version1)
    
    Lnew=A.^(1/2);
    L=reshape(Lnew,n*m,1);
    L=repmat(L,4);
    %New dx&dy
    cx = (L(:,2)+L(:,4))/2;
    cy = (L(:,2)+L(:,4))/2;
    
    dx = Lnew;
    dy = Lnew;
    
    %% Define Fluxes of forward difference step
    
    for i=1:n
        for j=1:m
            if i==n && j~=1 && j~=m
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==n&&j==1
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==n&&j==m
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j~=1 && j~=m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j==1
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==1&&j==m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif j==1&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif j==m&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            else
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            end
        end
    end
    RHSC = newRHSC;
    UB=Lnew.*JB; %Bottom
    UR=Lnew.*JR; %Right
    UT=Lnew.*JT; %Top
    UL=Lnew.*JL; %Left
    RHSW=-2*(UT+UB+UL+UR);
    
    
    %% Update Concentrations, Area, and Pressures
    %Update Concentrations, and Area
    
    C=RHSC*dt+C;
    W=RHSW*dt+W;
    A=W/(phi);
    
    
    
    %Update Pressures
    Q=zeros(n*m,1);
    %for i=1:n*m
        Q=kappa/4*((L(1:n*m,1)-Linit(1:n*m,1))./Linit(1:n*m,1)+...
            (L(1:n*m,2)-Linit(1:n*m,2))./Linit(1:n*m,2)+...
            (L(1:n*m,3)-Linit(1:n*m,3))./Linit(1:n*m,3)+...
            (L(1:n*m,4)-Linit(1:n*m,4))./Linit(1:n*m,4));
    %end
    for i=1:n
        for j=1:m
            P(i,j)=Q((j-1)*4+i);
        end
    end
    
    
    % Plot Cell Sizes
    Cell_Growth(A, n);
    % % Plot the cells
    %     for k=1:TotCells
    %         for r=1:4
    %             x1=Nodes(Cells(k,r),1); y1=Nodes(Cells(k,r),2);
    %             if r<4
    %                 x2=Nodes(Cells(k,r+1),1); y2=Nodes(Cells(k,r+1),2);
    %             else
    %                 x2=Nodes(Cells(k,1),1); y2=Nodes(Cells(k,1),2);
    %             end
    %             x=[x1 x2]; y=[y1 y2];
    %             plot(x,y);
    %         end
    %     end
end
delete(gcp('nocreate'));




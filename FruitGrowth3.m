%Fruit Growth Algorithm
%For each increment of time, compute dc(j)-> dw(k)->dA(k)
%Update network of cells to adjust cell area/volume
clear all
%Constants:
R=8.31; %Gas constant (std units)
T=293; %Temp in Kelvin
alpha=1;
D=10^(-14); %Diffusivity of sugar
phi=1;
Z=10^(-7);
kappa=1;
Ciso=.15; %isotonic sugar concentration
Wiso=.0001; %isotonic water mass
Cout=1.25*Ciso; %outer sugar concentration
Pout=0; %outer pressure
tfinal=10; %number of time steps
dt=5; %time step size

%Initialize Mesh
m=8; n=8;                   % of cells in each direction
TotCells=m*n;               % total number of cells
TotNodes=(m+1)*(n+1);       % total number of nodes
CellXWid=.01; CellYWid=.01; % side lengths for an individual cell
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
% for k=1:TotNodes
%     plot(Nodes(k,1),Nodes(k,2),'bO');
% end
 
 
% Set up Cells, i.e. specify the 4 global node indices for each cell vertex
for i=1:m
    for j=1:n
        k = (i-1)*m+j;
        Cells(k,1)=k+i-1; Cells(k,2)=k+1+i-1; Cells(k,3)=Cells(k,2)+n+1; Cells(k,4)=Cells(k,2)+n;
    end
end 
% Plot the cells
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

%Calculate Lengths
L=zeros(n*m,4);
for i=1:n*m
    for j=1:3
        L(i,j)=norm(Nodes(Cells(i,j),:)-Nodes(Cells(i,j+1),:));
    end
    L(i,4)=norm(Nodes(Cells(i,4),:)-Nodes(Cells(i,1),:));
end
%store initial lengths
Linit=L;
%Calculate Areas with shoelace formula
B=zeros(n*m,1);
for i=1:n*mswz
    B(i)=.5*(det(Nodes([Cells(i,1),Cells(i,2)],:))+det(Nodes([Cells(i,2),Cells(i,3)],:))+det(Nodes([Cells(i,3),Cells(i,4)],:))+det(Nodes([Cells(i,4),Cells(i,1)],:)));
end
A=zeros(n,m);
for i=1:n
    for j=1:m
        A(i,j)=B((j-1)*4+i);
    end
end

%Initialize dx and dy
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




%Water
% JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
% JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
% JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
% JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom

%Water Flows

% UB(i,j)=L(i,1)*JB(i,j); %Bottom
% UR(i,j)=L(i,2)*JR(i,j); %Right
% UT(i,j)=L(i,3)*JT(i,j); %Top
% UL(i,j)=L(i,4)*JL(i,j); %Left



%
%
%

%Right Hand Side
% 
% RHS_C(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
% RHS_W(i,j)=2*(UT(i,j)+UB(i,j)+UL(i,j)+UR(i,j));

%
%
%
%
%Update
Lnew=zeros(n,m);
for t=1:tfinal
    %Move Vertices
    %Assume squares (version1)
    
    for i=1:n
        for j=1:m
            Lnew(i,j)=(A(i,j))^(1/2);
        end
    end
    L=reshape(Lnew',n*m,1);
    L=repmat(L,4);
    %New dx&dy
    for i=1:n*m
        cx(i)=(L(i,2)+L(i,4))/2;
        cy(i)=(L(i,3)+L(i,1))/2;
    end
    for i=1:n
        for j=1:m
            dx(i,j)=Lnew(i,j);
            dy(i,j)=Lnew(i,j);
        end
    end
    %Define Fluxes and Right hand sides of forward difference step
    for i=1:n
        for j=1:m
            if i==n && j~=1 && j~=m
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==n&&j==1
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==n&&j==m
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j~=1 && j~=m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j==1
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==1&&j==m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif j==1&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif j==m&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            else            
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JL(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JR(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JT(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JB(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            end
            UB(i,j)=Lnew(i,j)*JB(i,j); %Bottom
            UR(i,j)=Lnew(i,j)*JR(i,j); %Right
            UT(i,j)=Lnew(i,j)*JT(i,j); %Top
            UL(i,j)=Lnew(i,j)*JL(i,j); %Left
            RHSW(i,j)=2*(UT(i,j)+UB(i,j)+UL(i,j)+UR(i,j));
        end
    end
    for i=1:n
        for j=1:m
            if i==n && j~=1 && j~=m
                RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==n&&j==1
                RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==n&&j==m
                RHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j~=1 && j~=m
                RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j==1
                RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==1&&j==m
                RHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif j==1&&i~=1&&i~=n
                RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif j==m&&i~=1&&i~=n
                RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            else 
                RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            end
        end
    end
    %Update Concentrations, and Area
    for i=1:n
        for j=1:m
            C(i,j)=RHSC(i,j)*dt+C(i,j);
            W(i,j)=RHSW(i,j)*dt+W(i,j);
            A(i,j)=W(i,j)/(phi);
        end
    end
    
    
    %Update Pressures
    Q=zeros(n*m,1);
    for i=1:n*m
        Q(i)=kappa/4*((L(i,1)-Linit(i,1))/Linit(i,1)+(L(i,2)-Linit(i,2))/Linit(i,2)+(L(i,3)-Linit(i,3))/Linit(i,3)+(L(i,4)-Linit(i,4))/Linit(i,4));
    end
    for i=1:n
        for j=1:m
            P(i,j)=Q((j-1)*4+i);
        end
    end
    
    %Plot Growth
     delta = 10;
    
    axis( [0 50 0 50] );
    for r = 1:n
        for c = 1:m
        
            s = sqrt( A(r,c)/Wiso );
            x = r*delta; y = c*delta;
            fill( [x x+s x+s x], [y y y+s y+s], 'b' ); hold on; axis( [0 10*(n+2) 0 10*(m+2)] );
        end
    end
    M(t)=getframe;

end
movie(M,5,1)
    



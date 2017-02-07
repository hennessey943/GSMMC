clear all
% Couldn't fix the problem, but I optimized a couple for loops

%Fruit Growth Algorithm
%For each increment of time, compute dc(j)-> dw(k)->dA(k)
%Update network of cells to adjust cell area/volume
fid = figure;
writerObj = VideoWriter('Fruit_Growth3.avi');
writerObj.FrameRate = 7;
open(writerObj);

%Constants:
R=8.31; %Gas constant std units
T=293; %Temp (kelvin)
alpha=1;
D=10^(-14)*10000; %Diffusivity of sugar
phi=1;
Z=10^(-11); % Z = 10^-7 <- testing 
kappa=1000/15;
Ciso=.15; %isotonic sugar concentration
Wiso=0.0001; %isotonic water mass
Cout=.75*Ciso;%.5*Ciso; %sugar concentration out of bounds
Pout=0; %Pressure out of bounds
tfinal=60; %number of iterations
dt=3600*12; %time step  dt=5 <-- testing
Cin=Cout; %Sugar brought in
Win=Wiso; %Water brought in
Pin=0.3;
%Initialize Mesh
m=8; n=8;                   % of cells in each direction
TotCells=m*n;               % total number of cells
TotNodes=(m+1)*(n+1);       % total number of nodes
CellXWid=.01; CellYWid=.01; % side lengths (in meters) for an individual cell
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
%Store initial lengths
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
                %newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==n&&j==1
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==n&&j==m
                Fx(i,j)=-D*(Cout-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %newRHSC(i,j)=-(-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j~=1 && j~=m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            elseif i==1&&j==1
                Fx(i,j)=-D*(C(i+1,j)-Cin)/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cin)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pin-alpha*R*T*(C(i,j)-Cin)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pin-alpha*R*T*(C(i,j)-Cin)); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif i==1&&j==m
                Fx(i,j)=-D*(C(i+1,j)-Cout)/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            elseif j==1&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-Cout)/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1))/(2*dy(i,j));
            elseif j==m&&i~=1&&i~=n
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(Cout-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-Pout-alpha*R*T*(C(i,j)-Cout)); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(-Fy(i,j-1))/(2*dy(i,j));
            else
                Fx(i,j)=-D*(C(i+1,j)-C(i-1,j))/(2*dx(i,j));
                Fy(i,j)=-D*(C(i,j+1)-C(i,j-1))/(2*dy(i,j));
                JB(i,j)=Z*(P(i,j)-P(i,j-1)-alpha*R*T*(C(i,j)-C(i,j-1))); %Left
                JT(i,j)=Z*(P(i,j)-P(i,j+1)-alpha*R*T*(C(i,j)-C(i,j+1))); %Right
                JR(i,j)=Z*(P(i,j)-P(i+1,j)-alpha*R*T*(C(i,j)-C(i+1,j))); %Top
                JL(i,j)=Z*(P(i,j)-P(i-1,j)-alpha*R*T*(C(i,j)-C(i-1,j))); %Bottom
                %newRHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            end
        end
    end
    
    for i=1:n
        for j=1:m
            if i==n && j~=1 && j~=m
                RHSC(i,j)=-(-D/(2*dx(n,j))*(Cout-C(n,j))-Fx(n-1,j))/(2*dx(n,j))...
                    - (Fy(n,j+1) - Fy(n,j-1))/(2*dy(n,j));
            elseif i==n&&j==1
                nFx1 = -D/(2*dx(n,1))*(Cout - C(n-1,1));
                nFx2 = -D/(2*dx(n-1,1))*(C(n,1)-C(n-2,1));
                nFy1 = -D/(2*dy(n,2))*(C(n,3)-C(n,1));
                nFy2 = -D/(2*dy(n,1))*(C(n,1)-Cout);
                RHSC(i,j)=-(nFx1-nFx2)/(2*dx(n,1)) ...
                    -(nFy1-nFy2)/(2*dy(n,1));
            elseif i==n&&j==m
                nFx1 = -D/(2*dx(n,m))*(Cout - C(n,m));
                nFx2 = -D/(2*dx(n-1,m))*(C(n,m)-C(n-2,m));
                nFy1 = -D/(2*dy(n,m))*(Cout-C(n,m));
                nFy2 = -D/(2*dy(n,m-1))*(C(n,m)-C(n,m-2));
                RHSC(i,j)=-(nFx1-nFx2)/(2*dx(n,m)) ...
                    -(nFy2-nFy1)/(2*dy(n,m));
            elseif i==1&&j~=1 && j~=m
                RHSC(i,j)=-(Fx(2,j)-Fx(1,j)+D/(2*dx(1,j))*(C(1,j)-Cout))/(2*dx(1,j))...
                    -(Fy(1,j+1)*Fy(1,j-1))/(2*dy(1,j));
            elseif i==1&&j==1
                RHSC(i,j)=-(Fx(2,1)+D/(2*dx(1,1))*(C(1,1)-Cout))/(2*dx(1,1))...
                    -(Fy(1,2)+D/(2*dy(1,1))*(C(1,1)-Cout))/(2*dy(1,1));
            elseif i==1&&j==m
                RHSC(i,j)=-(Fx(2,m)+D/(2*dx(1,m))*(C(1,m)-Cout))/(2*dx(1,m))...
                    -(-D/(2*dy(1,m))*(Cout-C(1,m))-Fy(1,m-1))/(2*dy(1,m));
            elseif j==1&&i~=1&&i~=n
                RHSC(i,j)=-(Fx(i+1,1)-Fx(i-1,1))/(2*dx(i,1))...
                    -(Fy(i,2) + D/(2*dy(i,1))*(C(i,1)-Cout))/(2*dy(i,1));
            elseif j==m&&i~=1&&i~=n
                RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(-D*(Cout-C(i,j-1))/(2*dy(i,j))-Fy(i,j-1))/(2*dy(i,j));
            else
                RHSC(i,j)=-(Fx(i+1,j)-Fx(i-1,j))/(2*dx(i,j))-(Fy(i,j+1)-Fy(i,j-1))/(2*dy(i,j));
            end
        end
    end
    
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
    for i=1:n*m
        Q(i)=kappa/4*((L(i,1)-Linit(i,1))./Linit(i,1)+...
        (L(i,2)-Linit(i,2))./Linit(i,2)+...
        (L(i,3)-Linit(i,3))./Linit(i,3)+...
        (L(i,4)-Linit(i,4))./Linit(i,4));
    end
    for i=1:n
        for j=1:m
            P(i,j)=Q((j-1)*4+i);
        end
    end
    %Plot Concentration
    
    
    % Plot Cell Sizes
    figure(fid);
    Cell_Growth(A/Wiso, n, t/2);
    
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
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
close(writerObj);



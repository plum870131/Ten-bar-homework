%fmincon parameter
x0=[0.5;0.5];
A=[];
b=[];
Aeq=[];  
beq=[];
ub=[0.5;0.5];   
lb=[0.001;0.001];   

%optimize
options = optimoptions('fmincon','Algorithm', 'sqp','Display','iter');
[x,fval] = fmincon(@objfun,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

%get required answer
[stress,disp,force] = answer(x)


% ¥Ø¼Ð¨ç¼Æ
function f=objfun(r)
    rho=7860;
    r1=r(1);
    r2=r(2);
    L=[1 1 1 1 1 1 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi)];
    L=L*9.14;
    A=[r1 r1 r1 r1 r1 r1 r2 r2 r2 r2];
    A=A.^2;
    A=A*pi;
    
    total_weight=0;
    for i=1:10
        total_weight = total_weight + L(i)*A(i)*rho;
    end
    f=total_weight;

end


%nonlinear boundary condition
function [c,ceq]=nonlcon(r)
    yield=250*10^6;
    E=200*10^9;
    r1=r(1);
    r2=r(2);

    node=[3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4;];
    
    A=[r1 r1 r1 r1 r1 r1 r2 r2 r2 r2];
    A=A.^2;
    A=A*pi;

    L=[1 1 1 1 1 1 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi)];
    L=L*9.14;

    Co=[-1 -1 -1 -1 0 0  -cos(1/4*pi) -cos(1/4*pi) -cos(1/4*pi) -cos(1/4*pi)];
    Si=[0 0 0 0 -1 -1 sin(1/4*pi) -sin(1/4*pi) sin(1/4*pi) -sin(1/4*pi)];
    
    %stiffness matrix of each element
    M = cell(10);
    for i=1:10
        co=Co(i);
        si=Si(i);
        M{i} = [co^2 co*si -co^2 -co*si; co*si si^2 -co*si -si^2; -co^2 -co*si co^2 co*si; -co*si -si^2 co*si si^2]*E*A(i)/L(i);
    end
    
    %combined stiffness matrix
    stiffness=zeros(12,12);
    for i =1:10
        node_a=node(i,1);
        node_b=node(i,2);
        stiffness_temp = zeros(4,12);
        stiffness_temp(:,2*node_a-1) = M{i}(:,1);
        stiffness_temp(:,2*node_a) =   M{i}(:,2);
        stiffness_temp(:,2*node_b-1) = M{i}(:,3);
        stiffness_temp(:,2*node_b) =   M{i}(:,4);

        stiffness_temp2 = zeros(12,12);
        stiffness_temp2(2*node_a-1,:) =   stiffness_temp(1,:);
        stiffness_temp2(2*node_a,:) = stiffness_temp(2,:);
        stiffness_temp2(2*node_b-1,:) =   stiffness_temp(3,:);
        stiffness_temp2(2*node_b,:) = stiffness_temp(4,:);

        stiffness = stiffness + stiffness_temp2;
    end 

    %calculate displacement Q
    F=[0 0 0 -10^7 0 0 0 -10^7 0 0 0 0];
    F_de = F(1:8);
    stiffness_de = stiffness(1:8, 1:8);
    %Q=K-1 Ft
    Q_de = inv(stiffness_de)*transpose(F_de);
    Q=[Q_de ;0;0;0;0];
    
    %calculate stress
    stress=size(10);
    for i =1:10
        stress(i) = E/L(i)*[-Co(i) -Si(i) Co(i) Si(i)]*transpose([Q(2*node(i,1)-1) Q(2*node(i,1)) Q(2*node(i,2)-1) Q(2*node(i,2))]);
    end
    
    c=size(11);
    for i =1:10
       c(i)=abs( stress(i))-yield;
    end
    %node=[3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4;];
    
    %node 2 displacement boundary equation
    c(11)=sqrt(Q(3)^2+Q(4)^2)-0.02;
    ceq=[];

end

%output answer
function [stress,disp,force]= answer(r)
    yield=250*10^6;
    E=200*10^9;
    r1=r(1);
    r2=r(2);

    node=[3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4;];
    
    A=[r1 r1 r1 r1 r1 r1 r2 r2 r2 r2];
    A=A.^2;
    A=A*pi;

    L=[1 1 1 1 1 1 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi) 1/cos(1/4*pi)];
    L=L*9.14;

    Co=[-1 -1 -1 -1 0 0  -cos(1/4*pi) -cos(1/4*pi) -cos(1/4*pi) -cos(1/4*pi)];
    Si=[0 0 0 0 -1 -1 sin(1/4*pi) -sin(1/4*pi) sin(1/4*pi) -sin(1/4*pi)];

    M = cell(10);
    for i=1:10
        co=Co(i);
        si=Si(i);
        M{i} = [co^2 co*si -co^2 -co*si; co*si si^2 -co*si -si^2; -co^2 -co*si co^2 co*si; -co*si -si^2 co*si si^2]*E*A(i)/L(i);
    end

    stiffness=zeros(12,12);

    for i =1:10
        node_a=node(i,1);
        node_b=node(i,2);
        stiffness_temp = zeros(4,12);
        stiffness_temp(:,2*node_a-1) = M{i}(:,1);
        stiffness_temp(:,2*node_a) =   M{i}(:,2);
        stiffness_temp(:,2*node_b-1) = M{i}(:,3);
        stiffness_temp(:,2*node_b) =   M{i}(:,4);

        stiffness_temp2 = zeros(12,12);
        stiffness_temp2(2*node_a-1,:) =   stiffness_temp(1,:);
        stiffness_temp2(2*node_a,:) = stiffness_temp(2,:);
        stiffness_temp2(2*node_b-1,:) =   stiffness_temp(3,:);
        stiffness_temp2(2*node_b,:) = stiffness_temp(4,:);

        stiffness = stiffness + stiffness_temp2;

    end 

    F=[0 0 0 -10^7 0 0 0 -10^7 0 0 0 0];
    F_de = F(1:8);
    stiffness_de = stiffness(1:8, 1:8);
    %Q=K-1 Ft
    Q_de = inv(stiffness_de)*transpose(F_de);
    Q=[Q_de ;0;0;0;0];
    disp=Q;
    
   
    
    
    stress=size(10);
    for i =1:10
        stress(i) = E/L(i)*[-Co(i) -Si(i) Co(i) Si(i)]*transpose([Q(2*node(i,1)-1) Q(2*node(i,1)) Q(2*node(i,2)-1) Q(2*node(i,2))]);
    end
    
    c=size(11);
    for i =1:10
       c(i)=abs( stress(i))-yield;
    end
    %node=[3 5; 1 3; 4 6; 2 4; 3 4; 1 2; 4 5; 3 6; 2 3; 1 4;];
    
  
    force=stiffness(9:12,1:12)*Q;
    

end
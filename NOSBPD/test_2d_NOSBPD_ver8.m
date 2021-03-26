clear all;  clc;

% Lattice
% Lb = 1.0/150.0; % Lattice spacing 
% tb = Lb; % Thickness of the plate
% dA = Lb^2; % area of a paricle
% dV = Lb^2*tb; % volume of a particle

H = 0.44;
HL = .16;
W = 0.48;
TD = .16;
Lb = 0.01;
tb = Lb;
dV = Lb^2*tb;

% [p, Np, Nx, Ny] = square_latttice(2,Lx,Ly,Lb); % reference coordinate
% [p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,Lb); % reference coordinate
[p, Np, Nx, Ny] = cooks_membrane(H,HL,W,TD,Lb);

figure(1)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

horizon = 2.015*Lb; % 3.015*L ----- horizon

% inputs

% Rubber
rho = 1100; % density
C = 1e6; % damping
nu = 0.4; % poisson's ratio
E = 1e6; % Young's modulus 
lambda = E*nu/((1+nu)*(1-2*nu)); % Lame parameters
mu = E/(2*(1+nu));
P_ext = 1e4;
F_ext = P_ext/Lb; % external force

dt = 0.0001; % time stepsize
Nt = 20000; % number of steps (iterations)
M = 1000; % the time when the full amount of force is acting on a material body

% % Steel
% rho = 7850; % density
% C = 4e8; % damping
% nu = 1/3; % poisson's ratio
% E = 200e9; % Young's modulus 
% lambda = E*nu/((1+nu)*(1-2*nu)); % Lame parameters
% mu = E/(2*(1+nu));
% P_ext = 200e6;
% F_ext = P_ext/Lb; % external force
% 
% dt = 0.0000001; % time stepsize
% Nt = 1000000; % number of steps (iterations)
% M = 5000; % the time when the full amount of force is acting on a material body

%% Constructing initial bonds and shape tensor

[xi, K, v, omega, num_fam, fam_ind] = initial(Np,p,Lb,dV,horizon);
xi = sparse(xi); v = sparse(v); omega = sparse(omega); fam_ind = sparse(fam_ind);

%% Simulation
Ib = eye(2);

p_ini = p;
u_ini = zeros(Np,2); % initial speed

p_current = p_ini;
u_current = u_ini;

A = [p_ini u_ini];
% fileID = fopen('initial_Data','w');
% fprintf(fileID,'%e %e %e %e\r\n',A');
% fclose(fileID);

for ii = 0:Nt-1
    Time = dt*(ii+1)
    %% Dynamics using the advanced trapezoidal rule
    
    % apply the max-load at M-th step
    if ii < M
        nnn = (ii+1)/M;
    end
    
    [p_new, u_new] = Advanced_traperzoidal_1(p_current,u_current,p_ini,Np,H,W,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn,C,dt,rho);

    u_current = u_new;
    p_current = p_new;
    
%     if mod(Time,0.001) == 0
%         A = [p_new u_new];
%         fname = sprintf ('%s%i%s%.3f', 'D_',ii+1 ,'_',Time);
%         fileID = fopen(fname,'w');
%         fprintf(fileID,'%e %e %e %e\r\n',A');
%         fclose(fileID);
%     end
    
%     u(:,2*(ii+1)+1:2*(ii+2)) = u_new;
    p(:,2*(ii+1)+1:2*(ii+2)) = p_new;
end

% plotting
% figure(1)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on

figure(2)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
scatter(p_new(:,1),p_new(:,2),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

figure(3)
scatter(p_new(:,1),p_new(:,2),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

% figure(2)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% hold on
% scatter(p(:,2*ii+1),p(:,2*(ii+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on

% figure(5)
% scatter(p(:,2*ii+1),p(:,2*(ii+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on
% 
% figure(6)
% movie(p,Nt,Lx,Ly)

% reconstructing position and velocity vectors
 
p1 = p(:,2*(ii+1)+1:2*(ii+2));

for i = 1:Nx+1
    for j = 1:Ny+1
        X(j,i) = p_ini(i+(Nx+1)*(j-1),1);
        Y(j,i) = p_ini(i+(Nx+1)*(j-1),2);
        U(j,i) = p_new(i+(Nx+1)*(j-1),1)- p_ini(i+(Nx+1)*(j-1),1);
        V(j,i) = p_new(i+(Nx+1)*(j-1),2)- p_ini(i+(Nx+1)*(j-1),2);
        
        pX(j,i) = p1(i+(Nx+1)*(j-1),1);
        pY(j,i) = p1(i+(Nx+1)*(j-1),2);
        pU(j,i) = p1(i+(Nx+1)*(j-1),1)- p_ini(i+(Nx+1)*(j-1),1);
        pV(j,i) = p1(i+(Nx+1)*(j-1),2)- p_ini(i+(Nx+1)*(j-1),2);
    end
end

figure(8)
pcolor(X,Y,U)
title('U_x')
shading flat
colorbar

figure(9)
pcolor(X,Y,V)
title('U_y')
shading flat
colorbar

figure(10)
pcolor(pX,pY,pU)
title('U_x')
shading flat
colorbar

figure(11)
pcolor(pX,pY,pV)
title('U_y')
shading flat
colorbar

% figure(12)
% plot(X(Ny/2+1,:),U(Ny/2+1,:),'-*b',X(Ny/2+1,:),(P_ext/E)*X(Ny/2,:),'--r')
% ylabel('U_x')
% xlabel('x')
% legend('PD_solution','Analytical')
% 
% figure(13)
% plot(Y(:,Nx/2+1),V(:,Nx/2+1),'-*b',Y(:,Nx/2+1),-nu*(P_ext/E)*Y(:,Nx/2+1),'--r')
% ylabel('U_y')
% xlabel('y')
% legend('PD_solution','Analytical')
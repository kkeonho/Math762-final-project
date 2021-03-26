clear all; close all; clc;

% Lattice
Lx = 1.; 
Ly = .4; % sqrt(3)*Ly
Lb = .1; % Lattice spacing 
tb = 0.01; % Thickness of the plate
dA = Lb^2; % area of a paricle
dV = Lb^2*tb; % volume of a particle

[p, Np, Nx, Ny] = square_latttice(2,Lx,Ly,Lb); % reference coordinate
% [p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,Lb); % reference coordinate
u_ini = zeros(Np,2); % initial speed

horizon = 3.015*Lb; % 3.015*L ----- horizon

% inputs

% Rubber
rho = 1100; % density
C = 1e6; % damping
nu = 0.4; % poisson's ratio
E = 1e6; % Young's modulus 
lambda = E*nu/((1+nu)*(1-2*nu)); % Lame parameters
mu = E/(2*(1+nu));
P_ext = 1e3;
F_ext = P_ext/Lb; % external force

dt = 0.001; % time stepsize
Nt = 1000; % number of steps (iterations)
M = 500; % the time when the full amount of force is acting on a material body


% Steel
% rho = 7850; % density
% C = 4e8; % damping
% nu = 1/3; % poisson's ratio
% E = 200e9; % Young's modulus 
% lambda = E*nu/((1+nu)*(1-2*nu)); % Lame parameters
% mu = E/(2*(1+nu));
% P_ext = 200e6;
% F_ext = P_ext/Lb; % external force

% dt = 0.00001; % time stepsize
% Nt = 1000; % number of steps (iterations)
% M = 500; % the time when the full amount of force is acting on a material body

%% Constructing initial bonds and shape tensor

[xi, K, v, omega, num_fam, fam_ind] = initial(Np,p,Lb,dV,horizon);
% [xi, K, v, omega, num_fam, fam_ind] = initial(Np,p,Lb,dA,horizon);


%% Simulation
Ib = eye(2);

p_ini = p;
p_current = p_ini;
u_current = u_ini;

% indi = 0; eps = 1e-9;

for ii = 0:Nt-1
    Time = dt*(ii+1)
    %% Dynamics using the advanced trapezoidal rule
    
    % apply the max-load at M-th step
    if ii < M
        nnn = (ii+1)/M;
    end
    
    % net force matrix
    F_net_current = deformed(p_ini,p_current,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
%     F_net_current = deformed(p_ini,p_current,Np,Lx,Ly,Lb,dA,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
    
    % predictor
    u_new = (1 - C*dt/rho)*u_current + dt/rho*F_net_current;
    p_new = p_current + dt*u_new;
    
    % net force matrix
    F_net_new = deformed(p_ini,p_new,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
%     F_net_new = deformed(p_ini,p_new,Np,Lx,Ly,Lb,dA,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
    
    % corrector
    u_new = (1 - C*dt/rho)*u_current + 0.5*dt/rho*(F_net_current + F_net_new);
    p_new = p_current + 0.5*dt*(u_new + u_current);

    u_current = u_new;
    p_current = p_new;
    
    
    
%     u(:,2*(ii+1)+1:2*(ii+2)) = u_new;
%     p(:,2*(ii+1)+1:2*(ii+2)) = p_new;
end

initial_dat = [p_ini u_ini];
deformed_dat = [p_new u_new];

fileID = fopen('initial_Data.txt','w');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\n',initial_dat);
fclose(fileID);

fileID = fopen('deformed_Data.txt','w');
fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f\n',deformed_dat);
fclose(fileID);



% plotting
% figure(1)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on
% 
% figure(2)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% hold on
% scatter(p_new(:,1),p_new(:,2),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on

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
for i = 1:Nx+1
    for j = 1:Ny+1
        Xp(j,i) = p_new(i+(Nx+1)*(j-1),1);
        Yp(j,i) = p_new(i+(Nx+1)*(j-1),2);
        Vx(j,i) = u_new(i+(Nx+1)*(j-1),1);
        Vy(j,i) = u_new(i+(Nx+1)*(j-1),2);
    end
end

% Numerical and analytic solutions
n1 = 0; n2 = 0;
for i = 1:Np
    if p_ini(i,2) == p_ini((Np+1)/2,2)
        n1 = n1+1;
        x(n1) = p_ini(i,1);
        uxa(n1) = p_ini(i,1)*(P_ext/E);
%         uxa1(n1) = p_new(i,1)*(Pressure/E);
        ux(n1) = p_new(i,1)- p_ini(i,1);
    elseif p_ini(i,1) == p_ini((Np+1)/2,1)
        n2 = n2+1;
        y(n2) = p_ini(i,2);
        uya(n2) = -nu*p_ini(i,2)*(P_ext/E);
%         uya1(n2) = -P*p_new(i,2)*(Pressure/E);
        uy(n2) = p_new(i,2) - p_ini(i,2);
    end
end

% plotting
figure(10)
plot(x,uxa,'--r',x,ux,'--b')
ylabel('U_x')
xlabel('x')
legend('Analytical','PD_solution')
xlim([p(1,1) p(Np,1)])

figure(11)
plot(y,uya,'--r',y,uy,'--b')
ylabel('U_y')
xlabel('y')
legend('Analytical','PD_solution')
xlim([p(1,2) p(Np,2)])
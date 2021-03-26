clear all; close all; clc;

% Lattice
Lx = 25/20; 
Ly = 0.5; % sqrt(3)*Ly
Lb = 1.0/20.0; % Lattice spacing 
tb = Lb; % Thickness of the plate
dA = Lb^2; % area of a paricle
dV = Lb^2*tb; % volume of a particle

[p, Np, Nx, Ny] = square_latttice(2,Lx,Ly,Lb); % reference coordinate
% [p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,Lb); % reference coordinate

horizon = 3.015*Lb; % 3.015*L ----- horizon

% inputs
% Rubber
rho = 1100; % density
C = 1e6; % damping
nu = 0.4; % poisson's ratio
E = 1e6; % Young's modulus 
lambda = E*nu/((1+nu)*(1-2*nu)); % Lame parameters
mu = E/(2*(1+nu));
K_bulk = lambda + 2*mu/3; % Bulk modulus
P_ext = 4e5;
F_ext = P_ext/Lb; % external force

dt = 0.0005; % time stepsize
Nt = 3000; % number of steps (iterations)
M = 500; % the time when the full amount of force is acting on a material body

%% Constructing initial bonds and shape tensor

n_crack = 6;
x_crack = 0;
y_crack = [p(1,2) p(1+(Nx+1)*n_crack/2,2) p((Nx+1)*(Ny-n_crack/2)+1,2) p((Nx+1)*Ny+1,2)];

[xi, v, omega, num_fam, fam_ind, fail] = initial(Np,Nx,p,Lb,horizon,x_crack,y_crack);
xi = sparse(xi); v = sparse(v); omega = sparse(omega); fam_ind = sparse(fam_ind);

figure(1)
clf
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on
hold on
for i = 1:Np
    for j = 1:num_fam(i)
        if fail(i,j) == 0
            x = [p(i,1) p(fam_ind(i,j),1)];
            y = [p(i,2) p(fam_ind(i,j),2)];
            figure(1)
            plot(x,y,'r--o','LineWidth',1.5)
            hold on
        end
    end
end

% figure(2)
% clf
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on
% hold on
% for i = 1:Np
%     for j = 1:num_fam(i)
%         if fail(i,j) == 1
%             x = [p(i,1) p(fam_ind(i,j),1)];
%             y = [p(i,2) p(fam_ind(i,j),2)];
%             figure(2)
%             plot(x,y,'b--o','LineWidth',1.5)
%             hold on
%         end
%     end
% end

%% Simulation
Ib = eye(2);

p_ini = p;
u_ini = zeros(Np,2); % initial speed

p_current = p_ini;
u_current = u_ini;

for ii = 0:Nt-1
    Time = dt*(ii+1)
    %% Dynamics using the advanced trapezoidal rule
    
    % apply the max-load at M-th step
    if ii < M
        nnn = (ii+1)/M;
    end
    
    [p_new, u_new, fail, sigma_vm] = Advanced_traperzoidal(p_current,u_current,p_ini,Np,Nx,Ny,Lx,Ly,Lb,dV,xi,num_fam,fam_ind,fail,v,omega,lambda,mu,P_ext,nnn,C,dt,rho);

    u_current = u_new;
    p_current = p_new;
    
%     u(:,2*(ii+1)+1:2*(ii+2)) = u_new;
    p(:,2*(ii+1)+1:2*(ii+2)) = p_new;
    
    VM_stress(:,ii+1) = sqrt(3.*sigma_vm);
end

figure(2)
scatter(p_new(:,1),p_new(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on
hold on

for i = 1:Np
    for j = 1:num_fam(i)
        if fail(i,j) == 0
            x = [p(i,1) p(fam_ind(i,j),1)];
            y = [p(i,2) p(fam_ind(i,j),2)];
            figure(2)
            plot(x,y,'r--o','LineWidth',1.5)
            hold on
        end
    end
end

% movie(p,ii,Lx,Ly)

N_step = 10000;
for i = 1:Nx+1
    for j = 1:Ny+1
        X(j,i) = p(i+(Nx+1)*(j-1),N_step*2+1);
        Y(j,i) = p(i+(Nx+1)*(j-1),N_step*2+2);
        X_ini(j,i) = p_ini(i+(Nx+1)*(j-1),1);
        Y_ini(j,i) = p_ini(i+(Nx+1)*(j-1),2);
        VM(j,i) = VM_stress(i+(Nx+1)*(j-1),N_step);
    end
end

figure(8)
pcolor(X,Y,VM)
title('Von Mises Stress at T = 1')
shading flat
colorbar

figure(9)
pcolor(X_ini,Y_ini,VM)
title('Von Mises Stress at T = 1')
shading flat
colorbar

figure(11)
plot(1:ii+1,VM_stress((Nx+1)*(Ny+1)/2+1,:),'b-','LineWidth',1.5)
xlabel('N_{steps}')
ylabel('\sigma_{VM}')
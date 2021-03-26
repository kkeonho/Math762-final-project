clear all; close all; clc;

%% 2d plane Lx*(sqrt(3)*Ly)

Lx = 20; 
Ly = 10; % sqrt(3)*Ly
L = 1; % Lattice spacing

[p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,L); % Generate a 2d hexagonal lattice
dp0 = 0*p; 

E = 1;
m = 1;
nu = 0.2;
dt = 0.0001;
tb = L*5;

F_ext  = 1;

K = K(E,L,tb,nu);

Max_Iter = 10000;

for i=1:Max_Iter
    
    if i<=1000
        F_ext = i/1000;
    end

%% Forward Euler

%     S = stretch_matrix_1(p(:,2*i-1:2*i),Np,Nx,Ny,L);
%     F = force_matrix(K,S);
%     
%     [F_net_x, F_net_y] = net_force(F,F_ext,Nx,Ny);
% 
%     if i == 1
%         p(:,1+2*i) = F_net_x'.*dt^2/m + dp0(:,1).*dt + p(:,1+2*(i-1));
%         p(:,2*(i+1)) = F_net_y'.*dt^2/m + dp0(:,2).*dt + p(:,2*i);
%     else
%         p(:,1+2*i) = F_net_x'*dt^2/m + 2*p(:,1+2*(i-1)) - p(:,1+2*(i-2));
%         p(:,2*(i+1)) = F_net_y'*dt^2/m + 2*p(:,2*i) - p(:,2*(i-1));
%     end
% 

%% Advance scheme

    if i == 1
        U_current = dp0;
        chi_current = p;
    end
    
    %% stage 1
    S = stretch_matrix_1(chi_current,Np,Nx,Ny,L);
    F_net_current = force_matrix_1(K,S,F_ext,Nx,Ny);
%     F = force_matrix(K,S);
%     
%     [F_net_x, F_net_y] = net_force(F,F_ext,Nx,Ny);
%     F_net_current = [F_net_x' F_net_y'];
    
    
    U_new = U_current + dt/m*F_net_current;
    chi_new = chi_current + dt*U_new;
    
    %% stage 2
    S = stretch_matrix_1(chi_new,Np,Nx,Ny,L);
    F = force_matrix(K,S);
    
    [F_net_x, F_net_y] = net_force(F,F_ext,Nx,Ny);
    F_net_new = [F_net_x' F_net_y'];
        
    U_new = U_current + 0.5*dt/m*(F_net_current + F_net_new);
    chi_new = chi_current + 0.5*dt*(U_new + U_current);
    
    p(:,2*i+1:2*(i+1)) = chi_new;

    U_current = U_new;
    chi_current = chi_new;
    
    

end

figure(2)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
scatter(p(:,1+2*i),p(:,2*(i+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

% movie(p,Max_Iter,Lx,Ly)




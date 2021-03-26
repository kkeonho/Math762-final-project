clear all; close all; clc;

%% 2d plane Lx*(sqrt(3)*Ly)

Lx = 20; 
Ly = 10; % sqrt(3)*Ly
L = 1; % Lattice spacing

[p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,L); % Generate a 2d hexagonal lattice
U_current = 0*p;
Chi_current = p;

E = 1;
m = 1;
nu = 0.3;
dt = 0.1;
tb = L*5;

F_ext  = 1;

K = K(E,L,tb,nu);

Max_Iter = 10000;

for i=1:Max_Iter
    
    if i<=5000
        F_ext = i/5000;
    end
    
    [Chi_new, U_new] = Advanced_traperzoidal(Chi_current,U_current,Np,Nx,Ny,L,K,F_ext,dt,m);
    
    p(:,2*i+1:2*(i+1)) = Chi_new;
    U_current = U_new;
    Chi_current = Chi_new;
end

figure(2)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
scatter(p(:,1+2*i),p(:,2*(i+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

figure(3)
scatter(p(:,1+2*i),p(:,2*(i+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

% movie(p,Max_Iter,Lx,Ly)




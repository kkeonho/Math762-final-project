clear all; clf; clc;

%% 2d plane Lx*(sqrt(3)*Ly)

Lx = 10; 
Ly = 10; % sqrt(3)*Ly
L = 1; % Lattice spacing
E = 1;
m = 1;
nu = 0.1;
dt = 0.01;
tb = L*5;
F_ext  = 1;

Max_Iter = 1000;

p = SPLM_2d(Lx,Ly,L,E,m,nu,dt,tb,F_ext,Max_Iter);
% 
% figure(1)
% scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)

figure(2)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
scatter(p(:,1+2*Max_Iter),p(:,2*(Max_Iter+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

% figure(3)
% scatter(p(:,1+2*Max_Iter),p(:,2*(Max_Iter+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
% grid on

% movie(p,Max_Iter,Lx,Ly)




clear all; close all; clc;

% Lattice
Lx = 10; 
Ly = 10; % sqrt(3)*Ly
Lb = 1; % Lattice spacing
tb = 0.1; % Thickness of the plate
dA = Lb^2; % area of a paricle
dV = dA*tb; % volume of a particle

[p, Np, Nx, Ny] = square_latttice(2,Lx,Ly,Lb); % reference coordinate
% [p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,Lb); % reference coordinate
u_current = zeros(Np,2); % initial speed

delta = 2.015*Lb; % 3.015*L ----- horizon

% inputs
rho = 1100; % density
C = 1; % damping
P = 1/3; % poisson's ratio
E = 1e6; % Young's modulus
lambda = E*P*((1+P)*(1-2*P)); % Lame parameters
mu = E/(2*(1+P));
F_ext = 1;

% t0 = 0; % initial time
% tEnd = .05; % end time
% Nt = 5000; % number of steps (iterations)
% dt = (tEnd - t0)/Nt; % time stepsize

dt = 0.001; % time stepsize
Nt = 1000; % number of steps (iterations)


%% Constructing initial bonds and shape tensor

[xi, K, num_fam, fam_ind] = initial_configuration(Np,p,Lb,dA);
for i = 1:Np
    num_fam(i) = 0;
    B = zeros(2,2);
    for j = 1:Np
        if i~=j
            ini_dist = sqrt((p(j,1) - p(i,1))^2 + (p(j,2) - p(i,2))^2); % initial distance
            w = influence_func(ini_dist,Lb);
            if w~=0
                omega(i,j) = w; % influence function
                v(i,j) = vol_frac(Lb,ini_dist,delta); % volume fraction
                
                num_fam(i) = num_fam(i) + 1; % number of particles in the horizon
                fam_ind(i,num_fam(i)) = j; % indeces of particles in the horizon
                
                xi(num_fam(i),2*i-1:2*i) = [(p(j,1)-p(i,1)) (p(j,2)-p(i,2))]; % reference bonds xi
                
                B = B + v(i,j)*dA*omega(i,j)*xi(num_fam(i),2*i-1:2*i)'*xi(num_fam(i),2*i-1:2*i);
            end
        end
    end
    
    K(2*i-1:2*i,:) = B; % Shape tensor at each particle
end

%% Simulation
Ib = eye(2);

for ii = 0:Nt-1
%     if ii <1000
%         nn = ii+1;
%     end

    
    %% Constructing deformed bonds,non-local deformation gradient tensor, strain tensor, and the first Piola-Kirchoff tensor
    for j = 1:Np
        FF = zeros(2,2);
        for k = 1:num_fam(j)
            % Deformed bond
            eta(k,2*j-1:2*j) =[(p(fam_ind(j,k),2*ii+1)-p(j,2*ii+1)) (p(fam_ind(j,k),2*(ii+1))-p(j,2*(ii+1)))]; % deformed bonds eta

            % Deformation gradient tensor
            FF = FF + v(j,fam_ind(j,k))*dA*omega(j,fam_ind(j,k))*eta(k,2*j-1:2*j)'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
        end
        
        % strain tensor
        EE = 0.5*(FF'*FF - Ib); 
        trEE = trace(EE); 
        
        % First Piola-Kirchoff tensor
        PPK1 = lambda * trEE * FF + 2 * mu * FF * EE; 
        
%         F(2*j-1:2*j,:) = FF;
        PK1(2*j-1:2*j,:) = PPK1;
    end
    
    %% Force acting on particles
    F_net = zeros(Np,2);
    for j  = 1:Np
        trac = zeros(2,1);
        for k = 1:num_fam(j)
            % compute the traction vector at each point
            trac = trac + omega(j,fam_ind(j,k))*[PK1(2*j-1:2*j,:)*(K(2*j-1:2*j,:))^-1 + PK1(2*fam_ind(j,k)-1:2*fam_ind(j,k),:)*(K(2*fam_ind(j,k)-1:2*fam_ind(j,k),:))^-1]*xi(k,2*j-1:2*j)';
        end
        
        % internal force
        F_net(j,:) = trac'*dA; 
        
        % compute the net-force at each particle; each side of the plane
        % is pulling 
%         if p(j,1) == -0.5 
%             F_net(j,1) = F_net(j,1) - F_ext;
%         elseif p(j,1) == 0.5 
%             F_net(j,1) = F_net(j,1) + F_ext;
%         end
        
        if p(j,1) == -Lx/2 
            F_net(j,1) = F_net(j,1) - F_ext;
        elseif p(j,1) == Lx/2 
            F_net(j,1) = F_net(j,1) + F_ext;
        end
    end

    %% Dynamics
    
    u_new = (1 - C*dt/rho)*u_current + dt/rho*F_net;
    p(:,2*(ii+1)+1:2*(ii+2)) = p(:,2*ii+1:2*ii+2) + 0.5*dt*(u_new + u_current);

    u_current = u_new;
end

figure(1)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

figure(2)
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
scatter(p(:,2*ii+1),p(:,2*(ii+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on

figure(5)
scatter(p(:,2*ii+1),p(:,2*(ii+1)),'MarkerEdgeColor','r','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
grid on
% 
% figure(6)
% movie(p,Nt,Lx,Ly)
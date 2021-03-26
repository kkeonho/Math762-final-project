clear all; close all; clc;

Lx = 1; 
Ly = 1; % sqrt(3)*Ly
Lb = .1; % Lattice spacing

delta = 2.015*Lb; % 3.015*L ----- horizon
 

[p, Np, Nx, Ny] = hexa_latttice(2,Lx,Ly,Lb);


%initial
 
for i = 1:Np
    num_fam(i) = 0;
    for j = 1:Np
        
        ini_dist = sqrt(power(p(j,1) - p(i,1),2) + power(p(j,2) - p(i,2),2)); % initial distance
        
        w = influence_func(ini_dist,Lb);
        if i~=j
            if w~=0
                omega(i,j) = w; % influence function
                v(i,j) = vol_frac(Lb,ini_dist,delta); % volume fraction
                
                num_fam(i) = num_fam(i) + 1; % number of particles in the horizon
                fam_ind(i,num_fam(i)) = j; % indeces of particles in the horizon
                
                xi(num_fam(i),2*i-1:2*i) = [(p(j,1)-p(i,1)) (p(j,2)-p(i,2))]; % reference bonds xi
            end
        end
    end
    
    B = zeros(2,2);
    for j = 1:num_fam(i)
        B = B + v(i,fam_ind(i,j))*omega(i,fam_ind(i,j))*xi(j,2*i-1:2*i)'*xi(j,2*i-1:2*i);
    end
    K(2*i-1:2*i,2*i-1:2*i) = B;
end



figure(1) % initial configuration
scatter(p(:,1),p(:,2),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)

lambda = 1; % lame parameters
mu = 1;
I = eye(Np*2);
Num_Iter = 1;
for i = 1:Num_Iter
    p(:,2*i+1) = p(:,2*i-1) + 2*p(:,2*i); % deformation
    p(:,2*(i+1)) = 3*p(:,2*i-1) + p(:,2*i);
    
%     p(:,2*i+1) = p(:,2*i-1).^3 + 2*p(:,2*i);
%     p(:,2*(i+1)) = 3*p(:,2*i-1) + p(:,2*i).^3;
    
    for j = 1:Np
        FF = zeros(2,2);
        for k = 1:num_fam(j)
            eta(k,2*j-1:2*j) =[(p(fam_ind(j,k),2*i+1)-p(j,2*i+1)) (p(fam_ind(j,k),2*(i+1))-p(j,2*(i+1)))]; % deformed bonds eta

            % Deformation gradient tensor
            FF = FF + v(j,fam_ind(j,k))*omega(j,fam_ind(j,k))*eta(k,2*j-1:2*j)'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,2*j-1:2*j))^-1;
        end
        F(2*j-1:2*j,2*j-1:2*j) = FF;
        
%         E = (FF'*FF - I)/2;
%         trE = trace(E);
%         
%         PK1 = lambda * trE * FF + 2 * mu * FF * E;
%         
%         trac = W*[ PK1(X) * B^-1(X) + PK1(X?) * B^-1(X?) ]  (X? - X)
    end    
    
    E = (F'*F - I)/2;
    trE = trace(E);
    figure(2)
    scatter(p(:,2*i+1),p(:,2*(i+1)),'MarkerEdgeColor','b','MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
end

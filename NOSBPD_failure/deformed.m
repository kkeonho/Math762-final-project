function [F_net, fail, sigma_vm] = deformed(p_ini,p,Np,Nx,Ny,Lx,Ly,Lb,dA,xi,num_fam,fam_ind,fail,v,omega,lambda,mu,P,nnn)

kbulk = lambda + 2/3*mu;
cstr = 3.1;

PK1 = zeros(2*Np,2);
K = zeros(2*Np,2);
sigma_vm = zeros(Np,1);

%% Updating the fail parameter
    for j = 1:Np
        for k = 1:num_fam(j)
            R = sqrt((p(fam_ind(j,k),1) - p(j,1))^2 + (p(fam_ind(j,k),2) - p(j,2))^2); % deformed length
            R0 = sqrt((p_ini(fam_ind(j,k),1) - p_ini(j,1))^2 + (p_ini(fam_ind(j,k),2) - p_ini(j,2))^2); % initial length
            str = R/R0;
            
            if (str > cstr) && (fail(j,k) == 1)
                fail(j,k) = 0;
                for l=1:num_fam(fam_ind(j,k))
                    if fam_ind(fam_ind(j,k),l) == j
                        fail(fam_ind(j,k),l) = 0;
                        break
                    end
                end
            end
        end
    end
    
%% Constructing deformed bonds,non-local deformation gradient tensor, strain tensor, and the first Piola-Kirchoff tensor    
    for i = 1:Np
        FF = zeros(2,2);
        B = zeros(2,2);
        PPK1 = zeros(2,2);
        
        for j = 1:num_fam(i)
            eta = zeros(2,2);
            eta =[(p(fam_ind(i,j),1)-p(i,1)) (p(fam_ind(i,j),2)-p(i,2))]; % deformed bond eta
            B = B + fail(i,j)*v(i,j)*dA*omega(i,j)*xi(j,2*i-1:2*i)'*xi(j,2*i-1:2*i); % reconstructing the shape tensor with the updated fail
            FF = FF + fail(i,j)*v(i,j)*dA*omega(i,j)*eta'*xi(j,2*i-1:2*i); % reconstructing the deformation gradient tensor with the updated fail
        end
        FF = FF*B^-1; % nonlocal deformation gradient tensor
        
        C = transpose(FF)*FF; % Rigth Cauchy-Green tensor
        TFF = transpose(FF); % traspose of the nonlocal deformation gradient tensor
        trC = trace(C); % trace
        J = det(FF); % Jacobian
        
        PPK1 = mu*(FF - 1/2*trC*TFF^-1)/J + kbulk/4*(J^2 - 1/J^2)*TFF^-1; %% First Piola-Kirchoff stress tensor
        
        sigma = 1/J*PPK1*TFF;
        sigma_vm(i) = 1/2*(trace(sigma^2) - 1/3*trace(sigma)^2);
        
        PK1(2*i-1:2*i,:) = PPK1;
        K(2*i-1:2*i,:) = B;
    end
    
%% Net force acting on particles
    
    F_net = net_force_NOSB(p_ini,Np,num_fam,fam_ind,fail,omega,v,PK1,K,xi,Lb,dA,Lx,Ly,Nx,Ny,P,nnn);

end


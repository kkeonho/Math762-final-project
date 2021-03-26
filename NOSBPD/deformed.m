function [F_net] = deformed(p_ini,p,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P,nnn)

PK1 = zeros(2*Np,2);
%% Constructing deformed bonds,non-local deformation gradient tensor, strain tensor, and the first Piola-Kirchoff tensor
    for j = 1:Np
        FF = zeros(2,2);
        eta = zeros(2,2);
        EE = zeros(2,2);
        PPK1 = zeros(2,2);
        
        for k = 1:num_fam(j)
            % Deformed bond
%             eta(k,2*j-1:2*j) =[(p(fam_ind(j,k),1)-p(j,1)) (p(fam_ind(j,k),2)-p(j,2))]; % deformed bonds eta
            eta =[(p(fam_ind(j,k),1)-p(j,1)) (p(fam_ind(j,k),2)-p(j,2))]; % deformed bonds eta

            % Deformation gradient tensor
%             FF = FF + v(j,fam_ind(j,k))*dA*omega(j,fam_ind(j,k))*eta(k,2*j-1:2*j)'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
            FF = FF + v(j,fam_ind(j,k))*dV*omega(j,fam_ind(j,k))*eta'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
            
            R = sqrt((p(k,1) - p(j,1))^2 + (p(k,2) - p(j,2))^2);
            R0 = sqrt((p_ini(k,1) - p_ini(j,1))^2 + (p_ini(k,2) - p_ini(j,2))^2);
        end

        % strain tensor
%         EE = 0.5*(FF'*FF - Ib); 
        EE = 0.5*(FF'+FF) - Ib; % linearized
        trEE = trace(EE); 
        
        % First Piola-Kirchoff tensor
%         PPK1 = lambda * trEE * FF + 2 * mu * FF * EE; 
        PPK1 = lambda * trEE * eye(2) + 2.0 * mu * EE; % linearized
        
        J = det(FF);
        
        sigma = 1/J*FF*transpose(PPK1);
        F(2*j-1:2*j,:) = FF;
        PK1(2*j-1:2*j,:) = PPK1;
    end
    
    %% Net force acting on particles
    
    F_net = net_force_NOSB(p_ini,Np,num_fam,fam_ind,omega,v,PK1,K,xi,Lb,dV,Lx,Ly,P,nnn);

end


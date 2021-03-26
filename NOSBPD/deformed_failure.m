function [F_net, fail] = deformed_failure(p_ini,p,Np,Lx,Ly,Lb,dA,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,kbulk,P,nnn,cstr,fail,horizon)


PK1 = zeros(2*Np,2);
%% Constructing deformed bonds,non-local deformation gradient tensor, strain tensor, and the first Piola-Kirchoff tensor
    for j = 1:Np
        FF = zeros(2,2);
        C = zeros(2,2);
        TFF = zeros(2,2);
        eta = zeros(2,2);
        EE = zeros(2,2);
        PPK1 = zeros(2,2);
        B = zeros(2,2);
        
        for k = 1:num_fam(j)
            % Deformed bond
%             eta(k,2*j-1:2*j) =[(p(fam_ind(j,k),1)-p(j,1)) (p(fam_ind(j,k),2)-p(j,2))]; % deformed bonds eta
            eta =[(p(fam_ind(j,k),1)-p(j,1)) (p(fam_ind(j,k),2)-p(j,2))]; % deformed bonds eta
            
            R = sqrt((p(k,1) - p(j,1))^2 + (p(k,2) - p(j,2))^2);
%             omega(j,k) = influence_func(R,Lb,horizon);
            R0 = sqrt((p_ini(k,1) - p_ini(j,1))^2 + (p_ini(k,2) - p_ini(j,2))^2);
            str = (R)/R0;
            
            if (str > cstr) && (fail(j,k) == 1)
                fail(j,k) = 0;
                for l=1:num_fam(fam_ind(j,k))
                    if fam_ind(fam_ind(j,k),l) == j
                        fail(fam_ind(j,k),l) = 0;
                        break
                    end
                end
            end
            
%             B = B + fail(j,k)*v(j,fam_ind(j,k))*dA*omega(j,fam_ind(j,k))*xi(k,2*j-1:2*j)'*xi(k,2*j-1:2*j);
            % Deformation gradient tensor
%             FF = FF + v(j,fam_ind(j,k))*dA*omega(j,fam_ind(j,k))*eta(k,2*j-1:2*j)'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
%             FF = FF + fail(j,k)*v(j,fam_ind(j,k))*dA*omega(j,fam_ind(j,k))*eta'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
            FF = FF + fail(j,k)*v(j,k)*dA*omega(j,k)*eta'*xi(k,2*j-1:2*j)*(K(2*j-1:2*j,:))^-1;
        end

%         % strain tensor
% %         EE = 0.5*(FF'*FF - Ib); 
%         EE = 0.5*(FF'+FF) - Ib; % linearized
%         trEE = trace(EE);
% 
% %         C = transpose(FF)*FF;
% %         trc = trace(C);
% %         J = det(FF);
% %         
% %         TFF = transpose(FF);
% %         
% %         PPK1 = mu*(FF - 1/2*trc*TFF^-1)/J + kbulk/4*(J^2 - 1/J^2)*TFF^-1;
%         
%         % First Piola-Kirchoff tensor
% %         PPK1 = lambda * trEE * FF + 2 * mu * FF * EE; 
%         PPK1 = lambda * trEE * eye(2) + 2 * mu * EE; % linearized

        C = transpose(FF)*FF;
        TFF = transpose(FF);
        trC = trace(C);
        J = det(FF);
        
        PPK1 = mu*(FF - 1/2*trC*TFF^-1)/J + kbulk/4*(J^2 - 1/J^2)*TFF^-1;
        
%         sigma = 1/J*FF*transpose(PPK1);
%         F(2*j-1:2*j,:) = FF;
        PK1(2*j-1:2*j,:) = PPK1;
%         K(2*j-1:2*j,:) = B;
    end
    
    %% Net force acting on particles
    
    F_net = net_force_NOSB_failure(p_ini,Np,num_fam,fam_ind,omega,v,PK1,K,xi,Lb,dA,Lx,Ly,P,nnn,fail);

end


function F_net = net_force_NOSB(p_ini,Np,num_fam,fam_ind,omega,v,PK1,K,xi,Lb,dV,Lx,Ly,P,nnn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F_net = zeros(Np,2);

N = [0;1];
    for j  = 1:Np
        
        
        trac = zeros(2,1);
        for k = 1:num_fam(j)
            % compute the traction vector at each point
            trac = trac + v(j,fam_ind(j,k))*omega(j,fam_ind(j,k))*(PK1(2*j-1:2*j,:)*(K(2*j-1:2*j,:))^-1 + PK1(2*fam_ind(j,k)-1:2*fam_ind(j,k),:)*(K(2*fam_ind(j,k)-1:2*fam_ind(j,k),:))^-1)*xi(k,2*j-1:2*j)';
%             trac = trac + omega(j,fam_ind(j,k))*(PK1(2*j-1:2*j,:)*(K(2*j-1:2*j,:))^-1 + PK1(2*fam_ind(j,k)-1:2*fam_ind(j,k),:)*(K(2*fam_ind(j,k)-1:2*fam_ind(j,k),:))^-1)*xi(k,2*j-1:2*j)';
        end
        
        % internal force
        F_net(j,:) = trac'*dV;
        
        %% plane strain
        % compute the net-force at each particle; each side of the plane is pulling 
        if p_ini(j,1) == -Lx/2 
            F_net(j,1) = F_net(j,1) - P/Lb*nnn;
        elseif p_ini(j,1) == Lx/2 
            F_net(j,1) = F_net(j,1) + P/Lb*nnn;
        end

        %% Shearing
        
%         if p_ini(j,1) == Ly 
%             F_net(j,2) = F_net(j,2) + P/Lb*nnn;
%         end

    end
end


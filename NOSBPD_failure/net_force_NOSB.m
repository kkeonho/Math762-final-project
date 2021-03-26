function F_net = net_force_NOSB(p_ini,Np,num_fam,fam_ind,fail,omega,v,PK1,K,xi,Lb,dA,Lx,Ly,Nx,Ny,P,nnn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
F_net = zeros(Np,2);

N = [0;1];

indi1 = 0;
indi2 = 0;
for i = 1:num_fam((Nx+1)*(Ny+1)/2)
    if (p_ini(fam_ind((Nx+1)*(Ny+1)/2,i),1)>0)&&(fail((Nx+1)*(Ny+1)/2,i)==0)
        indi1 = indi1 +1;
    end
end

for i = 1:num_fam((Nx+1)*(Ny+1)/2+1)
    if (p_ini(fam_ind((Nx+1)*(Ny+1)/2+1,i),1)<0)&&(fail((Nx+1)*(Ny+1)/2+1,i)==0)
        indi2 = indi2 +1;
    end
end

fail_ind = 1;
if (indi1 == 11)&&(indi2 == 11)
    fail_ind = 0;
end



    for j  = 1:Np
        
        %% plane strain
        trac = zeros(2,1);
        for k = 1:num_fam(j)
            % compute the traction vector at each point
            trac = trac + fail(j,k)*v(j,k)*omega(j,k)*(PK1(2*j-1:2*j,:)*(K(2*j-1:2*j,:))^-1 + PK1(2*fam_ind(j,k)-1:2*fam_ind(j,k),:)*(K(2*fam_ind(j,k)-1:2*fam_ind(j,k),:))^-1)*xi(k,2*j-1:2*j)';
        end
        
        % internal force
        F_net(j,:) = trac'*dA;
        
        % compute the net-force at each particle; each side of the plane is pulling 
        if p_ini(j,1) == -Lx/2 
            F_net(j,1) = F_net(j,1);
        elseif p_ini(j,1) == Lx/2 
            F_net(j,1) = F_net(j,1) + P/Lb*nnn*fail_ind;
        end

    end
end


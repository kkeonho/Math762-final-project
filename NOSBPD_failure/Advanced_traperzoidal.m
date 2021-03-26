function [p_new, u_new, fail_new, sigma_vm] = Advanced_traperzoidal(p_current,u_current,p_ini,Np,Nx,Ny,Lx,Ly,Lb,dV,xi,num_fam,fam_ind,fail,v,omega,lambda,mu,P_ext,nnn,C,dt,rho)

    %% stage 1
    [F_net_current, fail_current, sigma_vm] = deformed(p_ini,p_current,Np,Nx,Ny,Lx,Ly,Lb,dV,xi,num_fam,fam_ind,fail,v,omega,lambda,mu,P_ext,nnn);
    
    u_new = (1 - C*dt/rho)*u_current + dt/rho*F_net_current;
    for i = 1:Np
        if p_ini(i,1) == -Lx/2
            u_new(i,:) = 0;
        elseif p_ini(i,1) == Lx/2
            u_new(i,2) = 0;
        end
    end
    p_new = p_current + dt*u_new;
    
    %% stage 2
    [F_net_new, fail_new, sigma_vm] = deformed(p_ini,p_new,Np,Nx,Ny,Lx,Ly,Lb,dV,xi,num_fam,fam_ind,fail_current,v,omega,lambda,mu,P_ext,nnn);

    u_new = (1 - C*dt/rho)*u_current + 0.5*dt/rho*(F_net_current + F_net_new);
    for i = 1:Np
        if p_ini(i,1) == -Lx/2
            u_new(i,:) = 0;
        elseif p_ini(i,1) == Lx/2
            u_new(i,2) = 0;
        end
    end
    p_new = p_current + 0.5*dt*(u_new + u_current);
end
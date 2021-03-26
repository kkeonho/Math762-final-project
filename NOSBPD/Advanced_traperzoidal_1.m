function [p_new, u_new] = Advanced_traperzoidal_1(p_current,u_current,p_ini,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn,C,dt,rho)

    %% Midpoint
    F_net_current = deformed(p_ini,p_current,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
    
    u_new = (1 - C*dt/rho)*u_current + dt/rho*F_net_current;
    p_new = p_current + dt/2*(u_new+u_current);
    
    %% Trapezoidal 
%     %% stage 1
%     F_net_current = deformed(p_ini,p_current,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
%     
%     u_new = (1 - C*dt/rho)*u_current + dt/rho*F_net_current;
%     p_new = p_current + dt*u_new;
%     
%     %% stage 2
%     F_net_new = deformed(p_ini,p_new,Np,Lx,Ly,Lb,dV,xi,K,num_fam,fam_ind,v,omega,Ib,lambda,mu,P_ext,nnn);
% 
%     u_new = (1 - C*dt/rho)*u_current + 0.5*dt/rho*(F_net_current + F_net_new);
%     p_new = p_current + 0.5*dt*(u_new + u_current);
    
    %% plane shearing
%     for i = 1:Np
%         if p_ini(i,1) == 0
%             p_new(i,:) = p_ini(i,:);
%         end
%     end
end
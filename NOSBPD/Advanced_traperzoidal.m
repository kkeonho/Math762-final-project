function [Chi_new, U_new] = Advanced_traperzoidal(Chi_current,U_current,Np,Nx,Ny,L,K,F_ext,dt,m)


    %% stage 1
    S = stretch_matrix_1(Chi_current,Np,Nx,Ny,L);
    F_net_current = force_matrix_1(K,S,F_ext,Nx,Ny);
    
    U_new = U_current + dt/m*F_net_current;
    Chi_new = Chi_current + dt*U_new;
    
    %% stage 2
    S = stretch_matrix_1(Chi_new,Np,Nx,Ny,L);
    F_net_new = force_matrix_1(K,S,F_ext,Nx,Ny);
        
    U_new = U_current + 0.5*dt/m*(F_net_current + F_net_new);
    Chi_new = Chi_current + 0.5*dt*(U_new + U_current);
end


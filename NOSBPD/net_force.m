function [F_net_x, F_net_y] = net_force(F,F_ext,Nx,Ny)

    F_net_x = F(1,:)-F(2,:) + (F(3,:) - F(4,:) - F(5,:) + F(6,:))*1/2;
    F_net_y = (F(3,:) - F(4,:) + F(5,:) - F(6,:))*sqrt(3)/2;
    
    
    for j=1:Ny+1
        if j==Ny+1
            F_net_x(1+(2*Nx+1)*(j-1)) = F_net_x(1+(2*Nx+1)*(j-1)) - F_ext;
            
%             F_net_x(1+(2*Nx+1)*(j-1)) = 0;
%             F_net_y(1+(2*Nx+1)*(j-1)) = 0;
            
            F_net_x(Nx+1+(2*Nx+1)*(j-1)) = F_net_x(Nx+1+(2*Nx+1)*(j-1)) + F_ext;
            
%             F_net_y(Nx+1+(2*Nx+1)*(j-1)) = F_net_y(Nx+1+(2*Nx+1)*(j-1)) + F_ext;
        else
            F_net_x(1+(2*Nx+1)*(j-1)) = F_net_x(1+(2*Nx+1)*(j-1)) - F_ext;
            F_net_x(1+(2*Nx+1)*(j-1)+Nx+1) = F_net_x(1+(2*Nx+1)*(j-1)+Nx+1) - F_ext;

%             F_net_x(1+(2*Nx+1)*(j-1)) = 0;
%             F_net_x(1+(2*Nx+1)*(j-1)+Nx+1) = 0;
%             F_net_y(1+(2*Nx+1)*(j-1)) = 0;
%             F_net_y(1+(2*Nx+1)*(j-1)+Nx+1) = 0;
            
            F_net_x(Nx+1+(2*Nx+1)*(j-1)) = F_net_x(Nx+1+(2*Nx+1)*(j-1)) + F_ext;
            F_net_x((2*Nx+1)*(j-1)+2*Nx+1) = F_net_x((2*Nx+1)*(j-1)+2*Nx+1) + F_ext;
%             
%             F_net_y(Nx+1+(2*Nx+1)*(j-1)) = F_net_y(Nx+1+(2*Nx+1)*(j-1)) + F_ext;
%             F_net_y((2*Nx+1)*(j-1)+2*Nx+1) = F_net_y((2*Nx+1)*(j-1)+2*Nx+1) + F_ext;
        end
    end
    
%     F_net_y(1:Nx+1) = 0;
%     F_net_y(1+(2*Nx+1)*(Ny):1+(2*Nx+1)*(Ny)+Nx) = 0;
end



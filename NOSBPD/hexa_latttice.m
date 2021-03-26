function [p, Np, Nx, Ny] = hexa_latttice(dim,Lx,Ly,Lb)

Nx = Lx/Lb;
Ny = Ly/Lb;

Np = (2*Nx+1)*Ny + Nx+1;

% Np = (Nx+1)*(Ny+1)+Nx*Ny;

p = zeros(Np,2);

if dim==2
for j=1:Ny+1
    for i=1:Nx+1
%         p(i+(2*Nx+1)*(j-1),1) = Lb*(i-1);
%         p(i+(2*Nx+1)*(j-1),2) = sqrt(3)*Lb*(j-1);
        p(i+(2*Nx+1)*(j-1),1) = -Lx/2 + Lb*(i-1);
        p(i+(2*Nx+1)*(j-1),2) = -Ly*sqrt(3)/2 + sqrt(3)*Lb*(j-1);
        
        if (j<Ny+1) && (i<Nx+1)
%         p(i+(Nx+1)+(2*Nx+1)*(j-1),1) = Lb/2+Lb*(i-1);
%         p(i+(Nx+1)+(2*Nx+1)*(j-1),2) = sqrt(3)/2*Lb + sqrt(3)*Lb*(j-1);
        p(i+(Nx+1)+(2*Nx+1)*(j-1),1) = -Lx/2 + Lb/2+Lb*(i-1);
        p(i+(Nx+1)+(2*Nx+1)*(j-1),2) = -Ly*sqrt(3)/2 + sqrt(3)/2*Lb + sqrt(3)*Lb*(j-1);
        end
    end
end
elseif dim==3
end
end


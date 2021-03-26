function S = stretch(p,Np,Nx,Ny)

S = zeros(Np,6);

for j=1:Ny+1
    if j==1
        for i=1:2*Nx+1
            if i==1
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = 0;
                S2 = 0;
                S4 = 0;
                S6 = 0;
            elseif i==Nx+1
                S1 = 0;
                S3 = 0;
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = 0;
                S6 = 0;
            elseif (i<Nx+1)&&(i>1)
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = 0;
                S6 = 0;
            elseif i==Nx+2
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = 0;
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            elseif (i>Nx+2)&&(i<2*Nx+1)
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            elseif i==2*Nx+1
                S1 = 0;
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            end
            S(i+(2*Nx+1)*(j-1),:) = [S1 S2 S3 S4 S5 S6];
        end
    elseif j==Ny+1
        for i=1:Nx+1
            if i==1
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = 0;
                S5 = 0;
                S2 = 0;
                S4 = 0;
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            elseif i==Nx+1
                S1 = 0;
                S3 = 0;
                S5 = 0;
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = 0;
            else
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = 0;
                S5 = 0;
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            end
            S(i+(2*Nx+1)*(j-1),:) = [S1 S2 S3 S4 S5 S6];
        end
    else
        for i=1:2*Nx+1
            if i==1
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = 0;
                S2 = 0;
                S4 = 0;
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            elseif i==Nx+1
                S1 = 0;
                S3 = 0;
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = 0;
            elseif i==Nx+2
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = 0;
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            elseif i==2*Nx+1
                S1 = 0;
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            else
                S1 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+1,:));
                S3 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+(Nx+1),:));
                S5 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)+Nx,:));
                S2 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-1,:));
                S4 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx+1),:));
                S6 = norm(p(i+(2*Nx+1)*(j-1),:)-p(i+(2*Nx+1)*(j-1)-(Nx),:));
            end
            S(i+(2*Nx+1)*(j-1),:) = [S1 S2 S3 S4 S5 S6];
        end
    end
end
end


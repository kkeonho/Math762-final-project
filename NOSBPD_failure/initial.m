function [xi, v, omega, num_fam, fam_ind, fail] = initial(Np,Nx,p,Lb,horizon,x_crack,y_crack)

if horizon == 2.015*Lb
    max_pts = 12;
elseif horizon == 3.015*Lb
    max_pts = 28;
end

K = zeros(2*Np,2);
v = sparse(Np,max_pts);
fail  = zeros(Np,max_pts);
fam_ind = zeros(Np,max_pts);
xi = sparse(max_pts,2*Np);
omega = sparse(Np,max_pts);

%% Constructing initial bonds, weights, volume fraction and fail parameters
for i = 1:Np
    num_fam(i) = 0;
    B = zeros(2,2);
    for j = 1:Np
        if i~=j
            ini_dist = sqrt((p(j,1) - p(i,1))^2 + (p(j,2) - p(i,2))^2); % initial distance
            w = influence_func(ini_dist,Lb,horizon);
            if w~=0
                num_fam(i) = num_fam(i) + 1; % number of particles in the horizon
                fam_ind(i,num_fam(i)) = j; % indeces of particles in the horizon
                
                omega(i,num_fam(i)) = w; % influence function
                v(i,num_fam(i)) = vol_frac(Lb,ini_dist,horizon); % volume fraction
                xi(num_fam(i),2*i-1:2*i) = [(p(j,1)-p(i,1)) (p(j,2)-p(i,2))]; % reference bonds xi
                
                fail(i,num_fam(i)) = 1;
            end
        end
    end
end

%% Cut bonds across the cracks
for i = 1:Np
    B = zeros(2,2);
    for j = 1:num_fam(i)
        if ((p(i,1)<x_crack)&&(p(fam_ind(i,j),1)>x_crack))&&((p(i,2)<=y_crack(2))&&(p(fam_ind(i,j),2)<=y_crack(2))) || ((p(i,1)<x_crack)&&(p(fam_ind(i,j),1)>x_crack))&&((p(i,2)>=y_crack(3))&&(p(fam_ind(i,j),2)>=y_crack(3)))
            fail(i,j) = 0;
            
            for k = 1:num_fam(fam_ind(i,j))
                if fam_ind(fam_ind(i,j),k) == i
                    fail(fam_ind(i,j),k) = 0;
                    break
                end
            end
        end
        
        if horizon == 3.015*Lb
            if ((p(i,1)==p((Nx+1)/2,1))&&(p(fam_ind(i,j),1)>x_crack))&&((p(i,2)<y_crack(2))&&(p(fam_ind(i,j),2)<=(y_crack(2)+0.00001+Lb))) || ((p(i,1)==p((Nx+1)/2+1,1))&&(p(fam_ind(i,j),1)<x_crack))&&((p(i,2)<y_crack(2))&&(p(fam_ind(i,j),2)<=(y_crack(2)+0.00001+Lb))) || ((p(i,1)==p((Nx+1)/2,1))&&(p(fam_ind(i,j),1)>x_crack))&&((p(i,2)>y_crack(3))&&(p(fam_ind(i,j),2)>=(y_crack(3)-0.00001-Lb))) || ((p(i,1)==p((Nx+1)/2+1,1))&&(p(fam_ind(i,j),1)<x_crack))&&((p(i,2)>y_crack(3))&&(p(fam_ind(i,j),2)>=(y_crack(3)-0.00001-Lb)))
                fail(i,j) = 0;
                
                for k = 1:num_fam(fam_ind(i,j))
                    if fam_ind(fam_ind(i,j),k) == i
                        fail(fam_ind(i,j),k) = 0;
                        break
                    end
                end
            end
        end
    end
end

% %% construct the inital shape tensor
% for i = 1:Np
%     B = zeros(2,2);
%     for j = 1:num_fam(i)
%         B = B + fail(i,j)*v(i,j)*dA*omega(i,j)*xi(j,2*i-1:2*i)'*xi(j,2*i-1:2*i);
%     end
%     K(2*i-1:2*i,:) = B;
% end

end


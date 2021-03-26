function [xi, K, v, omega, num_fam, fam_ind] = initial(Np,p,Lb,dA,horizon)

K = zeros(2*Np,2);

for i = 1:Np
    num_fam(i) = 0;
    B = zeros(2,2);
    for j = 1:Np
        if i~=j
            ini_dist = sqrt((p(j,1) - p(i,1))^2 + (p(j,2) - p(i,2))^2); % initial distance
            w = influence_func(ini_dist,Lb,horizon);
            if (ini_dist<horizon)
                omega(i,j) = w; % influence function
                v(i,j) = vol_frac(Lb,ini_dist,horizon); % volume fraction
                
                num_fam(i) = num_fam(i) + 1; % number of particles in the horizon
                fam_ind(i,num_fam(i)) = j; % indeces of particles in the horizon
                
                xi(num_fam(i),2*i-1:2*i) = [(p(j,1)-p(i,1)) (p(j,2)-p(i,2))]; % reference bonds xi
                
                % Shape tensor
                B = B + v(i,j)*dA*omega(i,j)*xi(num_fam(i),2*i-1:2*i)'*xi(num_fam(i),2*i-1:2*i);
            end
        end
    end
    
    K(2*i-1:2*i,:) = B; % Shape tensor at each particle
    invK(2*i-1:2*i,:) = inv(B);
end

end


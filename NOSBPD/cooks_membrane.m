function [p, Np, N_W, N_H] = cooks_membrane(H,HL,W,TD,Lb)

N_H = ceil(H/Lb);
N_W = ceil(W/Lb);
N_HL = ceil(N_H);

dx = W/N_W;
dy1 = H/N_H;
dy2 = HL/N_H;

x = linspace(0,W,ceil(N_W+1));

for i = 1: N_H + 1
    y1 = (i-1)*dy1;
    y2 = H + TD - HL + dy2*(i-1);
    
    y = linspace(y1,y2,N_W+1);
    for j = 1:N_W + 1
        p(j+(N_W+1)*(i-1),1) = x(j);
        p(j+(N_W+1)*(i-1),2) = y(j);
    end
    
end

Np = (N_W + 1)*(N_H + 1);

end
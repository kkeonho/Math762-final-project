function w = influence_func(S0,Lb,horizon)

A = 2;
C = 60.0 / (7.0 * pi * A * A);

% s = S0/Lb; % S0/horizon // 2*S0/horizon

s = 2*S0/horizon;

if s<1
    w = C*(2.0 / 3.0 - s^2 + 0.5*s^3);
elseif s<=2
    w = C*(2.0 - s)^3/6.0;
else
    w = 0.0;
end

end


function v = vol_frac(Lb,S0,horizon) 

r = Lb/2;
    if (S0 <= (horizon - r))
        v = 1.0;
    elseif (S0 <= (horizon))
        v = (horizon + r - S0) / (2.0 * r);
    else
        v = 0.0;
    end
    

%     if (S0 <= (horizon - Lb))
%         v = 1.0;
%     elseif (S0 <= (horizon +Lb ))
%         v = (horizon + Lb - S0) / (2.0 * Lb);
%     else
%         v = 0.0;
%     end
end

%%%% Ref: Amani_etal_IJIE_2015_A_non_ordinary_state_based_peridynamics_formulation.pdf




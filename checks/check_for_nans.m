% Check for nans and infs

if strcmp(check_this,'inc')

    qqq = sum(1 - isfinite(inc_p));
    if qqq > 0
        disp('Non-finite inc_p')
    end
    qqq = sum(1 - isfinite(inc_w1));
    if qqq > 0
        disp('Non-finite inc_w1')
    end
    qqq = sum(1 - isfinite(inc_w2));
    if qqq > 0
        disp('Non-finite inc_w2')
    end
    qqq = sum(1 - isfinite(inc_m1));
    if qqq > 0
        disp('Non-finite inc_m1')
    end
    qqq = sum(1 - isfinite(inc_m2));
    if qqq > 0
        disp('Non-finite inc_m2')
    end
    qqq = sum(1 - isfinite(inc_eta1));
    if qqq > 0
        disp('Non-finite inc_eta1')
    end
    qqq = sum(1 - isfinite(inc_eta2));
    if qqq > 0
        disp('Non-finite inc_eta2')
    end
    qqq = sum(1 - isfinite(inc_q1));
    if qqq > 0
        disp('Non-finite inc_q1')
    end
    qqq = sum(1 - isfinite(inc_q2));
    if qqq > 0
        disp('Non-finite inc_q2')
    end
    qqq = sum(1 - isfinite(inc_T1));
    if qqq > 0
        disp('Non-finite inc_T1')
    end
    qqq = sum(1 - isfinite(inc_T2));
    if qqq > 0
        disp('Non-finite inc_T2')
    end
    qqq = sum(1 - isfinite(inc_Tw1));
    if qqq > 0
        disp('Non-finite inc_Tw1')
    end
    qqq = sum(1 - isfinite(inc_Tw2));
    if qqq > 0
        disp('Non-finite inc_Tw2')
    end
    qqq = sum(1 - isfinite(inc_u1));
    if qqq > 0
        disp('Non-finite inc_u1')
    end
    qqq = sum(1 - isfinite(inc_u2));
    if qqq > 0
        disp('Non-finite inc_u2')
    end
    qqq = sum(1 - isfinite(inc_v1));
    if qqq > 0
        disp('Non-finite inc_v1')
    end
    qqq = sum(1 - isfinite(inc_v2));
    if qqq > 0
        disp('Non-finite inc_v2')
    end
    qqq = sum(1 - isfinite(inc_tke1));
    if qqq > 0
        disp('Non-finite inc_tke1')
    end
    qqq = sum(1 - isfinite(inc_tke2));
    if qqq > 0
        disp('Non-finite inc_tke2')
    end
    qqq = sum(1 - isfinite(inc_vareta1));
    if qqq > 0
        disp('Non-finite inc_vareta1')
    end
    qqq = sum(1 - isfinite(inc_vareta2));
    if qqq > 0
        disp('Non-finite inc_vareta2')
    end
    qqq = sum(1 - isfinite(inc_varq1));
    if qqq > 0
        disp('Non-finite inc_varq1')
    end
    qqq = sum(1 - isfinite(inc_varq2));
    if qqq > 0
        disp('Non-finite inc_varq2')
    end

elseif strcmp(check_this,'tend')

    qqq = sum(1 - isfinite(tend.fluid(1).m.tot));
    if qqq > 0
        disp('Non-finite tend_m1')
    end
    qqq = sum(1 - isfinite(tend.fluid(1).meta.tot));
    if qqq > 0
        disp('Non-finite tend_meta1')
        qqq = sum(1 - isfinite(tend.fluid(1).meta.transport));
        if qqq > 0
            disp('-- transport')
        end
        qqq = sum(1 - isfinite(tend.fluid(1).meta.diffuse));
        if qqq > 0
            disp('-- diffuse')
        end
        qqq = sum(1 - isfinite(tend.fluid(1).meta.diffent));
        if qqq > 0
            disp('-- diffent')
        end
        qqq = sum(1 - isfinite(tend.fluid(1).meta.bflux));
        if qqq > 0
            disp('-- bflux')
        end
        qqq = sum(1 - isfinite(tend.fluid(1).meta.dissn));
        if qqq > 0
            disp('-- dissn')
        end
        qqq = sum(1 - isfinite(tend.fluid(1).meta.relabel));
        if qqq > 0
            disp('-- relabel')
        end
    end
    qqq = sum(1 - isfinite(tend.fluid(1).mq.tot));
    if qqq > 0
        disp('Non-finite tend_mq1')
    end
    qqq = sum(1 - isfinite(tend.fluid(1).mw.tot));
    if qqq > 0
        disp('Non-finite tend_mw1')
    end
    qqq = sum(1 - isfinite(tend.fluid(2).m.tot));
    if qqq > 0
        disp('Non-finite tend_m2')
    end
    qqq = sum(1 - isfinite(tend.fluid(2).meta.tot));
    if qqq > 0
        disp('Non-finite tend_meta2')
    end
    qqq = sum(1 - isfinite(tend.fluid(2).mq.tot));
    if qqq > 0
        disp('Non-finite tend_mq2')
    end
    qqq = sum(1 - isfinite(tend.fluid(2).mw.tot));
    if qqq > 0
        disp('Non-finite tend_mw2')
    end

elseif strcmp(check_this,'state')

    % Check for nans and infs
    qqq = sum(1 - isfinite(state_new.p));
    if qqq > 0
        disp('Non-finite p')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).m));
    if qqq > 0
        disp('Non-finite m1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).m));
    if qqq > 0
        disp('Non-finite m2')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).w));
    if qqq > 0
        disp('Non-finite w1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).w));
    if qqq > 0
        disp('Non-finite w2')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).eta));
    if qqq > 0
        disp('Non-finite eta1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).eta));
    if qqq > 0
        disp('Non-finite eta2')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).q));
    if qqq > 0
        disp('Non-finite q1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).q));
    if qqq > 0
        disp('Non-finite q2')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).T));
    if qqq > 0
        disp('Non-finite T1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).T));
    if qqq > 0
        disp('Non-finite T2')
    end
    qqq = sum(1 - isfinite(state_new.fluid(1).Tw));
    if qqq > 0
        disp('Non-finite Tw1')
    end
    qqq = sum(1 - isfinite(state_new.fluid(2).Tw));
    if qqq > 0
        disp('Non-finite Tw2')
    end

end
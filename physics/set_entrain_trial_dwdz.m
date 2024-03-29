% Dynamical entrainment/detrainment based on vertical velocity convergence

relabel.M12_dwdz = 0;
relabel.M21_dwdz = 0;

dM12dm1_dwdz = 0;
dM12dm2_dwdz = 0;
dM21dm1_dwdz = 0;
dM21dm2_dwdz = 0;

relabel.uhat12_dwdz = 0;
relabel.uhat21_dwdz = 0;
relabel.vhat12_dwdz = 0;
relabel.vhat21_dwdz = 0;
relabel.what12_dwdz = 0;
relabel.what21_dwdz = 0;
relabel.qhat12_dwdz = 0;
relabel.qhat21_dwdz = 0;
relabel.etahat12_dwdz = 0;
relabel.etahat21_dwdz = 0;

relabel.duhat12du1_dwdz = 0;
relabel.duhat12du2_dwdz = 0;
relabel.duhat21du1_dwdz = 0;
relabel.duhat21du2_dwdz = 0;
relabel.dvhat12dv1_dwdz = 0;
relabel.dvhat12dv2_dwdz = 0;
relabel.dvhat21dv1_dwdz = 0;
relabel.dvhat21dv2_dwdz = 0;
relabel.dwhat12dw1_dwdz = 0;
relabel.dwhat12dw2_dwdz = 0;
relabel.dwhat21dw1_dwdz = 0;
relabel.dwhat21dw2_dwdz = 0;
relabel.dqhat12dq1_dwdz = 0;
relabel.dqhat12dq2_dwdz = 0;
relabel.dqhat21dq1_dwdz = 0;
relabel.dqhat21dq2_dwdz = 0;
relabel.detahat12deta1_dwdz = 0;
relabel.detahat12deta2_dwdz = 0;
relabel.detahat21deta1_dwdz = 0;
relabel.detahat21deta2_dwdz = 0;

% Default transfer coefficients
bdetrainw_dwdz = dwdz.bdetrainw;
bentrainw_dwdz = dwdz.bentrainw;
bdetrainq_dwdz = dwdz.bdetrainq;
bentrainq_dwdz = dwdz.bentrainq;
bdetraint_dwdz = dwdz.bdetraint;
bentraint_dwdz = dwdz.bentraint;

if (dwdz.entrain | dwdz.detrain) & ischeme == 4
    % Vertical derivative of vertical velocity
    dw1dz = (w1(2:nzp) - w1(1:nz))./grid.dzp;
    dw2dz = (w2(2:nzp) - w2(1:nz))./grid.dzp;
    
    % Vertical derivative of vertical velocity variance
    ww2_below = circshift(ww2, 1);
    ww2(1) = 0;
    ww2_below(1) = 0;
    dww2dz = (ww2 - ww2_below)./grid.dzp;
    gauss_factor = exp(0.5*((w2(2:nzp) - w2(1:nz)).^2)./(ww2+ww2_below)) ./ sqrt(2*3.14159*(ww2+ww2_below));
    
    % Maximum area fraction for downdraft (w < 0) regions from LES
    sigmad = 0.5;
    % Fraction of fluid 1 that is a downdraft
    frac1_down = min(sigma1, sigmad);
    % Fraction of transfer that actually enters fluid 2
    frac1_transfer = 2*min(sigma2, 0.5);
    
    rate_sort12 = dwdz.detrain * min(max(0, -dw2dz * dwdz.detrain_factor), rdt);
    % rate_sort12 = dwdz.detrain * min(max(0, -(dw2dz - dww2dz.*gauss_factor) * dwdz.detrain_factor), rdt);
    rate_sort21 = dwdz.entrain * min(max(0, -dw1dz * dwdz.entrain_factor), rdt) .* frac1_down .* frac1_transfer;

    relabel.M12_dwdz = m2.*rate_sort12;
    relabel.M21_dwdz = m1.*rate_sort21;

    dM12dm1_dwdz = 0 * rate_sort12;
    dM12dm2_dwdz =     rate_sort12;

    dM21dm1_dwdz =     rate_sort21;
    dM21dm2_dwdz = 0 * rate_sort21;
    
    % Calculate the transfer (b) coefficients based on the transfer rate and the PDFs
    if constants.param.dwdz.use_pdf
        deltaw   = w2-w1 + 1e-8*(abs(w2-w1) == 0);
        deltaq   = q2-q1 + 1e-8*(abs(q2-q1) == 0);
        deltaeta = eta2-eta1 + 1e-8*(abs(eta2-eta1) == 0);
        
        % Detrainment
        if dwdz.detrain
            w2std   = sqrt(tke2*2/3);
            q2std   = sqrt(state.fluid(2).varq);
            eta2std = sqrt(state.fluid(2).vareta);
            
            % Avoid division by zero
            w2std   = w2std   + 1e-8*(w2std   == 0);
            q2std   = q2std   + 1e-8*(q2std   == 0);
            eta2std = eta2std + 1e-8*(eta2std == 0);
            
            % Standard deviation at w-levels
            w2stdw   = weight_to_w(grid, w2std);
            q2stdw   = weight_to_w(grid, q2std);
            eta2stdw = weight_to_w(grid, eta2std);
            
            % Quantity which described how much of a pdf is transferred based on the transfer rate
            % See McIntyre 2020 Thesis, chapter 4.2
            pdf_factor12 = sqrt(2)*erfinv(2*dt*rate_sort12 - 1);
            
            % Part of PDF transferred is -infinity to cutoff below
            wcutoff12   = w2   + weight_to_w(grid, w2std   .* pdf_factor12);
            qcutoff12   = q2   + weight_to_w(grid, q2std   .* pdf_factor12);
            etacutoff12 = eta2 + weight_to_w(grid, eta2std .* pdf_factor12);
            
            % Mean of part of PDF transferred
            what12_dwdz = w2 - w2stdw .* exp(-0.5 * ((w2-wcutoff12)./w2stdw).^2) / sqrt(2*pi);
            qhat12_dwdz = q2 - q2stdw .* exp(-0.5 * ((q2-qcutoff12)./q2stdw).^2) / sqrt(2*pi);
            etahat12_dwdz = eta2 - eta2stdw .* exp(-0.5 * ((eta2-etacutoff12)./eta2stdw).^2) / sqrt(2*pi);
            
            bdetrainw_dwdz = (what12_dwdz - w1)./deltaw;
            bdetrainq_dwdz = (qhat12_dwdz - q1)./deltaq;
            bdetraint_dwdz = (etahat12_dwdz - eta1)./deltaeta;
            
            % Precaution to prevent too-extreme values
            bdetrainw_dwdz = min(2, max(-1, bdetrainw_dwdz));
            bdetrainq_dwdz = min(2, max(-1, bdetrainq_dwdz));
            bdetraint_dwdz = min(2, max(-1, bdetraint_dwdz));
            
            % Remove noise from the b-coefficients in places where no transfer is occuring.
            % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
            % If this is not done, the noise/spikes can be interpolated onto regions where transfer
            % is occuring which can cause the simulation to crash.
            filter = weight_to_w(grid, relabel.M12_dwdz) > 1e-5;
            % filter = weight_to_w(grid, relabel.M12_dwdz/(relabel.M12_dwdz + 1e-8*(relabel.M12_dwdz == 0)));
            bdetrainw_dwdz = bdetrainw_dwdz .* filter + (1-filter);
            bdetrainq_dwdz = bdetrainq_dwdz .* filter + (1-filter);
            bdetraint_dwdz = bdetraint_dwdz .* filter + (1-filter);
        end
        
        % Entrainment
        if dwdz.entrain
            w1std   = sqrt(tke1*2/3);
            q1std   = sqrt(state.fluid(1).varq);
            eta1std = sqrt(state.fluid(1).vareta);
            
            % Avoid division by zero
            w1std   = w1std   + 1e-8*(w1std   == 0);
            q1std   = q1std   + 1e-8*(q1std   == 0);
            eta1std = eta1std + 1e-8*(eta1std == 0);
            
            % Standard deviation at w-levels
            w1stdw   = weight_to_w(grid, w1std);
            q1stdw   = weight_to_w(grid, q1std);
            eta1stdw = weight_to_w(grid, eta1std);
            
            % Quantity which described how much of a pdf is transferred based on the transfer rate
            % See McIntyre 2020 Thesis, chapter 4.2
            pdf_factor21 = sqrt(2)*erfinv(2*dt*rate_sort21 - 1);
            
            % Part of PDF transferred is -infinity to cutoff below
            wcutoff21   = w1   - weight_to_w(grid, w1std   .* pdf_factor21);
            qcutoff21   = q1   - weight_to_w(grid, q1std   .* pdf_factor21);
            etacutoff21 = eta1 - weight_to_w(grid, eta1std .* pdf_factor21);
            
            % Mean of part of PDF transferred
            what21_dwdz = w1 + w1stdw .* exp(-0.5 * ((w1-wcutoff21)./w1stdw).^2) / sqrt(2*pi);
            qhat21_dwdz = q1 + q1stdw .* exp(-0.5 * ((q1-qcutoff21)./q1stdw).^2) / sqrt(2*pi);
            etahat21_dwdz = eta1 + eta1stdw .* exp(-0.5 * ((eta1-etacutoff21)./eta1stdw).^2) / sqrt(2*pi);
            
            bentrainw_dwdz = (what21_dwdz - w2)./-deltaw;
            bentrainq_dwdz = (qhat21_dwdz - q2)./-deltaq;
            bentraint_dwdz = (etahat21_dwdz - eta2)./-deltaeta;
            
            % Precaution to prevent too-extreme values
            bentrainw_dwdz = min(2, max(-1, bentrainw_dwdz));
            bentrainq_dwdz = min(2, max(-1, bentrainq_dwdz));
            bentraint_dwdz = min(2, max(-1, bentraint_dwdz));
            
            % Remove noise from the b-coefficients in places where no transfer is occuring.
            % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
            % If this is not done, the noise/spikes can be interpolated onto regions where transfer
            % is occuring which can cause the simulation to crash.
            filter = weight_to_w(grid, relabel.M21_dwdz) > 1e-5;
            % filter = weight_to_w(grid, relabel.M12_dwdz/(relabel.M12_dwdz + 1e-8*(relabel.M12_dwdz == 0)));
            bentrainw_dwdz = bentrainw_dwdz .* filter + (1-filter);
            bentrainq_dwdz = bentrainq_dwdz .* filter + (1-filter);
            bentraint_dwdz = bentraint_dwdz .* filter + (1-filter);
        end
    end

    % Entrained and detrained values of eta
    [relabel.etahat12_dwdz,relabel.detahat12deta1_dwdz,relabel.detahat12deta2_dwdz,...
     relabel.etahat21_dwdz,relabel.detahat21deta1_dwdz,relabel.detahat21deta2_dwdz] ...
        = findqhat(eta1, eta2, bentraint_dwdz, bdetraint_dwdz);

    % Entrained and detrained values of q
    [relabel.qhat12_dwdz,relabel.dqhat12dq1_dwdz,relabel.dqhat12dq2_dwdz,...
     relabel.qhat21_dwdz,relabel.dqhat21dq1_dwdz,relabel.dqhat21dq2_dwdz] ...
        = findqhat(q1, q2, bentrainq_dwdz, bdetrainq_dwdz);

    % Entrained and detrained values of velocities u, v and w
    [relabel.uhat12_dwdz,relabel.duhat12du1_dwdz,relabel.duhat12du2_dwdz,...
     relabel.uhat21_dwdz,relabel.duhat21du1_dwdz,relabel.duhat21du2_dwdz] ...
        = findqhat(u1, u2, dwdz.bentrainu, dwdz.bdetrainu);

    [relabel.vhat12_dwdz,relabel.dvhat12dv1_dwdz,relabel.dvhat12dv2_dwdz,...
     relabel.vhat21_dwdz,relabel.dvhat21dv1_dwdz,relabel.dvhat21dv2_dwdz] ...
        = findqhat(v1, v2, dwdz.bentrainu, dwdz.bdetrainu);

    [relabel.what12_dwdz,relabel.dwhat12dw1_dwdz,relabel.dwhat12dw2_dwdz,...
     relabel.what21_dwdz,relabel.dwhat21dw1_dwdz,relabel.dwhat21dw2_dwdz] ...
        = findqhat(w1, w2, bentrainw_dwdz, bdetrainw_dwdz);



    % Experimental features

    % Prevent double counting from instab entrainment 
    % relabel.M21_dwdz = max(relabel.M21_dwdz-relabel.M21_instab, 0);
    % dM21dm1_dwdz = max(dM21dm1_dwdz-dM21dm1_instab, 0);

    % Force detrained velocity to be zero
    % relabel.what12_dwdz = 0*relabel.what12_dwdz;
    % relabel.dwhat12dw1_dwdz = 0*relabel.dwhat12dw1_dwdz;
    % relabel.dwhat12dw2_dwdz = 0*relabel.dwhat12dw2_dwdz;

    % Force entrained velocity to be zero
    % relabel.what21_dwdz = 0*relabel.what21_dwdz;
    % relabel.dwhat21dw1_dwdz = 0*relabel.dwhat21dw1_dwdz;
    % relabel.dwhat21dw2_dwdz = 0*relabel.dwhat21dw2_dwdz;

    % Force transferred velocities to have closer values to their new fluid
    % relabel.what12_dwdz = min(relabel.what12_dwdz, 0);
    % relabel.what21_dwdz = max(relabel.what21_dwdz, 0);
end
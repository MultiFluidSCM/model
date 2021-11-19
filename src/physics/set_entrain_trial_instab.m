% Instability source entrainment where stratification is unstable
% *** Put tuneable coefficients in constants.params ***
relabel.M12_instab = instab.detrain * instab.detrain_factor * m2 .* n2pos;
relabel.M21_instab = instab.entrain * instab.entrain_factor * m1 .* n1neg;

dM12dm1_instab = instab.detrain * 0.  * n2pos;
dM12dm2_instab = instab.detrain * instab.detrain_factor * n2pos;
dM21dm1_instab = instab.entrain * instab.entrain_factor * n1neg;
dM21dm2_instab = instab.entrain * 0.  * n1neg;

% Default transfer coefficients
bdetrainw_instab = instab.bdetrainw;
bentrainw_instab = instab.bentrainw;
bdetrainq_instab = instab.bdetrainq;
bentrainq_instab = instab.bentrainq;
bdetraint_instab = instab.bdetraint;
bentraint_instab = instab.bentraint;

% Calculate the transfer (b) coefficients based on the transfer rate and the PDFs
if constants.param.instab.use_pdf
    deltaw   = w2-w1 + 1e-8*(abs(w2-w1) == 0);
    deltaq   = q2-q1 + 1e-8*(abs(q2-q1) == 0);
    deltaeta = eta2-eta1 + 1e-8*(abs(eta2-eta1) == 0);
    
    % Detrainment
    if instab.detrain
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
        pdf_factor12 = sqrt(2)*erfinv(2*dt*dM12dm2_instab - 1);
        
        % Part of PDF transferred is -infinity to cutoff below
        wcutoff12   = w2   + weight_to_w(grid, w2std   .* pdf_factor12);
        qcutoff12   = q2   + weight_to_w(grid, q2std   .* pdf_factor12);
        etacutoff12 = eta2 + weight_to_w(grid, eta2std .* pdf_factor12);
        
        % Mean of part of PDF transferred
        what12_instab = w2 - w2stdw .* exp(-0.5 * ((w2-wcutoff12)./w2stdw).^2) / sqrt(2*pi);
        qhat12_instab = q2 - q2stdw .* exp(-0.5 * ((q2-qcutoff12)./q2stdw).^2) / sqrt(2*pi);
        etahat12_instab = eta2 - eta2stdw .* exp(-0.5 * ((eta2-etacutoff12)./eta2stdw).^2) / sqrt(2*pi);
        
        bdetrainw_instab = (what12_instab - w1)./deltaw;
        bdetrainq_instab = (qhat12_instab - q1)./deltaq;
        bdetraint_instab = (etahat12_instab - eta1)./deltaeta;
        
        % Precaution to prevent too-extreme values
        bdetrainw_instab = min(2, max(-1, bdetrainw_instab));
        bdetrainq_instab = min(2, max(-1, bdetrainq_instab));
        bdetraint_instab = min(2, max(-1, bdetraint_instab));
        
        % Remove noise from the b-coefficients in places where no transfer is occuring.
        % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
        % If this is not done, the noise/spikes can be interpolated onto regions where transfer
        % is occuring which can cause the simulation to crash.
        filter = weight_to_w(grid, relabel.M12_instab) > 1e-4;
        % filter = weight_to_w(grid, relabel.M12_instab/(relabel.M12_instab + 1e-8*(relabel.M12_instab == 0)));
        bdetrainw_instab = bdetrainw_instab .* filter + (1-filter);
        bdetrainq_instab = bdetrainq_instab .* filter + (1-filter);
        bdetraint_instab = bdetraint_instab .* filter + (1-filter);
    end
    
    % Entrainment
    if instab.entrain
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
        pdf_factor21 = sqrt(2)*erfinv(2*dt*dM21dm1_instab - 1);
        
        % Part of PDF transferred is -infinity to cutoff below
        wcutoff21   = w1   - weight_to_w(grid, w1std   .* pdf_factor21);
        qcutoff21   = q1   - weight_to_w(grid, q1std   .* pdf_factor21);
        etacutoff21 = eta1 - weight_to_w(grid, eta1std .* pdf_factor21);
        
        % Mean of part of PDF transferred
        what21_instab = w1 + w1stdw .* exp(-0.5 * ((w1-wcutoff21)./w1stdw).^2) / sqrt(2*pi);
        qhat21_instab = q1 + q1stdw .* exp(-0.5 * ((q1-qcutoff21)./q1stdw).^2) / sqrt(2*pi);
        etahat21_instab = eta1 + eta1stdw .* exp(-0.5 * ((eta1-etacutoff21)./eta1stdw).^2) / sqrt(2*pi);
        
        bentrainw_instab = (what21_instab - w2)./-deltaw;
        bentrainq_instab = (qhat21_instab - q2)./-deltaq;
        bentraint_instab = (etahat21_instab - eta2)./-deltaeta;
        
        % Precaution to prevent too-extreme values
        bentrainw_instab = min(2, max(-1, bentrainw_instab));
        bentrainq_instab = min(2, max(-1, bentrainq_instab));
        bentraint_instab = min(2, max(-1, bentraint_instab));
        
        % Remove noise from the b-coefficients in places where no transfer is occuring.
        % If the transfer from dw/dz detrainment is very small, set the coefficient to 1.
        % If this is not done, the noise/spikes can be interpolated onto regions where transfer
        % is occuring which can cause the simulation to crash.
        filter = weight_to_w(grid, relabel.M21_instab) > 1e-4;
        % filter = weight_to_w(grid, relabel.M12_instab/(relabel.M12_instab + 1e-8*(relabel.M12_instab == 0)));
        bentrainw_instab = bentrainw_instab .* filter + (1-filter);
        bentrainq_instab = bentrainq_instab .* filter + (1-filter);
        bentraint_instab = bentraint_instab .* filter + (1-filter);
    end
end

% Entrained and detrained values of eta
[relabel.etahat12_instab,relabel.detahat12deta1_instab,relabel.detahat12deta2_instab,...
 relabel.etahat21_instab,relabel.detahat21deta1_instab,relabel.detahat21deta2_instab] ...
    = findqhat(eta1, eta2, bentraint_instab, bdetraint_instab);

% Entrained and detrained values of q
[relabel.qhat12_instab,relabel.dqhat12dq1_instab,relabel.dqhat12dq2_instab,...
 relabel.qhat21_instab,relabel.dqhat21dq1_instab,relabel.dqhat21dq2_instab] ...
    = findqhat(q1, q2, bentrainq_instab, bdetrainq_instab);

% Entrained and detrained values of velocities u, v and w
[relabel.uhat12_instab,relabel.duhat12du1_instab,relabel.duhat12du2_instab,...
 relabel.uhat21_instab,relabel.duhat21du1_instab,relabel.duhat21du2_instab] ...
    = findqhat(u1, u2, instab.bentrainu, instab.bdetrainu);

[relabel.vhat12_instab,relabel.dvhat12dv1_instab,relabel.dvhat12dv2_instab,...
 relabel.vhat21_instab,relabel.dvhat21dv1_instab,relabel.dvhat21dv2_instab] ...
    = findqhat(v1, v2, instab.bentrainu, instab.bdetrainu);

[relabel.what12_instab,relabel.dwhat12dw1_instab,relabel.dwhat12dw2_instab,...
 relabel.what21_instab,relabel.dwhat21dw1_instab,relabel.dwhat21dw2_instab] ...
    = findqhat(w1, w2, bentrainw_instab, bdetrainw_instab);

% relabel.what12_instab = min(relabel.what21_instab, 0);
% relabel.what21_instab = max(relabel.what21_instab, 0);
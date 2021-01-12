function [chi_hat,deriv] = compute_chi_hat(x)
% Compute chi_hat and also the derivative of w_hat wrt sigmaw
% 

r2 = sqrt(2);
rr2 = 1/r2;
r2bypi = sqrt(2/pi);

if x < -5
    
    % Use asymptotic formula to avoid underflow errors
    rden = 1/(1 - 0.5/x^2);
    chi_hat = r2*x*rden;
    deriv = r2*rden*rden/x;
    
else
    
    % Use exact formula
    E = exp(-x^2);
    R = (1 + erf(x));
    chi_hat = -r2bypi*E/R;
    deriv = (-r2bypi*R*E - 2*r2bypi*(x^2)*R*E - (2*r2/pi)*x*E*E)/(R*R);

end


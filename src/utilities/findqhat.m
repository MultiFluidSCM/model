function [ qhat12, dqhat12dq1, dqhat12dq2, qhat21, dqhat21dq1, dqhat21dq2 ] ...
           = findqhat(q1,q2,bentrain,bdetrain)

% Determine entrained and detrained values of q
% Also compute terms needed for linearization

% Use a simple linear interpolation between plume and
% environment with prescribed coefficients

% Detrained values
qhat12 = (1 - bdetrain).*q1 + bdetrain.*q2;
dqhat12dq1 = ones(size(q1)).*(1 - bdetrain);
dqhat12dq2 = ones(size(q1)).*bdetrain;

% Entrained values
qhat21 = (1 - bentrain).*q2 + bentrain.*q1;
dqhat21dq2 = ones(size(q2)).*(1 - bentrain);
dqhat21dq1 = ones(size(q2)).*bentrain;

end


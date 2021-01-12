% Compare relabelling terms computed with old and new routines

% What's in the structures?
relabel
relabelx
pause

% E and D rates
disp('diff M21')
relabel.M21 - relabelx.M21
disp('diff M21')
relabel.M12 - relabelx.M12
pause

% Derivatives
disp('disp dM21dm1')
relabel.dM21dm1 - relabelx.dM21dm1
disp('disp dM12dm1')
relabel.dM12dm1 - relabelx.dM12dm1
pause
disp('disp dM21dm2')
relabel.dM21dm2 - relabelx.dM21dm2
disp('disp dM12dm2')
relabel.dM12dm2 - relabelx.dM12dm2
pause
disp('disp dM21dw1')
relabel.dM21dw1 - relabelx.dM21dw1
disp('disp dM12dw1')
relabel.dM12dw1 - relabelx.dM12dw1
pause
disp('disp dM21dw2')
relabel.dM21dw2 - relabelx.dM21dw2
disp('disp dM12dw2')
relabel.dM12dw2 - relabelx.dM12dw2
pause
disp('disp dM21deta1')
relabel.dM21deta1 - relabelx.dM21deta1
disp('disp dM12deta1')
relabel.dM12deta1 - relabelx.dM12deta1
pause
disp('disp dM21deta2')
relabel.dM21deta2 - relabelx.dM21deta2
disp('disp dM12deta2')
relabel.dM12deta2 - relabelx.dM12deta2
pause
disp('disp dM21dq1')
relabel.dM21dq1 - relabelx.dM21dq1
disp('disp dM12dq1')
relabel.dM12dq1 - relabelx.dM12dq1
pause
disp('disp dM21dq2')
relabel.dM21dq2 - relabelx.dM21dq2
disp('disp dM12dq2')
relabel.dM12dq2 - relabelx.dM12dq2
pause
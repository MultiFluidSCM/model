function [ drag, ddragdw1, ddragdw2 ] = find_drag( m2bar, w1, w2, zdrag )
%FIND_DRAG Determine drag on updrafts
% and terms needed for linearization

% Thermals are assumed to be of scale zdrag

dw = w2 - w1;
% disp('** No Drag **')
% dw = dw*0;
adw = abs(dw);
drag = m2bar.*dw.*adw/zdrag;
ddragdw1 = - 2*m2bar.*adw/zdrag;
ddragdw2 =   2*m2bar.*adw/zdrag;

end


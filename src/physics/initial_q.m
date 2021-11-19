function qv = initial_q( z, initial )
%INITIAL_Q Set initial water specific humidity
%   Provides a water vapour as a function of height z

for i=1:length(initial.z)-1
   if (initial.z(i) <= z && z < initial.z(i+1))
        qv = initial.qv(i)+(initial.qv(i+1)-initial.qv(i)) .* (z-initial.z(i))./(initial.z(i+1)-initial.z(i)) ;
   end
end
if z >= initial.z(end)
    qv = initial.qv(end);
end

end
function qv = initial_q( z, initial )
%INITIAL_Q Set initial water vapour mixing ratio
%   Provides a water vapour as a function of height z

initial.rv = initial.rv./(1 + initial.rv);

for i=1:length(initial.z)-1
   if (initial.z(i) <= z && z < initial.z(i+1))
        qv = initial.rv(i)+(initial.rv(i+1)-initial.rv(i)) .* (z-initial.z(i))./(initial.z(i+1)-initial.z(i)) ;
   end
end
if z >= initial.z(end)
    qv = initial.rv(end);
end

end
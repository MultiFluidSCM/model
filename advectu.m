function [ Au, dAudua, dAuduc, dAudub ] = find_advectu( dzw, w, u )
%   Calculate an upwind finite difference approximation to
%   the vertical advection of u or v
%   w du/dz , along with terms needed in linearization.
%   w is assumed to have been interpolated to u levels

nzp = numel(dzw);
nz  = nzp - 1;

% Now compute upwind gradient and derivatives
% wrt u at levels above, centred, and below
k = 1;
if w(k) > 0
    Au(k) = w(k)*(u(k  )         )/dzw(k  );
    dAudua(k) = 0;
    dAuduc(k) =   w(k)/dzw(k  );
    dAudub(k) = 0;
else
    Au(k) = w(k)*(u(k+1) - u(k  ))/dzw(k+1);
    dAudua(k) =   w(k)/dzw(k+1);
    dAuduc(k) = - w(k)/dzw(k+1);
    dAudub(k) = 0;
end
for k = 2:nz-1
    if w(k) > 0
        Au(k) = w(k)*(u(k  ) - u(k-1))/dzw(k  );
        dAudua(k) = 0;
        dAuduc(k) =   w(k)/dzw(k  );
        dAudub(k) = - w(k)/dzw(k  );
    else
        Au(k) = w(k)*(u(k+1) - u(k  ))/dzw(k+1);
        dAudua(k) =   w(k)/dzw(k+1);
        dAuduc(k) = - w(k)/dzw(k+1);
        dAudub(k) = 0;
    end
end
k = nz;
if w(k) > 0
    Au(k) = w(k)*(u(k  ) - u(k-1))/dzw(k  );
    dAudua(k) = 0;
    dAuduc(k) =   w(k)/dzw(k  );
    dAudub(k) = - w(k)/dzw(k  );
else
    Au(k) = 0;
    dAudua(k) = 0;
    dAuduc(k) = 0;
    dAudub(k) = 0;
end
    

end


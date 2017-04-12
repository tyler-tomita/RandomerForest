function theta = vector_angle(u,v,varargin)
    % VECTOR_ANGLE computes the angle between two vectors in the plane
    % spanned by the vectors

    if nargin < 3
        unit = 'radians';
    else
        unit = varargin(1);
    end
    
    if ~strcmp(unit,'degrees') && ~strcmp(unit,'radians')
        error('unit must be either ''degrees'' or ''radians''')
    end
    
    theta = real(acos(dot(u,v)/norm(u)/norm(v)));
%     theta = atan2d(norm(cross(u,v)),dot(u,v));
    
    if strcmp(unit,'degrees')
        theta = theta*180/pi;
    end

end
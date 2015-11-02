%generates uniformly random dxd rotation matrix
function R = random_rotation(d)
    [R,~] = qr(randn(d));
    if det(R) < 0
        if d < 3
            R = R(:,[2,1]);
        else
            R = R(:,[2,1,3:end]);
        end
    end
end
function Xr = rotate2d(X,theta)

    Xr = X*[cos(theta) sin(theta);-sin(theta) cos(theta)];
    
end
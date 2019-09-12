function [psi] = velocity2stream(u1, u2)
%VELOCITY2STREAM calculate the streamfunction field from the u, v velocity
% fields, given that the flow is incompressible
%   Math stolen from
%   https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2001/msg00303.html
%
%   u = d(Psi)/dy and v = -d(Psi)/dx.
%   u=d(Psi)/dy -> Psi = Int{y0:y}u dy + a(x)
%   d(Psi)/dx = Int{y0:y}du/dx dy + da/dx
%
%   Psi = Int{y0:y}u(x,y) dy - Int{x_0:x}v(x,y_0) dx

    
    X_length = size(u1, 1);
    Y_length = size(u1, 2);
    
    term_1 = cumsum(u1(1, :));
    term_2 = zeros(size(u1));
    for y = 1:Y_length
        term_2(:, y) = cumsum(u2(:, y));
    end
    
    psi = repmat(term_1, [X_length, 1]) + term_2;
    
end


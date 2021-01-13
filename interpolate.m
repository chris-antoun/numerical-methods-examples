% This script implements the Lagrange form of polynomial interpolation
% Define a function called interpolate
function y = interpolate(f, x_i, x)
    % f is the function we are trying to approximate
    % x_i is a vector containing the interpolating points
    % x is a vector of points that we want to interpolate over
    y = zeros(size(x)); % This will be used to write the outcome 
    for j = 1:length(x)
        for k = 1:length(x_i)
            L_k = 1; % This is the case where j = k
            for l = 1:length(x_i)
                if l == k
                    continue;
                end
                L_k = L_k * (x(j)-x_i(l))/(x_i(k) - x_i(l));
            end
            y(j) = y(j) + f(x_i(k))*L_k;
        end
    end
                    
                
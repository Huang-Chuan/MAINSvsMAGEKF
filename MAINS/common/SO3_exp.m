function R = SO3_exp(v)
% Copyright (C) 2024 by Chuan Huang
% SO3_EXP  Exponential map for SO(3)
%   R = SO3_EXP(v) computes the exponential map for SO(3) given a tangent
%   vector v in R^3. The output R is a 3x3 rotation matrix.
    theta = norm(v);
    if theta < 1e-8
        R = eye(3);
    else
        k = v / theta;
        K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
    end
end

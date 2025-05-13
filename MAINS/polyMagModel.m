% Copyright 2023, Chuan Huang
    
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at

%   http://www.apache.org/licenses/LICENSE-2.0

% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
classdef polyMagModel
    %   The class for handling polynomial model fitting and updates 

    properties
        A
        pinvA
        B
        J1J2
    end
    
    properties (Access = private)
        points
        theta
        dim
        order
        dpsi
        db_dpsi
        phi
        pinvphi
    end


    methods (Access = public)
        
        function obj = polyMagModel(order)
            obj.points = polyMagModel.get_points(order);
            obj.dim = order^2 + 4*order+3;
            obj.A    = [];
            % points
            for i = 1 : size(obj.points, 1)
                obj.A = [obj.A; polyMagModel.get_ABnull(obj.points(i, :), order)];
            end
            obj.pinvA = pinv(obj.A);
            obj.order = order;

            % choose db_dpsi function based on the order
            if order == 1
                obj.db_dpsi = @dbtheta_dpsi1;
            elseif order == 2
                obj.db_dpsi = @dbtheta_dpsi2;
            elseif order == 3
                obj.db_dpsi = @dbtheta_dpsi3;
            end
            obj.theta = zeros(obj.dim, 1);
        end

        function obj = set_theta(obj, theta)
            obj.theta = theta;
        end

        function obj = init_theta(obj, mag)
           [obj.theta, ~] =  obj.LS_coeff(mag);
        end

        function theta = get_theta(obj)
            theta = obj.theta;
        end 


        function obj = update_theta(obj, dpsi)
            % keep record of dpsi
            obj.dpsi = dpsi;
            % calculate B
            obj.B = obj.getB(dpsi);
            % calculate d(B*theta)/dpsi
            obj.J1J2 = zeros(3 * size(obj.points, 1), 6);
            for i = 1 : size(obj.points, 1)
                obj.J1J2(3 * (i - 1) + 1 : 3 * i, :) = obj.db_dpsi(dpsi, obj.points(i, :)', obj.theta);
            end
            % update theta
            obj.theta = obj.pinvA * obj.B * obj.theta;
        end


        function obj = set_phi(obj, sensor_pos)
            
            obj.phi = [];
            for i = 1 : size(sensor_pos, 2)
                obj.phi =  [obj.phi; polyMagModel.get_ABnull(sensor_pos(:, i), obj.order)];
            end
            obj.pinvphi = pinv(obj.phi);
        end



        function [coeff, var_res] = LS_coeff(obj, mag)
        % calculate LS estimate of theta based on magnetic field data
        % Input
        %       mag: magnetic field data   --- (dim of the magnetic field data) x 1 vector 
        % Output
        %       coeff:  estimated theta    --- (dim of theta) x 1 vector
        %       var_res: fitting residuals --- (dim of the magnetic field data) x 1 vector

            coeff = obj.pinvphi*mag;
            var_res = mean((mag - obj.phi * coeff).^2);
        end 



    end

    methods (Access = private)
        function B = getB(obj, dpsi)
            % dpsi: 3 x 1
            dp   = dpsi(1:3);
            dphi = dpsi(4:6);

            %R=eye(3)+vect2skew(dphi);        % small angle approximation
            R = SO3_exp(dphi);
            B=zeros(3*size(obj.points,1), obj.dim);
            
            for ii=1:size(obj.points, 1)
                r=R*obj.points(ii,:)'+dp;
                B(3*(ii-1)+1:3*ii,:) = R.'*polyMagModel.get_ABnull(r,obj.order);
            end
        end

    end



    methods (Static)
        
        function points = get_points(order)
        % Pre-defined point sets for determine theta's dynamics
        % Input
        %       order: the order of the polynomial model --- scalar
        % Output
        %       points: the selected points where the magnetic field is expressed --- ceil(kappa/3) x 3 vector (row)
            switch order
                case 1
                    points = [0 0 1;...
                             0 1 0;...
                             1 0 0;...
                            ];            
                case 2
                    points = [0 0 1;...
                             0 1 0;...
                             1 0 0;...
                             1 1 1;...
                             0 0 0];
                case 3
                    points =  [0 0 1;...
                             0 1 0;...
                             1 0 0;...
                             1 1 1;...
                             0 0 0;...
                             1 1 0;...
                             1 0 1;...
                            -1 0 1];
                otherwise
                    disp('Not supported!')
            end
        end



        function ABnull=get_ABnull(r, order)
        % calculate the basis for representing magnetic field
        %     r: location                          --- 1 x 3 or 3 x 1 vector
        % order: the order of the polynomial model --- scalar
            rx=r(1);
            ry=r(2);
            rz=r(3);
            
            switch order
            
                case 1
                    ABnull=[ 0, 0, 1,  0,     0, rz, ry,  2*rx;...
                        0, 1, 0, rz,  2*ry,  0, rx,     0;...
                        1, 0, 0, ry, -2*rz, rx,  0, -2*rz];
            
                case 2
                    ABnull=[ 0, 0, 1,  0,     0, rz, ry,  2*rx,           0,               0, ry*rz, ry^2 - rz^2,     2*rx*rz,     2*rx*ry, 3*rx^2 - 3*rz^2;...
                        0, 1, 0, rz,  2*ry,  0, rx,     0,     2*ry*rz, 3*ry^2 - 3*rz^2, rx*rz,     2*rx*ry,           0, rx^2 - rz^2,               0;...
                        1, 0, 0, ry, -2*rz, rx,  0, -2*rz, ry^2 - rz^2,        -6*ry*rz, rx*ry,    -2*rx*rz, rx^2 - rz^2,    -2*ry*rz,        -6*rx*rz];
            
                case 3
                    ABnull=[ 0, 0, 1,  0,     0, rz, ry,  2*rx,           0,               0, ry*rz, ry^2 - rz^2,     2*rx*rz,     2*rx*ry, 3*rx^2 - 3*rz^2,                0,                   0,  ry^2*rz - rz^3/3,      ry^3 - 3*ry*rz^2,        2*rx*ry*rz,                2*rx*ry^2 - 2*rx*rz^2, 3*rx^2*rz - rz^3, 3*ry*rx^2 - 3*ry*rz^2, 4*rx^3 - 12*rx*rz^2;...
                        0, 1, 0, rz,  2*ry,  0, rx,     0,     2*ry*rz, 3*ry^2 - 3*rz^2, rx*rz,     2*rx*ry,           0, rx^2 - rz^2,               0, 3*ry^2*rz - rz^3, 4*ry^3 - 12*ry*rz^2,        2*rx*ry*rz, 3*rx*ry^2 - 3*rx*rz^2,  rx^2*rz - rz^3/3,                2*ry*rx^2 - 2*ry*rz^2,                0,      rx^3 - 3*rx*rz^2,                   0;...
                        1, 0, 0, ry, -2*rz, rx,  0, -2*rz, ry^2 - rz^2,        -6*ry*rz, rx*ry,    -2*rx*rz, rx^2 - rz^2,    -2*ry*rz,        -6*rx*rz, ry^3 - 3*ry*rz^2, 4*rz^3 - 12*ry^2*rz, rx*ry^2 - rx*rz^2,           -6*rx*ry*rz, ry*rx^2 - ry*rz^2, - 2*rx^2*rz - 2*ry^2*rz + (4*rz^3)/3, rx^3 - 3*rx*rz^2,           -6*rx*ry*rz, 4*rz^3 - 12*rx^2*rz];
                case 4  
                    ABnull=[ 0, 0, 1,  0,     0, rz, ry,  2*rx,           0,               0, ry*rz, ry^2 - rz^2,     2*rx*rz,     2*rx*ry, 3*rx^2 - 3*rz^2,                0,                   0,  ry^2*rz - rz^3/3,      ry^3 - 3*ry*rz^2,        2*rx*ry*rz,                2*rx*ry^2 - 2*rx*rz^2, 3*rx^2*rz - rz^3, 3*ry*rx^2 - 3*ry*rz^2, 4*rx^3 - 12*rx*rz^2,                         0,                              0,      ry^3*rz - ry*rz^3,   ry^4 - 6*ry^2*rz^2 + rz^4,               2*rx*rz*ry^2 - (2*rx*rz^3)/3,                       2*rx*ry^3 - 6*rx*ry*rz^2, 3*ry*rx^2*rz - ry*rz^3, 3*rx^2*ry^2 - 3*rx^2*rz^2 - 3*ry^2*rz^2 + rz^4,     4*rx^3*rz - 4*rx*rz^3,   4*ry*rx^3 - 12*ry*rx*rz^2, 5*rx^4 - 30*rx^2*rz^2 + 5*rz^4;...
                           0, 1, 0, rz,  2*ry,  0, rx,     0,     2*ry*rz, 3*ry^2 - 3*rz^2, rx*rz,     2*rx*ry,           0, rx^2 - rz^2,               0, 3*ry^2*rz - rz^3, 4*ry^3 - 12*ry*rz^2,        2*rx*ry*rz, 3*rx*ry^2 - 3*rx*rz^2,  rx^2*rz - rz^3/3,                2*ry*rx^2 - 2*ry*rz^2,                0,      rx^3 - 3*rx*rz^2,                   0,     4*ry^3*rz - 4*ry*rz^3, 5*ry^4 - 30*ry^2*rz^2 + 5*rz^4, 3*rx*ry^2*rz - rx*rz^3,   4*rx*ry^3 - 12*rx*ry*rz^2,               2*ry*rz*rx^2 - (2*ry*rz^3)/3, 3*rx^2*ry^2 - 3*rx^2*rz^2 - 3*ry^2*rz^2 + rz^4,      rx^3*rz - rx*rz^3,                       2*ry*rx^3 - 6*ry*rx*rz^2,                         0,   rx^4 - 6*rx^2*rz^2 + rz^4,                              0;...
                         1, 0, 0, ry, -2*rz, rx,  0, -2*rz, ry^2 - rz^2,        -6*ry*rz, rx*ry,    -2*rx*rz, rx^2 - rz^2,    -2*ry*rz,        -6*rx*rz, ry^3 - 3*ry*rz^2, 4*rz^3 - 12*ry^2*rz, rx*ry^2 - rx*rz^2,           -6*rx*ry*rz, ry*rx^2 - ry*rz^2, - 2*rx^2*rz - 2*ry^2*rz + (4*rz^3)/3, rx^3 - 3*rx*rz^2,           -6*rx*ry*rz, 4*rz^3 - 12*rx^2*rz, ry^4 - 6*ry^2*rz^2 + rz^4,      - 20*ry^3*rz + 20*ry*rz^3, rx*ry^3 - 3*rx*ry*rz^2, - 12*rx*ry^2*rz + 4*rx*rz^3, rx^2*ry^2 - rx^2*rz^2 - ry^2*rz^2 + rz^4/3,         - 6*rx^2*ry*rz - 2*ry^3*rz + 4*ry*rz^3, ry*rx^3 - 3*ry*rx*rz^2,         - 2*rx^3*rz - 6*rx*ry^2*rz + 4*rx*rz^3, rx^4 - 6*rx^2*rz^2 + rz^4, - 12*ry*rx^2*rz + 4*ry*rz^3,      - 20*rx^3*rz + 20*rx*rz^3];
            
                otherwise
                    error('Only model orders 1-4 supported');
            end
        end
    end

end
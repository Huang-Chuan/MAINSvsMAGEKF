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
function rotm = q2r(q)
    % q2r: converts a quaternion to a rotation matrix
    %   q: quaternion         4 x 1 or 1 x 4
    %   rotm: rotation matrix 3 x 3
    q_w = q(1);
    q_x = q(2);
    q_y = q(3);
    q_z = q(4);


    rotm = [q_w^2+q_x^2-q_y^2-q_z^2,2*(q_x*q_y-q_w*q_z),2*(q_x*q_z+q_w*q_y);
           2*(q_x*q_y+q_w*q_z), q_w^2-q_x^2+q_y^2-q_z^2,2*(q_y*q_z-q_w*q_x);
           2*(q_x*q_z-q_w*q_y), 2*(q_y*q_z+q_w*q_x), q_w^2-q_x^2-q_y^2+q_z^2];


end
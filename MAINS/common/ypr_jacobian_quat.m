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
function J = ypr_jacobian_quat(quatRef)
    % calculate the Jacobian of the ypr angles with respect to the quaternion perturbation
    eulAngRef = quat2eul(quatRef, 'ZYX');
    J = zeros(3);

    dqx = [1, 1/2 * 1e-7, 0, 0];
    dqy = [1, 0, 1/2 * 1e-7, 0];
    dqz = [1, 0,   0,  1/2 * 1e-7];
    
    eulAng1 = quat2eul(quatmultiply(quatRef, dqx), 'ZYX');
    eulAng2 = quat2eul(quatmultiply(quatRef, dqy), 'ZYX');
    eulAng3 = quat2eul(quatmultiply(quatRef, dqz), 'ZYX');
    
    
    J(:, 1) = angdiff(eulAng1, eulAngRef).' / 1e-7;
    J(:, 2) = angdiff(eulAng2, eulAngRef).' / 1e-7;
    J(:, 3) = angdiff(eulAng3, eulAngRef).' / 1e-7;

end
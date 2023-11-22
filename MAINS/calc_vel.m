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
function [vel] = calc_vel(in_data)
% compute velocity in navigation frame

    timeVector = in_data.t(1:end);
            dt = mean(diff(timeVector));
          pos  = in_data.gt.pos(:, 1:end);

           vel = zeros(size(pos));
    for i = 3 : length(pos) - 2
        vel(:, i) = (-pos(:, i + 2) +  8 * pos(:, i + 1) - 8 * pos(:, i - 1) + pos(:, i - 2))/(12*dt);
    end
    vel(:, 1) = vel(:, 2);
    vel(:, end) = vel(:, end - 1);
end
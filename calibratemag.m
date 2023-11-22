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
function [calibratedMag] = calibratemag(MAG, mag_paras, N)
    calibratedMag = zeros(size(MAG));
    
    for i = 1 : N
        MAG_i = MAG(:, (i - 1) * 3 + 1 : 3 * i)';
        D_i = reshape(mag_paras(1:9, i), 3, 3);
        h_i = reshape(mag_paras(10:end, i), 3, 1);
        calibratedMag(:, (i - 1) * 3 + 1 : 3 * i) = (inv(D_i) * (MAG_i - h_i))';
    
    end
end
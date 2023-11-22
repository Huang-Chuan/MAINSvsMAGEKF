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
function r=PosMagArray()
% This function returns the position of the sensors in the magnetic array
% Output
%        r:  position of the sensors    --- 3 x (number of sensors) vector
    r = NaN(3, 30);
%   Comment out the below and replace r with the correponding positions 
%   if you want to use your own configuration
    dx = 0.05;
    dy = 0.05;
    kk = 0;

    for jj=1:5
        for ii=1:6
            kk=kk+1;
            if ii<4
                r(1, kk) = (ii-3.5)*dx-0.5*dx;
            else
                r(1, kk) = (ii-3.5)*dx+0.5*dx;
            end
            r(2, kk) = -(jj-3)*dy;
            r(3, kk) = 0;
        end
    end
end

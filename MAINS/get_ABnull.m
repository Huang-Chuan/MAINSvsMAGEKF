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
function ABnull=get_ABnull(r,order)
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

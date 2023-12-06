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
function y = vect2skew(v)
%This is a function to generate the skew symmectric form [v] of a vector v
%Usage:
%   vect2skew(v)
%       v: a 3x1 or 1x3 vector

    y=[0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];

end
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
function [stat] = compute_stats(in_data, out_data)
    % 
    stat.duration = in_data.t(end) - in_data.t(1);    
    stat.dist = sum(vecnorm(diff(in_data.gt.pos')'));
    
    

    speed_t = vecnorm(calc_vel(in_data));

    t_start = 60;
    idx = in_data.t > t_start;

    stat.duration_test = stat.duration - t_start;
    stat.dist_test = sum(vecnorm(diff(in_data.gt.pos(:, idx)')'));
    stat.height = mean(in_data.gt.pos(3, idx));
    
    posEst = out_data.x_h(1:3, idx);
    posTrue = in_data.gt.pos(:, idx);

    % horizontal error 
    stat.h_err = mean(vecnorm(posEst(1:2, :) - posTrue(1:2, :)));
    fprintf('Horizontal error: %.2f m.\n', stat.h_err);
    stat.h_err_end = mean(vecnorm(posEst(1:2, end) - posTrue(1:2, end)));
    fprintf('Horizontal error (end): %.2f m.\n',stat.h_err_end);
    % vertical error
    stat.v_err = mean(abs(posEst(3, :) - posTrue(3, :)));
    fprintf('Vertical error: %.2f m.\n', stat.v_err);
    stat.v_err_end = mean(vecnorm(posEst(3, end) - posTrue(3, end)));
    fprintf('Vertical error (end): %.2f m.\n',stat.v_err_end);

    
    % speed rmse
    speed_magOdometry = sqrt(sum(out_data.x_h(4:6, :).^2,1));
    speed_rmse = mean(abs(speed_t(idx) - speed_magOdometry(idx)));
    fprintf('The speed rmse is: %f.\n', speed_rmse);

    stat.speed_rmse = speed_rmse;
end
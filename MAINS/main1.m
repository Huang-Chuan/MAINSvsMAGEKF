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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           
% Main script for the MAINS system. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clearvars;
    
%% Load filter settings
disp('Loads settings')
settings = readConfig('config/config.json');

%% Load data
disp('Loads data')
load(settings.dataPath);


%% preprocess imu_array
data.u = data.imu_array.u;
data.t = data.gt.t;
data.numFrames = data.gt.numFrames;
data.sampleFreq = 100;
load(settings.calibrationFilePath)
data.mag_array.field = calibratemag(data.mag_array.field, theta_all(1:end-1,:), 30);
% interpolate
data.mag_array.field = interp1(data.mag_array.t, data.mag_array.field, data.t);
data.mag_array.field(end-5:end, :) = repmat(data.mag_array.field(end - 5, :), 6, 1); 
in_data = data;

%% Run the MAG-aided INS
B_b = permute(reshape(in_data.mag_array.field', 3, 30, []), [1,3,2]);
B_b = B_b(:, :,  settings.availableSensorIdx);
disp('Runs the MAG-aided INS')
out_data=magaidedINS(in_data.sampleFreq,...
                     in_data.u(:, 1:3)', ...
                     in_data.u(:, 4:6)', ...
                     B_b,...
                     settings.sensor_locs(:, settings.availableSensorIdx), ...
                     in_data.gt, ...
                     settings);

%% Compute the error
disp('Plot data')
viewresult(out_data, in_data, settings);

%% The following function works with in_data that has the field "gt"
stat=compute_stats(in_data, out_data);

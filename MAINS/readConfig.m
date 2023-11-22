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
function [settings] = readConfig(configFilePath)

    try
        configData = jsondecode(fileread(configFilePath));
    catch
        error('Error reading or parsing the JSON configuration file.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           LOAD DATA & Calibration       %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scenarioName = configData.ScenarioName;
    settings.scenarioName = scenarioName;
    
    settings.dataPath = fullfile(configData.DataFolderPath,[scenarioName, '.mat']);
    
    calibrationFileName = configData.CalibrationFileName;
    settings.calibrationFilePath = fullfile(configData.CalibrationFolderPath,[calibrationFileName, '.mat']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             MAG        &          POS Aiding  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings.magAiding   = configData.MagAiding.IsTurnedOn;
    settings.magAidingStartTime = configData.MagAiding.StartTime;
    settings.magAidingEndTime = configData.MagAiding.EndTime;
    
    
    settings.posAiding = configData.PosAiding.IsTurnedOn;
    settings.posAidingStartTime = configData.PosAiding.StartTime;
    settings.posAidingEndTime = configData.PosAiding.EndTime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             SENSOR PARAMETERS           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % sensor locations
    settings.sensor_locs = configData.SensorConfig.Locations';
    % active sensor number
    settings.availableSensorIdx = configData.SensorConfig.ActiveSensors';
    % number of sensors 
    settings.numSensors = length(settings.availableSensorIdx);
    % sampling frequnecy
    settings.fs = configData.SensorConfig.Frequency;
    % sampling interval
    settings.dT = 1 / settings.fs;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             POLYNOMIAL ORDER            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings.polyOrder = configData.PolynomialModel.Order;
    settings.dimTheta  = settings.polyOrder^2 + 4 * settings.polyOrder + 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             STATE MASKS                 %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    settings.magSensorBiasInclude = false;
    [settings.numErrorStates, settings.errorStateMask] = makeErrorStateMask(settings);
    [settings.numStates, settings.stateMask] = makeStateMask(settings);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             GRAVITY VECTOR              %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    settings.g = configData.Gravity;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    INIT  uncertainties （standard deviation)     %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    settings.init_pos_std = configData.InitUncertainty.Init_Pos_Std;                                        % Position [m]
    settings.init_vel_std = configData.InitUncertainty.Init_Vel_Std;                                        % Velocity [m/s]
    settings.init_q_std = (pi/180)*configData.InitUncertainty.Init_Ori_Std*ones(3,1);                       % Attitude (roll,pitch,yaw) [rad]
    settings.init_accBias_std = configData.InitUncertainty.Init_Acc_Bias_Std;                               % Accelerometer biases [m/s^2]
    settings.init_gyroBias_std = (pi/180)*configData.InitUncertainty.Init_Gyro_Bias_Std;                    % Gyro biases [rad/s]       
    settings.init_coeff_std = configData.InitUncertainty.Init_Theta_Std;                                    % Coefficients    
    settings.init_magBias_std = configData.InitUncertainty.Init_Mag_Bias_Std;                               % magnetometer biases [mu T]


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             FILTER PARAMETERS           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %-------------------------------------------- -------------------%
    %--------------- Process noise covariance (Q) -------------------%
    %-------------------------------------------- -------------------%
    % IMU measurement noise
    settings.sigma_acc_w = configData.FilterParameters.STD_Acc_Noise;                                            % unit: [m/s^2]
    settings.sigma_gyro_w = pi/180*configData.FilterParameters.STD_Gyro_Noise;                                   % unit: [rad/s]

    % IMU bias random walk noise 
    settings.sigma_acc_bias_rw = configData.FilterParameters.STD_Acc_Bias_Random_Walk;                           % unit: [m/s^(5/2)]
    settings.sigma_gyro_bias_rw = pi/180*configData.FilterParameters.STD_Gyro_Bias_Random_Walk;                  % unit: [rad/s^(3/2)]
    % MAG bias random walk noise
    settings.sigma_mag_bias_rw =  configData.FilterParameters.STD_Mag_Bias_Random_Walk;                          % unit: [mu T/s^(3/2)]
    
    % Theta random walk noise
    settings.sigma_coeff_w = configData.FilterParameters.STD_Theta_Random_Walk';
    if length(settings.sigma_coeff_w) ~= settings.dimTheta
        error('Dimension of theta random walk noise does not match the dimension of theta')
    end

    %-------------------------------------------- -------------------%
    % ---------- Measurement noise covariance (R) -------------------%
    %-------------------------------------------- -------------------%
    
    % Position aiding 
    settings.RPos = (configData.FilterParameters.STD_Pos)^2 * eye(3);                     % unit: [m^2]


end

function [numErrorStates, masks] = makeErrorStateMask(settings)
    
    numErrorStates = 15 + settings.dimTheta;
    
    pos = 1 : 3;
    vel = 4 : 6;
    epsilon = 7 : 9;
    acc_bias = 10 : 12;
    gyro_bias = 13 : 15;
    theta     = 16 : numErrorStates;

    masks.pos       = false(numErrorStates, 1);
    masks.vel       = false(numErrorStates, 1);
    masks.epsilon   = false(numErrorStates, 1);
    masks.acc_bias  = false(numErrorStates, 1);
    masks.gyro_bias = false(numErrorStates, 1);
    masks.theta     = false(numErrorStates, 1);

    masks.pos(pos) = true;
    masks.vel(vel) = true;
    masks.epsilon(epsilon) = true;
    masks.acc_bias(acc_bias) = true;
    masks.gyro_bias(gyro_bias) = true;
    masks.theta(theta)  = true;

    if (settings.magSensorBiasInclude)
        numMagBiasStates = (settings.numSensors - 1) * 3;
        % update error state dimentions
        numErrorStates = numErrorStates + numMagBiasStates;
        % update masks
        masks.pos = [masks.pos; false(numMagBiasStates, 1)];
        masks.vel = [masks.vel; false(numMagBiasStates, 1)];
        masks.epsilon  = [masks.epsilon; false(numMagBiasStates, 1)];
        masks.acc_bias  = [masks.acc_bias; false(numMagBiasStates, 1)];
        masks.gyro_bias = [masks.gyro_bias; false(numMagBiasStates, 1)];
        masks.theta     = [masks.theta; false(numMagBiasStates, 1)];
        masks.mag_bias = false(numErrorStates, 1);
        masks.mag_bias(end - (numMagBiasStates - 1) : end) = true;
    end

end



function [numStates, masks] = makeStateMask(settings)
    
    numStates = 16 + settings.dimTheta;

    pos = 1 : 3;
    vel = 4 : 6;
    q_nb = 7 : 10;
    acc_bias = 11 : 13;
    gyro_bias = 14 : 16;
    theta     = 17 : numStates;
    
    masks.pos = false(numStates, 1);
    masks.vel = false(numStates, 1);
    masks.q_nb = false(numStates, 1);
    masks.acc_bias = false(numStates, 1);
    masks.gyro_bias = false(numStates, 1);
    masks.theta     = false(numStates, 1);

    masks.pos(pos) = true;
    masks.vel(vel) = true;
    masks.q_nb(q_nb) = true;
    masks.acc_bias(acc_bias) = true;
    masks.gyro_bias(gyro_bias) = true;
    masks.theta(theta)  = true;

    if (settings.magSensorBiasInclude)
        numMagBiasStates = (settings.numSensors - 1) * 3;
        % update error state dimentions
        numStates = numStates + numMagBiasStates;
        % update masks
        masks.pos = [masks.pos; false(numMagBiasStates, 1)];
        masks.vel = [masks.vel; false(numMagBiasStates, 1)];
        masks.q_nb  = [masks.q_nb; false(numMagBiasStates, 1)];
        masks.acc_bias  = [masks.acc_bias; false(numMagBiasStates, 1)];
        masks.gyro_bias = [masks.gyro_bias; false(numMagBiasStates, 1)];
        masks.theta     = [masks.theta; false(numMagBiasStates, 1)];
        masks.mag_bias = false(numStates, 1);
        masks.mag_bias(end - (numMagBiasStates - 1) : end) = true;
    end
end


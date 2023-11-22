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
function out_data=magaidedINS(fs, acc_b, gyr_b, B_b, sensor_pos, gt, settings)

    % Copy data to variables with shorter name
    u=[acc_b; gyr_b]';
    mag=reshape(permute(B_b,[2,1,3]), size(u, 1), []);
    % overwrite 
    settings.sensor_locs = sensor_pos;
    settings.numSensors = size(settings.sensor_locs, 2);

    magEnable = settings.magAiding;
    magStartTime = settings.magAidingStartTime;
    magEndTime = settings.magAidingEndTime;


    posEnable = settings.posAiding;
    posStartTime = settings.posAidingStartTime;
    posEndTime = settings.posAidingEndTime;

    stateMask = settings.stateMask;


    % Get measurement equations
    [H, HPos, Phi]=getHs(settings);
   
    %% Initialization
    % Initialize the navigation state
    z = reshape(mag(1, :), 3, []);
    % initiate polynomial magnetic model
    m = polyMagModel(settings.polyOrder);
    m = m.set_phi(sensor_pos);
    m = m.init_theta(z(:));
    [x_h, P]=init_navigation_state(m,  u(1, 1:3), gt, settings);

    % Get process noise covariance and measurement noise covariance
    Q = getProcessNoiseCov(settings);
    settings.Q = Q;
    RPos = settings.RPos;

    % Allocate memory for the output data
    N=size(u,1);
    out_data.x_h=zeros(settings.numStates,N);
    out_data.x_h(:,1)=x_h;
    out_data.diag_P=zeros(settings.numErrorStates,N);
    out_data.diag_P(:,1)=diag(P);
    out_data.res=zeros(3*settings.numSensors, N);
    out_data.nis=zeros(1, N);
    out_data.debug.theta = zeros(settings.dimTheta, N);
    out_data.debug.theta_pred = zeros(settings.dimTheta, N);
    out_data.debug.theta(:, 1) = x_h(end-settings.dimTheta+1 : end);
    out_data.debug.theta_pred(:, 1) = x_h(end-settings.dimTheta+1 : end);
    out_data.debug.noisevar = zeros(N, 1);

    t=(0:N-1)/fs;
    %% Information fusion
    for k=2:N
        
        % Sampling period
        Ts=1/fs;

        % calculate psi
        psi = calc_psi(x_h, u(k - 1, :), Ts, settings);

        % update theta
        m = m.update_theta(psi);
        out_data.debug.theta_pred(:, k) = m.get_theta();
        % Get state space model matrices
        [F,G]=state_space_model(x_h, u(k - 1, :), Ts, m, settings);
        
        % Update the nominal state state
        [x_h]=nominalStateProp(x_h, u(k - 1, :), Ts, m, settings);
        
        % Time update of the Kalman filter state covariance.
        P=F * P * F'+ G * Q * G';
        
        % Position aiding update
        if posEnable && (t(k) >= posStartTime && t(k) <= posEndTime) && (~any(isnan(gt.pos(:, k)), 'all')) 
            K = P * HPos' / (HPos * P * HPos' + RPos); 
            delta_z = gt.pos(:, k) - x_h(stateMask.pos);
            P = (eye(size(P)) - K * HPos) * P * (eye(size(P)) - K * HPos)' + K * RPos * K';
            delta_x = K * delta_z;
            [x_h, P]=correctNominalState(x_h, P, delta_x, settings);
        end    


        % Mag aiding update
        if magEnable && (t(k) >= magStartTime && t(k) <= magEndTime)
            z = reshape(mag(k, :), 3, []);
            delta_z = z(:) - Phi * x_h(settings.stateMask.theta);
            if(settings.magSensorBiasInclude)
                delta_z = delta_z - [x_h(settings.stateMask.mag_bias); zeros(3, 1)];
            end
            out_data.res(:, k) = delta_z;
            [trueCoeff, resVar] = m.LS_coeff(z(:));
            out_data.debug.noisevar(k) = resVar;
            out_data.debug.theta(:, k) = trueCoeff;
            R = resVar * eye(3*settings.numSensors);
            K = P * H' / (H * P * H' + R);   
            P = (eye(size(P)) - K * H) * P * (eye(size(P)) - K * H)' + K * R * K';
            delta_x = K * delta_z;
            [x_h, P]=correctNominalState(x_h, P, delta_x, settings);
            % update theta in model 
            m = m.set_theta(x_h(end-settings.dimTheta + 1 : end));
        end

        % Save the data to the output data structure
        out_data.x_h(:,k)=x_h;
        out_data.diag_P(:,k)=diag(P);
        
    end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          SUB-FUNCTIONS                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Process Noise  Covariance     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q] = getProcessNoiseCov(settings)
    Q = blkdiag(settings.sigma_acc_w^2 * eye(3), ...
        settings.sigma_gyro_w^2 * eye(3), ...
        settings.sigma_acc_bias_rw^2 * eye(3), ...
        settings.sigma_gyro_bias_rw^2* eye(3), ...
        diag(settings.sigma_coeff_w.^2));
    if(settings.magSensorBiasInclude)
        Q = blkdiag(Q, settings.sigma_mag_bias_rw^2 * eye((settings.numSensors - 1) * 3));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Init navigation state     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_h, P]=init_navigation_state(m, u, gt, settings)
    % extract initial pose from gt
    if isfield(gt,'pos')
        init_pos = gt.pos(:, 1);
    else
        init_pos = [0;0;0];
    end

    if isfield(gt,'ori')
        init_q = rotm2quat(gt.ori(:, :, 1))';
    else
        init_q = initialize_pose(u);
    end

    % Initial nominal state vector (position, velocity, orientation, inertial sensor biases, theta, magnetometer bias) 
    if settings.magAiding
        x_h=[init_pos; 
            zeros(3, 1); 
            init_q; 
            zeros(6, 1); 
            m.get_theta()];
    else
        x_h=[init_pos; 
            zeros(3, 1); 
            init_q; 
            zeros(6, 1); 
            zeros(settings.dimTheta, 1)];
    end
    
    if(settings.magSensorBiasInclude)
        x_h = [x_h; zeros(3 * (settings.numSensors - 1), 1)];
    end

    % Initial error state uncertainties (position, velocity, orientation, inertial sensor biases, theta, magnetometer bias) 
    P = diag([settings.init_pos_std^2 * ones(3, 1); 
            settings.init_vel_std^2 * ones(3, 1);  
            settings.init_q_std.^2;
            settings.init_accBias_std^2 * ones(3, 1); 
            settings.init_gyroBias_std^2 * ones(3, 1);
            settings.init_coeff_std^2 * ones(settings.dimTheta, 1)]);
    if(settings.magSensorBiasInclude)
        P = blkdiag(P, settings.init_magBias_std^2 * eye(3 * (settings.numSensors - 1)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Calculate Psi       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = calc_psi(x_h, u, dt, settings)
    mask = settings.stateMask;
    
    vk = x_h(mask.vel);
    R_nb = q2r(x_h(mask.q_nb));
    acc_bias = x_h(mask.acc_bias);
    gyro_bias = x_h(mask.gyro_bias);
    
    %  parse u
    acc_m = u(1:3)';
    omega_m = u(4:end)';

    acc_nav = R_nb * (acc_m - acc_bias) + settings.g;
    dp = vk * dt + 1/2 * acc_nav * dt^2;
    dp_body = R_nb.' * dp;

    omega_h = omega_m - gyro_bias;
    rotvec = omega_h * dt;

    psi = [dp_body; rotvec];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  State transition matrix   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,G]=state_space_model(x_h,u,dt,m,settings)
    masks = settings.stateMask;
    numErrorStates = settings.numErrorStates;
    errorMasks = settings.errorStateMask;           

    % parse x_h
    vk = x_h(masks.vel);
    R_nb = q2r(x_h(masks.q_nb));
    acc_bias = x_h(masks.acc_bias);
    gyro_bias = x_h(masks.gyro_bias);
    %  parse u
    acc_m = u(1:3)';
    omega_m = u(4:end)';

    % remove bias from gyro readings
    omega_h = omega_m - gyro_bias;
    rotvec = omega_h * dt;

    J1J2 = m.J1J2;

    % M as in eq. 23d
    M = zeros(settings.dimTheta + 6, numErrorStates);                        
    M(1:settings.dimTheta, errorMasks.theta) = eye(settings.dimTheta);
    M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.vel)   = R_nb.' * dt;
    M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.epsilon)  = vect2skew(R_nb.' *  dt * (vk  + settings.g * dt / 2));
    M(settings.dimTheta + 1: settings.dimTheta + 3, errorMasks.acc_bias)   = -dt^2 / 2 * eye(3);
    M(settings.dimTheta + 4: settings.dimTheta + 6, errorMasks.gyro_bias)  = -eye(3)*dt;


    % construct F
    F = zeros(numErrorStates);
    F(errorMasks.pos, errorMasks.pos) = eye(3);
    F(errorMasks.pos, errorMasks.vel) = eye(3) * dt;
    F(errorMasks.vel, errorMasks.vel) = eye(3);
    F(errorMasks.vel, errorMasks.epsilon) = -R_nb * vect2skew(acc_m - acc_bias) * dt;
    F(errorMasks.vel, errorMasks.acc_bias) = -R_nb * dt;
    F(errorMasks.epsilon, errorMasks.epsilon) = eye(3) - vect2skew(rotvec);
    F(errorMasks.epsilon, errorMasks.gyro_bias) = -eye(3)*dt;
    F(errorMasks.acc_bias, errorMasks.acc_bias) = eye(3);
    F(errorMasks.gyro_bias, errorMasks.gyro_bias) = eye(3);    
    F(errorMasks.theta, :) = m.pinvA * ([m.B J1J2] * M) ;

    % construct G
    G = zeros(numErrorStates, size(settings.Q, 1));
    G(errorMasks.vel, 1:3) = eye(3) * dt;
    G(errorMasks.epsilon, 4:6) = eye(3) * dt;
    G(errorMasks.acc_bias, 7:9) = eye(3) * sqrt(dt);
    G(errorMasks.gyro_bias, 10:12) = eye(3) * sqrt(dt);
    G(errorMasks.theta, 1 : 3) = -m.pinvA * (J1J2(:, 1:3) * (dt^2/2));
    G(errorMasks.theta, 4 : 6) = -m.pinvA * (J1J2(:, end-2:end))*dt;
    G(errorMasks.theta, 12 + 1 : 12 + settings.dimTheta) = eye(settings.dimTheta);

    if(settings.magSensorBiasInclude)
        F(errorMasks.mag_bias, errorMasks.mag_bias) = eye(3 * (settings.numSensors - 1));
        G(errorMasks.mag_bias, end - 3 * (settings.numSensors - 1) + 1: end) = eye(3 * (settings.numSensors - 1)) * sqrt(dt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Nominal state propogation        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x] = nominalStateProp(xk, u, dt, m, settings)
    %   INPUT:
    %                 xk:  current state estimate at time k
    %                  u:  accelerometer reading and gyroscope reading
    %                 dt:  sampling interval
    %           settings:  struct setting 
    %   
    %   OUTPUT:
    %                  x:  predict nominal state at time k + 1            

        masks = settings.stateMask;        
        
        % parse xk
        pk = xk(masks.pos);
        vk = xk(masks.vel);
        q_nb = xk(masks.q_nb);
        R_nb = q2r(q_nb);
        acc_bias = xk(masks.acc_bias);
        gyro_bias = xk(masks.gyro_bias);
        %  parse u
        acc_m = u(1:3)';
        omega_m = u(4:end)';
    
        acc_nav = R_nb * (acc_m - acc_bias) + settings.g;
        dp = vk * dt + 1/2 * acc_nav * dt^2;
        
        % nominal state \hat{x}_k as in eq. 16
        x = zeros(size(xk));
        x(masks.pos) = pk + dp;
        x(masks.vel) = vk + acc_nav * dt;
        x(masks.acc_bias) = acc_bias;
        x(masks.gyro_bias) = gyro_bias;
        
        omega_h = omega_m - gyro_bias;
        rotvec = omega_h * dt;
        dq = rotvec2quat(rotvec');
        qw = dq(1);
        qv = dq(2:4)';
        x(masks.q_nb) = (qw*eye(4) + [0 -qv'; qv -vect2skew(qv)]) * q_nb;
        x(masks.theta) = m.get_theta();

        if(settings.magSensorBiasInclude)
            x(masks.mag_bias) = xk(masks.mag_bias);
        end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Measurement Equation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, HPos, Phi]=getHs(settings)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             Measurement Matrix          %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Phi = [];
    for i = 1 : length(settings.sensor_locs)
        Phi =  [Phi; get_ABnull(settings.sensor_locs(:, i), settings.polyOrder)];
    end
    
    H = zeros(settings.numSensors * 3, settings.numErrorStates);
    H(:, settings.errorStateMask.theta) = Phi;
    if(settings.magSensorBiasInclude)
        H(1:end-3, settings.errorStateMask.mag_bias) = eye(3 * (settings.numSensors - 1));
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%             Auxiliary aiding            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HPos = zeros(3, settings.numErrorStates);
    HPos(:, settings.errorStateMask.pos) = eye(3);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Nominal State Correction         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xh, P]=correctNominalState(xh, P, delta_x, settings)
    
    stateMask = settings.stateMask;
    errorStateMask =  settings.errorStateMask;

    xh(stateMask.pos) = xh(stateMask.pos) +  delta_x(errorStateMask.pos);
    xh(stateMask.vel) = xh(stateMask.vel) +  delta_x(errorStateMask.vel);
    qv = 1/2 * delta_x(errorStateMask.epsilon);
    xh(stateMask.q_nb) = (eye(4) + [0 -qv'; qv -vect2skew(qv)]) * xh(stateMask.q_nb);
    xh(stateMask.q_nb) = xh(stateMask.q_nb).' / norm(xh(stateMask.q_nb)); 
    xh(stateMask.acc_bias) = xh(stateMask.acc_bias) +  delta_x(errorStateMask.acc_bias);
    xh(stateMask.gyro_bias) = xh(stateMask.gyro_bias) +  delta_x(errorStateMask.gyro_bias);
    xh(stateMask.theta)     = xh(stateMask.theta) +  delta_x(errorStateMask.theta);

    if(settings.magSensorBiasInclude)
        xh(stateMask.mag_bias)     = xh(stateMask.mag_bias) +  delta_x(errorStateMask.mag_bias);
    end

end


function [q0] = initialize_pose(u)
    % initialize orientation based accelerometer
    yaw   = 0;      
    pitch = atan2(-u(1), sqrt(u(2)^2 + u(3)^2));
    roll  = atan2(u(2), u(3));
    q0 = eul2quat([yaw pitch roll])';
end
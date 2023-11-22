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
function [] = viewresult(out_data, in_data, settings)
    timeVector = in_data.t(1:end);
    X = out_data.x_h(:, 1:end)';
    Ps = out_data.diag_P(:, 1:end).';
    ref.rotations = in_data.gt.ori(:, :, 1:end);
    ref.positions = in_data.gt.pos(:, 1:end);
    stateMask = settings.stateMask;
    errorStateMask = settings.errorStateMask;
    fieldRes = out_data.res(:,(1:end));
    %nis      = out_data.nis(1:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%               position plot             %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    t = tiledlayout(3,2);
    % plot horizontal plane trajectory
    nexttile([2, 2]);
    plot(X(:,1),X(:,2),'k')
    hold on;
    plot(X(1,1),X(1,2),'rs')
    plot(X(end,1),X(end,2),'r*')
    error_ellipse(diag(Ps(end,1:2)),X(end,[1 2]),'conf',0.95,'style','b');

    if exist('ref','var')
        plot(ref.positions(1,:),ref.positions(2,:),'b-.');
    end

    % plot final heading
    quiver(X(end,1), X(end,2), ref.rotations(1, 1, end), ref.rotations(2, 1, end), 'r');
    quiver(X(end,1), X(end,2), ref.rotations(1, 2, end), ref.rotations(2, 2, end), 'b');
    rotm = quat2rotm(X(end, 7:10));
    quiver(X(end,1), X(end,2), rotm(1, 1), rotm(2, 1), 'r--');
    quiver(X(end,1), X(end,2), rotm(1, 2), rotm(2, 2), 'b--');
    title('2D Trajectory')
    legend('Trajectory','Start point','End point','95% conf.','Qualisys reference')
    xlabel('x [m]')
    ylabel('y [m]')
    
    grid on
    box on
    
    % plot vertical height
    nexttile([1, 2]);
    plot(timeVector, X(:, 3),'k')
    hold on
    plot(timeVector, ref.positions(3,:),'b-.');
    title('Heigth')
    legend('Estimated height','Qualisys reference')
    xlabel('time [s]')
    ylabel('z [m]')
    grid on
    box on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           height and speed plot         %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    subplot(2,1,1)
    plot(timeVector, X(:, 3),'r')
    hold on;
    curve1 =  X(:, 3)+sqrt(Ps(:, 3));
    curve2 =  X(:, 3)-sqrt(Ps(:, 3));
    h = fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5],'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
    %ylim([-2 2])
    title('Heigth')
    xlabel('time [s]')
    ylabel('z [m]')
    grid on
    box on

    subplot(2,1,2)
    plot(timeVector,sqrt(sum(X(:, 4:6).^2, 2)))
    title('Speed')
    xlabel('time [s]')
    ylabel('|v| [m/s]')
    grid on
    box on

    figure
    for ii=1:3
        subplot(3, 1, ii)
        plot(timeVector,X(:, 3 + ii),'r')
        if ii==1
            title('Velocity')
        end
        hold on;
        curve1 =  X(:, 3 + ii)+sqrt(Ps(:, 3 + ii));
        curve2 =  X(:, 3 + ii)-sqrt(Ps(:, 3 + ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel('Speed [m/s]')
        xlabel('time [s]')
        grid minor
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%              Orientation plot           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    labels = ['Yaw', 'Pitch', "Roll"];
    eul_deg = rad2deg(quat2eul(X(:, 7:10),'ZYX'));
    eulCov = zeros(length(timeVector), 3);
    for i = 1 : length(timeVector)
        J = ypr_jacobian_quat(X(i, 7:10));
        eulCov(i, :) = diag(J * diag(out_data.diag_P(7:9 , i)) * J');
    end
    figure
    for ii=1:3
        subplot(3,1,ii)
        plot(timeVector,eul_deg(:, ii),'r')
        hold on
        plot(timeVector,in_data.gt.rpys(ii, :), 'k')
        if ii==1
            title('Attitude')
        end
        hold on;
        % compute correct jacobian for each estimated orientation
        curve1 =  eul_deg(:, ii) + rad2deg(sqrt(eulCov(:, ii)));
        curve2 =  eul_deg(:, ii) - rad2deg(sqrt(eulCov(:, ii)));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel(strcat(labels(ii), '[deg]'))
        xlabel('time [s]')
        grid minor
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%          Acceleration bias plot         %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    acc_bias = X(:, 11:13);
    acc_bias_var = Ps(:, 10:12);
    figure;
    for ii = 1 : 3
        subplot(3, 1, ii)
        plot(timeVector,acc_bias(:, ii),'r')
        if ii==1
            title('Accelerometer bias')
        end
        hold on;
        curve1 =  acc_bias(:, ii) + sqrt(acc_bias_var(:, ii));
        curve2 =  acc_bias(:, ii) - sqrt(acc_bias_var(:, ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        hold on
        ylabel('Bias [m/s^2]')
        xlabel('time [s]')
        grid minor
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           Gyroscope bias plot           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gyro_bias = rad2deg(X(:, 14:16));
    gyro_bias_var = Ps(:, 13:15);
    figure;
    for ii = 1 : 3
        subplot(3, 1, ii)
        plot(timeVector, gyro_bias(:, ii),'r')
        if ii==1
            title('Gyro bias')
        end
        hold on;
        curve1 =  gyro_bias(:, ii) + rad2deg(sqrt(gyro_bias_var(:, ii)));
        curve2 =  gyro_bias(:, ii) - rad2deg(sqrt(gyro_bias_var(:, ii)));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        ylabel('Bias [deg/s]')
        xlabel('time [s]')
        grid minor
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%           Magnetometer bias plot        %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if settings.magSensorBiasInclude
        mag_bias = X(:, settings.stateMask.mag_bias);
        mag_bias_var = Ps(:, settings.errorStateMask.mag_bias);
        t = tiledlayout(5, 6);
        for i = 1 : 29
            nexttile,
            hold on;
            plot(timeVector, mag_bias(:, 3 * (i - 1) + 1), 'r', 'LineWidth', 2);
            curve1 =  mag_bias(:, 3 * (i - 1) + 1) + sqrt(mag_bias_var(:, 3 * (i - 1) + 1));
            curve2 =  mag_bias(:, 3 * (i - 1) + 1) - sqrt(mag_bias_var(:, 3 * (i - 1) + 1));
            fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)],  'r', 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
            
            plot(timeVector, mag_bias(:, 3 * (i - 1) + 2), 'g', 'LineWidth', 2);
            curve1 =  mag_bias(:, 3 * (i - 1) + 2) + sqrt(mag_bias_var(:, 3 * (i - 1) + 2));
            curve2 =  mag_bias(:, 3 * (i - 1) + 2) - sqrt(mag_bias_var(:, 3 * (i - 1) + 2));
            fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)],  'g', 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);

            plot(timeVector, mag_bias(:, 3 * (i - 1) + 3), 'b', 'LineWidth', 2);
            curve1 =  mag_bias(:, 3 * (i - 1) + 3) + sqrt(mag_bias_var(:, 3 * (i - 1) + 3));
            curve2 =  mag_bias(:, 3 * (i - 1) + 3) - sqrt(mag_bias_var(:, 3 * (i - 1) + 3));
            fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)],  'b', 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
        
            % title(strcat('$\theta_{', num2str(i), '}$'),'FontSize',12,'FontName','Times New Roman','interpreter','latex');
            grid minor;
            
            box on
        end

        title(t, 'Magnetometer bias');
        xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
        ylabel(t, ' ','FontSize',12,'FontName','Times New Roman')

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%        coefficient      plot            %%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coeff = X(:, stateMask.theta);
    coeff_var = Ps(:, errorStateMask.theta);
    f=figure;
    t = tiledlayout(3, 5);
    for i = 1 : min(15, size(coeff, 2))
        nexttile,
        hold on;
        plot(timeVector, coeff(:, i), 'r', 'LineWidth', 2);
        
        curve1 =  coeff(:, i) + sqrt(coeff_var(:, ii));
        curve2 =  coeff(:, i) - sqrt(coeff_var(:, ii));
        fill([timeVector fliplr(timeVector)], [curve1; flipud(curve2)], [.5 .5 .5], 'FaceAlpha',.2,'EdgeColor', [.5 .5 .5]);
    
        title(strcat('$\theta_{', num2str(i), '}$'),'FontSize',12,'FontName','Times New Roman','interpreter','latex');
        grid minor;
        
        box on
    end
    h  = axes(f, 'visible', 'off');
    title(t, 'Coefficient');
    xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel(t, ' ','FontSize',12,'FontName','Times New Roman')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%          residual       plot            %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    f=figure;
    t = tiledlayout(5, 6);
    for i = 1 : settings.numSensors
        nexttile,
        hold on;
        plot(timeVector, fieldRes(3*(i-1)+1, :), 'r', 'LineWidth', 2);
        plot(timeVector, fieldRes(3*(i-1)+2, :), 'g', 'LineWidth', 2);
        plot(timeVector, fieldRes(3*(i-1)+3, :), 'b', 'LineWidth', 2);
        grid minor;
        
        box on
    end
    h  = axes(f, 'visible', 'off');
    title(t, 'Residual');
    xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
    ylabel(t, ' ','FontSize',12,'FontName','Times New Roman')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%             NIS          plot           %% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure;
    % plot(timeVector,nis);
    % hold on;
    % plot(timeVector, ones(length(timeVector),1) * settings.threshold, '--');
    % title('NIS plot');
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%          signal variance plot           %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magSig = reshape(mag', 3, size(mag, 2)/3, []);
% sigVar = squeeze(mean((magSig-mean(magSig,2)).^2, [1 2]));
% figure;
% plot(timeVector, 10 * log10(sigVar),'k','LineWidth',2);
% xlabel(t,'time [s]','FontSize',12,'FontName','Times New Roman')
% ylabel('Signal variance [dB]','FontSize',12,'interpreter','latex','FontName','Times New Roman')
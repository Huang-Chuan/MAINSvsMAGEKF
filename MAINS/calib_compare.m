
clearvars

datalists={'bigsquare1.mat',...
'bigsquare2.mat',...
'bigsquare3.mat',...
'bigsquare1normalheight.mat',...
'bigsquare2normalheight.mat',...
'bigsquare3normalheight.mat',...
'bigsquare1tilted.mat',...
'bigsquare2tilted.mat',...
};

calibfiles={'calibJointMAP.mat',...
'calibEKF.mat',...
'calibML.mat',...
'nocalib.mat'};





%% Load filter settings
disp('Loads settings')
settings = readConfig('config/config.json');



%% show example LP-2





for j = 1 : length(calibfiles)
    
    disp('Loads data')
    settings.dataPath = ['data/', datalists{2}];
    load(settings.dataPath);
    
    settings.calibrationFilePath = ['calibration/', calibfiles{j}];    
    load(settings.calibrationFilePath);

    
    data.u = data.imu_array.u;
    data.t = data.gt.t;
    data.numFrames = data.gt.numFrames;
    data.sampleFreq = 100;

    data.mag_array.field = calibratemag(data.mag_array.field, theta_all(1:end,:), 30);
    % interpolate
    data.mag_array.field = interp1(data.mag_array.t, data.mag_array.field, data.t);
    data.mag_array.field(end-5:end, :) = repmat(data.mag_array.field(end - 5, :), 6, 1); 
    in_data = data;
    clear data

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
    if(strcmp(calibfiles{j}, 'calibJointMAP.mat'))
        traj.proposed = out_data.x_h(1:3,:);
    elseif(strcmp(calibfiles{j}, 'calibEKF.mat'))
        traj.wu = out_data.x_h(1:3,:);
    elseif (strcmp(calibfiles{j}, 'calibML.mat'))
        traj.kok = out_data.x_h(1:3,:);
    else
        traj.nocalib = out_data.x_h(1:3,:);
    end

end


% plot
bar_colors = [0 0.4470 0.7410;  % EKF (Blue)
              0.8500 0.3250 0.0980;  % ML (Orange)
              0.4660 0.6740 0.1880]; % Joint MAP (Green)

figure;
plot(traj.proposed(1,:),traj.proposed(2,:), 'Color',bar_colors(3,:), 'DisplayName', 'Proposed');
hold on;
plot(traj.kok(1,:),traj.kok(2,:), 'Color',bar_colors(2,:), 'DisplayName', 'Kok et al.');
plot(traj.wu(1,:),traj.wu(2,:), 'Color',bar_colors(1,:), 'DisplayName', 'Wu et al.');
plot(in_data.gt.pos(1, :),in_data.gt.pos(2, :), 'k', 'DisplayName', 'True trajectory');
grid on;
xlabel('x [m]','FontName','Times New Roman')
ylabel('y [m]','FontName','Times New Roman')

legend;








%% 
fp = fopen('results.txt', 'w');

for i = 1 : length(datalists)
    settings.dataPath = ['data/', datalists{i}];
    fprintf(fp, '%s\n', datalists{i});
    

    for j = 1 : length(calibfiles)
        disp('Loads data')
        load(settings.dataPath);

        
        settings.calibrationFilePath = ['calibration/', calibfiles{j}];    
        load(settings.calibrationFilePath);
        %% preprocess imu_array
        data.u = data.imu_array.u;
        data.t = data.gt.t;
        data.numFrames = data.gt.numFrames;
        data.sampleFreq = 100;

        data.mag_array.field = calibratemag(data.mag_array.field, theta_all(1:end,:), 30);
        % interpolate
        data.mag_array.field = interp1(data.mag_array.t, data.mag_array.field, data.t);
        data.mag_array.field(end-5:end, :) = repmat(data.mag_array.field(end - 5, :), 6, 1); 
        in_data = data;
        clear data

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
        %% The following function works with in_data that has the field "gt"
        stat=compute_stats(in_data, out_data);        

        fprintf(fp, '%s %f %f\n', calibfiles{j}, stat.h_err, stat.v_err);
    end


end


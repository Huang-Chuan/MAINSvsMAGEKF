function [output] = compare_MAGEKF_MAINS(DATAFILE, run_MAGEKF)
% Input:  DATAFILE: name of the data file
%         run_MAGEKF: boolean, if true, run MAGEKF
% Output: output: struct containing the results

data_name = DATAFILE;
% select sensors
sensor_sel = [3, 6, 27, 30, 16];
%sensor_sel = [1:30];


%% Load data
disp('Loads data')
load(fullfile('data/',data_name));

data.u = data.imu_array.u;
data.u(:, 4:6) = data.u(:, 4:6) - mean(data.u(1:200, 4:6));  % remove bias
data.t = data.gt.t;
data.numFrames = data.gt.numFrames;
data.sampleFreq = 100;
load('calibration/calibResults_magcal0601_EKF.mat')
data.mag_array.field = calibratemag(data.mag_array.field, theta_all(1:end-1,:), 30);
% interpolate
data.mag_array.field = interp1(data.mag_array.t, data.mag_array.field, data.t);
data.mag_array.field(end-5:end, :) = repmat(data.mag_array.field(end - 5, :), 6, 1); 
in_data = data;

% export data
out_data.mag = in_data.mag_array.field';
out_data.sensorPos = PosMagArray();
out_data.acc = data.u(:, 1:3)';
out_data.gyro = data.u(:, 4:6)';
out_data.gt_pos = in_data.gt.pos;
out_data.gt_ori = in_data.gt.ori;
[vel] = calc_vel(in_data);
out_data.gt_vel = vel;
%% 
clear data in_data setting theta_all

%% data conversion

% rotation matrix to convert data
Rot = eul2rotm([0 0 pi/2]);          % Ro2o' where o is the body frame in MAINS, o' is the body frame in MAGEKF
r = Rot * PosMagArray();             % Convert to the body frame in MAGEKF
sensor_pos = r(:, sensor_sel);
display_active_sensors(sensor_pos, r);
acc_b = Rot * out_data.acc;
gyr_b = Rot * out_data.gyro;
B_b   = zeros(size(out_data.mag));
for i = 1 : 30
    B_b((i-1)*3+1: 3*i, :) = Rot * out_data.mag((i-1)*3+1: 3*i, :); 
end
Rb2n = pagemtimes(out_data.gt_ori, Rot');
qtrue = rotm2quat(Rb2n);
qtrue = qtrue';

% reshape B_b
B_b = permute(reshape(B_b, 3, 30, []), [1,3,2]);
% select B_b
B_b = B_b(:, :, sensor_sel);

fs = 100;
T = 1/fs;
N=length(qtrue);



% Ground truth velocity
vntrue = out_data.gt_vel;
ptrue  = out_data.gt_pos;
gt.pos = ptrue;
gt.ori = Rb2n;
%% Hassen
startT = 60;
x0 = [qtrue(:, 1); quat2rotm(qtrue(:, 1)')'*vntrue(:, 1); B_b(:,1,5); zeros(5,1); ptrue(:,1)];
if run_MAGEKF
  % Start measuring computation time
  tic;
  Xestimated2g = EKF_only(fs, ...
      acc_b,...
      gyr_b,...
      B_b,...
      sensor_pos,...
      x0,...
      gt,...
      startT * fs);
  % Stop measuring computation time
  computation_time = toc;
  % Display the computation time
  disp(['MAGEKF Computation Time: ' num2str(computation_time) ' seconds']);
else
  Xestimated2g = zeros(18, length(qtrue));
  Xestimated2g(1:4,:)= repmat([1;0;0;0], 1, length(qtrue));
end

%% MAINS
cd MAINS
settings=readConfig('config/config_rotated_body_frame.json');

% Start measuring computation time
tic;
% Your existing code here
mains = magaidedINS(fs, acc_b, gyr_b, B_b, sensor_pos, gt, settings);
% Stop measuring computation time
computation_time = toc;

% Display the computation time
disp(['MAINS Computation Time: ' num2str(computation_time) ' seconds']);

settings.magAiding = false;
mains_no_aiding = magaidedINS(fs, acc_b, gyr_b, B_b, sensor_pos, gt, settings);
cd ..

%% Calculate postional errors
output.INS = struct('pos', [], 'err', []);
output.MAINS = struct('pos', [], 'err', []);
output.hassen = struct('pos', [], 'err', []);
output.GT = struct('pos', ptrue);
[h_err1, v_err1] = calc_pos_err(ptrue, mains_no_aiding.x_h(1:3, :));
output.INS.err = [h_err1; v_err1];
output.INS.pos = mains_no_aiding.x_h(1:3, :);
[h_err2, v_err2] = calc_pos_err(ptrue, mains.x_h(1:3, :));
output.MAINS.err = [h_err2; v_err2];
output.MAINS.pos = mains.x_h(1:3, :);
if run_MAGEKF
  [h_err3, v_err3] = calc_pos_err(ptrue, Xestimated2g(16:18, :));
  output.hassen.err = [h_err3; v_err3];
  output.hassen.pos = Xestimated2g(16:18, :);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CUT DATA        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptrue = ptrue(:, startT * fs : end);
vntrue = vntrue(:, startT * fs : end);
qtrue = qtrue(:, startT * fs : end);

mains.x_h = mains.x_h(1:10, startT * fs: end);
mains_no_aiding.x_h = mains_no_aiding.x_h(1:10, startT * fs: end);

if run_MAGEKF
  Xestimated2g = Xestimated2g(:, startT * fs : end);

  Velb_Hassen = Xestimated2g(5:7,:);
  q_Hassen = Xestimated2g(1:4,:);
  Veln_Hassen = zeros(size(Velb_Hassen));
  for i = 1 : length(q_Hassen)
    Veln_Hassen(:, i) =  quat2rotm(q_Hassen(:, i)') * Velb_Hassen(:, i);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  t=startT-T:T:(N-1)*T;
  h1=figure('NumberTitle', 'off', 'Name', 'Velocity Plot');
  subplot(3,1,1)
  plot(t,vntrue(1,:),'r',t,Veln_Hassen(1,:),'b',t,mains.x_h(4,:),'g');
  l=legend('$v^n_x$','$\hat{v}^n_x$ (MAGEKF)', '$\hat{v}^n_x$ (MAINS)');
  set(l, 'interpreter', 'latex')
  xlabel('Time [s]');
  ylabel('$v^n_x$ [m$s^{-1}$]', 'interpreter', 'latex')
  grid on
  subplot(3,1,2)
  plot(t,vntrue(2,:),'r',t,Veln_Hassen(2,:),'b',t,mains.x_h(5,:),'g');
  l=legend('$v^n_y$','$\hat{v}^n_y$ (MAGEKF)', '$\hat{v}^n_y$ (MAINS)');
  set(l, 'interpreter', 'latex')
  xlabel('Time [s]');
  ylabel('$v^n_y$ [m$s^{-1}$]', 'interpreter', 'latex')
  grid on
  subplot(3,1,3)
  plot(t,vntrue(3,:),'r',t,Veln_Hassen(3,:),'b',t,mains.x_h(6,:),'g');
  l=legend('$v^n_z$','$\hat{v}^n_z$ (MAGEKF)', '$\hat{v}^n_z$ (MAINS)');
  set(l, 'interpreter', 'latex')
  xlabel('Time [s]');
  ylabel('$v^n_z$ [m$s^{-1}$]', 'interpreter', 'latex')
  grid on

  h2=figure('NumberTitle', 'off', 'Name', 'Walked trajectory3d');
  plot3(ptrue(1,:), ptrue(2,:),ptrue(3,:),'r',...
        Xestimated2g(16,:), Xestimated2g(17,:),Xestimated2g(18,:),'b',...
        mains.x_h(1,:),mains.x_h(2,:),mains.x_h(3,:),'g');
  legend('True trajectory', 'EKFMAG trajectory', 'MAINS trajectory');
  xlabel('x [m]');
  ylabel('y [m]');
  zlabel('z [m]');
  grid on;

  folderPath = fullfile('results',DATAFILE(1:end-4));
  % Check if the folder exists
  if ~exist(folderPath, 'dir')
      % The folder doesn't exist, so create it
      mkdir(folderPath);
      disp(['Folder created: ' folderPath]);
  else
      disp(['Folder already exists: ' folderPath]);
  end
  saveas(h1, fullfile(folderPath, 'velocity comparison'),'jpg')
  saveas(h2, fullfile(folderPath, 'position comparison'),'jpg')
  % %%%%%%%%%%%%%%%%%%%%%% Velocity Error Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(fullfile(folderPath, 'output.txt'), 'w');
if run_MAGEKF
  PrintStats(fileID,[ptrue;vntrue;qtrue], [Xestimated2g(16:18, :); Veln_Hassen; Xestimated2g(1:4, :)], "Method: MAGEKF")
end
PrintStats(fileID,[ptrue; vntrue; qtrue], mains.x_h(1:10, :), "Method: MAINS")
PrintStats(fileID,[ptrue; vntrue; qtrue], mains_no_aiding.x_h(1:10, :), "Method: Stand-alone INS")
end

function PrintStats(fileID,gt, x, title)
  dx = x(1:3, :)-gt(1:3,:);
  dv = x(4:6, :)-gt(4:6,:);
  dq = 2*acos(abs(sum(gt(7:10, :).*x(7:10, :), 1)));
  fprintf(fileID,'%s\n', title);
  fprintf(fileID,'\tMean pos [m]:    [%.4f, %.4f, %.4f]''\n', mean(dx, 2));
  fprintf(fileID,'\tHorizontal RMSE [m]:    %.4f\n', sqrt(mean(sum(dx(1:2,:).^2, 1))));
  fprintf(fileID,'\tHorizontal error at the end [m]:    %.4f\n', sqrt(mean(sum(dx(1:2,end).^2, 1))));
  fprintf(fileID,'\tVertical RMSE [m]:    %.4f\n', sqrt(mean(sum(dx(3,:).^2, 1))));
  fprintf(fileID,'\tVertical error at the end [m]:    %.4f\n', sqrt(mean(sum(dx(3,end).^2, 1))));
  fprintf(fileID,'\tMean vel [m/s]:  [%.4f, %.4f, %.4f]''\n', mean(dv, 2));
  fprintf(fileID,'\tRMSE vel [m/s]:  %.4f\n', sqrt(mean(sum(dv.^2, 1))));
  fprintf(fileID,'\tMean ori [rad]:  %.4f\n', mean(dq, 2));
  fprintf(fileID,'\tRMSE ori [rad]:  %.4f\n', sqrt(mean(dq.^2)));
end

function [h_err, v_err] = calc_pos_err(gt, x)
  dx = x(1:3, :) - gt(1:3, :);
  h_err = vecnorm(dx(1:2, :));
  v_err = abs(dx(3, :));
end




function display_active_sensors(sensor_pos, r)
  figure;    
  % All sensors
      plot3(r(1,:), r(2,:), r(3,:),'s');
      xlabel('x')
      ylabel('y')
      zlabel('z')
      hold on;
      % Active sensors
      plot3(sensor_pos(1,:), sensor_pos(2,:), sensor_pos(3,:),'s','Color','r');
      xlabel('x')
      ylabel('y')
      zlabel('z')
      ylim([-0.2 0.2])
  end
  
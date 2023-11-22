function   Xestimated2g = EKF_only( ...
    fs,... Sampling rate
    ab,... acc body frame
    w,... gyro body frame
    B,... Mag field, body frame
    sensor_pos,... Position of sensors
    x0,... Initial guess
    gt,...  Position aiding
    startT ... Start time
    )

  T=1/fs;
  gravity=9.82;

  %%

  B1= B(:, :, 1);
  B2= B(:, :, 2);
  B3= B(:, :, 3);
  B4= B(:, :, 4);
  B5= B(:, :, 5);



  wx=w(1, :);
  wy=w(2, :);
  wz=w(3, :);

  %% Time

  N=length(B1);
  t=0:T:(N-1)*T;

  %% Gradient calculation

  pos1 = sensor_pos(:, 1);
  pos2 = sensor_pos(:, 2);
  pos3 = sensor_pos(:, 3);
  pos4 = sensor_pos(:, 4);
  pos5 = sensor_pos(:, 5);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  matrix_B=eye(3);
  matrix_G1=[0     pos1(3) pos1(2)  pos1(1)   0;
    pos1(3)   0     pos1(1)    0     pos1(2);
    pos1(2)   pos1(1)    0    -pos1(3)     -pos1(3)];
  matrix_T1=[ (pos1(1)^2 -pos1(3)^2)   2*(pos1(1)*pos1(2))  0  (pos1(2)^2 -pos1(3)^2)  2*(pos1(1)*pos1(3))  0  2*(pos1(2)*pos1(3));
    0  (pos1(1)^2 -pos1(3)^2)  (pos1(2)^2 -pos1(3)^2)   2*(pos1(1)*pos1(2))  0   2*(pos1(2)*pos1(3))  2*(pos1(1)*pos1(3));
    -2*pos1(1)*pos1(3)   -2*pos1(2)*pos1(3)   -2*pos1(2)*pos1(3)  -2*pos1(1)*pos1(3) (pos1(1)^2 -pos1(3)^2)  (pos1(2)^2 -pos1(3)^2) 2*(pos1(1)*pos1(2)) ];
  P1=[matrix_B matrix_G1 0.5*matrix_T1];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  matrix_G2=[0     pos2(3) pos2(2)  pos2(1)   0 ;
    pos2(3)   0     pos2(1)    0     pos2(2);
    pos2(2)   pos2(1)    0    -pos2(3)     -pos2(3)];


  matrix_T2=[ (pos2(1)^2 -pos2(3)^2)   2*(pos2(1)*pos2(2))  0  (pos2(2)^2 -pos2(3)^2)  2*(pos2(1)*pos2(3))  0  2*(pos2(2)*pos2(3));
    0  (pos2(1)^2 -pos2(3)^2)  (pos2(2)^2 -pos2(3)^2)   2*(pos2(1)*pos2(2))  0   2*(pos2(2)*pos2(3))  2*(pos2(1)*pos2(3));
    -2*pos2(1)*pos2(3)   -2*pos2(2)*pos2(3)   -2*pos2(2)*pos2(3)  -2*pos2(1)*pos2(3) (pos2(1)^2 -pos2(3)^2)  (pos2(2)^2 -pos2(3)^2) 2*(pos2(1)*pos2(2))];
  P2=[matrix_B matrix_G2 0.5*matrix_T2];


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  matrix_G3=[0     pos3(3) pos3(2)  pos3(1)   0;
    pos3(3)   0     pos3(1)    0     pos3(2);
    pos3(2)   pos3(1)    0    -pos3(3)     -pos3(3)];

  matrix_T3=[(pos3(1)^2 -pos3(3)^2)   2*(pos3(1)*pos3(2))  0  (pos3(2)^2 -pos3(3)^2)  2*(pos3(1)*pos3(3))  0  2*(pos3(2)*pos3(3));
    0  (pos3(1)^2 -pos3(3)^2)  (pos3(2)^2 -pos3(3)^2)   2*(pos3(1)*pos3(2))  0   2*(pos3(2)*pos3(3))  2*(pos3(1)*pos3(3));
    -2*pos3(1)*pos3(3)   -2*pos3(2)*pos3(3)   -2*pos3(2)*pos3(3)  -2*pos3(1)*pos3(3) (pos3(1)^2 -pos3(3)^2)  (pos3(2)^2 -pos3(3)^2) 2*(pos3(1)*pos3(2))];
  P3=[matrix_B matrix_G3 0.5*matrix_T3];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  matrix_G4= [0     pos4(3) pos4(2)  pos4(1)   0 ;
    pos4(3)   0     pos4(1)    0     pos4(2);
    pos4(2)   pos4(1)    0    -pos4(3)     -pos4(3)];


  matrix_T4=[(pos4(1)^2 -pos4(3)^2)   2*(pos4(1)*pos4(2))  0  (pos4(2)^2 -pos4(3)^2)  2*(pos4(1)*pos4(3))  0  2*(pos4(2)*pos4(3));
    0  (pos4(1)^2 -pos4(3)^2)  (pos4(2)^2 -pos4(3)^2)   2*(pos4(1)*pos4(2))  0   2*(pos4(2)*pos4(3))  2*(pos4(1)*pos4(3));
    -2*pos4(1)*pos4(3)   -2*pos4(2)*pos4(3)   -2*pos4(2)*pos4(3)  -2*pos4(1)*pos4(3) (pos4(1)^2 -pos4(3)^2)  (pos4(2)^2 -pos4(3)^2) 2*(pos4(1)*pos4(2))];
  P4=[matrix_B matrix_G4 0.5*matrix_T4];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% P5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  matrix_G5=zeros(3,5);

  matrix_T5=zeros(3,7);

  P5=[matrix_B matrix_G5 0.5*matrix_T5]; % should be [I 0]


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%% Gradient & Tensor Calculation %%%%%%%%%%%%%%%%%%%%%%%%

  P=[P1;P2;P3;P4;P5];
  Pcroix=pinv(P);

  for k= 1:N
    vf(:,k)=Pcroix*[B1(:,k);B2(:,k);B3(:,k);B4(:,k);B5(:,k)];
  end
  gradient_B = zeros(3,3,N);
  TensorB = zeros(3,9,N);
  for k= 1:N
    gradient_B(:,:,k)=[vf(7,k) vf(6,k) vf(5,k);
      vf(6,k) vf(8,k) vf(4,k);
      vf(5,k) vf(4,k) -vf(7,k)-vf(8,k)];

    TensorB(:,:,k)=[vf(9,k) vf(10,k) vf(13,k) vf(10,k) vf(12,k) vf(15,k) vf(13,k) vf(15,k) -vf(9,k)-vf(12,k);
      vf(10,k) vf(12,k) vf(15,k) vf(12,k) vf(11,k) vf(14,k) vf(15,k) vf(14,k) -vf(10,k)-vf(11,k);
      vf(13,k) vf(15,k) -vf(9,k)-vf(12,k) vf(15,k) vf(14,k) -vf(10,k)-vf(11,k) -vf(9,k)-vf(12,k) -vf(10,k)-vf(11,k) -vf(13,k)-vf(14,k)];
    %[vect_propre(:,:,k),valeur_propre(:,:,k)]=eig(gradient_B(:,:,k));
  end
  TensorB1=TensorB(:,1:3,:);
  TensorB2=TensorB(:,4:6,:);
  TensorB3=TensorB(:,7:9,:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%% Indexing data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % TMI

  gradient_B=gradient_B(:,:,1:N);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% Kalman + grad  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Previous state (initial guegg)
  X_chap2g = x0;

  % P is the estimation covariance
  P2g =eye(18);
  P2g(1,1) = 0; P2g(2,2) = 0; P2g(3,3) = 0; P2g(4,4)=0;
  P2g(5,5) = 0; P2g(6,6) =0; P2g(7,7)=0;
  P2g(8,8) = 10; P2g(9,9) = 10; P2g(10,10) = 10;
  P2g(16,16) = 0; P2g(17,17) =0; P2g(18,18)=0;

  % Q is the process noise covariance
  Q2g = zeros(18);
  Q2g(5,5)=0.1;Q2g(6,6)=0.1;Q2g(7,7)=0.1;
  %% ceinture & cheville
  Q2g(8,8)=0.001;Q2g(9,9)=0.001;Q2g(10,10)=0.001;
  %% poche & dos
  Q2g(11,11)=100;Q2g(12,12)=100;Q2g(13,13)=100; Q2g(14,14)=100;Q2g(15,15)=100;


  % R is the measurement noise covariance
  Rmeas2g =5*eye(8);
  Rmeas2pos = 0.01 * eye(3);

  Xestimated2g = zeros(size(x0, 1), length(w));
  for k=1:length(w)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% The state function %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Defining the state function
    f1g=@(t,x)[
      (0.5 * quaternProd([x(1) x(2) x(3) x(4)], [0 wx(k) wy(k) wz(k)]))';
      cross(-w(:,k),[x(5) x(6) x(7)])'+ab(:,k)-quatern2rotMat([x(1) x(2) x(3) x(4)])*[0;0;gravity];
      (cross(-w(:,k),[x(8) x(9) x(10)]))'+[x(11) x(12) x(13); x(12) x(14) x(15); x(13) x(15) -x(11)-x(14)]*[x(5) x(6) x(7)]';
      TensorB1(1,:,k)*x(5:7)+[x(11) x(12) x(13)]*[0; wz(k);-wy(k)]-[0 -wz(k) wy(k)]*[x(11); x(12); x(13)];
      TensorB2(1,:,k)*x(5:7)+[x(11) x(12) x(13)]*[-wz(k); 0; wx(k)]-[0 -wz(k) wy(k)]*[x(12); x(14); x(15)];
      TensorB3(1,:,k)*x(5:7)+[x(11) x(12) x(13)]*[wy(k); -wx(k); 0]-[0 -wz(k) wy(k)]*[x(13); x(15) ;-x(11)-x(14)];
      TensorB2(2,:,k)*x(5:7)+[x(12) x(14) x(15)]*[- wz(k); 0; wx(k)]-[ wz(k) 0 -wx(k)]*[x(12) ;x(14); x(15)];
      TensorB3(2,:,k)*x(5:7)+[x(12) x(14) x(15)]*[wy(k); -wx(k); 0]-[wz(k) 0 -wx(k)]*[x(13); x(15) ;-x(11)-x(14)]
      quatern2rotMat([x(1) x(2) x(3) x(4)])\x(5:7);
      ];
    %Discretization
    fgg=@(x)x + (1/6)*(f1g(t(k),x)+2*f1g(t(k)+0.5*T,x+0.5*T*f1g(t(k),x))+2*f1g((t(k)+0.5*T),(x+0.5*T*f1g(t(k)+0.5*T,x+0.5*T*f1g(t(k),x))))+f1g((t(k)+T),(x+f1g((t(k)+0.5*T),(x+0.5*T*f1g(t(k)+0.5*T,x+0.5*T*f1g(t(k),x))))*T)))*T;  % main equation;
    %Defining the Measurement function
    h2g=@(x)x(8:15);
    %Generating measurmentg
    Z2g = [B5(:,k);squeeze(gradient_B(1,1,k));squeeze(gradient_B(2,1,k));squeeze(gradient_B(3,1,k));squeeze(gradient_B(2,2,k));squeeze(gradient_B(3,2,k))];
    %Calculating the jacobian of the state function
    Xestimated2g(:,k)=X_chap2g;
    F2g=numjacobian(fgg,X_chap2g);
    %Applying the equations of the filter
    X_chap2g=fgg(X_chap2g);
    P2g=F2g*P2g*F2g'+Q2g;
    if k <= startT
        H2g  = [numjacobian(h2g,X_chap2g); [zeros(3, 15) eye(3)]];
        Y_tilt2g=[Z2g-h2g(X_chap2g); gt.pos(:, k) - X_chap2g(end-2:end,:)];
        S2g=H2g*P2g*H2g'+blkdiag(Rmeas2g, Rmeas2pos);
    else
        H2g  = numjacobian(h2g,X_chap2g);
        Y_tilt2g=Z2g-h2g(X_chap2g);
        S2g=H2g*P2g*H2g'+Rmeas2g;
    end
    K2g=P2g*H2g'/S2g;
    X_chap2g=X_chap2g+K2g*Y_tilt2g;
    P2g=(eye(18)-K2g*H2g)*P2g;
  end


end

function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2022
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Enrico Mion, Bhavya Sukhija, Jin Cheng
% enmion@ethz.ch
% sukhijab@ethz.ch
% jicheng@student.ethz.ch
%
% --
% Revision history
% [24.04.18, MH]    first version by Matthias & Carlo
% [07.05.19, CS]    updated for 2019
% [29.03.21, TB]    updated for 2021
% [21.04.22 Noah Baumann] updated for 2022

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % Replace the following:
    posEst = [0,0]; % initial pos means are zero
    linVelEst = [0,0]; % initial vel means are zero
    oriEst = 0; % initial rot mean is zero
    windEst = 0; % initial wind direction is zero
    driftEst = 0; % initial drift is zero
    
    % initial position variance is calculated analytically
    posVar = [1/4*estConst.StartRadiusBound^2,1/4*estConst.StartRadiusBound^2];
    linVelVar = [0,0]; % initial velocity variance is zero
    % initial orientation variance is calculated analytically
    oriVar = 1/3*estConst.RotationStartBound^2;
    windVar = 1/3*estConst.WindAngleStartBound^2;
    driftVar = 0; % initial drift variance is zero
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar, linVelVar, oriVar, driftVar, windVar]);
    % estimator state (pos, rot, drift)
    estState.xm = [posEst, linVelEst, oriEst, driftEst, windEst]';
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.
% x = [p_x p_y s_x s_y phi b rho]' in R^7
% v = [v_d, v_r, v_b v_rho]' in R^4
% w = [w_a, w_b, w_c, w_g, w_n]' in R^5

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% process noise Q
Q = diag([estConst.DragNoise, estConst.RudderNoise, estConst.GyroDriftNoise, estConst.WindAngleNoise]);

%% prior update
% use ODE to push mean and variance forward:
[~,xsim] = ode45(@dynamics,[0 dt], [estState.xm; estState.Pm(:)]);

xp = xsim(end,1:7)';
Pp = reshape(xsim(end,8:end)',7,7);

%% measurement update
% sensor constants:
x_a = estConst.pos_radioA(1);
y_a = estConst.pos_radioA(2);
x_b = estConst.pos_radioB(1);
y_b = estConst.pos_radioB(2);
x_c = estConst.pos_radioC(1);
y_c = estConst.pos_radioC(2);

H = [];
M = [];
R = [];
z = [];

% Build measurement update matrices (assume distance C meas. is available)
% z_a:
H = [(xp(1)-x_a)/((xp(1)-x_a)^2+(xp(2)-y_a)^2)^(1/2),...
     (xp(2)-y_a)/((xp(1)-x_a)^2+(xp(2)-y_a)^2)^(1/2) 0 0 0 0 0];
% z_b:    
H = [H;
    (xp(1)-x_b)/((xp(1)-x_b)^2+(xp(2)-y_b)^2)^(1/2),...
    (xp(2)-y_b)/((xp(1)-x_b)^2+(xp(2)-y_b)^2)^(1/2) 0 0 0 0 0]; 
% z_c:
H = [H;
    (xp(1)-x_c)/((xp(1)-x_c)^2+(xp(2)-y_c)^2)^(1/2),...
    (xp(2)-y_c)/((xp(1)-x_c)^2+(xp(2)-y_c)^2)^(1/2) 0 0 0 0 0];
% z_g:
H = [H;
     0 0 0 0 1 1 0]; 
% z_n:
H = [H;
     0 0 0 0 1 0 0]; 
% M:
M = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
% R:
R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);
% z:
z = sense(1:5)';

% Nonlinear measurement prediction
h_k = [((xp(1)-x_a)^2+(xp(2)-y_a)^2)^(1/2);
       ((xp(1)-x_b)^2+(xp(2)-y_b)^2)^(1/2);
       ((xp(1)-x_c)^2+(xp(2)-y_c)^2)^(1/2);
       xp(5)+xp(6);
       xp(5)];

% If distance A,B & C meas. are not available, reduce H, M, R, z, h_k matrices
if(isinf(sense(3)))
    H = H([4,5],:);
    M = M(4:5,4:5); 
    R = R(4:5,4:5);   
    z = z([4,5]);
    h_k = h_k([4,5]);
end



if(~isempty(z))
    % do update only if there was a measurement, otherwise, just keep 
    % prior update
    K = Pp*H'/(H*Pp*H' + M*R*M');
    estState.xm = xp + K*(z - h_k);
    % Joseph form for EKF
    estState.Pm = (eye(size(Pp))-K*H)*Pp*(eye(size(Pp))-K*H)'+K*M*R*M'*K';
else
    % if there was no meas. update, use prior update data
    estState.Pm = Pp;
    estState.xm = xp;
end

% Get resulting estiamtes and variances
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
driftEst = estState.xm(6);
windEst = estState.xm(7);

posVar = [estState.Pm(1,1), estState.Pm(2,2)];
linVelVar = [estState.Pm(3,3), estState.Pm(4,4)];
oriVar = estState.Pm(5,5);
driftVar = estState.Pm(6,6);
windVar = estState.Pm(7, 7);
% ----------------------------------------------------------------------- %
% helper
    function [dx] = dynamics(t,x)       
        % helper vars:
        thrust = tanh(actuate(1));
        c_dh = estConst.dragCoefficientHydr;
        c_da = estConst.dragCoefficientAir;
        c_r = estConst.rudderCoefficient;
        c_w = estConst.windVel;
       
        % state ODE
        dx = zeros(56,1);
        dx(1) = x(3);                                     % dp_x/dt
        dx(2) = x(4);                                     % dp_y/dt
        dx(3) = cos(x(5))*(thrust-c_dh*(x(3)^2+x(4)^2))-c_da*(x(3)-c_w*cos(x(7)))*sqrt((x(3)-c_w*cos(x(7)))^2+(x(4)-c_w*sin(x(7)))^2);   % ds_x/dt
        dx(4) = sin(x(5))*(thrust-c_dh*(x(3)^2+x(4)^2))-c_da*(x(4)-c_w*sin(x(7)))*sqrt((x(3)-c_w*cos(x(7)))^2+(x(4)-c_w*sin(x(7)))^2);   % ds_y/dt
        dx(5) = c_r*actuate(2);                           % dphi/dt
        dx(6) = 0;                                        % db/dt
        dx(7) = 0;                                        % drho/dt
        
        % linearize for variance update
        A = ...
            [ 0, 0,                                                                                                                                                                                                       1,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0; ...
              0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       1,                                            0, 0,                                                                                                                                                                                                 0; ...
              0, 0, - c_da*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2) - 2*c_dh*x(3)*cos(x(5)) - (c_da*(x(3) - c_w*cos(x(7)))*(2*x(3) - 2*c_w*cos(x(7))))/(2*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2)),                                                                - 2*c_dh*x(4)*cos(x(5)) - (c_da*(x(3) - c_w*cos(x(7)))*(2*x(4) - 2*c_w*sin(x(7))))/(2*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2)), -sin(x(5))*(thrust - c_dh*(x(3)^2 + x(4)^2)), 0, (c_da*c_w*(x(3) - c_w*cos(x(7)))*(x(4)*cos(x(7)) - x(3)*sin(x(7))))/((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2) - c_da*c_w*sin(x(7))*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2); ...
              0, 0,                                                                - 2*c_dh*x(3)*sin(x(5)) - (c_da*(x(4) - c_w*sin(x(7)))*(2*x(3) - 2*c_w*cos(x(7))))/(2*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2)), - c_da*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2) - 2*c_dh*x(4)*sin(x(5)) - (c_da*(x(4) - c_w*sin(x(7)))*(2*x(4) - 2*c_w*sin(x(7))))/(2*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2)),  cos(x(5))*(thrust - c_dh*(x(3)^2 + x(4)^2)), 0, c_da*c_w*cos(x(7))*((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2) + (c_da*c_w*(x(4) - c_w*sin(x(7)))*(x(4)*cos(x(7)) - x(3)*sin(x(7))))/((x(3) - c_w*cos(x(7)))^2 + (x(4) - c_w*sin(x(7)))^2)^(1/2); ...
              0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0; ...
              0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0; ...
              0, 0,                                                                                                                                                                                                       0,                                                                                                                                                                                                       0,                                            0, 0,                                                                                                                                                                                                 0];

        L = [0, 0, 0, 0;
             0, 0, 0, 0;
             -cos(x(5))*c_dh*(x(3)^2+x(4)^2), 0, 0, 0;
             -sin(x(5))*c_dh*(x(3)^2+x(4)^2), 0, 0, 0;
             0, dx(5), 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1];
        
        P = reshape(x(8:end),7,7);
        
        % variance ODE
        dP = A*P + P*A' + L*Q*L';
        % combine
        dx(8:end,1) = dP(:);
    end
end
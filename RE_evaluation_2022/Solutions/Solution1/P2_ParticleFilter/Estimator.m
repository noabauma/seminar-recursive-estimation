function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==1, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index, scalar
%                       corresponds to continous time t = k*Ts
%                       If tm==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
%
% --
% Revision history
% [21.04.19, CS]    first version, template

% Set number of particles:
N_particles = 4096; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    
    % ################################################################### %
    % These particles are the posterior particles at discrete time k = 0
    % which will be fed into your estimator again at k = 1
    % Replace the following:
   
    corners = [estConst.pA;estConst.pB];
    [postParticles.x_r, postParticles.y_r] = UniformCircle(estConst.d,corners,N_particles);
    postParticles.phi = UniformN(estConst.phi_0,N_particles);
    postParticles.rho = UniformN(estConst.m,N_particles);
    postParticles.kappa = UniformN(estConst.l,N_particles);

    % ################################################################### %
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!
% ####################################################################### %

% Copy particles into single array:
state = [prevPostParticles.x_r;
         prevPostParticles.y_r;
         prevPostParticles.phi;
         prevPostParticles.rho;
         prevPostParticles.kappa];

%% Prior Update:
% Push particles through system dynamics:
for p = 1:N_particles
    u0 = act(1) + Uniform(estConst.sigma_f/2);
    state(1:2, p) = propagatePosition(state(1:2,p)', state(3,p),u0)';
    state(3, p) = myWrapToPi(state(3,p) + act(2) + Uniform(estConst.sigma_phi/2));
end

%% Posterior Update:
% If there is a measurement, calculate the weights:
betas = NaN(N_particles,1);

% Calculate expected measurement from prior particles
for p = 1:N_particles
    contour_ = estConst.contour;
    contour_(1, 2) = state(4,p);
    contour_(2, 2) = state(4,p);
    contour_(8, 1) = state(5,p);
    contour_(9, 1) = state(5,p);
    dist2Wall_p = getDist2Wall(state(1,p), state(2,p), state(3,p), contour_);
    betas(p) = triTriangular(sens - dist2Wall_p, estConst);
end

% Redraw robot samples:
cdf = cumsum(betas);
stateTemp = zeros(5,N_particles);
if(cdf(end) == 0)
    [stateTemp] = reInitEstimator(N_particles, estConst);
else
    for d = 1:N_particles
        ind = find(cdf > rand*cdf(end),1,'first');
        stateTemp(:,d) = state(:,ind);
    end
end
state = stateTemp;

% Add roughening:
Kr = 0.05;
D = 4;
% Calculate dmax, maximum inter-sample variability for each dimension:
Emax = zeros(3,1);
Emax(1) = max(state(1,:)) - min(state(1,:));
Emax(2) = max(state(2,:)) - min(state(2,:));
Emax(3) = maxVariabilityAngles(state(3,:));
Emax(4) = max(state(4,:)) - min(state(4,:));
Emax(5) = max(state(5,:)) - min(state(5,:));

for dim = 1:5
    state(dim,:) = state(dim,:) + Kr * Emax(dim) * N_particles^(-1/D) * randn(1, N_particles);
end

% Ensure all particles are in valid range:
contour_outside = estConst.contour;
contour_outside(1, 2) = -estConst.m;
contour_outside(2, 2) = -estConst.m;
contour_outside(8, 1) = -estConst.l;
contour_outside(9, 1) = -estConst.l;
[inside,~] = inpolygon(state(1,:), state(2,:), contour_outside(:,1), ...
    contour_outside(:,2)); 
for idx = 1:N_particles
    if(~inside(idx))
        [state(1,idx), state(2,idx)] = ...
            findClosestPointOnContour(state(1,idx),state(2,idx),...
            contour_outside);
    end
end
state(3,:) = myWrapToPi(state(3,:));

postParticles.x_r = state(1,:);
postParticles.y_r = state(2,:);
postParticles.phi = state(3,:);
postParticles.rho = state(4,:);
postParticles.kappa = state(5,:);
    
end % end estimator

% ----------------------------------------------------------------------- %
% --------------------------helper functions----------------------------- %
% ----------------------------------------------------------------------- % 

% find closest point on contour
function [xC, yC] = findClosestPointOnContour(x,y,contour)    
    dMin = inf;
    for i = 1:size(contour,1)
        if i < size(contour,1)
            dist = lineToSegmentDistance(contour(i,1),contour(i+1,1),x,contour(i,2),contour(i+1,2),y);
        else
            dist = lineToSegmentDistance(contour(i,1),contour(1,1),x,contour(i,2),contour(1,2),y);
        end
        if dist < dMin
            contIdx = i;
            dMin = dist;
        end
    end
    % find coordinates onclosest contour
    if contIdx < size(contour,1) 
        [~, xC, yC] = lineToSegmentDistance(contour(contIdx,1),...
            contour(contIdx+1,1), x, contour(contIdx,2), contour(contIdx+1,2), y);
    else
        [~, xC, yC] = lineToSegmentDistance(contour(contIdx,1),...
            contour(1,1), x, contour(contIdx,2), contour(1,2), y);
    end   
end

% ------------------------------------------------------------------- %
function [state] = reInitEstimator(N_particles, estConst)   
   state = [rand(1,N_particles)*2;
            rand(1,N_particles)*2;
            -pi + 2*pi*rand(1,N_particles);
            UniformN(estConst.m,N_particles);
            UniformN(estConst.l,N_particles)];            
end

% ------------------------------------------------------------------- %
function E = maxVariabilityAngles(angles)
    angles = myWrapToPi(angles);
    positive_angles = angles(angles>=0);
    negative_angles = angles(angles<0);
    
    diff_positive = max(positive_angles) - min(positive_angles);
    diff_negative = max(negative_angles) - min(negative_angles);
    
    mixed = sortrows([negative_angles,positive_angles-pi;
                  zeros(size(negative_angles)),...
                  ones(size(positive_angles))]')';
          
    temp = diff(mixed,[],2);
    diff_mixed = min(abs(temp(1,temp(2,:)~=0)));       
    E = max([diff_positive,diff_negative,diff_mixed]);     
end

% ------------------------------------------------------------------- %
function lHood = triTriangular(xBar, estConst)
    eps = estConst.epsilon;
    if((0 <= abs(xBar)) && (abs(xBar) <= 2*eps))
        lHood = 2/(5*eps) - 1/(5*eps^2)*abs(xBar);            
    elseif((2*eps < abs(xBar)) && (abs(xBar) <= 2.5*eps))
        lHood = 2/(5*eps^2)*(abs(xBar)-2*eps);
    elseif((2.5*eps < abs(xBar)) && (abs(xBar)<= 3*eps))
        lHood = 1/(5*eps) - 2/(5*eps^2)*(abs(xBar)-5*eps/2); 
    else
        lHood = 0;
    end
end

% ------------------------------------------------------------------- %
function rnext = propagatePosition(r,phi,u)
    rnext = r + u*[cos(phi),sin(phi)];
end

% ------------------------------------------------------------------- %
function distance = getDist2Wall(xR,yR,phi,contour)
    Nwalls = size(contour,1);
    wallsX = [contour(:,1),[contour(2:end,1);contour(1,1)]];
    wallsY = [contour(:,2),[contour(2:end,2);contour(1,2)]];
    distance2walls = inf*ones(Nwalls,1);

    for i = 1:Nwalls
        distance2walls(i) = computeRayLineIntersection(wallsX(i,1),...
            wallsX(i,2),xR,wallsY(i,1),wallsY(i,2),yR,phi);
    end

    distance = min(distance2walls);
end

% ------------------------------------------------------------------- %
function distance = computeRayLineIntersection(x1, x2, xR, y1, y2, yR, phiR)
    w = [x1-xR; y1-yR];
    v = [cos(phiR); sin(phiR)];
    r = [x2-x1; y2-y1];
    % check if ray and line are parallel
    if (abs(v(1) * r(2) - v(2) * r(1)) < 1e-10)
        distance = inf;
        return;
    end
    % t is the ray parameter; s is the line parameter
    t = (r(1) * w(2) - r(2) * w(1)) / (r(1) * v(2) - r(2) * v(1));
    s = (v(2) * w(1) - v(1) * w(2)) / (v(1) * r(2) - v(2) * r(1));
    % check if intersection is within line and on positive side of ray
    if t < 0 || s < 0 || s > 1
        distance = inf;
    else
        distance = t;
    end
end

% ------------------------------------------------------------------- %
function [dist,x,y] = lineToSegmentDistance(x1,x2,x3,y1,y2,y3)
    % (x1,y1) and (x2,y2) define the line, (x3,y3) is the point
    px = x2-x1;
    py = y2-y1;
    norm = px*px + py*py;
    u =  ((x3 - x1) * px + (y3 - y1) * py) / norm;
    if u > 1
        u = 1;
    elseif u < 0
        u = 0;
    end
    x = x1 + u * px;
    y = y1 + u * py;
    dx = x - x3;
    dy = y - y3;
    dist = sqrt(dx * dx + dy * dy); % sqrt is not necessary for comparisons
end
% ------------------------------------------------------------------- %

function lambda = myWrapToPi(lambda)
    q = (lambda < -pi) | (pi < lambda);
    lambda(q) = mod(lambda(q) + pi, 2*pi) - pi;
end

% ------------------------------------------------------------------- %
function [x, y] = UniformCircle(l, corners, N)
    % corners rows: [x0,y0]
    % number of rows: 2
    samples = rand(3,N); % columns 1: corner, 2: radius, 3: angle
    r = sqrt(samples(2,:)) * l;
    angle = samples(3,:) * 2 * pi;
    corner1 = samples(1,:) > 0.5;
    x = corner1 * corners(1,1) + (1-corner1) * corners(2,1) + r .* cos(angle);
    y = corner1 * corners(1,2) + (1-corner1) * corners(2,2) + r .* sin(angle);
end

% ----------------------------------------------------------------------- %
function val = Uniform(mx)
    val = -mx + 2 * mx * rand;
end

% ----------------------------------------------------------------------- %
function val = UniformN(mx,N)
    val = -mx + 2*mx*rand(1,N);
end

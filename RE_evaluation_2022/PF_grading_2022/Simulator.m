function [km, state, input, sense] = Simulator( simConst )
%% Simulation: Generate data
% Simulation of the system to generate: inputs, measurements, true states.
% Inputs and measurements will be visible to the estimator.
%
simConst.alpha = 0.01; % Constant that defines the magnitude of the input
simConst.d_safety = 0.1;  % Safety distance from the wall

%% Initialize
% Number of samples of the simulation
N = simConst.N;

% The robot position
r = zeros(N,2); %[x_r,y_r]

% The robot orientation
phi = zeros(N,1);

% Initialize uncertain points in contour
rho = Uniform(simConst.m) * ones(N, 1);
simConst.contour(1, 2) = simConst.contour(1, 2) + rho(1);
simConst.contour(2, 2) = simConst.contour(2, 2) + rho(1);

kappa = Uniform(simConst.l) * ones(N, 1);
simConst.contour(8, 1) = simConst.contour(8, 1) + kappa(1);
simConst.contour(9, 1) = simConst.contour(9, 1) + kappa(1);

% The initial position and orientation
corners = [simConst.pA; simConst.pB];
[x0,y0] = UniformCircle(simConst.d, corners, 1);
r(1,:) = [x0, y0];
phi(1,1) = Uniform(simConst.phi_0);

% The input
u = zeros(N-1,2); %[u_f,u_phi]

% The measurements
distance = zeros(N-1,1);

% Sign of last input
positive = true;

for n = 1:N-1
    % Store control input (known to estimator)
    u(n,1) = computeInput(r(n,:), phi(n,1), positive, simConst);
    u(n,2) = 0.01;
    
    % Update the last input sign
    positive = (u(n,1) > 0);
    
    % Add noise to the input
    u0 = u(n,1) + Uniform(simConst.sigma_f/2);
    
    % Update position and orientation
    r(n+1,:) = propagatePosition(r(n,:), phi(n,1), u0);
    phi(n+1,1) = wrapToPi(phi(n,1) + u(n,2) + Uniform(simConst.sigma_phi/2));
   
    % Store distance measurement    
    distance(n) = getDistMeas(r(n+1,1), r(n+1,2), phi(n+1,1), simConst.contour, simConst.epsilon);
end

km = 0:simConst.N;
state = [r, phi, rho, kappa];
input = u;
sense = [Inf; distance];

end


% ----------------------------------------------------------------------- %
% ------------------------component functions---------------------------- %
% ----------------------------------------------------------------------- % 

function rnext = propagatePosition(r, phi, u)
    rnext = r + u * [cos(phi), sin(phi)];
end

% ------------------------------------------------------------------- %
function u = computeInput(r, phi, positive, simConst)
    u_plus = UniformMinMax(0, simConst.alpha);
    u_minus = - u_plus;
    if positive
        r_plus = propagatePosition(r, phi, u_plus + simConst.sigma_f/2);
        d_plus = computeMinDistance(r_plus, simConst.contour);
        u = u_plus;
        if d_plus < simConst.d_safety
            r_minus = propagatePosition(r, phi, u_minus - simConst.sigma_f/2);
            d_minus = computeMinDistance(r_minus, simConst.contour);
            if d_minus > d_plus
               u = u_minus;
            elseif d_minus == d_plus
               u = 0;
            end
        end
    else
        r_minus = propagatePosition(r, phi, u_minus - simConst.sigma_f/2);
        d_minus = computeMinDistance(r_minus, simConst.contour);
        u = u_minus;
        if d_minus < simConst.d_safety
            r_plus = propagatePosition(r, phi, u_plus + simConst.sigma_f/2);
            d_plus = computeMinDistance(r_plus, simConst.contour);
            if d_plus > d_minus
               u = u_plus;
            elseif d_minus == d_plus
               u = 0;
            end
        end
    end
end

% ----------------------------------------------------------------------- %
function distance = getDistMeas(xR, yR, phiR, contour, epsilon)
    % Compute distance measurement based on current state
    Nwalls = size(contour,1);
    wallXpairs = [contour(:,1), [contour(2:end, 1); contour(1,1)]];
    wallYpairs = [contour(:,2), [contour(2:end, 2); contour(1,2)]];
    distance2walls = inf * ones(Nwalls,1);
    for idx = 1:Nwalls
        distance2walls(idx) = computeRayLineIntersection( ...
            wallXpairs(idx,1), wallXpairs(idx,2), xR, ...
            wallYpairs(idx,1), wallYpairs(idx,2), yR, ...
            phiR);
    end
    distance = min(distance2walls);
    % Add noise
    uSamp = rand(1);
    if (uSamp <= 1/10)
        w_bar = inverseCumTrian(-3 * epsilon, -2 * epsilon, 1 / (5 * epsilon), uSamp);
    elseif (1/10 <= uSamp) && (uSamp <= 9/10)
        uSamp = uSamp-1/10;
        w_bar = inverseCumTrian(-2 * epsilon, 2 * epsilon, 2 / (5 * epsilon), uSamp);
    elseif (9/10 <= uSamp)
        uSamp = uSamp-9/10;
        w_bar = inverseCumTrian(2 * epsilon, 3 * epsilon, 1 / (5 * epsilon), uSamp);
    end
    distance = distance + w_bar;
end


% ----------------------------------------------------------------------- %
% --------------------------helper functions----------------------------- %
% ----------------------------------------------------------------------- % 

function x = inverseCumTrian(a,b,c,y)
    if (0<=y) && (y <= c*(b-a)/4)
        x = a+sqrt(y*(b-a)/c);
    elseif (c*(b-a)/4 <= y) && (y <= c*(b-a)/2)
        x = b-sqrt((b-a)^2/2-(b-a)/c*y);
    end
end

% ----------------------------------------------------------------------- %
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

% ----------------------------------------------------------------------- %
function d = computeMinDistance(r, contour)
    d = inf;
    for i = 1:size(contour,1)
        if i < size(contour,1)
            dist = lineToSegmentDistance(contour(i,1), contour(i+1,1), r(1), contour(i,2), contour(i+1,2), r(2));
        else
            dist = lineToSegmentDistance(contour(i,1), contour(1,1), r(1), contour(i,2), contour(1,2), r(2));
        end
        if dist < d
            d = dist;
        end
    end
end

% ----------------------------------------------------------------------- %
function dist = lineToSegmentDistance(x1,x2,x3,y1,y2,y3)
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

% ----------------------------------------------------------------------- %
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
function val = UniformMinMax(mn, mx)
    val = mn + (mx - mn) * rand;
end

function [t, r, v, stats] = Encke(a_disturb, tspan, rt0, vt0, muC, options) 
% Following Battin, pages 448-449

    %% Initialize

    % Check I/O argument
    error(nargchk(5,6,nargin,'struct'));
    error(nargoutchk(0,4,nargout,'struct'));    
    assert(isa(a_disturb, 'function_handle'), ...
        'Encke:ioarg_error',...
        'Argument ''a'' must be a valid function handle.');
    assert(isnumeric(tspan) && isvector(tspan) && numel(tspan) >= 2 && ...
           all(isfinite(tspan)) && all(isreal(tspan)) && ...
           issorted(tspan) && numel(tspan) == numel(unique(tspan)), ...
        'Encke:ioarg_error',...
        'The time span vector (''tspan'') must be given as real, sorted values.');
    assert(isnumeric(rt0) && isvector(rt0) && numel(rt0) == 3 && all(isreal(rt0)) && all(isfinite(rt0)) && ...
           isnumeric(vt0) && isvector(vt0) && numel(vt0) == 3 && all(isreal(vt0)) && all(isfinite(vt0)),...
        'Encke:ioarg_error',...
        'The initial position and velocity (''rt0'', ''vt0'') must both be real 3-element vectors.');
     assert(isnumeric(muC) && isscalar(muC) && muC > 0 && isreal(muC),...
        'Encke:ioarg_error',...
        ['The standard gravitational parameter of the central body (''muC'') must be givan ',...
        'as a real positive scalar.']);
    
    % Initialize constants 
    stats.rectifications = 0;
    stats.fevals = 1; % NOTE: one eval is done to check the user-provided function
    done = false;
    hinitial = (tspan(end)-tspan(1))/100;
    t = [];
    r = [];
    v = [];
    
    % Initialize options
    if nargin < 6 || isempty(options)
        options = odeset; end
    
    assert(isstruct(options),...
        'Encke:invalid_options',...
        ['Argument ''options'' must be a structure similar to the one generated by odeset(), ',...
        'optionally containing the fields ''handle'' and ''velocity_required''. Note that ',...
        'those fields might be removed when odeset() or odeget() is called after the fields ',...
        'have been set.']);    

    if ~isfield(options, 'handle')
        % Defaults to "best" integrator in standard MATLAB
        ode_options = options;
        solver      = @ode45;
        solverType  = 0;
    else
        % Custom handle
        ode_options = rmfield(options, {'handle' 'velocity_required'});
        solver      = options.handle;
        solverType  = any(strcmpi(options.velocity_required, 'no')) || ...
            (islogical(options.velocity_required) && options.velocity_required);
    end
    
    % Make sure the disturbing acceleration function works
    try
        if solverType==0
            [~] = a_disturb(tspan(1),rt0,vt0);
        else
            [~] = a_disturb(tspan(1),rt0);
        end
        
    catch ME
        ME2 = MException(...
            'Encke:a_disturb_failure',...
            'Failed to evaluate disturbing acceleration function for initial values.');
        throw(addCause(ME,ME2));
    end

    % Set/get defaults. Make sure the rectification check is carried
    % out at each integration step
    % NOTE: An output function is used to detect an event. Although an
    % EvenFcn would seem better suited, in this case, we're not 
    % interested in *exactly* where the rectification should take place
    % (standard behavior in almost all ODE solvers), but rather, whether
    % it makes sense to do it in the next step. Therefore, to save on
    % function evaluations, the terminal property of the OutputFunction is
    % used, which makes sure the condition is only checked *once* per
    % iteration. 
    solver_options = odeset(ode_options, ...
        'AbsTol'     , min(odeget(ode_options,'AbsTol',1e-8), 1e-8),...
        'RelTol'     , min(odeget(ode_options,'RelTol',1e-8), 1e-8),...
        'InitialStep', min(odeget(ode_options,'InitialStep',hinitial), hinitial),...
        'OutputFcn'  , @(varargin)...
        check_rectify(solverType, varargin{:}) || ...
        feval(odeget(ode_options, 'OutputFcn', @(varargin)false), varargin{:})...
    );


    % Main loop:
    % - Let solver integrate until rectification is needed
    % - Break off integration and collect data 
    % - Apply recitifcation
    % - Repeat from the start, until the final time is reached.
    while ~done
        t0 = tspan(1);
        
        % We can in all likelihood keep using the step size from the
        % previous integration
        solver_options = odeset(solver_options,...
            'InitialStep', hinitial(end));
        
        % The 2 sovler types have 2 different calls
        if solverType == 0            
            D = solver(...
                @(t,y) de1(t,y, a_disturb,t0,rt0,vt0,muC), ...
                tspan, ...
                zeros(6,1), solver_options);
            
            delta    = D.y(1:3,:).';
            deltadot = D.y(4:6,:).';
            
        else            
            D = solver(...
                @(t,y) de2(t,y, a_disturb,t0,rt0,vt0,muC), ...
                tspan, ...
                zeros(3,1),zeros(3,1), solver_options);
            
            delta    = D.y;
            deltadot = D.yp;
            
        end
        
        tslice = D.x(:);
        stats.fevals = stats.fevals + D.stats.nfevals;
        
        % Save values for next iteration
        hinitial = tslice(end)-tslice(end-1);
        tspan    = [tslice(end) tspan(end)];
        done     = tslice(end)==tspan(end);
        
        % Insert data into output arrays  
        [rosc, vosc] = check_rectify();
        t = [t;  tslice(1:end-~done)];                                   %#ok<AGROW>
        r = [r   [rt0 rosc(:,1:end-~done)] + delta(1:end-~done,:).'   ]; %#ok<AGROW>
        v = [v   [vt0 vosc(:,1:end-~done)] + deltadot(1:end-~done,:).']; %#ok<AGROW>
        
        % Apply rectification 
        rt0 = rosc(:,end) + delta(end,:).';
        vt0 = vosc(:,end) + deltadot(end,:).';
        stats.rectifications = stats.rectifications + ~done;
        
    end
    
end


% Little wrapper function; checks if rectification is needed
function varargout = check_rectify(solverType, varargin)
    
    persistent rosc vosc
    if nargin == 0
        varargout{1} = rosc;
        varargout{2} = vosc;
        return;
    end
    
    terminate = false;
    flag = varargin{end};
    if isempty(flag)
        if solverType == 0
            [terminate, rosc(:,end+1), vosc(:,end+1)] = de1();
        else
            [terminate, rosc(:,end+1), vosc(:,end+1)] = de2();
        end
              
    elseif strcmp(flag, 'init')
        rosc = [];
        vosc = [];
    end
    
    varargout{1} = terminate;
    
end


% Compute first derivative of 
%
%    y(t) = [
%       delta(t)     = r(t) - rosc(t)
%       delta_dot(t) = v(t) - vosc(t)
%    ]
% 
% This DE should be used if the disturbing acceleration depends on the
% velocity, or if the solver is designed to solve systems of the type
%
%       dy/dt = f(t,y)
% s.t.  y(t0) = y0
%
function varargout = de1(t,y, a_disturb,t0,rt0,vt0,muC)
    
    persistent a_main a_perturb rosc vosc
    if nargin==0
        % Condition to trigger rectification
        varargout{1} = norm(a_main) > 5*norm(a_perturb);
        % Osculating orbit data
        varargout{2} = rosc;
        varargout{3} = vosc;        
        return;
    end
    
    delta     = y(1:3);
    delta_dot = y(4:6);
      
    [rosc, vosc] = osculating_orbit(rt0,vt0, t-t0, muC);
    r  = delta + rosc; 
    v  = delta_dot + vosc;
    q  = delta'*(delta-2*r) / (r'*r); 
    fq = 1-(1+q)^(3/2);
    
    a_main    = muC/norm(rosc)^3 * (fq*r - delta);
    a_perturb = a_disturb(t,r,v);
    a_total   = a_main + a_perturb;
    
    % Generate output
    varargout{1} = [
        delta_dot  % d (delta)/dt 
        a_total    % d�(delta)/dt� 
    ];
    
end


% Compute second derivative of 
%
%    delta(t) = r(t) - rosc(t)
%
% This DE should be used if the disturbing acceleration is independt of 
% the velocity. If this is the case, the solver must ba designed to solve 
% systems of type
%
%     d�y/dt� = f(t,y)
% s.t.  y(t0) = y0
%       dy/dt = v0
%
function varargout = de2(t,delta, a_disturb,t0,rt0,vt0,muC)
       
    % Condition to trigger rectification
    persistent a_main a_perturb rosc vosc
    if nargin==0
       % Condition to trigger rectification
        varargout{1} = norm(a_main) > 5*norm(a_perturb);
        % Osculating orbit data
        varargout{2} = rosc;
        varargout{3} = vosc;  
        return;
    end
    
    [rosc,vosc] = osculating_orbit(rt0,vt0, t-t0, muC);
    r  = delta + rosc;       
    q  = delta'*(delta-2*r) / (r'*r);  
    fq = 1-(1+q)^(3/2);
            
    % Generate output
    a_main       = muC/norm(rosc)^3 * (fq*r - delta);
    a_perturb    = a_disturb(t,r);
    varargout{1} = a_main + a_perturb;
    
end


% Compute position & velocity at a certain time in osculating orbit
%
% References: 
% [1] S.W. Shepperd, "Universal Keplerian State Transition Matrix".
% Celestial Mechanics 35(1985) pp. 129--144, DOI: 0008-8714/85.15
function [rosc, vosc] = osculating_orbit(rt,vt, dt, muC) 
% %#eml
% eml.extrinsic('warning', 'error');
    
    % Quick exit for trivial case    
    rosc = rt;
    vosc = vt;
    if (dt == 0)
        return; end
    
    % Intitialize
    maxiters = 50;
    r1   = rt(:).';
    v1   = vt(:).';     
    nu0  = r1*v1.';
    r1m  = norm(r1);
    beta = 2*muC/r1m - v1*v1.';
    ulim = 1/sqrt(abs(beta));
    
    % Take care of period effects
    DeltaU = 0;
    if (beta > 0)
        P = 2*pi*muC*beta^(-3/2);
        n = floor((dt + P/2 - 2*nu0/beta)/P);
        DeltaU = 2*pi*n*beta^(-5/2);
    end
    
    % Loop until convergence of the time step
    iter   = 0;      t = 0;   regulafalsi_count = 0;
    deltaT = t-dt;   u = 0; 
    while abs(deltaT) > 1e2*max(eps([t dt])) % accuracy in seconds                                              
        iter = iter + 1;                     % (two last digits are allowed to differ)
        
        % Compute q
        % NOTE: [q] may not exceed 1/2. In principle, this will never
        % occur, but the iterative nature of the procedure can bring
        % it above 1/2 for some iterations.
        bu = beta*u*u;
        q  = bu/(1 + bu);
                
        % Escape clause;
        % The value for [q] will almost always stabilize to a value less
        % than 1/2 after a few iterations, but NOT always. 
        if (iter > maxiters) 
            if (q > 0.5)
                error(...
                    'lagrange_coefs:q_unstable',...
                    'Could not find solution in %d iterations.', maxiters);
            else
                warning(...
                    'lagrange_coefs:time_accuracy_not_met',...
                    ['Could not find Lagrange coefficients with time accuracy better ',...
                    'than %e. Assuming convergence...'], 1e2*max(eps([t dt])));
                break
            end
        end
        
        % Evaluate continued fraction         
        if q < 1 % convergence
            A =  1;   B = 1;   G = 1;   n = 0;      
            k = -9;   d = 15;  l = 3;   Gprev = inf;
            while abs(G-Gprev) > 1e2*max(eps([G Gprev]))
                k = -k;                 l = l + 2;
                d = d + 4*l;            n = n + (1+k)*l;
                A = d/(d - n*A*q);      B = (A-1)*B;
                Gprev = G;              G = G + B;
            end % continued fraction evaluation
        else
            error(...
                'lagrange_coefs:continued_fraction_diverges',...
                'Continued fraction diverges for q = %e', q);
        end
        
        % Kepler iteration
        U0w2   = 1 - 2*q;
        U1w2   = 2*(1-q)*u;
        U      = 16/15*U1w2^5*G + DeltaU;
        U0     = 2*U0w2^2-1;
        U1     = 2*U0w2*U1w2;
        U2     = 2*U1w2^2;
        U3     = beta*U + U1*U2/3;
        r      = r1m*U0 + nu0*U1 + muC*U2;
        t      = r1m*U1 + nu0*U2 + muC*U3;
        deltaT = t - dt;
        
        % Newton-Raphson method works most of the time, but is
        % not too stable; the method fails far too often for my
        % liking...
        % u = u - deltaT/4/(1-q)/r;
        % Halley's method is much better in that respect:
        uPrev = u;
        u = u - deltaT/((1-q)*(4*r + deltaT*beta*u));
        
        % But that too might fail. If that is the case, try to rescue 
        % by making a single Regula-Falsi step 
        if abs(u) > ulim             
            u = (uPrev - sign(deltaT)*ulim)/2;
            regulafalsi_count = regulafalsi_count + 1;
            if regulafalsi_count > maxiters/2
                error(...
                    'lagrange_coefs:u_unstable',...
                    'Could not locate zero of deltaT(u).');
            end
        end
        
    end % time loop
    
    if iter ~= 0
        % Produce Lagrange coefficients
        f = 1 - muC/r1m*U2;     F = -muC*U1/r/r1m;
        g = r1m*U1 + nu0*U2;    G = 1 - muC/r*U2;
        
        % and the new orbital position and velocity
        rosc = f*rt + g*vt;
        vosc = F*rt + G*vt;
    end
    
end 


% TODO: test if progress_orbit is faster with coordinate transformations or
%       MEX function 
% TODO: dense output?
% TODO: single output argument is an object, with a plotting method
%       (and methods similar to DEVAL)


function [t,... % time
          r,...
          v,...
          stats] = Encke(solver,...
                         perturbation,...                         
                         tspan,...
                         rt0,...
                         vt0,...
                         muC,...
                         varargin)
                     
    % Please report bugs and inquiries to:
    %
    % Name       : Rody P.S. Oldenhuis
    % E-mail     : oldenhuis@gmail.com    (personal)
    %              oldenhuis@luxspace.lu  (professional)
    % Affiliation: LuxSpace sarl
    % Licence    : BSD


    % If you find this work useful, please consider a donation:
    % https://www.paypal.me/RodyO/3.5


% Following Battin, pages 448-449
    
    %% Initialize Encke
    
    check_input(nargin, nargout);
    encke_options = get_options(varargin{:});
    
    if isempty(solver)
        % Defaults to best integrator in standard MATLAB
        solver = ODESolver();
    end

    % Make sure the disturbing acceleration function works
    try
        switch solver.order
            case 1
                [~] = perturbation(tspan(1),rt0,vt0);
            case 2
                [~] = perturbation(tspan(1),rt0);
        end

    catch ME
        ME2 = MException([mfilename ':perturbation_failure'], [...
                         'Failed to evaluate disturbing acceleration function ',...
                         'for initial values.']);
        throw(addCause(ME,ME2));
    end

    % Deal with existing output functions
    existing_outputFcn = odeget(solver.options, ...
                                'OutputFcn', @(varargin)false );
    
    % Initialize Encke's outputs
    stats = struct('solver'              , solver,...                   
                   'delta'               , [],...
                   'delta_dot'           , [],...  
                   'rectification_times' , [],...
                   'function_evaluations', 1); % NOTE: one eval is done to check 
                                               % the user-provided function    

    % Initially, the osculating orbit equals the initial values
    tosc = tspan(1);
    rosc = rt0;
    vosc = vt0;

    %% Run solver

    % When integrating differences (instead of position/velocity directly),
    % the absolute tolerance should remain unchanged. The relative
    % tolerance should however be adjusted, because the basis against which 
    % the relative tolerance is measured is of course much smaller

    % Set default AbsTol
    % RelTol: rescale to 10% of initial position -- this agrees with the 
    % criterion for rectification
    reltol_rescale = 1/(encke_options.scale_factor);
    solver.options = odeset(solver.options, ...
                            'AbsTol',                  odeget(solver.options, 'AbsTol', 1e-8),...
                            'RelTol', reltol_rescale * odeget(solver.options, 'RelTol', 1e-6) );

    % Initilaize output variables
    t_out = [];
    r_out = [];
    v_out = [];

    % These must be visible to the inner functions ("locally global)"    
    rosc_slice = [];
    vosc_slice = [];

    % Rectification is carried out via OutputFcn.
    %
    % NOTE: An output function is used to detect an event. Although an
    % EvenFcn would seem better suited, in this case, we're not
    % interested in *exactly* where the rectification should take place
    % (standard behavior in almost all ODE solvers), but rather, whether
    % it makes sense to do it in the next step. Therefore, to save on
    % function evaluations, the terminal property of the OutputFunction is
    % used, which makes sure the condition is only checked *once* per
    % iteration.
    outputFcn_wrapper = @(varargin) check_rectify(varargin{:}) || ...
                                    existing_outputFcn_wrapper(varargin{:});

    solver.options = odeset(solver.options, ...
                            'OutputFcn', outputFcn_wrapper);

    % Main loop:
    % - Let solver integrate until rectification is needed
    % - Break off integration and collect data
    % - Apply recitifcation
    % - Repeat from the start, until the final time is reached.
    done = false;
    while ~done

        % The 2 solver types have 2 different calls
        switch solver.order
            case 1
                D = solver.funfcn(@de1, ...
                                  tspan, ...
                                  zeros(6,1),...
                                  solver.options);

                delta    = D.y(1:3,:).';
                deltadot = D.y(4:6,:).';

            case 2 
                D = solver.funfcn(@de2, ...
                                  tspan, ...
                                  zeros(3,1),...
                                  zeros(3,1),...
                                  solver.options);

                delta    = D.y;
                deltadot = D.yp;
                                
        end

        % Append solver outputs to stats output
        tslice = D.x(:);
        stats.function_evaluations = stats.function_evaluations + ...
                                     D.stats.nfevals;

        stats.delta     = [stats.delta;        delta(1:end-~done,:)];
        stats.delta_dot = [stats.delta_dot; deltadot(1:end-~done,:)];

        % Save values for next iteration
        hinitial = tslice(end) - tslice(end-1);
        if numel(tslice) > 2
            hinitial = tslice(end-1) - tslice(end-2); end

        tspan = [tslice(end) tspan(end)];
        done  = tslice(end)==tspan(end);

        % We can in all likelihood keep using the step size from the
        % previous integration
        solver.options = odeset(solver.options,...
                                'InitialStep', 0.8*hinitial);

        % Insert data into output arrays
        t_out = [t_out;  tslice(1:end-~done)];                                         %#ok<AGROW>
        r_out = [r_out   [rt0 rosc_slice(:,1:end-~done)] +    delta(1:end-~done,:).']; %#ok<AGROW>
        v_out = [v_out   [vt0 vosc_slice(:,1:end-~done)] + deltadot(1:end-~done,:).']; %#ok<AGROW>

        % Apply rectification
        rt0 = rosc;
        vt0 = vosc;

    end
    
    %% Finalize

    % Rename and transpose outputs for consistency with
    % outputs from ODE suite
    t = t_out;
    r = r_out.';
    v = v_out.';

    %% Nssted functions
    
    % Check function inputs
    function check_input(argc, argo)

        % Check I/O argument counts
        error(   nargchk(6,inf,argc,'struct')); %#ok<*NCHKN> (narginchk() without error does not work in this case; 
        error(nargoutchk(0,4,  argo,'struct')); %#ok<*NCHKE>  we're not checking THIS function's I/O counts)
        
        assert(isa(solver, 'ODESolver'), ...
              [mfilename ':ioarg_error'],...
              'Argument ''solver'' must be a valid ODESolver object.');

        assert(isa(perturbation, 'function_handle'), ...
              [mfilename ':ioarg_error'],...
              'Argument ''perturbation'' must be a valid function handle.');

        assert(any(solver.order == [1 2]),...
               [mfilename ':unsupported_solver_order'],...
               'Solver order must be equal to 1 or 2.');

        assert(isnumeric(tspan) && isvector(tspan) && numel(tspan) >= 2 && ...
               all(isfinite(tspan)) && all(isreal(tspan)) && ...
               issorted(tspan) && numel(tspan) == numel(unique(tspan)), ...
               [mfilename ':ioarg_error'], [...
               'The time span vector (''tspan'') must be given as real, sorted ',...
               'values.']);

        assert(isnumeric(rt0) && isvector(rt0) && numel(rt0) == 3 && all(isreal(rt0)) && all(isfinite(rt0)) && ...
               isnumeric(vt0) && isvector(vt0) && numel(vt0) == 3 && all(isreal(vt0)) && all(isfinite(vt0)),...
               [mfilename ':ioarg_error'], [...
               'The initial position and velocity (''rt0'', ''vt0'') must both be ',...
               'real 3-element vectors.']);

        assert(isnumeric(muC) && isscalar(muC) && muC > 0 && isreal(muC),...
               [mfilename ':ioarg_error'], [...
               'The standard gravitational parameter of the central body (''muC'') ',...
               'must be given as a real positive scalar.']);
           
    end
    
    function encke_options = get_options(varargin)
        
        % Default options
        encke_options = struct('scale_factor', 1.0 / 100); % Scale criterion for rectifications; 1.0%
        
        if nargin == 0
            return; end        
                
        % Non-defaults should be passed via parameter/value pairs
        assert(mod(nargin, 2)==0,...
               [mfilename ':pvpairs_expected'], ...
               'Parameter value pairs expected.');
        
        parameters = varargin(1:2:end);
        values     = varargin(2:2:end);
         
        for ii = 1:numel(parameters)
                
            parameter = parameters{ii};
            value     = values{ii};
            
            switch lower(parameter)
                
                case {'scale' 'scale_factor'}
                    % TODO: assertions
                    encke_options.scale_factor = value;   
                    
                otherwise
                    warning([mfilename ':unsupported_parameter'], ...
                            'Unsupported parameter: "%s"; ignoring...',...
                            parameter);
                    continue;
            end
            
        end
        
    end

    %{
    Wrapper for any existing output function. For the chained method, 
    if the existing output function asks the solver to terminate, also
    the outer loop should be broken
    %}
    function terminate = existing_outputFcn_wrapper(varargin)
        terminate = feval(existing_outputFcn, varargin{:});
        if (terminate)
            done = true; end        
    end


    % Little wrapper function; checks if rectification is needed
    function terminate = check_rectify(t, y, varargin)

        terminate = false;

        % solver order == 1: (t,y,flag)
        % solver order == 2: (t,y,dy,flag)
        flag = varargin{2 - (solver.order == 1)};

        if strcmp(flag, 'init')
            rosc_slice = [];
            vosc_slice = [];

        elseif isempty(flag)

            t_new = t(end);
            if solver.order == 1
                % (t,y,flag)                    
                delta_new    = y(1:3,end);
                deltadot_new = y(4:6,end);
            else
                % (t,y,dy,flag)                    
                delta_new    = y;
                deltadot_new = varargin{1};
            end

            [rosc_t, vosc_t] = advance_orbit(rosc,...
                                             vosc,...
                                             t_new - tosc,...
                                             muC);

            rosc_slice = [rosc_slice  rosc_t];
            vosc_slice = [vosc_slice  vosc_t];

            if norm(delta_new) > encke_options.scale_factor * norm(rosc_t)

                terminate = true;                                

                % Do rectification
                rosc = rosc_t + delta_new;
                vosc = vosc_t + deltadot_new;
                tosc = t_new;

                stats.rectification_times = [stats.rectification_times; t_new];

            end
        end
    end


    %{
    Compute first derivative of

          y(t) = [delta(t)      % == r(t) - rosc(t)
                  delta_dot(t)] % == v(t) - vosc(t)              

    This DE should be used if the disturbing acceleration depends on the
    velocity, or if the solver is designed to solve systems of the type

            dy/dt = f(t,y)

    s.t.  
            y(t0) = y0

    %}
    function dydt = de1(t,y)

        delta     = y(1:3);
        delta_dot = y(4:6);

        [rosc_t, vosc_t] = advance_orbit(rosc,...
                                         vosc,...
                                         t - tosc,...
                                         muC);

        r_t  = delta     + rosc_t;
        v_t  = delta_dot + vosc_t;
        rt2  = r_t.'*r_t;
        q_t  = delta.' * (delta - 2*r_t) / rt2;
        fq   = 1 - (1 + q_t)^(3/2);

        a_main    = muC/norm(rosc_t)^3 * (fq*r_t - delta);
        a_perturb = perturbation(t, r_t, v_t);
        a_total   = a_main + a_perturb;

        % Generate output
        dydt = [delta_dot  % d (delta)/dt
                a_total];  % d²(delta)/dt²

        % Check whether current scale still justifies the use of Encke 
        check_scale(a_perturb, rt2);

    end


    %{
    Compute second derivative of

          delta(t) = r(t) - rosc(t)

    This DE should be used if the disturbing acceleration is independent of
    the velocity. If this is the case, the solver must ba designed to solve
    systems of type

          d²y/dt² = f(t,y)

    s.t.
          y(t0) = y0
          dy/dt = v0
    %}
    function d2deltadt2 = de2(t, delta)

        rosc_t = advance_orbit(rosc,...
                               vosc,...
                               t - tosc,...
                               muC);

        r_t = delta + rosc_t;
        rt2  = r_t.'*r_t;
        q_t = delta.' * (delta - 2*r_t) / rt2;
        fq  = 1 - (1 + q_t)^(3/2);

        % Generate output
        a_main     = muC/norm(rosc_t)^3 * (fq*r_t - delta);
        a_perturb  = perturbation(t, r_t);
        d2deltadt2 = a_main + a_perturb;

        % Check whether current scale still justifies the use of Encke 
        check_scale(a_perturb, rt2);

    end
    
    % Check whether Encke's method still makes sense
    function check_scale(a_perturb, rt2)
        
        persistent scale_warning_issued 
        if isempty(scale_warning_issued)
            scale_warning_issued = false; end
        
        if ~scale_warning_issued && encke_options.scale_factor*muC/rt2 < norm(a_perturb)
            scale_warning_issued = true;
            warning([mfilename ':scale_issues'], [...
                    'Perturbative acceleration is over %.1f%% of the main ',...
                    'acceleration; Encke''s method may be inefficient for ',...
                    'this problem.'],...
                    encke_options.scale_factor*100);
        end
        
    end

end



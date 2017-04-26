% Compute position & velocity after a certain time step
%
% References:
% [1] S.W. Shepperd, "Universal Keplerian State Transition Matrix".
% Celestial Mechanics 35(1985) pp. 129--144, DOI: 0008-8714/85.15
function [r_new, v_new] = advance_orbit(rt, vt, dt, muC)%#eml
 
    % EML/CODER inits
    eml.extrinsic('mfilename', 'warning', 'error');
    r  = 0;
    U1 = 0;
    U2 = 0;
    
    % Quick exit for trivial case
    r_new = rt;
    v_new = vt;
    if dt == 0
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
                error([mfilename 'lagrange_coefs_q_unstable'],...
                      'Could not find solution in %d iterations.', maxiters);
            else
                warning([mfilename 'lagrange_coefs_time_accuracy_not_met'],...
                        ['Could not find Lagrange coefficients with time accuracy better ',...
                        'than %e. Assuming convergence...'], ...
                        1e2*max(eps([t dt])));
                break
            end
        end

        % Evaluate continued fraction
        A =  1;   B = 1;   G = 1;   n = 0;
        k = -9;   d = 15;  l = 3;   Gprev = inf;
        
        if q < 1 % convergence            
            while abs(G-Gprev) > 1e4*max(eps([G Gprev]))
                k = -k;                 l = l + 2;
                d = d + 4*l;            n = n + (1+k)*l;
                A = d/(d - n*A*q);      B = (A-1)*B;
                Gprev = G;              G = G + B;
            end % continued fraction evaluation

        else
            error([mfilename 'lagrange_coefs_continued_fraction_diverges'],...
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
        uPrev  = u;

        % Newton-Raphson method works most of the time, but is
        % not too stable; Halley's method is much better:
        % u = u - deltaT/4/(1-q)/r;
        u = u - deltaT/((1-q)*(4*r + deltaT*beta*u));

        % But that too might fail. If that is the case, try to rescue
        % by making a single Regula-Falsi step
        if abs(u) > ulim
            u = (uPrev - sign(deltaT)*ulim)/2;
            regulafalsi_count = regulafalsi_count + 1;
            if regulafalsi_count > maxiters/2
                error([mfilename 'lagrange_coefs_u_unstable'],...
                      'Could not locate zero of deltaT(u).');
            end
        end

    end % time loop

    if iter ~= 0
        % Produce Lagrange coefficients
        f = 1 - muC/r1m*U2;     F = -muC*U1/r/r1m;
        g = r1m*U1 + nu0*U2;    G = 1 - muC/r*U2;

        % and the new orbital position and velocity
        r_new = f*rt + g*vt;
        v_new = F*rt + G*vt;
    end

end


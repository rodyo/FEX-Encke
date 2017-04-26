function Encke_Demo
    
    %emlmex -eg {[0 0 0]', [0 0 0]', 0, 0} advance_orbit.m

    % Please report bugs and inquiries to:
    %
    % Name       : Rody P.S. Oldenhuis
    % E-mail     : oldenhuis@gmail.com    (personal)
    %              oldenhuis@luxspace.lu  (professional)
    % Affiliation: LuxSpace sarl
    % Licence    : BSD


    % If you find this work useful, please consider a small donation:
    % https://www.paypal.me/RodyO/3.5


    clc

    % initial values for test orbit
    RE     = 6378;
    muC    = 398600.44;
    rt0    = [RE 0 0]';
    vt0    = [0 1.0*sqrt(muC/RE) 0]';
    tspan0 = [0 2.50*86400];
    tspan  = tspan0;
    
 
    
    % Solver options
    options = odeset('Abstol', 1e-6,... % in km, so accuracy per step is < 1mm
                     'reltol', 1e-6);   % in km/km


    %% Encke's method

   
    solver = ODESolver('rkn1210', options); 
    %solver = ODESolver('ode113', options);   
    %solver = ODESolver('ode113', options);    
   

    tic
    [t, r, v, stats] = Encke(solver, @a_disturb, tspan, rt0, vt0, muC);
    
    
    fprintf(1, ...
            '\nEncke''s method required %f seconds (%d function evaluations).\n',...
            toc, stats.function_evaluations);
        

    %% Cowell's method
    
    function dydt = cowell1(tt,y, muC)
        
        drdt   = y(4:6);
        d2rdt2 = -muC * y(1:3) ./ norm(y(1:3)).^3 + ...
                a_disturb1(tt, y(1:3), y(4:6));

        dydt = [drdt
                d2rdt2];
    end

    function dy2dt2 = cowell2(t,y, muC)
        dy2dt2 = -muC*y/norm(y)^3 + a_disturb(t,y); 
    end
    
                 
    % Compare accuracy
    % NOTE: we need to have values at equal times to compute this,
    % hence this clumbsiness    
    switch solver.get_name
        
        case 'rkn1210'
            [~,rC,vC] = rkn1210(@(t,y)cowell2(t,y, muC), t, rt0,vt0);    
            
        case 'ode113'
            [~,rvC] = ode113(@(tt,y)cowell1(tt,y, muC), t, [rt0; vt0], options);
         
            rC = rvC(:,1:3);
            vC = rvC(:,4:6);    
    end
    
    % Plot comparison        
    r_ERR = sqrt(sum((rC-r).^2,2));
    v_ERR = sqrt(sum((vC-v).^2,2));

    figure(2), clf
    plot(t, r_ERR * 1e3)
    hold on, grid on, title('Position difference [m]')

    figure(3), clf
    plot(t, v_ERR * 1e3)
    hold on, grid on, title('Velocity difference [m/s]')

    % Compare CPU time and fevals
    % NOTE: to have a fair comparison of CPU time requirements, we
    % have to re-do the integration, this time with integrator settings
    % optimised for CPU time
    tic
    switch solver.get_name
        case 'rkn1210'    
            [~,rC,~,~,output] = rkn1210(@(t,y)cowell2(t,y, muC), tspan0, rt0,vt0, options);
            fevals = output.fevals;
    
        case 'ode113'
            PP = ode113(@(t,y)cowell1(t,y, muC), tspan0, [rt0; vt0], options);
            fevals = PP.stats.nfevals;
            rC = PP.y(1:3,:).';
    end
    
    fprintf(1,...
            'Cowell''s method required %f seconds (%d function evaluations).\n', ...
            toc, fevals);
    
    % Plot position results
    figure(1), clf
    plot( r(:,1),  r(:,2), 'b', ...
         rC(:,1), rC(:,2), 'r')
    axis equal
    legend('Encke''s method',...
           'Cowell''s method')

end

function a_perturb = a_disturb(~,r, varargin)

    a_perturb = [-r(2); r(1); 0];
    a_perturb = a_perturb/norm(a_perturb) * 9.2e-7; % NOTE: km/s²
    
end

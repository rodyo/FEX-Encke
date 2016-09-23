function Encke_Demo
    
    % Please report bugs and inquiries to:
    %
    % Name       : Rody P.S. Oldenhuis
    % E-mail     : oldenhuis@gmail.com    (personal)
    %              oldenhuis@luxspace.lu  (professional)
    % Affiliation: LuxSpace s�rl
    % Licence    : BSD
    
    
    % If you find this work useful, please consider a donation:
    % https://www.paypal.me/RodyO/3.5

   
    clc
    
    % initial values for test orbit
    RE  = 6378;
    muC = 398600.44;
    rt0 = [RE 0 0]';
    vt0 = [0 1.0*sqrt(muC/RE) 0]';
    tspan0 = [0 0.5*86400];
    tspan = tspan0;
    
 
    %% Encke's method
    
options = odeset('Abstol', 1e-14, 'reltol', 1e-7);
    
    options.handle = @rkn1210;
    options.velocity_required = 'no';
    a_disturb = @a_disturb2;
    
    %{
    options.handle = @ode113;
    options.velocity_required = 'yes';
    a_disturb = @a_disturb1;
    %}

    tic
    [t, r, v, stats] = Encke(a_disturb, tspan, rt0, vt0, muC, options) ;
    
    fprintf(1, '\nEncke''s method required %f seconds (%d function evaluations).\n',...
    toc, stats.fevals);




    %% Cowell's method

    % Compare accuracy 
    % NOTE: we need to have values at equal times to compute this, 
    % however, doing it like this is slow in MATLAB...
    % {
    [~,rC,vC] = rkn1210(@(t,y)cowell(t,y, muC), t, rt0,vt0);
    function dy2dt2 = cowell(t,y, muC)
        dy2dt2 = -muC*y/norm(y)^3 + a_disturb2(t,y); end
    %}
    
    %{
    [~,rvC] = ode113(@(t,y)cowell1(t,y, muC), t, [rt0 vt0], options);
    function dydt = cowell1(t,y, muC)
        dydt = [
            y(4:6)
            -muC*y(1:3)/norm(y(1:3))^3 + a_disturb2(t,y)];
    end
    rC = rvC(:,1:3);
    vC = rvC(:,4:6);
    %}
    
    
    % Plot accuracy comparison
    r_ERR = sqrt(sum((rC-r').^2,2));
    v_ERR = sqrt(sum((vC-v').^2,2));
    
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
    % {
    tic
    [~,rC,~,~,output] = rkn1210(@(t,y)cowell(t,y, muC), tspan0, rt0,vt0); 
    fprintf(1, '\nCowell''s method required %f seconds (%d function evaluations).\n', ...
        toc, output.fevals);
    %}
    %{
    tic
    PP = ode113(@(t,y)cowell1(t,y, muC), tspan0, [rt0 vt0]); 
    fprintf(1, '\nCowell''s method required %f seconds (%d function evaluations).\n', ...
        toc, PP.stats.nfevals);
    r = PP.y(:,1:3);
    v = PP.y(:,4:6);
    %}
    
    % Plot position results 
    figure(1), clf
    plot(...
         r(1,:),  r(2,:), 'b', ...
        rC(:,1), rC(:,2), 'r')    
    axis equal, legend('Encke''s method', 'Cowell''s method')
        
end


function a_perturb = a_disturb1(~,r,v)
    
    a_perturb = 0;    
    return
    
    a_perturb = [-r(2) r(1) 0]';
    a_perturb = a_perturb/norm(a_perturb) * 1e-6*6378/norm(r(1:3));
    
    % artifical delay to simulate CPU time requirements for more 
    % elaborate force models (0.01s in 1/10 cases averages 1ms per feval)
    %if rand < 0.1
    %    pause(0.01); end
end

function a_perturb = a_disturb2(~,r)
    
    %a_perturb = 0;    
    %return
    
    a_perturb = [-r(2) r(1) 0]';
    a_perturb = a_perturb/norm(a_perturb) * 1e-6*6378/norm(r(1:3));    
    
    % artifical delay to simulate CPU time requirements for more 
    % elaborate force models (0.01s in 1/10 cases averages 1ms per feval)
    %if rand < 0.1
    %    pause(0.01); end
end

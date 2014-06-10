function Encke_Demo
   
    clc
    
RE  = 6378;
muC = 398600.44;  
rt0 = [RE 0 0]';
vt0 = [-0.1 1.0*sqrt(muC/RE) 0.01]';
tspan0 = [0 1e6];
tspan = tspan0;

    
    
options.handle = @rkn1210;
options.velocity_required = 'no';
a_disturb = @a_disturb2;
% 
% options.handle = @ode45;
% options.velocity_required = 'yes';
% a_disturb = @a_disturb1;

[t, r, v, stats] = Encke(a_disturb, tspan, rt0, vt0, muC, options) ;
stats
  


    figure(1), clf, hold on
    plot(r(1,:), r(2,:))
    axis equal
    
    
    % Compare accuracy
    
    rt0 = [RE 0 0]';
    vt0 = [-0.1 1.0*sqrt(muC/RE) 0.01]';
    [~,rE,vE] =  rkn1210(@(t,y)cowell(t,y, muC), t, rt0,vt0);
    function dy2dt2 = cowell(t,y, muC)
        dy2dt2 = -muC*y/norm(y)^3 + a_disturb2(t,y); end
    
    figure(2)
    r_ERR = sqrt(sum((rE-r').^2,2));
    plot(t,r_ERR)
    hold on
    v_ERR=sqrt(sum((vE-v').^2,2));
    plot(t,100*v_ERR, 'r')
    
    % Compare fevals
    rt0 = [RE 0 0]';
    vt0 = [-0.1 1.0*sqrt(muC/RE) 0.01]';
    [t,r,v,~,output] =  rkn1210(@(t,y)cowell(t,y, muC), tspan0, rt0,vt0);
    output.fevals
    
    figure(1)
    plot(r(:,1), r(:,2), 'r')
    axis equal
        
end


function a_perturb = a_disturb1(~,r,v)
    a_perturb = [-r(2) r(1) 0]';
    a_perturb = a_perturb/norm(a_perturb) * 1e-4;
end

function a_perturb = a_disturb2(~,r)
    a_perturb = [-r(2) r(1) 0]';
    a_perturb = a_perturb/norm(a_perturb) * 1e-4;
end

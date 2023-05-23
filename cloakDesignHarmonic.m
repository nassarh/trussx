%% Reference medium
% properties
Eo = 1; 
Rho = 1;

% dimensions
W = 20;

% fiber direction
Fiber = [1,0];
    

%% Discretization of the reference medium
% "number of elements"
nb = 5;
N = 40*2^(nb+1);

% step sizes
dX = W/N;
dY = dX; % arbitrary

% bodies
% dimension, radius, mass, center of mass, moment of inertia, ports
myBody = rigidBody(2,dX/5,Rho*dX*dY,[0 0],Rho*dX*dY*(dX/5)^2,zeros(2,2));

% build reference medium
ref = trussx.fiber(N, dX, myBody, Eo*dY/dX);

% assemble
ref = ref.assemble();

%% cloaking frequency and reference solution
% boundary conditions (empty), number of eigenfrequencies, frequency level
[wref,uref,~] = ref.evSolve([], 1, 1);

% cloak's target frequency
wclk = wref;

%% Cloak geometry
% initialize cloak
cloak = trussx.fiber(N, dX, myBody, Eo*dY/dX);

% defect radius
a = 1;

% complexify, conformal transformation
X = cloak.x(:,1) + 1j*cloak.x(:,2) +1j*a/10 - W/2;
x = X/2 + X.*sqrt(1 - 4*a^2./X.^2)/2;

% decomplexify
cloak.x = [real(x)+W/2 imag(x)];

% transformation gradient (dx/dX)
F = x.^2./(x.^2-a^2);
stretch = abs(F);

% transformation Hessian
Hessian = -2*a^2*x./(x.^2-a^2).^2.*F;

%% Loop to stabilize cloaking frequency (by Fixed point iterations)
for iter = 1:1
    %% Cloak's properties
    % the bodies
    for i = 1:cloak.N
        % steps
        dx = dX*stretch(i);
        dy = dY*stretch(i);
        
        % theta
        theta = dx/10;
        
        % Young's
        E = 2*Eo*stretch(i)^2;
        
        % moment of inertia
        J = -E*theta^2*dx*dy/wclk^2;
        
        % center of mass
        cm = theta*E/wclk^2/Rho/stretch(i)^2*Hessian(i)*(Fiber(1)+1j*Fiber(2))^2;
        cm = [imag(cm), -real(cm)];
        
        % ports positions
        
        % transformed fiber
        fiber = F(i)*(Fiber(1)+1j*Fiber(2));
        fiber = [real(fiber), imag(fiber)];
        fiber = fiber/norm(fiber);
    
        fiber_n = [-fiber(2), fiber(1)];
        
        port = dx*theta*fiber_n/2;
        ports = [port; -port];
        
        % assign bodies
        cloak.bodies(i).M = Rho*dx*dy;
        cloak.bodies(i).cm = cm;
        cloak.bodies(i).J = J;
        cloak.bodies(i).x = ports;
    end
    
    for i = 1:cloak.E
        % the stiffness k
        p = cloak.edges(i,1);
        q = cloak.edges(i,2);
    
        sp = stretch(p);
        sq = stretch(q);
        spq = (sp+sq)/2;
    
        cloak.k(i) = 2*spq^2*ref.k(i);
    end

    %% assemble cloak
    cloak = cloak.assemble();
    
    %% get solution in cloak
    [wclk, u, ~] = cloak.evSolve([], 1, wclk^2);
end

%% visualize solutions

% pull back cloak displacement
Auref = (u(:,1)+1j*u(:,2)).*conj(F);
Auref = [real(Auref), imag(Auref)];

figure(1)

Fontsize = 15;

subplot(4,1,1)

scatter(ref.x(:,1),ref.x(:,2),40,vecnorm(uref(:,1:2),2,2)/abs(uref(N,1)),'filled')
axis equal off tight
hold on
scatter(cloak.x(:,1),cloak.x(:,2)+a,40,vecnorm(u(:,1:2),2,2)/abs(u(N,1)),'filled')
hold off

subplot(4,1,2)

plot(ref.x(:,1),Auref(:,1)/Auref(N,1), 'c', LineWidth=1)
hold on
plot(ref.x(:,1),u(:,1)/u(N,1),'k--',linewidth=2)
plot(ref.x(:,1),uref(:,1)/uref(N,1),'r',linewidth=1)
hold off
axis tight
grid on

xlabel('position','Interpreter','latex')
ylabel('horizontal displacement','Interpreter','latex')

xAX = get(gca,'XAxis'); 
set(xAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')
yAX = get(gca,'YAxis'); 
set(yAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')

subplot(4,1,3)

plot(ref.x(:,1),Auref(:,2)/Auref(N,1), 'c', LineWidth=1)
hold on
plot(ref.x(:,1),u(:,2)/u(N,1),'k--',linewidth=2)
plot(ref.x(:,1),uref(:,2)/uref(N,1),'r',linewidth=1)
hold off
axis tight
grid on

xlabel('position','Interpreter','latex')
ylabel('vertical displacement','Interpreter','latex')

xAX = get(gca,'XAxis'); 
set(xAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')
yAX = get(gca,'YAxis'); 
set(yAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')

lgd = legend('reference disp.','cloak disp.','cloak disp. pulled back');
lgd.Interpreter = 'latex';
lgd.Location = 'southwest';

subplot(4,1,4)

plot(ref.x(:,1),u(:,3)*dX/u(N,1), 'k--', LineWidth=2)
axis tight
grid on

xlabel('position','Interpreter','latex')
ylabel('rotation','Interpreter','latex')

xAX = get(gca,'XAxis'); 
set(xAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')
yAX = get(gca,'YAxis'); 
set(yAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')

% %% convergence analysis
% % error on horizontal displacement
% err_hd(nb) = norm(Auref(:,1)/Auref(N,1)-uref(:,1)/uref(N,1))/N;
% 
% % error on vertical displacement
% err_vd(nb) = norm(Auref(:,2)/Auref(N,1))/N;
% 
% % error on eigenfrequency
% err_w(nb) = abs(wref-wclk);
% 
% %% plot errors
% loglog(40*2.^(2:11),err_hd,'-o','MarkerSize',10,'LineWidth',2)
% hold on
% loglog(40*2.^(2:11),err_vd,'-x','MarkerSize',10,'LineWidth',2)
% loglog(40*2.^(2:11),err_w,'-^','MarkerSize',10,'LineWidth',2)
% 
% xlabel('number of nodes','Interpreter','latex')
% ylabel('error','Interpreter','latex')
% 
% Fontsize = 15;
% 
% xAX = get(gca,'XAxis'); 
% set(xAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')
% yAX = get(gca,'YAxis'); 
% set(yAX,'FontSize', Fontsize, 'TickLabelInterpreter','latex')
% title("error analysis, no fixed-point iteration")
% grid on
%
% %% bonus point with fixed point iteration
% bonus = abs(wref-wclk);
% scatter(40*2^11,bonus,90,'^','filled')
% 
% lgd = legend('quad. err. horiz. disp.','quad. err. vert. disp.','err. eigenfreq.', 'err. eigenfreq. + 10 FP');
% lgd.Interpreter = 'latex';
% lgd.Location = 'southwest';
% model parameters
m = 5;
MH = 15;
MT = 10;
r = 1;
l = 0.5;
g = 9.81;

parameters = [
    m, MH, MT, r, l, g
];

% control parameters
theta1d = pi/8;

theta3d = pi/6;
epsilon = 0.1;
alpha = 0.9;

controller_parameters = [
    alpha, epsilon, theta3d
];

% simulation parameters
T = 10;
dt = 1e-4;
N = floor(T/dt);

x0 = [
    -0.35; % TODO
    0; % TODO
    0;
    %pi/8; % TODO
    1;
    0;
    1.8;
];
% x0 = [
%     0.01;
%     0;
%     0;
%     0;
%     0;
%     0;
% ];
xsim = zeros(6,N);
xsim(:,1) = x0;
tsim = 0:dt:T;
tsim = tsim(1:end-1);
small_number = 1e-3;

% simulation using Forward Euler
for i=1:N-1
    t = tsim(i);
    x = xsim(:,i);

    theta1 = x(1);

    if theta1d - small_number < theta1 && theta1 < theta1d + small_number
        % impact
        fprintf('Impact at %.2fs\n', t)
        xnext = ImpactModel(x, parameters);
    else
        % continuous model
        u = Controller(x, parameters, controller_parameters);
        xnext = x + dt * MechanicalModel(x, u, parameters);
    end
    xsim(:,i+1) = xnext;
end

figure(10); clf;
plot(xsim(1,:), xsim(6,:));

figure(11); clf;
plot(tsim, xsim(1,:));

figure(12); clf;
plot3(xsim(1,:), xsim(4,:), xsim(6,:));
grid on;
ylim([0,2]);
xlim([-.4, .4]);
zlim([-0.5,2.5]);
xlabel("\(\theta_1\)","Interpreter","latex");
ylabel("\(\omega_1\)","Interpreter","latex");
zlabel("\(\omega_3\)","Interpreter","latex");


% animation
figure(13); clf;
pause(1)
for k=1:floor(N/200)
    i = k*200;
    theta1 = xsim(1,i);
    theta2 = xsim(2,i);
    theta3 = xsim(3,i);

    % connections
    plot([0, r*sin(theta1)], [0, r*cos(theta1)], 'k-','LineWidth',1); % stance leg
    hold on; grid on;
    plot([r*sin(theta1),r*sin(theta1)-r*sin(theta2)],[r*cos(theta1),0], 'k-','LineWidth',1); % swing leg
    plot([r*sin(theta1),r*sin(theta1)+l*sin(theta3)], [r*cos(theta1),r*cos(theta1)+l*cos(theta3)], 'k-','LineWidth',1); % torso
    
    % masses
    plot(0.5*r*sin(theta1), 0.5*r*cos(theta1), 'k.','MarkerSize',25); % stance leg
    plot(r*sin(theta1)-0.5*r*sin(theta2),0.5*r*cos(theta2), 'k.','MarkerSize',25); % swing leg
    plot(r*sin(theta1), r*cos(theta1), 'k.','MarkerSize',25); % hip
    plot(r*sin(theta1)+l*sin(theta3), r*cos(theta1)+l*cos(theta3), 'k.','MarkerSize',25); % torso
    hold off;


    xlim([-1.5, 1.5]);
    ylim([-0.5,2.5]);
    pbaspect([1,1,1]);
    title(sprintf("Time: %.2fs", tsim(i)));

    pause(5/1000)
end

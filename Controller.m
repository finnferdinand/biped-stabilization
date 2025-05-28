function u = Controller(x, system_parameters, controller_parameters)
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    m = system_parameters(1);
    MH = system_parameters(2);
    MT = system_parameters(3);
    r = system_parameters(4);
    l = system_parameters(5);
    g = system_parameters(6);

    alpha = controller_parameters(1);
    epsilon = controller_parameters(2);
    theta3d = controller_parameters(3);

    s12 = sin(theta1 - theta2);
    s13 = sin(theta1 - theta3);
    c12 = cos(theta1 - theta2);
    c13 = cos(theta1 - theta3);

    D = [
        (5/4*m + MH + MT) * r^2, -1/2*m*r^2*c12, MT*r*l*c13;
        -1/2*m*r^2*c12, 1/4*m*r^2, 0;
        MT*r*l*c13, 0, MT*l^2;
    ];
    C = [
        0, -1/2*m*r^2*s12*omega2, MT*r*l*s13*omega3;
        1/2*m*r^2*s12*omega1, 0, 0;
        -MT*r*l*s13*omega1, 0, 0;
    ];
    G = [
        -1/2*g*(2*MH + 3*m + 2*MT)*r*sin(theta1);
        1/2*g*m*r*sin(theta2);
        -g*MT*l*sin(theta3);
    ];

    f = D \ (-C * [omega1; omega2; omega3] - G);
    Lf2h = [
        f(3);
        f(2) + f(1);
    ];

    B = [
        -1, 0;
        0, -1;
        1, 1;
    ];
    g = D \ B;
    decoupling_matrix = [
        0, 0, 1;
        1, 1, 0;
    ] * g;

    % The decoupling matrix as presented in the paper.
%     detD = 1/4*m*MT*r^4*l^2*(5/4*m + MH + MT - m*c12^2 - MT*c13^2);
% 
%     R11 = m*r^3/4 * (5/4*m*r + MH*r + MT*r - m*r*c12^2 + MT*l*c13);
%     R12 = m*r^3/4 * (5/4*m*r + MH*r + MT*r - m*r*c12^2 + 2*MT*l*c12*c13);
%     R21 = -m*MT*l*r^2/4*(1 + 2*c12)*(r*c13 + l);
%     R22 = -MT*l*r^2/4*(5*m*l + 4*MH*l + 4*MT*l + m*r*c13 + 2*m*r*c12*c13 - 4*MT*l*c13^2 + 2*m*l*c12);
%     
%     decoupling_matrix = 1/detD * [
%         R11, R12;
%         R21, R22;
%     ];

    y1 = theta3 - theta3d;
    y2 = theta2 + theta1;
    y1dot = omega3;
    y2dot = omega2 + omega1;

    % The double integrator controller described in the paper. Possibly incorrect parentheses?
%     phia1 = y1 + (1/2-alpha)*sign(epsilon*y1dot)*abs(epsilon*y1dot)^(2-alpha);
%     psia1 = -sign(epsilon*y1dot)*abs(epsilon*y1dot)^alpha - sign(phia1)*abs(phia1)^(alpha/2-alpha);
%     
%     phia2 = y2 + (1/2-alpha)*sign(epsilon*y2dot)*abs(epsilon*y2dot)^(2-alpha);
%     psia2 = -sign(epsilon*y2dot)*abs(epsilon*y2dot)^alpha - sign(phia2)*abs(phia2)^(alpha/2-alpha);

    % The controller described by Bhat & Bernstein.
    phia1 = y1 + (1/(2-alpha))*sign(epsilon*y1dot)*abs(epsilon*y1dot)^(2-alpha);
    psia1 = -sign(epsilon*y1dot)*abs(epsilon*y1dot)^alpha - sign(phia1)*abs(phia1)^(alpha/(2-alpha));
    
    phia2 = y2 + (1/(2-alpha))*sign(epsilon*y2dot)*abs(epsilon*y2dot)^(2-alpha);
    psia2 = -sign(epsilon*y2dot)*abs(epsilon*y2dot)^alpha - sign(phia2)*abs(phia2)^(alpha/(2-alpha));


    Psi = [
        1/epsilon^2 * psia1;
        1/epsilon^2 * psia2;
    ];

    u = decoupling_matrix \ (Psi - Lf2h);
end

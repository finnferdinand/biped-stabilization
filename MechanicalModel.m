function xdot = MechanicalModel(x, u, parameters)
    m = parameters(1);
    MH = parameters(2);
    MT = parameters(3);
    r = parameters(4);
    l = parameters(5);
    g = parameters(6);

    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    omega1 = x(4);
    omega2 = x(5);
    omega3 = x(6);

    omega = [
        omega1;
        omega2;
        omega3;
    ];

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
    B = [
        -1, 0;
        0, -1;
        1, 1;
    ];
    %D,C,G,B

    xdot = [
        omega;
        D \ (-C * omega - G + B * u);
    ];
end

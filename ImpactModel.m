function xplus = ImpactModel(xminus, parameters)
    m = parameters(1);
    MH = parameters(2);
    MT = parameters(3);
    r = parameters(4);
    l = parameters(5);
    g = parameters(6);

    theta1minus = xminus(1);
    theta2minus = xminus(2);
    theta3minus = xminus(3);
    omega1minus = xminus(4);
    omega2minus = xminus(5);
    omega3minus = xminus(6);

    qedotminus = [
        omega1minus;
        omega2minus;
        omega3minus;
        0;
        0;
    ];

    % s12 = sin(theta1minus - theta2minus);
    % s13 = sin(theta1minus - theta3minus);
    c12 = cos(theta1minus - theta2minus);
    c13 = cos(theta1minus - theta3minus);

    De11 = (5/4*m + MH + MT)*r^2;
    De12 = -1/2*m*r^2*c12;
    De13 = MT*r*l*c13;
    De14 = (3/2*m + MH + MT)*r*cos(theta1minus);
    De15 = -(3/2*m + MH + MT)*r*sin(theta1minus);
    De22 = 1/4*m*r^2;
    De23 = 0;
    De24 = -1/2*m*r*cos(theta2minus);
    De25 = 1/2*m*r*sin(theta2minus);
    De33 = MT*l^2;
    De34 = MT*l*cos(theta3minus);
    De35 = -MT*l*sin(theta3minus);
    De44 = 2*m + MH + MT;
    De45 = 0;
    De55 = 2*m + MH + MT;
    
    De = [
        De11, De12, De13, De14, De15;
        De12, De22, De23, De24, De25;
        De13, De23, De33, De34, De35;
        De14, De24, De34, De44, De45;
        De15, De25, De35, De45, De55;
    ];
    
    E = [
        r*cos(theta1minus), -r*cos(theta2minus), 0, 1, 0;
        -r*sin(theta1minus), r*sin(theta2minus), 0, 0, 1;
    ];

    A = [
        De, -E';
        E, zeros(2,2);
    ] \ [
        De * qedotminus;
        zeros(2,1);
    ];

    omega1plus = A(1);
    omega2plus = A(2);
    omega3plus = A(3);

    xplus = [
        theta2minus;
        theta1minus;
        theta3minus;
        omega2plus;
        omega1plus;
        omega3plus;
    ];
end
function [A_, U_, Ux, L, R, D, x] = static_solver(l, cN)
close all;

vX = [0,l];
[L, R, D, x]=GetIM(cN,vX);

static_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux_1 = - eye(n)*B11;
wx_1 = zeros(n);
phx_1 = - NU31 * B11 * L;
gax_1 = zeros(n);

T11_1 = one_col;
ph0_1 = - one_col * NU31 * B11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux_2 = - R * B33 * NU13;
wx_2 = zeros(n);
phx_2 = - (D13*eye(n)  + R * B33 * L);
gax_2 = zeros(n);

R1_2 = one_col;
ph0_2 =  - B33 * (R * one_col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux_3 = zeros(n);
wx_3 = - R * B13;
phx_3 = zeros(n);
gax_3 = - (R*B13*L + D11*eye(n));

m11_3 = one_col;
ga_3 = - B13 * (R * one_col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ux_4 = zeros(n);
wx_4 = - B13 * eye(n);
phx_4 = zeros(n);
gax_4 = - B13 * L;

q1_4 = one_col;
ga_4 = - B13 * one_col;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [ux_1, wx_1, phx_1,gax_1;%zer_col,zer_col,ph0_1,zer_col, T11_1, zer_col, zer_col, zer_col;
    ux_4, wx_4, phx_4,gax_4;
    ux_2, wx_2, phx_2,gax_2;%zer_col,zer_col,ph0_2,zer_col, zer_col, R1_2, zer_col, zer_col;
    ux_3, wx_3, phx_3,gax_3;%zer_col,zer_col,zer_col,ga_3, zer_col, zer_col, m11_3, zer_col;
    ];%zer_col,zer_col,zer_col,ga_4, zer_col, zer_col, zer_col, q1_4;
    %zeros(8, 4*n), eye(8)];

Right = [P1;P4;P2;P3;];%u0;w0;ph0;ga0;T11_l;R1_l;M11_l;Q1_l];

U = A\Right;
u = L*U(1:n);% + U(4*n+1); % L*U(1:n) + U(4*n+1);
w = L*U(1+n:2*n);% + U(4*n+2);
phi = L*U(1+2*n:3*n);% + U(4*n+3);
ga = L*U(1+3*n:4*n);% + U(4*n+4);

% sum(abs((A(1:4*cN,1:4*cN)*U(1:4*cN)) - (Right(1:4*cN))))
% sum(abs((A(1:4*cN,1:4*cN)*[u;w;phi;ga]) - (Right(1:4*cN))))

sum(abs(A*U - Right))

U_ = [u;w;phi;ga];
A_ = A;
Ux = U;

end
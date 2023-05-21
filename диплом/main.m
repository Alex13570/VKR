%--------------------------------------------------------------------------
%mIM1(cN,cN)-интегрирующая матрица 1-го рода;
%mIM2(cN,cN)-интегрирующая матрица 2-го рода;
%cN-число сечений на всём интервале интегрирования;
%vX -вектор с элементами равными координатам граничных точек
%             подинтервалов, включая начало и конец интервала
%             интегрирования;
%--------------------------------------------------------------------------
clear variables;clc; close all;

cN = 10;
l = 1;

static_data
rho = 2700;

lambda = 3.52;
nhu = (lambda/(l*l))*(E1*h^2/(12*rho))^0.5;
% nhu = 10000;


cT = 100;% 11.1*pi*2;
t = 0;
dt =  1/cT/100;


[A, U, Ux, L, R, D, x] = static_solver(l, cN);
n = (length(U))/4;
Ux_1 = Ux;
% figure(1)
% plot(x, L*Ux(n+1:2*n))
% xlabel('x');
% ylabel('w');
% title('перемещение вдоль оси y. t = 0');
% saveas (gcf, 'u', 'jpeg')
% hold off
% 
% figure(2)
% plot(x,B13*(L*Ux(3*n+1:4*n) + Ux(1*n+1:2*n)))
% xlabel('x');
% ylabel('T13');
% title('усилия T13. t = 0');
% saveas (gcf, 'T13', 'jpeg')

mas = kron(diag([rho*h,rho*h,rho*h^3/12, rho*h^3/12]),(R*L))/(dt*dt);

for t = 0:dt:1000
[U1] = next_step(t, A, Ux, Ux_1, mas, cT);
Ux_1 = Ux;
Ux = U1;
figure(2)
plot(x, L*U1(1*n+1:2*n))
hold off
end

figure;
p = plot(x, u);
xlabel('x');
ylabel('u');
titlestr = 'перемещения вдоль оси х';
title(titlestr);
p(1).LineWidth = 2;
saveas (gcf, 'u', 'jpeg')





% xlswrite('res.xlsx',answertable, 'Table1')

% B = [12;12;12];
% C = [A,A, B,B;
%     A,B,B,A]
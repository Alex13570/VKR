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


cT = 10;% 11.1*pi*2;
t = 0;
dt =  1/cT/100;


[A, U, Ux, L, R, D, x] = static_solver(l, cN);
n = (length(U))/4;
Ux_1 = Ux;

mas = kron(diag([rho*h,rho*h,rho*h^3/12, rho*h^3/12]),(R*L))/(dt*dt);

figure(2)
filename = 'test.gif';

p = plot(x, L*Ux(1*n+1:2*n));
ylim([-4.1*1e-05 4.1*1e-05])
xlabel('x');
ylabel('w');
drawnow
frame = getframe(2);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);

p.XData = x;

for t = 0:dt:0.5
    [U1] = next_step(t, A, Ux, Ux_1, mas, cT);
    Ux_1 = Ux;
    Ux = U1;
    figure(2)
    p.YData = L*U1(1*n+1:2*n);
    title("перемещение вдоль оси y, t = "+t)
    drawnow
    frame = getframe(2);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append',"DelayTime",0.01);
end

% figure;
% p = plot(x, u);
% xlabel('x');
% ylabel('u');
% titlestr = 'перемещения вдоль оси х';
% title(titlestr);
% p(1).LineWidth = 2;
% saveas (gcf, 'u', 'jpeg')





% xlswrite('res.xlsx',answertable, 'Table1')

% B = [12;12;12];
% C = [A,A, B,B;
%     A,B,B,A]
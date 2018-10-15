% MAIN DRIVER PROGRAM
% WRITTEN BY: DRAGOS STRUGAR, B17-05
% 15 OCT 2018 ONWARDS
% DIFFERENTIAL EQUATIONS COURSE
%
%
% TASK -- use numerical methods to approximate dy/dx = 3xy + xy^2
% TASK -- euler method, improved euler method, and runge kutta method
% TASK -- plot the solutions
% TASK -- give error graphs
%
%

% clears working space
close all;
clc;
clear;

% gives IVP - initial value problem, with step = 0.1
a = 0; b = 5.5; ya = 3; h = 0.1;
% as the problem has an asymptote, divide the plot, and problem into two
% parts: before and after the asymptote
asymptote = sqrt(2/3) * sqrt(log(2));

% t1 - vector, from initial x0 to asymptote, evenly distributed with step h
% t1 - vector, from asymptote to ending X, evenly distributed with step h
t1 = a:h:asymptote;
t2 = (asymptote+h):h:b;

% USE FILES:
% EulerMethod.m
% ImprovedEulerMethod.m
% RungeKutta.m
% Problem12.m specifies the equation dy/dx = 3xy + xy^2
% Can be included as ''

% EULER METHOD, before the asymptote and after
[x1eul,y1eul] = EulerMethod('Problem12', a, asymptote, ya, h);
[x2eul,y2eul] = EulerMethod('Problem12', asymptote+h, b, -15.2575, h); % -15.2575 = sol. at asymptote + h

% IMPROVED EULER METHOD, before the asymptote and after
[x1improved,y1improved] = ImprovedEulerMethod('Problem12', a, asymptote, ya, h);
[x2improved,y2improved] = ImprovedEulerMethod('Problem12', asymptote+h, b, -15.2575, h);

% RUNGE-KUTTA METHOD, before the asymptote and after
[x1rk, y1rk] = RungeKuttaMethod('Problem12', a, asymptote, ya, h);
[x2rk, y2rk] = RungeKuttaMethod('Problem12', asymptote+h, b, -15.2575, h);

% EXACT SOLUTION OF THE PROBLEM, before the asymptote and after
sq1 = x1improved.^2;
sq2 = x2improved.^2;
yi1 = -1*(3*exp(3/2*sq1))./(exp(3/2*sq1)-2);
yi2 = -1*(3*exp(3/2*sq2))./(exp(3/2*sq2)-2);

% PLOT RESULTS
hold on;

% plot euler
plot(t1, y1eul, 'ro', 'LineWidth', 2);
plot(t2, y2eul, 'ro', 'LineWidth', 2);

% plot improved euler
plot(t1, y1improved, 'bx', 'LineWidth', 2);
plot(t2, y2improved, 'bx', 'LineWidth', 2);

% plot rk
plot(t1, y1rk, 'y*', 'LineWidth', 2);
plot(t2, y2rk, 'y*', 'LineWidth', 2);

% plot exact solution
plot(t1, yi1, 'g--', 'LineWidth', 2);
plot(t2, yi2, 'g--', 'LineWidth', 2);

hold off;
% ADD A BIT OF EXTRA INFO TO THE PLOT
title('Comparison between Euler, Improved Euler and Runge-Kutta for dy/dx=3xy+xy^2 on [0,5.5]');
legend('euler', 'euler', 'improved euler', 'improved euler','runge-kutta','runge-kutta', 'exact', 'exact');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
box on;

saveas(gcf, 'results.png')

% PLOT ALL ERRORS
% MAKE another graph

% euler method errors
eul1err = abs(yi1 - y1eul);
eul2err = abs(yi2 - y2eul);

% TOTAL EULER ERROR
eul_err = sum(eul1err) + sum(eul2err);

% euler method errors
imp1err = abs(yi1 - y1improved);
imp2err = abs(yi2 - y2improved);

% TOTAL IMPROVED EULER ERROR
imp_err = sum(imp1err) + sum(imp2err);

% runge kutta errors
rk1err = abs(yi1 - y1rk);
rk2err = abs(yi2 - y2rk);

% TOTAL RUNGE KUTTA ERROR
rk_err = sum(rk1err) + sum(rk2err);

% PLOT ERROR RESULTS
figure();
hold on;
% euler
plot(t1, eul1err, 'ro', 'LineWidth', 2);
plot(t2, eul2err, 'ro', 'LineWidth', 2);

% improved euler
plot(t1, imp1err, 'bx', 'LineWidth', 2);
plot(t2, imp2err, 'bx', 'LineWidth', 2);

% runge kutta
plot(t1, rk1err, 'y*', 'LineWidth', 2);
plot(t2, rk2err, 'y*', 'LineWidth', 2);

% add more information
hold off;
title('Comparison between Euler, Improved Euler and Runge-Kutta ERRORS for dy/dx=3xy+xy^2 on [0,5.5]');
legend('euler erorr', 'euler erorr', 'improved euler erorr', 'improved euler erorr','runge-kutta erorr','runge-kutta erorr');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
box on;

% save results in error_results.png
saveas(gcf, 'error_results.png')
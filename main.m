% MAIN DRIVER PROGRAM
% WRITTEN BY: DRAGOS STRUGAR, B17-05
% 15 OCT 2018 ONWARDS
% DIFFERENTIAL EQUATIONS COURSE
close all;
clc;
clear;

a = input('Enter x0: ');
ya = input('Enter y(x0): ');
b = 5.5;
h = input('Enter step h: ');

% as problem has asymptote, divide
% find constant of integration and asymptote
c = (log(ya/(ya+3))/3) - (0.5*a^2);
if c == -0.5 * a^2, asymptote = NaN; else asymptote = [sqrt(-2*c), -sqrt(-2*c)];
end

% as function is even, has two symmetrical asymptotes
asymptote_exists = 0;
if ~isnan(asymptote), as1 = asymptote(2); as2 = asymptote(1); asymptote_exists=1;
end

% storage of steps to determine delta-x
t = []; t1 = []; t2 = []; t3 = []; params = []; case_no = 0;

% determine # of sections
if (b<=as1 || (a>=as1 && b<=as1) || a>=as2 || ~asymptote_exists)
    t = a:h:b; params = [params, [a,b,ya]];
    case_no = 1;
elseif (a<as1 && b<=as2)
    % lower asymptote
    t1 = a:h:as1; t2 = as1+h:h:as2;
    params = [[a,as1,ya];[as1+h,b,DESolution(as1+h,c)]];
    case_no = 2;
elseif (a>=as1 && b>as2)
    % only upper
    t1 = a:h:as2; t2 = as2+h:h:b;
    params = [[a,as2,ya];[as2+h,b,DESolution(as2+h,c)]];
    case_no = 2;
elseif (a<as1 && b>as2)
    % both
    t1 = a:h:as1; t2 = as1+h:h:as2; t3 = as2+h:h:b;
    params = [[a,as1,ya];[as1+h,as2,DESolution(as1+h,c)];[as2+h,b,DESolution(as2+h,c)]];
    case_no = 3;
end


hold on;
if case_no == 1 % 1 interval only, no discontinuities
    [~, re]  = EulerMethod('Problem12', params(1,1), params(1,2), params(1,3), h);
    [~, rie] = ImprovedEulerMethod('Problem12', params(1,1), params(1,2), params(1,3), h);
    [~, rke] = RungeKuttaMethod('Problem12', params(1,1), params(1,2), params(1,3), h);
    sq = t.^2;
    y1 = -(3*exp(3*c+3/2*sq))./(exp(3*c+3/2*sq)-1); % compute exact solution
    
    % compute w/ all numerical methods
    err_eul = abs(y1 - re);
    err_ieu = abs(y1 - rie);
    err_rke = abs(y1 - rke);
    
    plot(t, re, 'ro', 'LineWidth', 2);
    plot(t, rie, 'bx', 'LineWidth', 2);
    plot(t, rke, 'm*', 'LineWidth', 2);
    plot(t, y1, 'g--', 'LineWidth', 2);
else
    sq1 = t1.^2;
    sq2 = t2.^2;
    
    y1 = -(3*exp(3*c+3/2*sq1))./(exp(3*c+3/2*sq1)-1);
    y2 = -(3*exp(3*c+3/2*sq2))./(exp(3*c+3/2*sq2)-1);
    
    plot(t1, y1, 'g--', 'LineWidth', 2);
    plot(t2, y2, 'g--', 'LineWidth', 2);
    
    if case_no == 2 % 2 intervals         
        for i = 1:2 % compute w/ all numerical methods
            [xre, re]  = EulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            [xrie, rie] = ImprovedEulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            [xrke, rke] = RungeKuttaMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            
            if i == 1, err_1_eul = abs(y1 - re); err_1_ieu = abs(y1 - rie); err_1_rke = abs(y1 - rke); else err_2_eul = abs(y2 - re); err_2_ieu = abs(y2 - rie); err_2_rke = abs(y2 - rke); end
             
            plot(xre, re, 'ro', 'LineWidth', 2);
            plot(xrie, rie, 'bx', 'LineWidth', 2);
            plot(xrke, rke, 'm*', 'LineWidth', 2);            
        end
         
    else % 3 intervals
        sq3 = t3.^2;
        y3 = -(3*exp(3*c+3/2*sq3))./(exp(3*c+3/2*sq3)-1); 
        plot(t3, y3, 'g--', 'LineWidth', 2);
        
        for i = 1:3 % compute w/ all numerical methods
            [xre, re]  = EulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            [xrie, rie] = ImprovedEulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            [xrke, rke] = RungeKuttaMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
            
            if i == 1, err_1_eul = abs(y1 - re); err_1_ieu = abs(y1 - rie); err_1_rke = abs(y1 - rke);
            elseif i == 2, err_2_eul = abs(y2 - re); err_2_ieu = abs(y2 - rie); err_2_rke = abs(y2 - rke); 
            else err_3_eul = abs(y3 - re); err_3_ieu = abs(y3 - rie); err_3_rke = abs(y3 - rke); 
            end
            
            plot(xre, re, 'ro', 'LineWidth', 2);
            plot(xrie, rie, 'bx', 'LineWidth', 2);
            plot(xrke, rke, 'm*', 'LineWidth', 2);            
        end
        plot(t3, y3, 'g--', 'LineWidth', 2);
    end
end

hold off;

title('Comparison between Euler, Improved Euler and Runge-Kutta for dy/dx=3xy+xy^2 on [x0,5.5]');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
box on;

saveas(gcf, 'results.png')

figure();
hold on;
% euler
if case_no == 1
    plot(t, err_eul, 'ro', 'LineWidth', 2);
    plot(t, err_ieu, 'bx', 'LineWidth', 2);
    plot(t, err_rke, 'm*', 'LineWidth', 2);
else
    plot(t1, err_1_eul, 'ro', 'LineWidth', 2);
    plot(t1, err_1_ieu, 'bx', 'LineWidth', 2);
    plot(t1, err_1_rke, 'm*', 'LineWidth', 2);
    
    plot(t2, err_2_eul, 'ro', 'LineWidth', 2);
    plot(t2, err_2_ieu, 'bx', 'LineWidth', 2);
    plot(t2, err_2_rke, 'm*', 'LineWidth', 2);
    
    if case_no == 3
        plot(t3, err_3_eul, 'ro', 'LineWidth', 2);
        plot(t3, err_3_ieu, 'bx', 'LineWidth', 2);
        plot(t3, err_3_rke, 'm*', 'LineWidth', 2);
    end
end

hold off;
title('Comparison between Euler, Improved Euler and Runge-Kutta ERRORS for dy/dx=3xy+xy^2 on [x0,5.5]');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
box on;
% save results in error_results.png
saveas(gcf, 'error_results.png')
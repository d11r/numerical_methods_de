% MAIN DRIVER PROGRAM
% WRITTEN BY: DRAGOS STRUGAR, B17-05
% 15 OCT 2018 ONWARDS
% DIFFERENTIAL EQUATIONS COURSE
close all;
clc;
clear;

a = input('Enter x0: ');
ya = input('Enter y(x0): ');
b = input('Enter X: ');
h = input('Enter step h: ');

% as problem has asymptote, divide
% find constant of integration and asymptote
c = (log(ya/(ya+3))/3) - (0.5*a^2);
if c == -0.5 * a^2, asymptote = NaN; else asymptote = [sqrt(-2*c), -sqrt(-2*c)]; end

% as function is even, has two symmetrical asymptotes
asymptote_exists = 0;
if ~isnan(asymptote), as1 = asymptote(2); as2 = asymptote(1); asymptote_exists=1; end

% storage of steps to determine delta-x
[t, t1, t2, t3, params, case_no] = NumberOfSections(a, as1, as2, b, h, c, asymptote_exists, ya);

hold on;
if case_no == 1 % 1 interval only, no discontinuities
    [~, ~, ~, re, rie, rke] = GetAllSolutionsWithI(params, h, 1);
    
    sq = t.^2;
    y1 = -(3*exp(3*c+3/2*sq))./(exp(3*c+3/2*sq)-1);
  
    % compute w/ all numerical methods
    err_eul = abs(y1 - re);
    err_ieu = abs(y1 - rie);
    err_rke = abs(y1 - rke);
    
    plot(t, re, 'ro', 'LineWidth', 2);
    plot(t, rie, 'bx', 'LineWidth', 2);
    plot(t, rke, 'm*', 'LineWidth', 2);
    plot(t, y1, 'g--', 'LineWidth', 2);
else
    sq1 = t1.^2; sq2 = t2.^2;
    
    y1 = -(3*exp(3*c+3/2*sq1))./(exp(3*c+3/2*sq1)-1);
    y2 = -(3*exp(3*c+3/2*sq2))./(exp(3*c+3/2*sq2)-1);
    
    plot(t1, y1, 'g--', 'LineWidth', 2);
    plot(t2, y2, 'g--', 'LineWidth', 2);
    
    if case_no == 2 % 2 intervals           
        for i = 1:2 % compute w/ all numerical methods
            [xre, xrie, xrke, re, rie, rke] = GetAllSolutionsWithI(params, h, i);
            
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
            [xre, xrie, xrke, re, rie, rke] = GetAllSolutionsWithI(params, h, i);
            
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
saveas(gcf, 'error_results.png');

figure();
hold on;

% third graph
st = 0.001:0.001:0.1;
f_err_eul = zeros(1, length(st)); f_err_ie = zeros(1, length(st)); f_err_rk = zeros(1, length(st));

for it = 1:length(st);
    h = st(it);
    [t_tmp, t1_tmp, t2_tmp, t3_tmp, params_tmp, case_no_tmp] = NumberOfSections(a, as1, as2, b, h, c, asymptote_exists, ya);
    
    sq1_tmp = t1_tmp.^2; sq2_tmp = t2_tmp.^2;
    
    y1_tmp = -(3*exp(3*c+3/2*sq1_tmp))./(exp(3*c+3/2*sq1_tmp)-1);
    y2_tmp = -(3*exp(3*c+3/2*sq2_tmp))./(exp(3*c+3/2*sq2_tmp)-1);
           
    for i = 1:2
        [~, ~, ~, re_tmp, rie_tmp, rke_tmp] = GetAllSolutionsWithI(params_tmp, h, i);
        
        if i == 1
             f_err_eul(it) = f_err_eul(it) + sum(abs(y1_tmp - re_tmp));
             f_err_ie(it)  = f_err_ie(it)  + sum(abs(y1_tmp - rie_tmp));
             f_err_rk(it)  = f_err_rk(it)  + sum(abs(y1_tmp - rke_tmp));
        elseif i == 2
            f_err_eul(it) = f_err_eul(it) + sum(abs(y2_tmp - re_tmp));
            f_err_ie(it)  = f_err_ie(it)  + sum(abs(y2_tmp - rie_tmp));
            f_err_rk(it)  = f_err_rk(it)  + sum(abs(y2_tmp - rke_tmp));
        end
    end 
end

N = floor((b-a)/h);
plot(st, f_err_eul, 'r--', 'LineWidth', 2);
plot(st, f_err_ie, 'b-.', 'LineWidth', 2);
plot(st, f_err_rk, 'm-', 'LineWidth', 2);

hold off;
title('Error as f-on of N');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
box on;
saveas(gcf, 'error_results_N.png');
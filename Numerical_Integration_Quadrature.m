%% Numerical Integration: Gauss Quadrature vs. Clenshaw-Curtis
% This project compares Gauss quadrature and Clenshaw-Curtis quadrature on
% sis test integrands over [-1,1]. For each function, this project computes the
% absolute integration error as the number of nodes increases, plots
% the error curves, and estimates convergence patterns from fitted slopes.

%% Problem Setup
functions = {struct('f',@(x)x.^20,'label','$\int_{-1}^{1} x^{20}\,dx$','Interpreter','latex')
    struct('f',@(x) exp(x),'label','$\int_{-1}^{1} e^{x}\,dx$','Interpreter','latex')
    struct('f',@(x) exp(-(x.^2)),'label','$\int_{-1}^{1} e^{-x^2}\,dx$','Interpreter','latex')
    struct('f',@(x) 1./(1+16*x.^2),'label','$\int_{-1}^{1} \frac{1}{1+16x^2}\,dx$','Interpreter','latex')
    struct('f',@(x) exp(-(1./(x.^2))),'label','$\int_{-1}^{1} e^{-\frac{1}{x^2}}\,dx$','Interpreter','latex')
    struct('f',@(x) (abs(x)).^3,'label','$\int_{-1}^{1} |x|^3\,dx$','Interpreter','latex')};

n = 1:30;
a = -1;
b = 1;
poly = zeros(length(functions),2);

%% Error Computation and Convergence Plots
%%
% Followings are the fitted slopes that summarize the observed error's
% decay rate for each test function. These values are presented in the
% command window and can be compared across methods.
%% Test Case 1: polynomial
% Function 1: f(x) = x^20
err_f = zeros(30,2);
f = functions{1}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig1 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(n(1:9),log(abs(err_f(1:9,1))),1);
x = linspace(1,9,64);
y = exp(p_g(2) + p_g(1) .* x);
poly(1,1) = p_g(1);
semilogy(x,y,'r--','LineWidth',1);

p_c = polyfit(n(1:19),log(abs(err_f(1:19,2))),1);
x = linspace(1,19,64);
y = exp(p_c(2) + p_c(1) .* x);
poly(1,2) = p_c(1);
semilogy(x,y,'b--','LineWidth',1);

xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{1}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(1,1),poly(1,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);

%% 
% Since f(x) = x^20 is an analytical polynomial of degree 20, we may utilize the two
% methods' degree of precisions to evaluate from which degree this
% polynomial will be exactly integrated.
%
% Gauss Quadrature:
%
% Its degree of precision is 2n-1, and therefore we have 2n - 1 >= 20.
% Since n is an integer, n is greater than or equal to 11.
%
% Clenshaw-Curtis:
%
% Its degree of precision is n+1, and therefore we have n + 1 >= 20. So, we
% have that n is greater than or equal to 19.
%
% As observed from the plot, these 2 results appear as a sharp drop toward 
% machine epsilon at n = 9 for Gauss' and n = 19 for Clenshaw-Curtis'
%% Test case 2 and 3: exponential
% Function 2: f(x) = e^x
err_f = zeros(30,2);
f = functions{2}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig2 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(n(1:6),log(abs(err_f(1:6,1))),1);
x = linspace(1,6,64);
y = exp(p_g(2) + p_g(1) .* x);
poly(2,1) = p_g(1);
semilogy(x,y,'r--','LineWidth',1);

p_c = polyfit(n(1:12),log(abs(err_f(1:12,2))),1);
x = linspace(1,12,64);
y = exp(p_c(2) + p_c(1) .* x);
poly(2,2) = p_c(1);
semilogy(x,y,'b--','LineWidth',1);

xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{2}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(2,1),poly(2,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);
%%
% For reference, the convergence rates for the second function are:
%
% Gauss Quadrature: -5.7596
%
% Clenshaw-Curtis: -3.2287
%% 
% Function 3: f(x) = e^(-x^2)
err_f = zeros(30,2);
f = functions{3}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig3 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(n(1:11),log(abs(err_f(1:11,1))),1);
x = linspace(1,11,64);
y = exp(p_g(2) + p_g(1) .* x);
poly(3,1) = p_g(1);
semilogy(x,y,'r--','LineWidth',1);

p_c = polyfit(n(1:19),log(abs(err_f(1:19,2))),1);
x = linspace(1,19,64);
y = exp(p_c(2) + p_c(1) .* x);
poly(3,2) = p_c(1);
semilogy(x,y,'b--','LineWidth',1);

xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{3}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(3,1),poly(3,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);
%%
% For reference, the convergence rates for the third function are:
%
% Gauss Quadrature: -3.253
%
% Clenshaw-Curtis: -1.9565
%%
% Function 2 and 3: e^x and e^(-x^2):
%
% Since both functions are smooth and analytic on [-1,1], both quadrature
% methods present rapid error decay as n increases.
%
% Meanwhile, in the semilogy plot, the error curves are approximately linear
% at the convergence phase, which is consistent with these two functions'
% exponential convergence behavior.
% 
% Moreover, as observed, Gauss quadrature generally exhibits smaller errors
% than Clenshaw does for the same n, which is consistent with Gauss' higher
% degree of precision when n is greater than 2.

%% Test case 4 and 5: functions w/ holes
% Function 4: f(x) = 1/(16+x^2)
err_f = zeros(30,2);
f = functions{4}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig4 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(n(1:end),log(abs(err_f(1:end,1))),1);
x = linspace(1,30,64);
y = exp(p_g(2) + p_g(1) .* x);
poly(4,1) = p_g(1);
semilogy(x,y,'r--','LineWidth',1);

p_c = polyfit(n(1:end),log(abs(err_f(1:end,2))),1);
x = linspace(1,30,64);
y = exp(p_c(2) + p_c(1) .* x);
poly(4,2) = p_c(1);
semilogy(x,y,'b--','LineWidth',1);

xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{4}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(4,1),poly(4,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);
%%
% For reference, the convergence rates for the fourth function are:
%
% Gauss Quadrature: -0.4937
%
% Clenshaw-Curtis: -0.49247
%%
% Function 5: f(x) = e^(-1/x^2)
err_f = zeros(30,2);
f = functions{5}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig5 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(n(1:end),log(abs(err_f(1:end,1))),1);
x = linspace(1,30,64);
y = exp(p_g(2) + p_g(1) .* x);
poly(5,1) = p_g(1);
semilogy(x,y,'r--','LineWidth',1);

p_c = polyfit(n(1:end),log(abs(err_f(1:end,2))),1);
x = linspace(1,30,64);
y = exp(p_c(2) + p_c(1) .* x);
poly(5,2) = p_c(1);
semilogy(x,y,'b--','LineWidth',1);

xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{5}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(5,1),poly(5,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);
%%
% For reference, the convergence rates for the fifth functions are:
%
% Gauss Quadrature: -0.46199
%
% Clenshaw-Curtis: -0.48751
%%
% Function 4 and 5: 1/(1+16x^2) and e^(-1/x^2)
% 
% Since the poles of rational function 1/(1+16x^2) lie outside the interval, 
% it is analytic on [-1,1]. Therefore, both Gauss and Clenshaw methods show
% rapid and smooth linear error decay in the semilog plot as we should
% expect from its bounded condition O(rho^(-n)), where rho is e in this
% case.
% 
% In comparison, although the function e^(-1/x^2) is also smooth on [-1,1],
% it is not analytic at x=0 since x^2 is then 0 otherwise. Therefore,
% unlike purely analytic test functions on [-1,1], even though its
% quadrature error still decay rapidly, the convergence may not appear as a
% perfect straight line over a longer semilogy interval.
%% Test case 6: finitely-differentiable function
% Function 6: f(x) = |x|^3
err_f = zeros(30,2);
f = functions{6}.f;
I_True = integral(f,a,b);

for i = 1:30
    I_Gauss = Gaussian_Quadrature(f,i);
    I_Clenshaw = ClenshawCurtis_Quadrature(f,i);
    err_f1_Gauss = abs(I_True - I_Gauss) + eps;
    err_f1_Clenshaw = abs(I_True - I_Clenshaw) + eps;
    err_f(i,1) = err_f1_Gauss;
    err_f(i,2) = err_f1_Clenshaw;
end

fig6 = figure;
semilogy(n,err_f(:,1),'r.','MarkerSize',12);
hold on;
semilogy(n,err_f(:,2),'bo','MarkerSize',5);

p_g = polyfit(log(n(1:end)),log(abs(err_f(1:end,1))),1);
x = linspace(1,30,64);
y = exp(p_g(2)) .* x .^(p_g(1));
poly(6,1) = p_g(1);
loglog(x,y,'r--','LineWidth',1);

p_c = polyfit(log(n(1:end)),log(abs(err_f(1:end,2))),1);
x = linspace(1,30,64);
y = exp(p_c(2)) .* x .^(p_c(1));
poly(6,2) = p_c(1);
loglog(x,y,'b--','LineWidth',1);


xlabel('n');
ylabel('|I - In|');
title(['Absolute Errors against n for ', ...
    functions{6}.label],'Interpreter','latex');
legend('Gauss Quadrature Abs. Errors','Clenshaw-Curtis Abs.Errors','Convergence Rate of Gauss','Convergence Rate of Clenshaw');
hold off;

convergeData = table(poly(6,1),poly(6,2),'VariableNames',{'Slope by Gauss','Slope by Clenshaw-Curtis'});
disp(convergeData);
%%
% For reference, the convergence rates for the fifth function are:
%
% Gauss Quadrature: -3.3939 
%
% Clenshaw-Curtis: -4.1556
%%
% Function 6: |x|^3
%
% Since |x|^3 is second differentiable, we may expect O(n^(-2)) bounds the 
% errors. Accordingly, the error decays more slowly as n increases in the
% current log-log plot.
%% Fitting bimodality to the DNA Content.

% This does not check if the distribution is indeed bimdodal - there are
% ways to check bimodality of a distribution. You may want to ememd the
% code to your taste.

% clear all;
% close all;
% clc;

function [G1_mean, G2_mean, std1, std2, p] = bimodefit(filename)

A = load(filename);
x = A(:,1);
bimode = @(x,p,mu1,mu2,sigma1,sigma2) p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);

pStart = .5;
muStart = quantile(x,[.25 .75]);
sigmaStart = sqrt(var(x) - .25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];


lb = [0 0 0 0 0];
ub = [1 Inf Inf Inf Inf];
                      
statset('mlecustom')

options = statset('MaxIter',20000, 'MaxFunEvals',40000);
paramEsts = mle(x, 'pdf',bimode, 'start',start, 'lower',lb, 'upper',ub, 'options',options);
                      

%% Added on April 20, 2017.

mu1 = paramEsts(1,2);
sigma1 = paramEsts(1,4);

bimode = @(x,p,mu2,sigma2) p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);

pStart = paramEsts(1,1);
mu2start = 2*mu1;           % This is to be set by hand until you get the fit of your choice.
sigmaStart = sigma1;
start = [pStart mu2start sigmaStart];
lb = [0 0.99*mu2start 0];
ub = [1 1.3*mu2start ((sigma1/mu1)*mu2start)];  % Constraining CV for G2/M peak. (May22, 2017)
% ub = [1 1.3*mu2start Inf];
                      
statset('mlecustom')

options = statset('MaxIter',20000, 'MaxFunEvals',40000);
paramEsts = mle(x, 'pdf',bimode, 'start',start, 'lower',lb, 'upper',ub, 'options',options);

p = paramEsts(1,1);
G1_mean = mu1;
std1 = sigma1;
G2_mean = paramEsts(1,2);
std2 = paramEsts(1,3);

%% Comment the previous section and un-comment the following if G2 != 2*G1.

% p = paramEsts(1,1);
% G1_mean = paramEsts(1,2);
% std1 = paramEsts(1,4);
% G2_mean = paramEsts(1,3);
% std2 = paramEsts(1,5);
    

    f = fopen('G1_peak.dat','w');  
    fprintf(f,'%f\t%f\t%f\t%f\t%f', G1_mean, G2_mean, std1, std2, p);

end

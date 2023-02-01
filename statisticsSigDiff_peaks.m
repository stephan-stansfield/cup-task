% Calculate statistical significance of difference between model mean
% values.

clear all
close all
clc

%% Looking at mean of subject **median** VAFs
C = readcell('data analysis/peak ratio results internal models.xlsx');      % simulations for SfN 2022 poster

% Read all columns:
x = C(1,2:end);                     % x-axis labels
u = 100*cell2mat(C(2,2:end));       % means of median peak ratios
s = 100*cell2mat(C(3,2:end));       % standard deviations

% % Original input shaping
% x1 = 0.501; % Mean of subject medians
% s1 = 0.182; % Standard deviation of subject medians
% 
% % Multi-mode input shaping
% x2 = 0.679;
% s2 = 0.148;
% 
% % Rigid body internal model (B|[7,20], K|[0,250])
% x3 = 0.777;
% s3 = 0.0856;
% 
% % Submovements & impedance direct fitting
% x4 = 0.959;
% s4 = 0.010;
% 
% n1 = 11;
% n2 = 11;
% n3 = 11;
% n4 = 11;
n = length(u); % same number of subjects for all simulations

% Compare experimental mean of medians to 1
[h5,p5,ci,stats] = ttest(cell2mat(medianBySubject(2:12,simNum+1)),1,...
        'Tail',tail)
disp('Original IS versus multi-mode IS w/ feedforward')
sigdiff(u(1),u(7),s(1),s(7),n,n)

% Compare original IS versus each simplified internal model IS with
% feedforward
disp('Original IS versus slow-mode internal model w/ feedforward')
sigdiff(u(1),u(8),s(1),s(8),n,n)

disp('Original IS versus fast-mode internal model w/ feedforward')
sigdiff(u(1),u(9),s(1),s(9),n,n)

disp('Original IS versus no-impedance internal model w/ feedforward')
sigdiff(u(1),u(10),s(1),s(10),n,n)

disp('Original IS versus rigid-body internal model w/ feedforward')
sigdiff(u(1),u(11),s(1),s(11),n,n)

% Analyze significance between all simulation combinations and show in a
% table
sigArray = strings(n,n);

for i = 1:n
        sigArray(i) = sigdiff(u(i),u(j),s(i),s(j),n,n);
    end
end

%% Looking at mean of subject **median** VAFs
%{
% Numbers for RQE presentation. Last updated 2022-01-14.

% Original input shaping
x1 = 0.469; % Mean of subject medians
s1 = 0.210; % Standard deviation of subject medians

% Multi-mode input shaping
x2 = 0.609;
s2 = 0.200;

% Rigid body internal model
x3 = 0.734;
s3 = 0.061;

% Submovements & impedance direct fitting
x4 = 0.959;
s4 = 0.010;

n1 = 11;
n2 = 11;
n3 = 11;
n4 = 11;

disp('Original IS versus multi-mode IS')
sigdiff(x1,x2,s1,s2,n1,n2)

disp('Original IS versus simplified internal model IS')
sigdiff(x1,x3,s1,s3,n1,n3)

disp('Original IS versus direct submovement fitting')
sigdiff(x1,x4,s1,s4,n1,n4)

disp('Multi-mode IS versus simplified internal model IS')
sigdiff(x2,x3,s2,s3,n2,n3)

disp('Multi-mode IS versus direct submovement fitting')
sigdiff(x2,x4,s2,s4,n2,n4)

disp('Simplified internal model IS versus direct submovement fitting')
sigdiff(x3,x4,s3,s4,n3,n4)
%}


%{ 
%% Looking at mean of subject **mean** VAFs

% Original input shaping
x1 = 0.441; % Mean of subject means
s1 = 0.136; % Standard deviation of subject means

% Multi-mode input shaping
x2 = 0.562;
s2 = 0.144;

% Rigid body internal model
x3 = 0.692;
s3 = 0.060;

% Submovements & impedance direct fitting
x4 = 0.951;
s4 = 0.011;

n1 = 11;
n2 = 11;
n3 = 11;
n4 = 11;

disp('Original IS versus multi-mode IS')
sigdiff(x1,x2,s1,s2,n1,n2)

disp('Original IS versus simplified internal model IS')
sigdiff(x1,x3,s1,s3,n1,n3)

disp('Original IS versus direct submovement fitting')
sigdiff(x1,x4,s1,s4,n1,n4)

disp('Multi-mode IS versus simplified internal model IS')
sigdiff(x2,x3,s2,s3,n2,n3)

disp('Multi-mode IS versus direct submovement fitting')
sigdiff(x2,x4,s2,s4,n2,n4)

disp('Simplified internal model IS versus direct submovement fitting')
sigdiff(x3,x4,s3,s4,n3,n4)
%}

%% Previous numbers. Unknown date.

%{
% 4-impulse versus 2-impulse with impedance. Feedforward F. RMSE
%%%%
disp('4-impulse vs. 2-impulse w/ impedance. No timing or amplitude variation, Feedforward F')

x4 = 0.145;
x2 = 0.151;

s4 = 0.066;
s2 = 0.043;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. Timing variation, no amplitude variation, Feedforward F')

x4 = 0.085;
x2 = 0.106;

s4 = 0.043;
s2 = 0.032;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. No timing variation, amplitude variation, Feedforward F')

x4 = 0.127;
x2 = 0.157;

s4 = 0.056;
s2 = 0.051;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. Timing and amplitude variation, Feedforward F')

x4 = 0.076;
x2 = 0.104;

s4 = 0.040;
s2 = 0.041;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%% 4-impulse vs. 2-impulse w/ impedance. No Feedforward F. RMSE.
%%%%
disp('4-impulse vs. 2-impulse w/ impedance. No timing or amplitude variation, No Feedforward F')

x4 = 0.209;
x2 = 0.245;

s4 = 0.080;
s2 = 0.036;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. Timing variation, no amplitude variation, No Feedforward F')

x4 = 0.115;
x2 = 0.185;

s4 = 0.052;
s2 = 0.041;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. No timing variation, amplitude variation, No Feedforward F')

x4 = 0.179;
x2 = 0.193;

s4 = 0.060;
s2 = 0.047;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('4-impulse vs. 2-impulse w/ impedance. Timing and amplitude variation, No Feedforward F')

x4 = 0.100;
x2 = 0.168;

s4 = 0.047;
s2 = 0.044;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%% 4-impulse versus 2-impulse with impedance. Feedforward F. STIFFNESS (K)
%%%%
disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. No timing or amplitude variation, Feedforward F')

x4 = 520.827;
x2 = 132.232;

s4 = 172.640;
s2 = 108.011;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. Timing variation, no amplitude variation, Feedforward F')

x4 = 585.543;
x2 = 152.745;

s4 = 158.744;
s2 = 82.191;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. No timing variation, amplitude variation, Feedforward F')

x4 = 514.975;
x2 = 171.915;

s4 = 169.067;
s2 = 113.015;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. Timing and amplitude variation, Feedforward F')

x4 = 558.195;
x2 = 174.126;

s4 = 160.319;
s2 = 87.448;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%% 4-impulse vs. 2-impulse w/ impedance. No Feedforward F. STIFFNESS (K)
%%%%
disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. No timing or amplitude variation, No Feedforward F')

x4 = 678.881;
x2 = 452.861;

s4 = 142.650;
s2 = 286.674;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. Timing variation, no amplitude variation, No Feedforward F')

x4 = 639.431;
x2 = 436.829;

s4 = 182.957;
s2 = 180.116;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. No timing variation, amplitude variation, No Feedforward F')

x4 = 555.906;
x2 = 321.048;

s4 = 175.419;
s2 = 267.917;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%%%%

disp('Stiffness. 4-impulse vs. 2-impulse w/ impedance. Timing and amplitude variation, No Feedforward F')

x4 = 618.621;
x2 = 370.354;

s4 = 159.688;
s2 = 196.851;

n4 = 799;
n2 = 799;

sigdiff(x2,x4,s2,s4,n2,n4)

%% Feedforward vs. no feedfoward; 2-impulse w/ impedance model

%%%%
disp('Force comparison. No timing or amplitude variation.')

xnof = 0.245;
xf = 0.151;

snof = 0.036;
sf = 0.043;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. Timing variation, no amplitude variation.')

xnof = 0.185;
xf = 0.106;

snof = 0.041;
sf = 0.032;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. No timing variation, amplitude variation.')

xnof = 0.193;
xf = 0.157;

snof = 0.047;
sf = 0.051;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. Timing and amplitude variation.')

xnof = 0.168;
xf = 0.104;

snof = 0.044;
sf = 0.041;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%% Feedforward vs. no feedfoward; 4-impulse model

%%%%
disp('Force comparison. No timing or amplitude variation.')

xnof = 0.209;
xf = 0.145;

snof = 0.080;
sf = 0.066;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. Timing variation, no amplitude variation.')

xnof = 0.115;
xf = 0.085;

snof = 0.052;
sf = 0.043;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. No timing variation, amplitude variation.')

xnof = 0.179;
xf = 0.127;

snof = 0.060;
sf = 0.056;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)

%%%%
disp('Force comparison. Timing and amplitude variation.')

xnof = 0.100;
xf = 0.076;

snof = 0.047;
sf = 0.040;

nnof = 799;
nf = 799;

sigdiff(xnof,xf,snof,sf,nnof,nf)
%}

% function sigdiff(x1,x2,s1,s2,n1,n2)
%     
%     Z = abs((x1 - x2)/sqrt((s1^2)/n1 + (s2^2)/n2));
%     if Z > C1p
%         disp('1% Significant Difference')
%     elseif Z > C5p
%         disp('5% Significant Difference')
%     else
%         disp('Difference not significant')
%     end
% end

function significance = sigdiff(x1,x2,s1,s2,n1,n2)
    % Critical values (two-tailed, large sample size)
    C1p = 2.576;
    C2p = 2.326;
    C5p = 1.96;

    Z = abs((x1 - x2)/sqrt((s1^2)/n1 + (s2^2)/n2));
    if Z > C1p
%         disp('1% Significant Difference')
        significance = '**';
    elseif Z > C5p
%         disp('5% Significant Difference')
        significance = '*';
    else
%         disp('Difference not significant')
        significance = '-';
    end

end
% plot cumulative FPT for first passage on a 1D interval
ancfpt = dlmread('test.out');
%
semilogx(ancfpt(:,2),ancfpt(:,3))

%% plot empirical cumulative distrib from BD simulations
bddata = dlmread('test.bd.out');
ancfpt = dlmread('test.out');
% look only at those particles that absorbed to the left
bddata = bddata(find(bddata(:,2)==1),:);

[bdcfpt,tvals] = ecdf(bddata(:,3));

semilogx(tvals,bdcfpt,ancfpt(:,2),ancfpt(:,3))

xlabel('time')
ylabel('cumulative distribution function of exit times')

%% Compare BD sims with sampling from analytic function
andata = dlmread('test.out');
bddata = dlmread('test.bd.out');


[bdcfpt,bdtvals] = ecdf(bddata(:,3));
[ancfpt,antvals] = ecdf(andata(:,3));

semilogx(antvals,ancfpt,bdtvals,bdcfpt)
title('Compare BD and sampling from analytic function (X0=0.1)')
xlabel('Time')
ylabel('cumulative FPT')
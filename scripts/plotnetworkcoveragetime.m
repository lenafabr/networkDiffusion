%% load in mean time to cover all edges (averaged among particles)
dataWT = dlmread('results_mito/20180316_timeCoverEdges_WT.out');
dataM = dlmread('results_mito/20180316_timeCoverEdges_mutant.out');


%% 
loglog(dataWT(:,2),dataWT(:,3),'b.','MarkerSize',20)
hold all
loglog(dataM(:,2),dataM(:,3),'r.','MarkerSize',20)
xlist = logspace(0.5,2);
loglog(xlist,xlist.^2*7/12,'k','LineWidth',2)
%loglog(xlist,xlist.^2*6/12,'g','LineWidth',2)
hold off
set(gca,'FontSize',16)
xlabel('total network edge length')
ylabel('avg time to cover full network')
leg=legend('wild-type','mutant', '$t = \frac{7x^2}{12D}$')
xlim([4,100])
ylim([9,3000])
set(leg,'Interpreter','latex')

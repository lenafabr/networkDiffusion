
% load network structure
[nodepos,edgenodes] = loadnetworkstruct('../examples/example1.net');


% get FPT
data = dlmread('../examples/example1.out');
fptdata = data(:,2:end);
mfpt = mean(fptdata(:,2:end),1);

%% plot network, with nodes color-coded by FPT

%cmat = [0,0,1;1,0,0];
cmat = jet(4); % color matrix;
% fpt values corresponding to color matrix rows
fptcol = linspace(log10(min(mfpt(2:end))),log10(max(mfpt(2:end))),size(cmat,1));


% plot edges
for ec = 1:size(edgenodes,1)
    plot(nodepos(edgenodes(ec,:),1),nodepos(edgenodes(ec,:),2),'k')
    hold all
end

for nc = 1:size(nodepos,1)
    if (mfpt(nc)==0)
        cval = [0 0 0];
    else
        cval = interp1(fptcol, cmat, log10(mfpt(nc)));
    end
    plot(nodepos(nc,1),nodepos(nc,2),'.','MarkerSize',25,'Color',cval)
    hold all
end
axis equal
hold off

%% plot network over time, nodes colored by fraction of particles hit
cmat = jet(4); % color matrix;

% times to visualize
nt = 100;
tlist = linspace(0,200,nt);

clear M
for tc = 1:nt
    % plot edges
    for ec = 1:size(edgenodes,1)
        plot(nodepos(edgenodes(ec,:),1),nodepos(edgenodes(ec,:),2),'k')
        hold all
    end
    
    for nc = 1:size(nodepos,1)   
        % fraction of particles hit
        fhit = nnz(fptdata(:,nc)<= tlist(tc))/size(fptdata,1);
        cval = interp1(linspace(0,1,size(cmat,1)), cmat, fhit);
        
        plot(nodepos(nc,1),nodepos(nc,2),'.','MarkerSize',25,'Color',cval)
        hold all
    end
    axis equal
    hold off
    set(gca,'Visible','off')
    set(gcf,'Color','w')
    text(0,-1,sprintf('Time: %0.0f', tlist(tc)),'FontSize',16)
    drawnow
    M(tc) = getframe(gcf);
end

%% make movie of particles hitting network nodes
animation2movie(M,'../examples/example1_nodehit.avi',5)

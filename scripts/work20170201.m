
% load network structure
[nodepos,edgenodes] = loadnetworkstruct('../example1.net');


% get FPT
fptdata = dlmread('../test.out');
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

%% just plot network
for ec = 1:size(edgenodes,1)
    plot(nodepos(edgenodes(ec,:),1),nodepos(edgenodes(ec,:),2),'k')
    hold all
end

plot(nodepos(:,1),nodepos(:,2),'.','MarkerSize',25)
axis equal
hold off
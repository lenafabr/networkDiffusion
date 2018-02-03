function [nodepos,edgenodes] = loadnetworkstruct(fname)
% load network structure from file fname

%%
fid = fopen(fname);

tline = fgetl(fid);
clear nodepos
while (ischar(tline))
    % ignore comments and blank lines
    if (length(tline)==0 || tline(1)=='#'); tline = fgetl(fid); continue; end
    
    %disp(tline);
    
    % split line by whitespace
    words= strsplit(tline);
    
    if(strcmpi(words(1),'NODE'))
        % read in node position       
        dim = nnz(cellfun(@(i) length(i), words(2:end))>0)-1;
        nums= cellfun(@(i) str2num(i), words(2:dim+2));
        nodepos(nums(1),:) = nums(2:end);
    elseif(strcmpi(words(1),'EDGE'))
        nums= cellfun(@(i) str2num(i), words(2:4));
        edgenodes(nums(1),:) = nums(2:3);
    end
    tline = fgetl(fid);
end

fclose(fid);
end
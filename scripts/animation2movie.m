function animation2movie(M,mfile,fps)
%% convert animation M to an avi movie saved in mfile
if (exist('VideoWriter','file'))
    writerObj = VideoWriter(mfile);
    writerObj.FrameRate=fps;
    open(writerObj);
    for tc = 1:length(M)
        tc
        writeVideo(writerObj,M(tc));
    end
    close(writerObj);
else
    % use older avi recorder
    aviobj = avifile(mfile,'fps',fps);
    for tc = 1:length(M)
        aviobj = addframe(aviobj,M(tc));
    end
    aviobj = close(aviobj)
end

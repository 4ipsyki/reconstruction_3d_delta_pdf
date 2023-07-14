function savepng(fullpath,name,res)

    if nargin == 3
        res = ['-r',num2str(res)];
    else
        res = '-r500';
    end
    
    if ~strcmpi(fullpath(end),'/')
        if ~strcmpi(name(1),'/')
            fullpath = [fullpath,'/'];
        end
    else
        if strcmpi(name(1),'/')
            name = name(2:end);
        end
    end
    
    if strcmpi(name(end-3:end),'.png')
        name = [name,'.png'];
    end

    set(gcf,'InvertHardcopy','on','PaperPositionMode','auto')
    
    print([fullpath,name],'-opengl','-dpng',res)

end
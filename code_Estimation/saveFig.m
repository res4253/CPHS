function saveFig(ax,ext,options)%% ax=gca; ext='pdf','jpg','png',...
arguments
    ax
    ext
    options.name = [];
end

if isempty(options.name)
user_input = input('filename : ','s');
else
    user_input = options.name;
end
filename = sprintf('%s.%s',user_input,ext);

exportgraphics(ax,filename);
end

% 
% for i=1:5
%     figure(i)
%     ax=gca;
%     saveFig(ax,'pdf')
% end





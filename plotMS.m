function plotMS(x,y)
% plot MS deployment

for i = 1:length(x);
    xTemp = x(i,1);
    yTemp = y(i,1);
    plot(xTemp, yTemp, 'kx');  % Plot only MS points
%     text(xTemp,yTemp,[' ' num2str(i) ' '],'Color','k','Margin',0.1,...
%         'FontSize',8,'HorizontalAlignment','center',...
%         'BackgroundColor','w');  % Plot MS labels
%     axis equal, axis([0.15,0.85,0.15,0.85]), zoom on
    hold on   
end
end
function plotMS(x,y)
% Plot MS deployment

for i = 1:length(x);
    xTemp = x(i,1);
    yTemp = y(i,1);
    %plot(xTemp, yTemp, 'kx');  % Plot only MS points
    text(xTemp,yTemp,[' ' num2str(i) ' '],'Color','k','Margin',0.1,...
        'FontSize',9,'HorizontalAlignment','center',...
        'BackgroundColor','w');  % Plot MS labels
    axis equal, axis([0.1,0.9,0.1,0.9]), zoom on
    hold on   
end
end
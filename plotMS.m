function plotMS(x,y)
% Plot MS deployment

for i = 1:length(x);
    xTemp = x(i,1);
    yTemp = y(i,1);
<<<<<<< HEAD
    plot(xTemp, yTemp, 'kx');  % Plot only MS points
%     text(xTemp,yTemp,[' ' num2str(i) ' '],'Color','k','Margin',0.1,...
%         'FontSize',8,'HorizontalAlignment','center',...
%         'BackgroundColor','w');  % Plot MS labels
%     axis equal, axis([0.15,0.85,0.15,0.85]), zoom on
=======
    %plot(xTemp, yTemp, 'kx');  % Plot only MS points
    text(xTemp,yTemp,[' ' num2str(i) ' '],'Color','k','Margin',0.1,...
        'FontSize',9,'HorizontalAlignment','center',...
        'BackgroundColor','w');  % Plot MS labels
    axis equal, axis([0.1,0.9,0.1,0.9]), zoom on
>>>>>>> b6a16daeedc2cf94d6090c9faf1b4019fc36d8d3
    hold on   
end
end
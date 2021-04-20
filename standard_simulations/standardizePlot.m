function [] = standardizePlot(mygcf,mygca,printName)
    % This function standardizes the figure panels which only contain x and y axes
    scaleFactor = 5;

    ax = mygca;
    ax.XLabel.FontSize = 7.5*scaleFactor;
    ax.YLabel.FontSize = 7.5*scaleFactor;
    ax.FontSize = 7.5*scaleFactor;
    set(mygca,'LineWidth',3,'FontName','Helvetica');
    set(mygcf,'units','centimeters');
    posgcf = mygcf.Position;
    
    % % square
    set(gcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 5 5].*scaleFactor);
    set(gca,'Position',[0.15 0.17 0.7 0.7]);
    set(mygcf,'renderer','Painters')
    ax.Color = 'None';
    ax.XColor = 'black'; ax.YColor = 'black'; ax.ZColor = 'black';
    set(gcf,'PaperPositionMode','auto');
    
    if(nargin==3)
        print('-depsc','-r300',printName,'-loose','-painters');
        saveas(gcf,strcat(printName,'.png'));
    end

end


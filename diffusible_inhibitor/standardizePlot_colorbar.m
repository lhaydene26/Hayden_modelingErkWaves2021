function [] = standardizePlot_colorbar(mygcf,mygca,printName)
    % This function standardizes the figure panels which contain plots with colorbars

    scaleFactor = 5;

    ax = mygca;
    ax.XLabel.FontSize = 7.5*scaleFactor;
    ax.YLabel.FontSize = 7.5*scaleFactor;
    ax.FontSize = 7.5*scaleFactor;
    set(mygca,'LineWidth',3,'FontName','Helvetica');
    set(mygcf,'units','centimeters');
    posgcf = mygcf.Position;
    posgca = mygca.Position;

    set(gcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 5 5].*scaleFactor);
    set(gca,'Position',[0.15 0.17 0.6 0.6]);
    set(mygcf,'renderer','Painters')
    ax.Color = 'None';
    ax.XColor = 'black'; ax.YColor = 'black'; ax.ZColor = 'black';
    set(gcf,'PaperPositionMode','auto');
    
    if(nargin==3)
        print('-depsc','-r300',printName,'-loose','-painters');
        saveas(gcf,strcat(printName,'.png'));
    end

end


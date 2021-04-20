function [] = standardizePlot(mygcf,mygca,printName)
    scaleFactor = 5;

    ax = mygca;
    ax.XLabel.FontSize = 7.5*scaleFactor;
    ax.YLabel.FontSize = 7.5*scaleFactor;
    ax.FontSize = 7.5*scaleFactor;
    set(mygca,'LineWidth',3,'FontName','Helvetica');
    set(mygcf,'units','centimeters');
    posgcf = mygcf.Position;
    posgca = mygca.Position;
    % rectangle
    %set(mygcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 4.4 3.8].*scaleFactor);
    %set(mygca,'Position',[0.13 0.17 0.85 0.75]);
    % 
    % % square
    set(gcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 5 5].*scaleFactor);
    set(gca,'Position',[0.18 0.16 0.7 0.7]);
    set(mygcf,'renderer','Painters')
    ax.Color = 'None';
    ax.XColor = 'black'; ax.YColor = 'black'; ax.ZColor = 'black';
%     legend boxoff;
    set(gcf,'PaperPositionMode','auto');
    
    if(nargin==3)
        print('-depsc','-r300',printName,'-loose','-painters');
        saveas(gcf,strcat(printName,'.png'));
    end

end


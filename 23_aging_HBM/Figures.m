clear 
clc

%%   Age profile
clear g
data_path='F:\lifespan';
cd(data_path);
sub_info=xlsread('sub_information_deleted.xlsx');
age_info=sub_info(:,3);
max_age=max(age_info);
min_age=min(age_info);


g=gramm('x',age_info);
g.stat_bin('nbins',8,'edges',10:10:80);
g.set_names('x','Age(year)','y','Num');
% g.set_title('Age profile','FontSize',16);
g.set_color_options('map',[0.5 0.5 0.5]);
g.set_text_options('font','arial','base_size',8,'label_scaling',1.3);
g.axe_property('XLim',[5,83],'YLim',[0,130]);
g.draw();
export_path='E:\GS Coherence\result\figure';
cd(export_path);
name='Age profile';
saveas(gca,[name,'.fig']);
g.export('file_name',name,'export_path',export_path,'file_type','eps',...
         'width',8,'height',5,'units','centimeters');
g.export('file_name',name,'export_path',export_path,'file_type','pdf',...
         'width',8,'height',5,'units','centimeters');
g.export('file_name',name,'export_path',export_path,'file_type','jpg',...
         'width',8,'height',5,'units','centimeters');
     

%%   Fig.1a  Heatmap of Cohere
clear
clc

% colormap
path = 'E:\GS Coherence\Colormap\jet.csv';
cmap = csvread(path);
cmap = cmap(129:end,:);

data_path='E:\GS Coherence\result\data';
cd(data_path);
load('GS_topo_Cohere.mat');
a=mean(ROI_Cohere,3);
a=a';

figure;
imagesc(a);

colormap(cmap);
c=colorbar;
set(c,'ticks',[0.2,0.3,0.4,0.5,0.6],'Fontname','arial','Fontsize',8);

Y_tick=[];
set(gca,'YTick',Y_tick);
% Y_ticklabl={'Frontal lobe','Temporal lobe','Parietal lobe','Insular lobe', ...
%             'Limbic lobe','Occipital lobe','Subcortical nuclei' };
% set(gca,'YTickLabel',Y_ticklabl);

%   load xlabel:freq point
xlabel_path='E:\GS Coherence\result\cluster2.0\Cohere_freq';
cd(xlabel_path);
load('cluster_frequency.mat');
AAA=zeros(1,4);
for i=1:length(AAA)
    AAA(i)=sum(AA(1:i)); 
end

%
X_tick=[1,AAA];
set(gca,'XTick',X_tick);
xlabel('Frq(Hz)');
digits(5);      %   Control calculation accuracy (significant number)
BB=AAA*0.25/257;
BB=round(BB,3);

X_ticklabl=[0,BB];
set(gca,'XTickLabel',X_ticklabl);

set(gca,'TickLength',[0.01 0]);
set(gca,'tickdir','in');
% set(gca,'XGrid','on','GridColor',[0 0 0],'GridLineStyle','-.', ...
%          'GridAlpha',1);
% set(gca,'YGrid','on','GridColor',[0 0 0],'GridLineStyle','-.', ...
%          'GridAlpha',1);
set(gca,'XTickLabelRotation',45);
% box off
% title('fig 1a','FontSize',16,'FontWeight','bold' )

save_path='E:\GS Coherence\result\figure';
cd(save_path);
name='fig_1a';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);

%   line diagram above the heat chart
x=1:257;
b=mean(a,1);
figure
plot(x,b,'linewidth',2,'color',[0.29 0.66 0.86]);
set(gca,'ylim',[0.32,0.45],'ytick',(0.3:0.05:0.45));
set(gca,'linewidth',1.8);

set(gca,'XTick',X_tick);
set(gca,'XTickLabel',X_ticklabl);


box off

save_path='E:\GS Coherence\result\figure';
cd(save_path);
name='fig_1a_up';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);

%%   Fig.1d  The line chart 
clear 
clc

xlabel_path='E:\GS Coherence\result\allsubcohere';
cd(xlabel_path);
load('cluster_region.mat');

clear g
X=1:257;
color=cell(2,1);
color{1,1}='areas 1';
color{2,1}='areas 2';

g=gramm('x',X,'y',C,'color',color);
g.geom_line();
g.set_names('x','Freq(Hz)','y','Cohere','color','color legend','linestyle','# Cyl');
g.set_title('The line chart');

color_map=[0 180 225;234 234 0]/255;
g.set_color_options('map',color_map);
g.set_line_options('base_size',1.5,'styles', {'-' ':' '--' '-.'});
g.axe_property('xtick',[0,50,100,150,200,250], ...
               'xticklabel',char('0','0.05','0.10','0.15','0.20','0.25'), ...
               'yticklabel',0:0.1:0.8,...
               'Xlim',[0,250],...
               'ylim',[0.2,0.5]);
g.set_text_options('font','arial',...
                   'base_size',10,...
                   'label_scaling',1,...
                   'legend_scaling',1);
g.draw();

set(gcf,'PaperUnits','centimeters',...
        'PaperPosition',[0 0 8 5],...
        'PaperPositionMode','auto');

export_path='E:\GS Coherence\result\figure';
cd(export_path);
name='Fig 1d';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);


%%   Fig.2a  Heatmap of Correlation
clear
clc

data_path='E:\GS Coherence\result\data';
cd(data_path);
load('AGECORR.mat');
path = 'E:\GS Coherence\Colormap\jet.csv';
cmap = csvread(path);

figure;
imagesc(Z');
colormap(cmap);
c=colorbar;
set(c,'Fontname','arial','Fontsize',8);
set(c,'ticks',[-0.4:0.1:0.3]);

% mycmap = get(gcf,'Colormap');
% save('MyColormaps','mycmap');

Y_tick=[];
set(gca,'YTick',Y_tick);

%   load xlabel:freq point
xlabel_path='E:\GS Coherence\result\cluster2.0\CORR_freq';
cd(xlabel_path);
load('cluster_frequency.mat');
AAA=zeros(1,4);
for i=1:length(AAA)
    AAA(i)=sum(AA(1:i)); 
end
%
X_tick=[1,AAA];
set(gca,'XTick',X_tick);
xlabel('Frq(Hz)');
digits(5);
BB=AAA*0.25/257;
BB=round(BB,3);

X_ticklabl=[0,BB];
set(gca,'XTickLabel',X_ticklabl);

set(gca,'TickLength',[0.01 0]);
set(gca,'tickdir','in');
% set(gca,'XGrid','on','GridColor',[0 0 0],'GridLineStyle','-.', ...
%          'GridAlpha',1);
% set(gca,'YGrid','on','GridColor',[0 0 0],'GridLineStyle','-.', ...
%         'GridAlpha',1);
set(gca,'XTickLabelRotation',45);

% box off
% title('fig 1a','FontSize',16,'FontWeight','bold' )

save_path='E:\GS Coherence\result\figure';
cd(save_path);
name='fig_2a_Z';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);

%   line diagram above the heat chart
x=1:257;
b=mean(Z,2);
figure
plot(x,b,'linewidth',2,'color',[0.29 0.66 0.86]);
set(gca,'ylim',[-0.15,0.1],'ytick',(-0.15:0.05:0.1));
set(gca,'linewidth',1.8);
box off

save_path='E:\GS Coherence\result\figure';
cd(save_path);
name='fig_2a_up';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);

%%   Fig.2d  The line chart 
clear 
clc

xlabel_path='D:\yangchengxiao\DATA_lifespan2.0\result\cluster\No_Zscore\264\CORR\region';
cd(xlabel_path);
load('cluster_region.mat');

clear g
X=1:257;
color=cell(3,1);
color{1,1}='areas 1';
color{2,1}='areas 2';
color{3,1}='areas 3';

g=gramm('x',X,'y',C,'color',color);
g.geom_line();
g.set_names('x','Freq(Hz)','y','Cohere','color','color legend','linestyle','# Cyl');
g.set_title('The line chart');

color_map=[0 204 225;255 68 1;225 225 0]/255;
g.set_color_options('map',color_map);
g.set_line_options('base_size',2.5,'styles', {'-' ':' '--' '-.'});
g.axe_property('xtick',[0,50,100,150,200,250], ...
               'xticklabel',char('0','0.05','0.10','0.15','0.20','0.25'), ...
               'ytick',-0.3:0.1:8,...
               'yticklabel',-0.3:0.1:0.3,...
               'Xlim',[0,250],...
               'ylim',[-0.3,0.2]);
g.set_text_options('font','arial',...
                   'base_size',10,...
                   'label_scaling',1,...
                   'legend_scaling',1);
g.draw();

set(gcf,'PaperUnits','centimeters',...
        'PaperPosition',[0 0 8 5],...
        'PaperPositionMode','auto');

export_path='C:\Users\Lenovo\Desktop\lifespan2.0\draw';
cd(export_path);
name='Fig 2d';
saveas(gca,[name,'.fig']);
saveas(gca,[name,'.eps']);
saveas(gca,[name,'.pdf']);
saveas(gca,[name,'.jpg']);

%% bar figure for neurosynth
clear g
path = 'E:\GS Coherence\result\neurosynth\Select_function.xlsx';
function_r = xlsread(path);
function_r(:,[2,4,6,8]) = [];

label = {'auditory';'speech';'language';'music';'painful';
    'comprehension';'secondary somatosensory';'auditory visual'};

g=gramm('x',label,'y',function_r(:,1));
g.stat_summary('geom',{'bar'},'dodge',0);
g.draw();
g.stat_bin('nbins',8);


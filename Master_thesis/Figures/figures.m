clc
clear
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

%% intensity
for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean'))
    GS_topo_mean = GS_topo_mean';
    GS_topo_mean = [GS_topo_mean(end-14:end,:);GS_topo_mean(1:end-15,:)];
    plot(GS_topo_mean);
    ylabel({'Fisher''s Z'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['Intensity_slow-',num2str(i)],'png');
end

%% correlation
for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean'))
    GS_topo_mean = [GS_topo_mean(end-14:end,:);GS_topo_mean(1:end-15,:)];
    CVC = corr(GS_topo_mean,GS_topo_mean,'type','pearson');
    imagesc(CVC);
    ylabel({'Cosine phase angle (degree)'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'YTick',[0 90 180 270 360],'YTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    ylim(gca,[0 360]);
    colorbar(gca)
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['Correlation_slow-',num2str(i)],'png');
end



%% brain regions in -90 and 90

for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean'))
    GS_topo_mean = GS_topo_mean';
    GS_topo_mean = [GS_topo_mean(end-14:end,:);GS_topo_mean(1:end-15,:)];
    for iroi = 1:246
        [v,pos] = min(GS_topo_mean(60:120,iroi));
        pos_diff(iroi) = 90-(pos+59);
    end
    for iroi = 1:246
        [v,pos] = min(GS_topo_mean(240:300,iroi));
        pos_diff2(iroi) = 270-(pos+239);
    end
    cd('E:\HCP\Results_PL\figures')
    a_vec2nii(pos_diff,['posdiff90_slow-',num2str(i)]);
    a_vec2nii(pos_diff2,['posdiff270_slow-',num2str(i)]);
end

%% intensity phase
for i = 1:5
    load(fullfile(sub_paths{i},'intensity_phase.mat'))
    strength_phase = [strength_phase(end-14:end,:);strength_phase(1:end-15,:)];
    phase_difference = [phase_difference(end-14:end,:);phase_difference(1:end-15,:)];
    figure
    plot(strength_phase)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'intensity phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['intensity_phase_slow-',num2str(i)],'png');
    
    plot(phase_difference)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'intensity phase angle difference (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['phase_difference_slow-',num2str(i)],'png');
    
    CVC_phase = corr(strength_phase',strength_phase','type','pearson');
    imagesc(CVC_phase);
    ylabel({'Cosine phase angle (degree)'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'YTick',[0 90 180 270 360],'YTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    ylim(gca,[0 360]);
    colorbar(gca)
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['CVC_phase-',num2str(i)],'png');
    
end

%% brain map
sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

[BrainNetViewerPath,fileN,extn]=fileparts(which('BrainNet.m'));
SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];
cd('E:\HCP\Results_PL\figures\States')
for i=1:5
    
    vol_file = dir(fullfile(sub_paths{i},'state_anova'));
    vol_file(1:2) = [];
    F_path = fullfile(sub_paths{i},'state_anova',vol_file(1).name);
    CfgFile='E:\HCP\Results_PL\figures\States\Cfg_f.mat';    
    H_BrainNet=BrainNet_MapCfg(SurfFileName,F_path,CfgFile);
    JpgFile=strcat('F','_slow-',num2str(i));       
    eval(['print -r500 -dtiff -noui ''',JpgFile,''';']);
    for ivol = 2:7
        vol_path = fullfile(sub_paths{i},'state_anova',vol_file(ivol).name);
        CfgFile='E:\HCP\Results_PL\figures\States\Cfg_t.mat';       
        H_BrainNet=BrainNet_MapCfg(SurfFileName,vol_path,CfgFile);
        name = vol_file(ivol).name;
        name = name(3:4);
        JpgFile=strcat('T_',name,'_slow-',num2str(i));       
        eval(['print -r500 -dtiff -noui ''',JpgFile,''';']);
        
        close all
    end
close all
end

%% efficiency
for i = 1:5
    load(fullfile(sub_paths{i},'Eg.mat'))
    Eg = [Eg(end-14:end),Eg(1:end-15)];
    plot(Eg)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'global efficiency'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['Eg_slow-',num2str(i)],'png');
    
    load(fullfile(sub_paths{i},'El.mat'))
    El = El';
    El = [El(end-14:end,:);El(1:end-15,:)];
    plot(El)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'local efficiency'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['El_slow-',num2str(i)],'png');
end

%% fcmat
for i = 1:5
    load(fullfile(sub_paths{i},'fcmat_mean.mat'))
    fcmat_mean(fcmat_mean == Inf) = 1;
    fcmat_mean = fcmat_mean(1:128,1:128,:);
    fcmat_mean2d = reshape(fcmat_mean,[],360);
    CVC = corr(fcmat_mean2d,fcmat_mean2d,'type','pearson');
    
    imagesc(CVC);
    ylabel({'phase angle (degree)'});
    xlabel({'phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'YTick',[0 90 180 270 360],'yTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    ylim(gca,[0 360]);
    colorbar(gca)
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_PL\figures')
    saveas(gca,['Correlationfcmat_slow-',num2str(i)],'png');
end

%% task intensity
clear
clc
sub_paths{1} = 'E:\HCP\Results_task\slow-1';
sub_paths{2} = 'E:\HCP\Results_task\slow-2';
sub_paths{3} = 'E:\HCP\Results_task\slow-3';
sub_paths{4} = 'E:\HCP\Results_task\slow-4';
sub_paths{5} = 'E:\HCP\Results_task\slow-5';
for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean'))
    GS_topo_mean = GS_topo_mean';
    GS_topo_mean = [GS_topo_mean(end-14:end,:);GS_topo_mean(1:end-15,:)];
    plot(GS_topo_mean);
    ylabel({'Fisher''s Z'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    set(gca,'FontSize',15);
    mkdir('E:\HCP\Results_task\figures')
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['Intensity_slow-',num2str(i)],'png');
end
close all
%% correlation
for i = 1:5
    load(fullfile(sub_paths{i},'GS_topo_mean'))
    GS_topo_mean = [GS_topo_mean(end-14:end,:);GS_topo_mean(1:end-15,:)];
    CVC = corr(GS_topo_mean,GS_topo_mean,'type','pearson');
    imagesc(CVC);
    ylabel({'Cosine phase angle (degree)'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'YTick',[0 90 180 270 360],'YTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    ylim(gca,[0 360]);
    colorbar(gca)
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['Correlation_slow-',num2str(i)],'png');
end
close all

%% intensity phase
for i = 1:5
    load(fullfile(sub_paths{i},'intensity_phase.mat'))
    strength_phase = [strength_phase(end-14:end,:);strength_phase(1:end-15,:)];
    phase_difference = [phase_difference(end-14:end,:);phase_difference(1:end-15,:)];
    figure
    plot(strength_phase)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'intensity phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['intensity_phase_slow-',num2str(i)],'png');
    
    plot(phase_difference)
    xlabel({'GS phase angle (degree)'});
    xlim(gca,[0 360]);
    ylabel({'intensity phase angle difference (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['phase_difference_slow-',num2str(i)],'png');
    
    CVC_phase = corr(strength_phase',strength_phase','type','pearson');
    imagesc(CVC_phase);
    ylabel({'Cosine phase angle (degree)'});
    xlabel({'Cosine phase angle (degree)'});
    set(gca,'XTick',[0 90 180 270 360],'XTickLabel',...
        {'-180','-90','0','90','180'});
    set(gca,'YTick',[0 90 180 270 360],'YTickLabel',...
        {'-180','-90','0','90','180'});
    xlim(gca,[0 360]);
    ylim(gca,[0 360]);
    colorbar(gca)
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['CVC_phase-',num2str(i)],'png');    
end

close all

%% RT ANOVA
for i = 1:5
    plot(p_phase(:,i))
    ylabel({'P'});
    xlabel({'time lag(sec)'});
    set(gca,'XTick',[0:2:21],'XTickLabel',...
        {'0','1','2','3','4','5','6','7','8','9','10'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['phase_F',num2str(i)],'png');  
    
    plot(p_task(:,i))
    ylabel({'P'});
    xlabel({'time lag(sec)'});
    set(gca,'XTick',[0:2:21],'XTickLabel',...
        {'0','1','2','3','4','5','6','7','8','9','10'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['task_F',num2str(i)],'png');  
    
    plot(p_inter(:,i))
    ylabel({'P'});
    xlabel({'time lag(sec)'});
    set(gca,'XTick',[0:2:21],'XTickLabel',...
        {'0','1','2','3','4','5','6','7','8','9','10'});
    set(gca,'FontSize',15);
    cd('E:\HCP\Results_task\figures')
    saveas(gca,['inter_F',num2str(i)],'png');  
end

%% 4 states RT
clear
clc

load('E:\HCP\Results_task\slow-3\MRT0.mat')
bar(MRT0([1,3,5,7,9,11,13,15,17,19,21],:),'DisplayName','MRT0([1,3,5,7,9,11,13,15,17,19,21],:)')
ylabel({'RT (milisec)'});
xlabel({'time lag(sec)'});
ylim(gca,[700 900]);
set(gca,'XTick',[1:12],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10'});
set(gca,'FontSize',15);
legend({'Trough','Rise','Peak','Fall'},'NumColumns',2,'FontSize',10);
set(gcf,'unit','centimeters','position',[0 0 20 8]);
cd('E:\HCP\Results_task\figures')
print('back-0_slow-3','-dpng','-r0')


load('E:\HCP\Results_task\slow-3\MRT2.mat')
bar(MRT2([1,3,5,7,9,11,13,15,17,19,21],:),'DisplayName','MRT2([1,3,5,7,9,11,13,15,17,19,21],:)')
ylabel({'RT (milisec)'});
xlabel({'time lag(sec)'});
ylim(gca,[900 1100]);
set(gca,'XTick',[1:12],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10'});
set(gca,'FontSize',15);
legend({'Trough','Rise','Peak','Fall'},'NumColumns',2,'FontSize',10);
set(gcf,'unit','centimeters','position',[0 0 20 8]);
cd('E:\HCP\Results_task\figures')
print('back-2_slow-3','-dpng','-r0') 

load('E:\HCP\Results_task\slow-4\MRT0.mat')
bar(MRT0([1,3,5,7,9,11,13,15,17,19,21],:),'DisplayName','MRT0([1,3,5,7,9,11,13,15,17,19,21],:)')
ylabel({'RT (milisec)'});
xlabel({'time lag(sec)'});
ylim(gca,[700 900]);
set(gca,'XTick',[1:12],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10'});
set(gca,'FontSize',15);
legend({'Trough','Rise','Peak','Fall'},'NumColumns',2,'FontSize',10);
set(gcf,'unit','centimeters','position',[0 0 20 8]);
cd('E:\HCP\Results_task\figures')
print('back-0_slow-4','-dpng','-r0')

load('E:\HCP\Results_task\slow-4\MRT2.mat')
bar(MRT2([1,3,5,7,9,11,13,15,17,19,21],:),'DisplayName','MRT2([1,3,5,7,9,11,13,15,17,19,21],:)')
ylabel({'RT (milisec)'});
xlabel({'time lag(sec)'});
ylim(gca,[900 1100]);
set(gca,'XTick',[1:12],'XTickLabel',...
    {'0','1','2','3','4','5','6','7','8','9','10'});
set(gca,'FontSize',15);
legend({'Trough','Rise','Peak','Fall'},'NumColumns',2,'FontSize',10);
set(gcf,'unit','centimeters','position',[0 0 20 8]);
cd('E:\HCP\Results_task\figures')
print('back-2_slow-4','-dpng','-r0')

%% local efficiency brain map

sub_paths{1} = 'E:\HCP\Results_PL\slow-1';
sub_paths{2} = 'E:\HCP\Results_PL\slow-2';
sub_paths{3} = 'E:\HCP\Results_PL\slow-3';
sub_paths{4} = 'E:\HCP\Results_PL\slow-4';
sub_paths{5} = 'E:\HCP\Results_PL\slow-5';

[BrainNetViewerPath,fileN,extn]=fileparts(which('BrainNet.m'));
SurfFileName=[BrainNetViewerPath,filesep,'Data',filesep,'SurfTemplate',filesep,'BrainMesh_ICBM152_smoothed.nv'];

for i = 1:5
    cd(sub_paths{i});
    vol_file = dir(fullfile(sub_paths{i},'EL_*'));
    
    for ivol = 1:4
        vol_path = fullfile(sub_paths{i},vol_file(ivol).name);
        CfgFile='E:\HCP\Results_PL\figures\local_efficiency\Cfg.mat';       
        H_BrainNet=BrainNet_MapCfg(SurfFileName,vol_path,CfgFile);
        name = vol_file(ivol).name;
        name = name(1:7);
        JpgFile=strcat(name,'_slow-',num2str(i));     
        cd('E:\HCP\Results_PL\figures\local_efficiency')
        eval(['print -r500 -dtiff -noui ''',JpgFile,''';']);      
        close all
    end
close all
end


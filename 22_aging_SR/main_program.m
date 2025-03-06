
%% compute GS frac and oci
clear
clc
out_path = 'E:\lifespan\time_diff\original';
sub_path = 'F:\lifespan\FunRawARWSC';
mask_path = 'D:\program\matlab\Common\Mask\AAL90_erzhi';
age_path = 'F:\lifespan\sub_information_deleted.xlsx';
age = xlsread(age_path);
age = age(:,3);

cd(out_path)
load('GS.mat');
GSV = std(GS,0,2);
[r,p] = corrcoef(GSV,age);
GSM = mean(GS,2);
[r2,p2] = corrcoef(GSM,age);
range = (max(GS')-min(GS'));
[r3,p3] = corrcoef(range,age);

% pwelch para
TR=2;%TR
nfft=512;
window=hamming(32);
overlap=16;
%  len=(nfft/2+1)*2;
fs=1/TR;

for isub = 1:322        
    GS_mean = mean(GS(isub,:), 2);
    GS_base = bsxfun(@minus,GS(isub,:),GS_mean); % correct base line
    [GS_spec(:,isub),f] = pwelch(GS_base,window,overlap,nfft,fs);
end

% detrend
GS_spec2 = GS_spec((8:end-1),:); % delete 0.0017Hz and 0.2Hz
f2 = f(8:end-1);
% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
for isub = 1:322
    isub
    y1 = GS_spec2(:,isub);
    [xData, yData] = prepareCurveData(f2, y1 );
    [fitresult, gof] = fit( xData, yData, ft, opts );
    y2 = fitresult(f2);
    GS_a(isub) = fitresult.a; % cofficient a y = a.*x^b
    GS_b(isub) = fitresult.b; % cofficient b y = a.*x^b
    GS_trend(:,isub) = y2;
    GS_osci(:,isub) = y1-y2;
end

GS_spec3 = reshape(GS_spec2,[],1);
f3=[];
for i = 1:322
    f3 = [f3;f2];
end

mkdir(out_path);
cd(out_path);
save spec.mat GS_spec GS_spec2 GS_spec3 GS_trend GS_osci GS_a GS_b f f2 f3

%% display GS picture and correlation
clc
clear
age_path = 'F:\lifespan\sub_information_deleted.xlsx';
age = xlsread(age_path);
age = age(:,3);
load('E:\lifespan\time_diff\original\spec.mat');

% correlate the correlation between sd/range and age
sd = std(GS_osci,0,1);
[r,p] = corrcoef(sd,age);
range = max(GS_osci)-min(GS_osci);
[r2,p2] = corrcoef(range,age);


mean_spec = plot_age_gs(age,GS_spec,f);
figure
mean_trend = plot_age_gs(age,GS_trend,f2);
figure
mean_osci = plot_age_gs(age,GS_osci,f2);


mean_spec(1:8,:) = [];
mean_spec(end,:) = []; % delete unreliable frequency
for i = 1:249
    [r,p] = corrcoef(GS_osci(i,:),age); % correlation of age and oscillation
    osci_r(i) = r(2);
    osci_p(i) = p(2);
end

% fdr correct
P_fdr = mafdr(osci_p,'BHFDR', true);
FDR_pos = find(P_fdr > 0.05);
osci_r(FDR_pos) = 0;

cd('E:\lifespan\time_diff\original\figure');
save mean_spec.mat mean_spec mean_osci mean_trend osci_r

%% a and b age corr
clc
clear
path = 'E:\lifespan\time_diff\original\spec.mat';
load(path);
age_path = 'F:\lifespan\sub_information_deleted.xlsx';
age = xlsread(age_path);
age = age(:,3);
out_path = 'E:\lifespan\time_diff\original';


[r1,p1] = corrcoef(GS_a,age);
[r2,p2] = corrcoef(GS_b,age);

cd('E:\lifespan\time_diff\original\figure')
save GS_coef.mat GS_a GS_b

%% compute max value
clc
clear
age_path = 'F:\lifespan\sub_information_deleted.xlsx';
age = xlsread(age_path);
age = age(:,3);
load('E:\lifespan\time_diff\original\spec.mat');
out_path = 'E:\lifespan\time_diff\original';
mean_osci = mean(GS_osci,2);

range = max(GS_osci)-min(GS_osci);
[r,p] = corrcoef(range,age);

[pks,locs]= findpeaks(-mean_osci);
mean_spec = mean(GS_spec(2:500,:),2);
deriv = diff(mean_spec,2);
[pks,locs]= findpeaks(deriv);

range1 = [1,46];
range2 = [47,147];

normalized_data = mapminmax(GS_osci', 0, 1);
for i = 1:322
    x = normalized_data(i,:);
    [max1(i),locs1(i)] = max(x(range1(1):range1(2)));
    [max2(i),locs2(i)] = max(x(range2(1):range2(2)));
    SC1(i) = a_SC((x(range1(1):range1(2))),f2(range1(1):range1(2)));
    SC2(i) = a_SC((x(range2(1):range2(2))),f2(range2(1):range2(2)));
end


[r1,p1] = corrcoef(SC1,age);
[r2,p2] = corrcoef(SC2,age);
cd('E:\lifespan\time_diff\original\figure')
save SC.mat SC1 SC2 age

%% pair t-test for SC
clear
clc
cd('E:\lifespan\time_diff\original\figure')
load('SC_deconv.mat')
SC1_deconv = SC1;
SC2_deconv = SC2;
load('SC_nonC.mat')
SC1_nonC = SC1;
SC2_nonC = SC2;
load('SC.mat')

[h,p,ci,stats] = ttest(SC1,SC1_deconv);
[h2,p2,ci2,stats2] = ttest(SC2,SC2_deconv);
ES1 = stats.tstat/sqrt(321);
ES2 = stats2.tstat/sqrt(321);

[h,p,ci,stats] = ttest(SC1,SC1_nonC);
[h2,p2,ci2,stats2] = ttest(SC2,SC2_nonC);
ES1 = stats.tstat/sqrt(321);
ES2 = stats2.tstat/sqrt(321);

%% correlate GS and FD
clc
clear
sub_path = 'F:\lifespan\FunRawARWSC';
FD_path = 'F:\lifespan\RealignParameter';
sub_file = dir(sub_path);
age_path = 'F:\lifespan\sub_information_deleted.xlsx';
age = xlsread(age_path);
age = age(:,3);

sub_file(1:2) = [];
load('E:\GS Coherence\result\data\GS.mat');

for i = 1:322
    cd(fullfile(FD_path,sub_file(i).name));
    name = dir('FD_Power*');
    FD(i,:) = importdata(name.name);
    name = dir('rp*');
    hm(:,:,i) = importdata(name.name);
end

% pwelch para
TR=2;
nfft= 512;
window=hamming(32);
overlap=16;
fs=1/TR;


for isub = 1:322        
    [FD_spec(:,isub),f] = pwelch(FD(isub,:),window,overlap,nfft,fs);
end
FD_spec2 = FD_spec((8:end-1),:); % delete 0.0017Hz and 0.2Hz
f2 = f(8:end-1);
mean_spec = plot_age_gs(age,FD_spec2,f2);

for ifre = 1:249
    [r,p] = corrcoef(FD_spec(ifre,:),age);
    R_FD(ifre) = r(2);
    P_FD(ifre) = p(2);
end

FDR = a_FDR(R_FD,P_FD);

cd('E:\lifespan\time_diff\original\figure')
save meanFD_spec mean_spec
save R_FD R_FD P_FD FD_spec

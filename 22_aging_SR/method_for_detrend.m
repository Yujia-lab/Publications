load('E:\lifespan8.30\time_diff\original\osci.mat');

%% wavelet package method(MP)
s=B;ls=length(s);
[c,l]=wavedec(s,6,'db4');%С���ֽ⺯��������ָ����С������ָ����ͨ�ʹ�ͨ�˲������зֽ⡣����ֵCС���ֽ���������Ӧ�ļ�¼����L��
subplot(7,1,1);plot(s);
title('ԭʼ�źż������ع��ź�');
ylabel('s');
for i=1:6
  decmp=wrcoef('a',c,l,'db4',7-i);%�ع���Ϊ�ź�approximation ('type' = 'a') or detail ('type' = 'd') coefficients are reconstructed.
subplot(7,1,1+i);
plot(decmp);
ylabel(['a',num2str(7-i)]);
end

%% smoothness priors approach (SPA)

n = 0:1/40:(1024-1)/40;
data = sin(2*pi*0.1*n)+sin(2*pi*5*n);
data = data(:);
N = length(data);
lambda = 100; % parameter
I = speye(N);
D2 = spdiags(ones(N-2,1)*[1 -2 1],[0 1 2],N-2,N);
trend = inv(I+lambda^2*D2'*D2)*data;
detrenddata = data-trend;
subplot(211);
plot(n,data,'r',n,trend,'g');
title('the orginal data and trend');
legend('the orginal data','the trend');
xlim([0 5]);
subplot(212);
plot(n,detrenddata,'m')
title('the data after detrenging');
xlim([0 5]);

%% emd
[imf,residual] = emd(X);
%��Ԥ�ڰ汾��trialʱ����1-1.4��֮�������ƽ������1.2�룬��������trial������֣�����10%��
function gradCPT_v6l

clear;clc;
addpath('materials');
addpath('func');
rng('Shuffle');
global allimg img_rect w time0 time1 time2 time3 xc yc  mountainindex cityindex endmask_texid imgR hz bgc
global sti_rect sti_R
global response_type stitime
%%
%ͼ������

bgc=155; %������ɫ
imgR=250; %ͼƬ���ֵİ뾶
sti_R = 20;%�̼����ֵİ뾶
response_type = 1; %1�Գ��з�Ӧ��2��ɽ��Ӧ

%ʱ������
time0=1; %��ʼͼƬͣ��ʱ��
time1=1.5; %�׶�1�е�ͼƬ����ʱ����Ҳ����Ӧʱ��
time2=1;  %�׶�1�ķ���ʱ��
time3=0.5; %ʵʱ����ʱ�����ڶ��׶��г��֣�
para.change_time=1.2; %�����л�ʱ��  �����ٶȣ� ����ϰ��3�׶Σ���ʽʵ��
stitime = 0.5;      %�̼�����ʱ��
para.reswindow=[0.7,0.4];
%��Ӧ����  �� 0.7-1-0.4��Ҳ����ͼƬ��0.7͸�������ӵ�1�ټ��ٵ�0.4����һ��ʱ����
para.response_key=51;       %��'3#'��Ӧ
%���Եİ����ᱻʶ��Ϊ�Ե�ǰͼƬ�ķ�Ӧ������һ���Է�Ӧ�����������ж�

structure_time=360;%�ṹ��ʱ��
rest_time = 480;%��Ϣ̬ʱ��

part = 1;  
address = hex2dec('D100'); %���ڵ�ַ


%%
%��д������Ϣ
defultip={'','','','','0','0'};
sub=inputdlg({'���','����','�Ա�(����1��Ů��0)', ...
    '��������(��:2000-01-01)','����(��������1����������0)', ...
    'ʵ��׶�'},'',1, ...
    defultip);
subID=sub{1};
if sub{6}~='0'
    part=str2double(sub{6});
end

%��ⱻ�Ա���Ƿ��ظ�
expname=[pwd,'/record/subID-',subID] ;
if exist([expname,'.mat'],'file')
    error('The data file exists! Please enter a different subject ID.');
end

%��������
[w,wrect,hz,xc,yc]=init_screen(bgc,1);
fliptime=round(para.change_time*hz);

%�Դ���������
%��Ԥ��
trialpara.num=400;%�Դ���
trialpara.targetpe_percentage=0.1;%Ŀ����ָ���
trialpara.targetpe_rand=1;%Ŀ���Ƿ�������� 1��0��
trialpara.rand_duration=1;%trial����ʱ���Ƿ���� 1��0��
trialpara.duration=para.change_time;
trialpara.rand_range=0.2;%�����Χ
%��Ԥ��
% trialpara.num=400;%�Դ���
% trialpara.targetpe_percentage=0.1;%Ŀ����ָ���
% trialpara.targetpe_rand=0;%Ŀ���Ƿ�������� 1��0��
% trialpara.rand_duration=0;%trial����ʱ���Ƿ���� 1��0��
% trialpara.duration=para.change_time;
% trialpara.rand_range=0;%�����Χ

list{1}=createList(trialpara);
list{2}=createList(trialpara);
%���ֱ�ʾͼƬ��
%1-10��ʾ����
%11-20��ʾɽ
%û��������ͬͼƬ

%%

allimg=readallimg;

endmask=ones(imgR*2,imgR*2,3)*bgc;
endmask_texid=Screen('MakeTexture',w,endmask);

cityindex=find(strcmp({allimg.type},'city'));
mountainindex=find(strcmp({allimg.type},'mountain'));
img_rect=[xc-imgR,yc-imgR,xc+imgR,yc+imgR];
sti_rect=[xc-imgR-sti_R,yc-imgR-sti_R,xc+imgR+sti_R,yc+imgR+sti_R];

%%

if sub{6}=='0' || sub{6}=='1'
    %��Ϣ̬1
reststate(rest_time);
end

%% ��ϰ&�ṹ
if sub{6}=='0'|| sub{6}=='2'
  train(para,structure_time);
end

%%
%��ʽʵ��
if sub{6}=='0' || sub{6}=='3'
    for i=1
        %����1
        
        disp_intro_mri(w,'task1.png');
        
        disp_info=one_step(0,list{i},0,para);
        
        result{i}.disp_info=disp_info;
        result{i}.alllist=list{i};
        result{i}.subinfo=sub;
        final=resultana(result{i},fliptime,0);
        result{i}.final=final;
        save([expname,'.mat'],'result');
    end
end

if sub{6}=='0' || sub{6}=='4'
    %��Ϣ̬2
  reststate(rest_time);
end

if sub{6}=='0' || sub{6}=='5'
    for i=2
        %����2
        
        disp_intro_mri(w,'task2.png');
        
        disp_info=one_step(0,list{i},0,para);
        
        result{i}.disp_info=disp_info;
        result{i}.alllist=list{i};
        result{i}.frame_num=list{i}.frame_num;
        result{i}.subinfo=sub;
        final=resultana(result{i},fliptime,0);
        result{i}.final=final;
        save([expname,'.mat'],'result');
    end
end

%%
ListenChar(0);
ShowCursor;
sca;
return

%�������
%�Դ�
%alllist�������ֵ
%ͼƬ���
%��㣬���㣬�յ�ʱ��
%��Ӧʱ��  �����˷�Ӧ�����¼��Ӧʱ��  ��Ӧʱ��ͼƬ���Ϊ��
%��Ӧ������
%����ǰӦ�÷�Ӧ�ҷ�Ӧ�ˣ� ���¼1����ʾ����
%����ǰ���÷�Ӧ��û��Ӧ�����¼2����ȷ�ܳ�
%����ǰ���÷�Ӧ�ҷ�Ӧ�ˣ� ���¼3����ʾ�鱨
%����ǰ�÷�Ӧ��û��Ӧ�� ���¼4����ʾ©��
%��Ӧ�����ȷ�ԣ�1��ȷ 0����
end

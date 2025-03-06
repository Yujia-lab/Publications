%低预期版本，trial时间在1-1.4秒之间随机，平均还是1.2秒，不按键的trial随机出现，概率10%；
function gradCPT_v6l

clear;clc;
addpath('materials');
addpath('func');
rng('Shuffle');
global allimg img_rect w time0 time1 time2 time3 xc yc  mountainindex cityindex endmask_texid imgR hz bgc
global sti_rect sti_R
global response_type stitime
%%
%图像设置

bgc=155; %背景颜色
imgR=250; %图片呈现的半径
sti_R = 20;%刺激呈现的半径
response_type = 1; %1对城市反应，2对山反应

%时长设置
time0=1; %初始图片停留时长
time1=1.5; %阶段1中的图片呈现时长，也即反应时长
time2=1;  %阶段1的反馈时长
time3=0.5; %实时反馈时长（第二阶段中出现）
para.change_time=1.2; %单次切换时长  正常速度， 即练习第3阶段，正式实验
stitime = 0.5;      %刺激呈现时长
para.reswindow=[0.7,0.4];
%反应窗，  即 0.7-1-0.4，也就是图片从0.7透明度增加到1再减少到0.4的这一段时间内
para.response_key=51;       %按'3#'反应
%被试的按键会被识别为对当前图片的反应。并进一步对反应正误性做出判断

structure_time=360;%结构像时长
rest_time = 480;%静息态时长

part = 1;  
address = hex2dec('D100'); %并口地址


%%
%填写被试信息
defultip={'','','','','0','0'};
sub=inputdlg({'编号','姓名','性别(男填1，女填0)', ...
    '出生日期(例:2000-01-01)','利手(左利手填1，右利手填0)', ...
    '实验阶段'},'',1, ...
    defultip);
subID=sub{1};
if sub{6}~='0'
    part=str2double(sub{6});
end

%检测被试编号是否重复
expname=[pwd,'/record/subID-',subID] ;
if exist([expname,'.mat'],'file')
    error('The data file exists! Please enter a different subject ID.');
end

%启动窗口
[w,wrect,hz,xc,yc]=init_screen(bgc,1);
fliptime=round(para.change_time*hz);

%试次序列设置
%低预期
trialpara.num=400;%试次数
trialpara.targetpe_percentage=0.1;%目标出现概率
trialpara.targetpe_rand=1;%目标是否随机出现 1是0否
trialpara.rand_duration=1;%trial持续时间是否随机 1是0否
trialpara.duration=para.change_time;
trialpara.rand_range=0.2;%随机范围
%高预期
% trialpara.num=400;%试次数
% trialpara.targetpe_percentage=0.1;%目标出现概率
% trialpara.targetpe_rand=0;%目标是否随机出现 1是0否
% trialpara.rand_duration=0;%trial持续时间是否随机 1是0否
% trialpara.duration=para.change_time;
% trialpara.rand_range=0;%随机范围

list{1}=createList(trialpara);
list{2}=createList(trialpara);
%数字表示图片，
%1-10表示城市
%11-20表示山
%没有连续相同图片

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
    %静息态1
reststate(rest_time);
end

%% 练习&结构
if sub{6}=='0'|| sub{6}=='2'
  train(para,structure_time);
end

%%
%正式实验
if sub{6}=='0' || sub{6}=='3'
    for i=1
        %任务1
        
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
    %静息态2
  reststate(rest_time);
end

if sub{6}=='0' || sub{6}=='5'
    for i=2
        %任务2
        
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

%结果含义
%试次
%alllist里的索引值
%图片类别
%起点，最大点，终点时刻
%反应时间  若做了反应，则记录反应时，  反应时以图片起点为记
%反应结果类别
%若当前应该反应且反应了， 则记录1，表示击中
%若当前不该反应且没反应，则记录2，正确拒斥
%若当前不该反应且反应了， 则记录3，表示虚报
%若当前该反应且没反应， 则记录4，表示漏报
%反应结果正确性，1正确 0错误
end

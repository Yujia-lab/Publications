function final=resultana(result,fliptime,stilist)
global response_type

disp_info=result.disp_info;
temp1=find([disp_info(:,8)]==0); %寻找中途退出迹象
if ~isempty(temp1)
    if temp1(end)==1  %虚假退出迹象 ，不做处理
    else %否则
        disp_info(temp1(end):end,:)=[]; %删掉未完成的后半段数据
        temp1=find([disp_info(:,4)]==0);  %寻找最近一个完成的整段数据
        disp_info(temp1(end)+1:end,:)=[]; %整段可分析
    end
end

a=1:fliptime:size(disp_info,1);
a(end)=[];
b=(fliptime*2):fliptime:size(disp_info,1);
readlist=[a;b];

final=[];
for trial=1:size(readlist,2)
    final(trial).trial=trial;
    thistrial=disp_info(readlist(1,trial):readlist(2,trial),:);
    final(trial).thistrial=thistrial;
    
    temp1=thistrial(thistrial(:,6)==trial,:); %反应区域
    
    index_alllist=temp1(1,7); %in alllist
    final(trial).index_alllist=index_alllist;
    
    if index_alllist<=10
        final(trial).type='city';
        thisans=response_type; 
    else
        final(trial).type='mountain';
        thisans=~response_type;  
    end
    
    final(trial).startT=thistrial(1,8);
    final(trial).maxT=thistrial(fliptime,8);
    final(trial).endT=thistrial(end,8);
    
    if any(temp1(:,9)) %第9列出现了响应，表示被试反应了
        if thisans==1  %如果应当反应，则
            signal_judge=1;  %击中
            judge=1;
        else
            signal_judge=3; %虚报
            judge=0;
        end
        RT=temp1(find(temp1(:,9)),8)-thistrial(1,8);
        
    else  %被试没反应
        if thisans==1 %应当反应
            signal_judge=4; %漏报
            judge=0;
        else
            signal_judge=2; %正确拒斥
            judge=1;
        end
        RT=[];
    end
    
    final(trial).RT=RT;
    final(trial).signal_judge=signal_judge;
    final(trial).judge=judge;
    if stilist==0
    else
    final(trial).stilist=stilist(trial);
    end
end

%记录，试次
%该试次图片的类别
%呈现起点，大点，终点时刻
%反应情况  1表示反应，0表示未反应
%反应时间  若反应，则记录反应时，  反应时以图片起点为记
%反应结果类别
%若当前应该反应且反应了， 则记录1，表示击中
%若当前不该反应且没反应，则记录2，正确拒斥
%若当前不该反应且反应了， 则记录3，表示虚报
%若当前该反应且没反应， 则记录4，表示漏报



%第一列，先绘制的图片句柄
%第二列，后绘制的图片句柄     变量的初始数行的第二列，其实是init的句柄
%第三列，先绘制的图片的透明度
%第四列，后绘制的图片的透明度
%第五列，0表示未到反应窗口，1表示到达反应窗口
%第六列，数字表示当前处于第几个图片试次
%第7列， 表示当前图片试次的图片索引值  也即alllist的信息，数字大小表示了图片的类别是山或者城市

%第8列，在呈现时候记录了每次呈现的时刻
%第9列， 反应窗内的初次按键，1表示按键
%第10列， 原始按键信息

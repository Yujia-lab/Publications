function final=resultana(result,fliptime,stilist)
global response_type

disp_info=result.disp_info;
temp1=find([disp_info(:,8)]==0); %Ѱ����;�˳�����
if ~isempty(temp1)
    if temp1(end)==1  %����˳����� ����������
    else %����
        disp_info(temp1(end):end,:)=[]; %ɾ��δ��ɵĺ�������
        temp1=find([disp_info(:,4)]==0);  %Ѱ�����һ����ɵ���������
        disp_info(temp1(end)+1:end,:)=[]; %���οɷ���
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
    
    temp1=thistrial(thistrial(:,6)==trial,:); %��Ӧ����
    
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
    
    if any(temp1(:,9)) %��9�г�������Ӧ����ʾ���Է�Ӧ��
        if thisans==1  %���Ӧ����Ӧ����
            signal_judge=1;  %����
            judge=1;
        else
            signal_judge=3; %�鱨
            judge=0;
        end
        RT=temp1(find(temp1(:,9)),8)-thistrial(1,8);
        
    else  %����û��Ӧ
        if thisans==1 %Ӧ����Ӧ
            signal_judge=4; %©��
            judge=0;
        else
            signal_judge=2; %��ȷ�ܳ�
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

%��¼���Դ�
%���Դ�ͼƬ�����
%������㣬��㣬�յ�ʱ��
%��Ӧ���  1��ʾ��Ӧ��0��ʾδ��Ӧ
%��Ӧʱ��  ����Ӧ�����¼��Ӧʱ��  ��Ӧʱ��ͼƬ���Ϊ��
%��Ӧ������
%����ǰӦ�÷�Ӧ�ҷ�Ӧ�ˣ� ���¼1����ʾ����
%����ǰ���÷�Ӧ��û��Ӧ�����¼2����ȷ�ܳ�
%����ǰ���÷�Ӧ�ҷ�Ӧ�ˣ� ���¼3����ʾ�鱨
%����ǰ�÷�Ӧ��û��Ӧ�� ���¼4����ʾ©��



%��һ�У��Ȼ��Ƶ�ͼƬ���
%�ڶ��У�����Ƶ�ͼƬ���     �����ĳ�ʼ���еĵڶ��У���ʵ��init�ľ��
%�����У��Ȼ��Ƶ�ͼƬ��͸����
%�����У�����Ƶ�ͼƬ��͸����
%�����У�0��ʾδ����Ӧ���ڣ�1��ʾ���ﷴӦ����
%�����У����ֱ�ʾ��ǰ���ڵڼ���ͼƬ�Դ�
%��7�У� ��ʾ��ǰͼƬ�Դε�ͼƬ����ֵ  Ҳ��alllist����Ϣ�����ִ�С��ʾ��ͼƬ�������ɽ���߳���

%��8�У��ڳ���ʱ���¼��ÿ�γ��ֵ�ʱ��
%��9�У� ��Ӧ���ڵĳ��ΰ�����1��ʾ����
%��10�У� ԭʼ������Ϣ

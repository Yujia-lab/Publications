%��������
clear,clc
%�Դ���������
trail_num = 480/1.2;

% ����ɽ����
m_num=0.1*trail_num;
m_id_rand = floor(4*rand(m_num,1))';
m_id_rand(0.5*m_num:end)=-m_id_rand(0.5*m_num:end);
m_idrand_s=sum(m_id_rand);
m_id_rand=Shuffle(m_id_rand);
m_id=(5:10:trail_num);
m_id=m_id+m_id_rand;
% ����ʵ������
list{1}=ceil(rand(1,trail_num)*10);
list{1}=getalllist(list{1});
list{1}(m_id)=list{1}(m_id)+10;
list{2}=list{1};
save triallist list
%���ֱ�ʾͼƬ��
%1-10��ʾ����
%11-20��ʾɽ
%û��������ͬͼƬ
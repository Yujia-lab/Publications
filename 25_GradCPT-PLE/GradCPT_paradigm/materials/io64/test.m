config_io;
global cogent;
if( cogent.io.status ~= 0 )
   error('inp/outp installation failed');
end

address = hex2dec('0FF8');
%��������Ϊ���ڵ�ַ�� ��Ҫ���ݸ��˵��Ե�ʵ�ʵ�ַ��ӡ� ���豸���������ҵ����ڣ��鿴���ԡ���Դ ���ɿ�����

%�·��������ÿ�뷢��һ��mark�� markֵ��1��ʼ������ӵ�100
for i=1:100;
WaitSecs(1);

outp(address,i);
WaitSecs(0.004);
outp(address,0);

disp(i);

end


%�ű����к� ��ctrl+c���жϳ���

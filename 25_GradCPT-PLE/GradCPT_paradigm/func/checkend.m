function checkend

[~,~,kc]=KbCheck;
if kc(27)
    ListenChar(0);
    ShowCursor;
    sca;
    error('�û���ֹ����');
end
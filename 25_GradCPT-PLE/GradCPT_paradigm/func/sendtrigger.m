function sendtrigger(mark)
global address 
outp(address,mark);
WaitSecs(0.004);
outp(address,0);

function alllist=getalllist(alllist)


while 1
    alllist=Shuffle(alllist);
    
    while 1
        needrep=0;
        indexrep=find(diff(alllist)==0);
        if isempty(indexrep)
            break
        end
        list1=alllist(1:indexrep(1));
        list2=alllist(indexrep(1)+1:end);      
        
        repcount=0;
        while 1
            repcount=repcount+1;
            if repcount>100
                needrep=1;
                break;
            end
            list2=Shuffle(list2);
            if list2(1)~=list1(end)
                break;
            end
        end
        if needrep==1
            break;
        end
        alllist=[list1,list2];
        
    end
    if needrep==0
        break;
    end
end
%alllist ,数字表示图片， 1-10表示城市的1-10
%11-20表示山的1-10
%没有连续相同图片

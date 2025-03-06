function mat = a_cell2mat(cell)

for i = 1:length(cell)
    data = cell{i};
    num_pos = isstrprop(data,'digit');
    num_str = data(num_pos);
    data_num = str2num(num_str);
    if isempty(data_num) == 1;
        data_num = -1;
    end
    mat(i,:) = data_num;
end
    
    
function data_array = readdata(filename)
% three dimension data pointsnum*dimension*group
fileID = fopen(filename);

num_read = regexp(filename,'\d*','Match');
data_dimension = zeros(1, 3);
for m = 1:3
    data_dimension(m) = str2double(num_read{m});
end
data_array = zeros(data_dimension(1), data_dimension(2), data_dimension(3));
tline = fgetl(fileID);
m = 0;
while ischar(tline)
    
    if tline == '#'
        m = m+1;
        n = 1;
    else
        data_array(n, :, m) = str2num(tline);
        n = n+1;
    end
    tline = fgetl(fileID);
end
    
end
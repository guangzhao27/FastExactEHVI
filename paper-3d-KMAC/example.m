clear all; close all;
pf = readdata('../data/ran.1000pts.3d.10');
ref = [0 0 0];
mu = [10 10 10];
s = 2.5*[1 1 1];
temp = [mu s];

for i=1:10
    main(pf(:, :, i),ref,temp)
% a = main_v10(pf,ref,temp)
end


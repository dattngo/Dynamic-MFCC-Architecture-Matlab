clc;
N=4000;%Number of Sample 
input_data = random('Normal',0,1,N,1);
save('C:\Users\G551\Desktop\MFCC\Input/input_data.mat','input_data');

% % % % Create INPUT_MFCC by floating point format % % % %

h = fopen('C:\Users\G551\Desktop\MFCC\Input/input_data.txt','w');
in_data = fpt2bin(input_data);
[p,q] =size(in_data)
for k=1:(p-1)
    fprintf(h,'%s\n',in_data(k,:));
end
fprintf(h,'%s',in_data(p,:));
fclose(h);
display('ok');




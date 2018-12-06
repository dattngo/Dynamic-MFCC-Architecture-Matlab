clc;
h = fopen('C:\Users\G551\Desktop\MFCC\Output_Hardware/output_file_log.txt','r'); %% output_file_log chinh la energy, chu khong phai ket qua LOG sau mel.
nextline = '';
str='';
while ischar(nextline)
    nextline = fgetl(h);
    if ischar(nextline)
        str = [str;nextline];
    end
end
[p,q] = size(str);
mem_size = p;
% hex >>>> bin
b_hang = '';
b_cot  = '';
for i=1: mem_size
    for j=1:8
        switch str(i,j)
            case '0'
                b_hang = [b_hang,'0','0','0','0'];
            case '1'
                b_hang = [b_hang,'0','0','0','1'];
            case '2'
                b_hang = [b_hang,'0','0','1','0'];
            case '3'
                b_hang = [b_hang,'0','0','1','1'];
            case '4'
                b_hang = [b_hang,'0','1','0','0'];
            case '5'
                b_hang = [b_hang,'0','1','0','1'];
            case '6'
                b_hang = [b_hang,'0','1','1','0'];
            case '7'
                b_hang = [b_hang,'0','1','1','1'];
            case '8'
                b_hang = [b_hang,'1','0','0','0'];
            case '9'
                b_hang = [b_hang,'1','0','0','1'];
            case 'a'
                b_hang = [b_hang,'1','0','1','0'];
            case 'b'
                b_hang = [b_hang,'1','0','1','1'];
            case 'c'
                b_hang = [b_hang,'1','1','0','0'];
            case 'd'
                b_hang = [b_hang,'1','1','0','1'];
            case 'e'
                b_hang = [b_hang,'1','1','1','0'];
            case 'f'
                b_hang = [b_hang,'1','1','1','1'];
        end
    end
    b_cot = [b_cot;b_hang];
    b_hang = '';
end
b = b_cot;
%%%%%%%%%%%%%% bin >>> floating_point ( ftp)
ftp = [];
for i=1:mem_size
    s = bin2dec(b(i,1))*(-2) + 1;
    e = bin2dec(b(i,2:9)) - 127;
    f = 0;
    for j=10:32
      f = bin2dec(b(i,j)) * 2^(9-j) + f;
    end
    f = f + 1;
    f_p = s * f * 2^e;
    ftp = [ftp;f_p];
end
Energy_hardware = ftp;

% % % % % % % % % % % % % % % % % % % % 
path = 'C:\Users\G551\Desktop\MFCC\Output_Matlab/energy_out_data.mat';
load(path);

Error_Energy = Energy_hardware (1 : size (energy_out_data)) - energy_out_data;
% [Max_Error_Energy, I] = max (abs(Error_Energy))
% Saisotuongdoi = abs (Error_Energy (I)/ energy_out_data (I)*100)
Max_Error_Energy = max (abs(Error_Energy))
Max_Energy = max (abs(energy_out_data))
Saisotuongdoi = abs ((Max_Error_Energy/Max_Energy)*100)
fprintf('==========Relative error of Energy result is %f %% ',Saisotuongdoi);


%==========================================================================

h = fopen('C:\Users\G551\Desktop\MFCC\Output_Hardware/output_file_cep.txt','r');
nextline = '';
str='';
while ischar(nextline)
    nextline = fgetl(h);
    if ischar(nextline)
        str = [str;nextline];
    end
end
[p,q] = size(str);
mem_size = p;
% hex >>>> bin
b_hang = '';
b_cot  = '';
for i=1: mem_size
    for j=1:8
        switch str(i,j)
            case '0'
                b_hang = [b_hang,'0','0','0','0'];
            case '1'
                b_hang = [b_hang,'0','0','0','1'];
            case '2'
                b_hang = [b_hang,'0','0','1','0'];
            case '3'
                b_hang = [b_hang,'0','0','1','1'];
            case '4'
                b_hang = [b_hang,'0','1','0','0'];
            case '5'
                b_hang = [b_hang,'0','1','0','1'];
            case '6'
                b_hang = [b_hang,'0','1','1','0'];
            case '7'
                b_hang = [b_hang,'0','1','1','1'];
            case '8'
                b_hang = [b_hang,'1','0','0','0'];
            case '9'
                b_hang = [b_hang,'1','0','0','1'];
            case 'a'
                b_hang = [b_hang,'1','0','1','0'];
            case 'b'
                b_hang = [b_hang,'1','0','1','1'];
            case 'c'
                b_hang = [b_hang,'1','1','0','0'];
            case 'd'
                b_hang = [b_hang,'1','1','0','1'];
            case 'e'
                b_hang = [b_hang,'1','1','1','0'];
            case 'f'
                b_hang = [b_hang,'1','1','1','1'];
        end
    end
    b_cot = [b_cot;b_hang];
    b_hang = '';
end
b = b_cot;
%%%%%%%%%%%%%% bin >>> floating_point ( ftp)
ftp = [];
for i=1:mem_size
    s = bin2dec(b(i,1))*(-2) + 1;
    e = bin2dec(b(i,2:9)) - 127;
    f = 0;
    for j=10:32
      f = bin2dec(b(i,j)) * 2^(9-j) + f;
    end
    f = f + 1;
    f_p = s * f * 2^e;
    ftp = [ftp;f_p];
end
cepstral_data_hardware = ftp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'C:\Users\G551\Desktop\MFCC\Output_Matlab/cepstral_out_data.mat';
load(path);
Error_Cepstral = cepstral_out_data - cepstral_data_hardware (1 : size (cepstral_out_data));
% [Max_Error_Cepstral, I] = max (abs(Error_Cepstral))
% Saisotuongdoi = Error_Cepstral (I)/ cepstral_out_data (I)*100

Max_Error_Cepstral = max (abs(Error_Cepstral))
Max_Cepstral = max (abs(cepstral_out_data))
Saisotuongdoi = abs ((Max_Error_Cepstral/Max_Cepstral)*100)
fprintf('==========Relative error of Cepstral result is %f %% ',Saisotuongdoi);

%==========================================================================

h = fopen('C:\Users\G551\Desktop\MFCC\Output_Hardware/output_file_delta.txt','r');
nextline = '';
str='';
while ischar(nextline)
    nextline = fgetl(h);
    if ischar(nextline)
        str = [str;nextline];
    end
end
[p,q] = size(str);
mem_size = p;
% hex >>>> bin
b_hang = '';
b_cot  = '';
for i=1: mem_size
    for j=1:8
        switch str(i,j)
            case '0'
                b_hang = [b_hang,'0','0','0','0'];
            case '1'
                b_hang = [b_hang,'0','0','0','1'];
            case '2'
                b_hang = [b_hang,'0','0','1','0'];
            case '3'
                b_hang = [b_hang,'0','0','1','1'];
            case '4'
                b_hang = [b_hang,'0','1','0','0'];
            case '5'
                b_hang = [b_hang,'0','1','0','1'];
            case '6'
                b_hang = [b_hang,'0','1','1','0'];
            case '7'
                b_hang = [b_hang,'0','1','1','1'];
            case '8'
                b_hang = [b_hang,'1','0','0','0'];
            case '9'
                b_hang = [b_hang,'1','0','0','1'];
            case 'a'
                b_hang = [b_hang,'1','0','1','0'];
            case 'b'
                b_hang = [b_hang,'1','0','1','1'];
            case 'c'
                b_hang = [b_hang,'1','1','0','0'];
            case 'd'
                b_hang = [b_hang,'1','1','0','1'];
            case 'e'
                b_hang = [b_hang,'1','1','1','0'];
            case 'f'
                b_hang = [b_hang,'1','1','1','1'];
        end
    end
    b_cot = [b_cot;b_hang];
    b_hang = '';
end
b = b_cot;
%%%%%%%%%%%%%% bin >>> floating_point ( ftp)
ftp = [];
for i=1:mem_size
    s = bin2dec(b(i,1))*(-2) + 1;
    e = bin2dec(b(i,2:9)) - 127;
    f = 0;
    for j=10:32
      f = bin2dec(b(i,j)) * 2^(9-j) + f;
    end
    f = f + 1;
    f_p = s * f * 2^e;
    ftp = [ftp;f_p];
end
delta_data_hardware = ftp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'C:\Users\G551\Desktop\MFCC\Output_Matlab/delta_out_data.mat';
load(path);
Error_Delta = delta_out_data - delta_data_hardware(1 : size (delta_out_data));
% [Max_Error_Delta, I] = max (abs(Error_Delta))
% Saisotuongdoi = abs (Error_Delta (I)/ delta_out_data (I)*100)

Max_Error_Delta = max (abs(Error_Delta))
Max_Delta = max (abs(delta_out_data))
Saisotuongdoi = abs ((Max_Error_Delta/Max_Delta)*100)
fprintf('==========Relative error of Delta result is %f %% ',Saisotuongdoi);

%==========================================================================

h = fopen('C:\Users\G551\Desktop\MFCC\Output_Hardware/output_file_delta_2nd.txt','r');
nextline = '';
str='';
while ischar(nextline)
    nextline = fgetl(h);
    if ischar(nextline)
        str = [str;nextline];
    end
end
[p,q] = size(str);
mem_size = p;
% hex >>>> bin
b_hang = '';
b_cot  = '';
for i=1: mem_size
    for j=1:8
        switch str(i,j)
            case '0'
                b_hang = [b_hang,'0','0','0','0'];
            case '1'
                b_hang = [b_hang,'0','0','0','1'];
            case '2'
                b_hang = [b_hang,'0','0','1','0'];
            case '3'
                b_hang = [b_hang,'0','0','1','1'];
            case '4'
                b_hang = [b_hang,'0','1','0','0'];
            case '5'
                b_hang = [b_hang,'0','1','0','1'];
            case '6'
                b_hang = [b_hang,'0','1','1','0'];
            case '7'
                b_hang = [b_hang,'0','1','1','1'];
            case '8'
                b_hang = [b_hang,'1','0','0','0'];
            case '9'
                b_hang = [b_hang,'1','0','0','1'];
            case 'a'
                b_hang = [b_hang,'1','0','1','0'];
            case 'b'
                b_hang = [b_hang,'1','0','1','1'];
            case 'c'
                b_hang = [b_hang,'1','1','0','0'];
            case 'd'
                b_hang = [b_hang,'1','1','0','1'];
            case 'e'
                b_hang = [b_hang,'1','1','1','0'];
            case 'f'
                b_hang = [b_hang,'1','1','1','1'];
        end
    end
    b_cot = [b_cot;b_hang];
    b_hang = '';
end
b = b_cot;
%%%%%%%%%%%%%% bin >>> floating_point ( ftp)
ftp = [];
for i=1:mem_size
    s = bin2dec(b(i,1))*(-2) + 1;
    e = bin2dec(b(i,2:9)) - 127;
    f = 0;
    for j=10:32
      f = bin2dec(b(i,j)) * 2^(9-j) + f;
    end
    f = f + 1;
    f_p = s * f * 2^e;
    ftp = [ftp;f_p];
end
delta_2nd_data_hardware = ftp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'C:\Users\G551\Desktop\MFCC\Output_Matlab/delta_2nd_out_data.mat';
load(path);
Error_Delta_2nd = delta_2nd_out_data - delta_2nd_data_hardware(1 : size (delta_2nd_out_data));
% [Max_Error_Delta, I] = max (abs(Error_Delta))
% Saisotuongdoi = abs (Error_Delta (I)/ delta_out_data (I)*100)

Max_Error_Delta_2nd = max (abs(Error_Delta_2nd))
Max_Delta_2nd = max (abs(delta_2nd_out_data))
Saisotuongdoi = abs ((Max_Error_Delta_2nd/Max_Delta_2nd)*100)
fprintf('==========Relative error of Delta-Delta result is %f %% ',Saisotuongdoi);



%==========================================================================

for i=3:Frame_Num_Blocking-2
        error_frame (j,i-2)= Error_Energy (i);
        j = j + 1;
    for c=1:Cepstral_Num
        error_frame (j,i-2)= Error_Cepstral ((i-1)*Cepstral_Num + c);
        j = j + 1;
    end
    
    for d=1:Cepstral_Num+1
        error_frame (j,i-2)= Error_Delta ((i-3)*Cepstral_Num + d);
        j = j + 1;
    end
  
    j = 1;
end
for i= 1:Frame_Num_Blocking-4
    saisotb (i) = sum (abs (error_frame (:,i))) / ((Cepstral_Num + 1) * 2);
end
saisotb = saisotb';
plot(saisotb);
xlabel('Frame');
ylabel('Average Error');
title(' AVERAGE ERROR FOR EACH FRAME ANALYSIS');
%saisotb = (round (saisotb * 10000))/10000;
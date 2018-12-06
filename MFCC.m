% ====================== Parameters ======================== %
Sample_Num = 4000;
Sample_In_Frame = 320;
Overlap_Value = 0.5; %% 50 %%
FFT_Num = 512;
FFT_Stage = 9;
Mel_Filter_Num = 63;
Fs = 8000;
Cepstral_Num = 31;

% ====================== Internal Parameters =============== %
Frame_Num = (Sample_Num /Sample_In_Frame)
Frame_Num_Blocking = 1 + ( (Sample_Num - Sample_In_Frame)/((1 - Overlap_Value)*(Sample_In_Frame)) );
Sample_Num_Blocking= Sample_In_Frame * (Frame_Num_Blocking);


% ====================== Pre-emphasis ======================= %

load('C:\Users\G551\Desktop\MFCC\Input/input_data.mat','input_data');
x=input_data;
preem_out_data=input_data;

for i_pre=2:Sample_Num 
      preem_out_data(i_pre) = x(i_pre) - 0.97* x(i_pre-1);
end
% save('C:\Users\G551\Desktop\MFCC\Output_Matlab/preem_out_data.mat','preem_out_data');

% Cau hinh theo dinh dang MEM cho Pre-emphasis %
x_1 = preem_out_data;
for i_pre_b = 1:(Frame_Num_Blocking)
        for j_pre_b = 1:Sample_In_Frame %%From 1 to 160
            preem_data_matrix(j_pre_b,i_pre_b) = x_1(j_pre_b + (Sample_In_Frame * (1 - Overlap_Value))*(i_pre_b-1)); %% Chu y (1 - Overlap_Value)
        end
end
preem_blocking_out_data = preem_data_matrix(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/preem_blocking_out_data.mat','preem_blocking_out_data');

% ====================== Window ======================= %

% Create WIN_COF_BLOCKING = wcb % 
for j_wcb = 1:(Frame_Num_Blocking)
    for i_wcb = 1:Sample_In_Frame
    window_cof_matrix(i_wcb,j_wcb) = 0.54 - 0.46 * cos( (2 * pi * (i_wcb-1)) / (Sample_In_Frame - 1) );
    window_cof_blocking = window_cof_matrix(:);
    end
end
% save('C:\Users\G551\Desktop\MFCC\Coefficient/window_cof_matrix.mat','window_cof_matrix');
% save('C:\Users\G551\Desktop\MFCC\Coefficient/window_cof_blocking.mat','window_cof_blocking');

% Create WIN_COF by floating point format = wc % 
for i_wc = 1:Sample_In_Frame
    a(i_wc) = 0.54 - 0.46 * cos( (2 * pi * (i_wc-1)) / (Sample_In_Frame - 1) );
    window_cof = a';
end
% save('C:\Users\G551\Desktop\MFCC\Coefficient/window_cof.mat','window_cof');
h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/window_cof.txt','w');
in_data = fpt2bin(window_cof);
[p,q] =size(in_data)
for k=1:(p-1)
    fprintf(h,'%s\n',in_data(k,:));
end
fprintf(h,'%s',in_data(p,:));
fclose(h);
display('ok');

% Window = w
for i_w = 1:Sample_Num_Blocking
window_out_data_temp(i_w) = preem_blocking_out_data(i_w) * window_cof_blocking(i_w);
window_out_data = window_out_data_temp';
end
% save('C:\Users\G551\Desktop\MFCC\Output_Matlab/window_out_data.mat','window_out_data');

% ====================== FFT = fft ======================= %
N = FFT_Num;
m = FFT_Stage;
x_window = window_cof_matrix .* preem_data_matrix;
%mo rong ra theo so diem FFT
x_fft_init=x_window;
[l_x_window, w_x_window]=size(x_window);
number_of_frame = w_x_window
if (length(x_fft_init)<N)
  x_fft_init(length(x_fft_init)+1:N,:)=0;
end
window_out_data_morong = x_fft_init(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/window_out_data_morong.mat','window_out_data_morong');

%tinh ma tran W
%san 
for i_fft=1:N
  W(i_fft)=exp(-1i*2*pi*(i_fft-1)/N);
  w_thuc(i_fft)= real(W(i_fft));
  w_ao(i_fft) =  imag(W(i_fft));

end

%Tao ma tran W dang MEM
input_data_1 = w_thuc';
h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_w_real.txt','w');
w_real = fpt2bin(input_data_1);
[p,q] =size(w_real)
for k=1:(p-1)
    fprintf(h,'%s\n',w_real(k,:));
end
fprintf(h,'%s',w_real(p,:));
fclose(h);

input_data_2 = w_ao';
h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_w_image.txt','w');
w_image = fpt2bin(input_data_2);
[p,q] =size(w_image)
for k=1:(p-1)
    fprintf(h,'%s\n',w_image(k,:));
end
fprintf(h,'%s',w_image(p,:));
fclose(h);

%sap xep lai de chuan bi nhan bang pp dao bit
for i_fft=1:N
  %doi i ra bit
  i_bis=i_fft-1;
  binary=0;
  for j_fft=1:m
    %binary(1)=LSB; binary(m)=MSB  
    binary(j_fft)=mod(i_bis,2);
    i_bis=floor(i_bis/2);
  end
  %dao bit, luc nay: binary(1)=MSB; binary(m)=LSB
  %tinh lai vi tri moi
  i_new=1;
  for j_fft=1:m
    i_new=i_new + binary(m+1-j_fft)*2^(j_fft-1);    
  end
  %gan lai cac gia tri trong x_fft_init vao x_fft
  x_fft_out_data(i_fft,:)=x_fft_init(i_new,:);
end
%tinh fft
for i_fft=1:m
  %de tinh fft cho moi stage_ung voi moi i
  %1.tinh buoc nhay
  %2.tinh xem k dang o nhom nao o stage nao
  %3.tinh vi tri cua k o trong nhom
  jump_factor=2^(i_fft-1); % tinh buoc nhay
  x_fft_temp(1:2^m,1:number_of_frame)=0;
  for k=1:N
    a=ceil((k-1)/(2^i_fft)); % a= nhom cua k
    b=mod(k-1,2^i_fft); %b= vi tri cua k trong nhom thu a
    if (b<=2^(i_fft-1)-1)
      x_fft_temp(k,:)=x_fft_out_data(k,:) + x_fft_out_data(k+jump_factor,:)*W(b*2^(m-i_fft)+1);  
    elseif (b>2^(i_fft-1)-1)
      x_fft_temp(k,:)=x_fft_out_data(k-jump_factor,:) + x_fft_out_data(k,:)*W(b*2^(m-i_fft)+1);   
    end
  end
  x_fft_out_data(:,:)=x_fft_temp;
end
% save('C:\Users\G551\Desktop\MFCC\Output_Matlab/x_fft_out_data','x_fft_out_data');

% ====================== Amplitude = amp ======================= %
real_x_fft = real(x_fft_out_data(:,:));
image_x_fft = imag(x_fft_out_data(:,:));
real_x_fft_abs(:,:)  = abs(real_x_fft(:,:));
image_x_fft_abs(:,:)  = abs(image_x_fft(:,:));
for i_amp = 1:(Frame_Num_Blocking)
  for j_amp = 1:(FFT_Num/2)  
    if(real_x_fft_abs(j_amp,i_amp) >= image_x_fft_abs(j_amp,i_amp)) 
      x_magnitude(j_amp,i_amp) = real_x_fft_abs(j_amp,i_amp) + 0.25 * image_x_fft_abs(j_amp,i_amp);
    else
      x_magnitude(j_amp,i_amp) = image_x_fft_abs(j_amp,i_amp)  + 0.25 * real_x_fft_abs(j_amp,i_amp);
  end
  end
end

amplitude_out_data = x_magnitude(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/amplitude_out_data','amplitude_out_data');

% ===================== Mel = mel ============================ %
% CREATE Mel_coefficient
temp1 = melbankm(Mel_Filter_Num,FFT_Num,Fs,0,0.5,'b');
temp1(:,((FFT_Num/2)+ 1))=[];
temp2 = full(temp1);
mel_cof_matrix = temp2;
temp3 = temp2';
mel_cof = temp3(:);
mel_cof_temp = mel_cof;
% save('C:\Users\G551\Desktop\MFCC\Coefficient/mel_cof.mat','mel_cof');
% save('C:\Users\G551\Desktop\MFCC\Coefficient/mel_cof_matrix.mat','mel_cof_matrix');

% CREATE Mel_coefficient for 8 Memory
for i_mel = 1:(((FFT_Num)/2)*Mel_Filter_Num)
   if(i_mel<= 4096)
       para =1;
       if (length(mel_cof)<4097)
       mel_cof(length(mel_cof)+1:4096,:)=0;
       end    
   elseif((i_mel>4097)&(i_mel<=8192))
       para =2;
       if (length(mel_cof)<8193)
       mel_cof(length(mel_cof)+1:8192,:)=0;
       end  
   elseif((i_mel>8193)&(i_mel<=12288))
       para =3; 
       if (length(mel_cof)<12289)
       mel_cof(length(mel_cof)+1:12288,:)=0;
       end  
   elseif((i_mel>12289)&(i_mel<=16384))
       para =4;
       if (length(mel_cof)<16385)
       mel_cof(length(mel_cof)+1:16384,:)=0;
       end  
   elseif((i_mel>16385)&(i_mel<=20480))
       para =5;
       if (length(mel_cof)<20481)
       mel_cof(length(mel_cof)+1:20480,:)=0;
       end
   elseif((i_mel>20481)&(i_mel<=24576))
       para =6;
       if (length(mel_cof)<24577)
       mel_cof(length(mel_cof)+1:24576,:)=0;
       end
   elseif((i_mel>24577)&(i_mel<=28672))
       para =7;
       if (length(mel_cof)<28673)
       mel_cof(length(mel_cof)+1:28672,:)=0;
       end
   elseif((i_mel>28673)&(i_mel<=32768))
       para =8;
       if (length(mel_cof)<32769)
       mel_cof(length(mel_cof)+1:32768,:)=0;
       end
   end
end
switch para
    case{1}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %     
    case{2}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
    case{3}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
    case{4}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
        mel_cof_3 = mel_cof(12289:16384);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_3.txt','w');
        in_data = fpt2bin(mel_cof_3);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % % 
    case{5}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_3 = mel_cof(12289:16384);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_3.txt','w');
        in_data = fpt2bin(mel_cof_3);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_4 = mel_cof(16385:20480);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_4.txt','w');
        in_data = fpt2bin(mel_cof_4);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
    case{6}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_3 = mel_cof(12289:16384);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_3.txt','w');
        in_data = fpt2bin(mel_cof_3);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_4 = mel_cof(16385:20480);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_4.txt','w');
        in_data = fpt2bin(mel_cof_4);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_5 = mel_cof(20481:24576);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_5.txt','w');
        in_data = fpt2bin(mel_cof_5);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
    case{7}
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_3 = mel_cof(12289:16384);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_3.txt','w');
        in_data = fpt2bin(mel_cof_3);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_4 = mel_cof(16385:20480);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_4.txt','w');
        in_data = fpt2bin(mel_cof_4);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_5 = mel_cof(20481:24576);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_5.txt','w');
        in_data = fpt2bin(mel_cof_5);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_6 = mel_cof(24577:28672);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_6.txt','w');
        in_data = fpt2bin(mel_cof_6);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
    otherwise
        mel_cof_0 = mel_cof(1:4096);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_0.txt','w');
        in_data = fpt2bin(mel_cof_0);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_1 = mel_cof(4097:8192);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_1.txt','w');
        in_data = fpt2bin(mel_cof_1);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_2 = mel_cof(8193:12288);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_2.txt','w');
        in_data = fpt2bin(mel_cof_2);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_3 = mel_cof(12289:16384);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_3.txt','w');
        in_data = fpt2bin(mel_cof_3);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_4 = mel_cof(16385:20480);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_4.txt','w');
        in_data = fpt2bin(mel_cof_4);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_5 = mel_cof(20481:24576);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_5.txt','w');
        in_data = fpt2bin(mel_cof_5);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_6 = mel_cof(24577:28672);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_6.txt','w');
        in_data = fpt2bin(mel_cof_6);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
        mel_cof_7 = mel_cof(28673:32768);
        % % % % % % % Create MEL_COF_0 by floating point format % % % % % %
        h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/mem_cof_7.txt','w');
        in_data = fpt2bin(mel_cof_7);
        [p,q] =size(in_data)
        for k=1:(p-1)
        fprintf(h,'%s\n',in_data(k,:));
        end
        fprintf(h,'%s',in_data(p,:));
        fclose(h);
        display('ok');
        % % % % % % %                   End               % % % % % % %
end

% Mel = mel
mel_out_data_matrix = mel_cof_matrix * x_magnitude;
mel_out_data = mel_out_data_matrix(:);
% save('C:\Users\G551\Desktop\MFCC\Output_Matlab/mel_out_data','mel_out_data');

% ===================== Log = log ============================ %
temp_log = abs(mel_out_data_matrix);
log_out_data_matrix = log10(temp_log);
log_out_data = log_out_data_matrix(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/log_out_data','log_out_data'); 
% Log_out_data chinh la Mel ben hardware vi trong khoi Mel o hardware da chen luon khoi LOG.

% ===================== Cepstral = cep ============================ %
%CREATE Cepstral coefficient
for i_cep=1:Cepstral_Num
  for j_cep=1:Mel_Filter_Num
    cos_factor(i_cep,j_cep)= cos((j_cep-0.5)*i_cep*pi/Mel_Filter_Num); 
  end
end
cepstral_cof_matrix = cos_factor;
temp_cep_cof = cos_factor';
cepstral_cof = temp_cep_cof(:);
% save ('C:\Users\G551\Desktop\MFCC\Coefficient/cepstral_cof.mat','cepstral_cof');

%CREATE Cepstral_coefficient by floating point format
h = fopen('C:\Users\G551\Desktop\MFCC\Coefficient/cepstral_cof.txt','w');
in_data = fpt2bin(cepstral_cof);
[p,q] =size(in_data)
for k=1:(p-1)
    fprintf(h,'%s\n',in_data(k,:));
end
fprintf(h,'%s',in_data(p,:));
fclose(h);
display('ok');

%Cepstral
cepstral_out_data_matrix = cepstral_cof_matrix * log_out_data_matrix;
cepstral_out_data = cepstral_out_data_matrix(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/cepstral_out_data','cepstral_out_data');

% ===================== Energy = E ============================ %
x_E = input_data;
for i_E = 1:(Frame_Num_Blocking)
        for j_E = 1:Sample_In_Frame 
            x_E_Blocking_matrix(j_E,i_E) = x_E(j_E + (Sample_In_Frame *(1 - Overlap_Value))*(i_E-1)); %% Chu y (1 - Overlap_Value)
        end
end
for i_E=1:(Frame_Num_Blocking)
  E(1,i_E)=0;  
  for j_E=1:Sample_In_Frame  
    E(1,i_E)=E(1,i_E)+ x_E_Blocking_matrix(j_E,i_E)*x_E_Blocking_matrix(j_E,i_E);
  end
end

E=log10(E);     % NOTE
energy_out_data = E';
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/energy_out_data','energy_out_data');

% ===================== Delta = D ============================ %
% x_mfcc_init_test=[cepstral_out_data_matrix;E]; /// NOTE 18/01/2016
x_mfcc_init=[E;cepstral_out_data_matrix];
x_delta=zeros(Cepstral_Num+1,Frame_Num_Blocking-4);
for i_D=3:Frame_Num_Blocking-2
  x_delta(:,i_D)=(2*(x_mfcc_init(:,i_D+2)-x_mfcc_init(:,i_D-2))+(x_mfcc_init(:,i_D+1)-x_mfcc_init(:,i_D-1)));
end

final= x_delta;
final(:,1)=[];
final(:,1)=[];
delta_out_data = final(:);

% ===================== Delta_Second Order = D_2nd ============================ %
x_delta_2nd=zeros(Cepstral_Num+1,Frame_Num_Blocking-8);
for i_D_2nd=3:Frame_Num_Blocking-6 
  x_delta_2nd(:,i_D_2nd)=(2*(final(:,i_D_2nd+2)-final(:,i_D_2nd-2))+(final(:,i_D_2nd+1)-final(:,i_D_2nd-1)));
  %x_delta_2nd(:,i_D_2nd)=(2*(x_delta(:,i_D_2nd+2)-x_delta(:,i_D_2nd-2))+(x_delta(:,i_D_2nd+1)-x_delta(:,i_D_2nd-1)));
end

final_2nd= x_delta_2nd;
final_2nd(:,1)=[];
final_2nd(:,1)=[];
delta_2nd_out_data = final_2nd(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/delta_2nd_out_data','delta_2nd_out_data');

% ===================== Result ============================ %
x_mfcc=[x_mfcc_init(:,3:Frame_Num_Blocking-2);x_delta(:,3:Frame_Num_Blocking-2)];
% delta_out_data = x_mfcc(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/delta_out_data','delta_out_data');
Result_matrix = [ E(:,3:number_of_frame-2);cepstral_out_data_matrix(:,3:number_of_frame-2);x_delta(:,3:number_of_frame-2) ];
Result = Result_matrix(:);
save('C:\Users\G551\Desktop\MFCC\Output_Matlab/Result','Result');




% ========================= END =========================== %







function bin_create = fpt2bin(in_fpt)
% % % % % % % % 
[mem_size,z] = size(in_fpt);
bin_result = [];
for i=1:mem_size
    if (in_fpt(i) == 0)
        fpt_str = '00000000000000000000000000000000 ';
    elseif (in_fpt(i) > 0) % Positive number
        k = 0;
        s = '0';
        if (in_fpt(i) >= 2)
            while (in_fpt(i) >= 2)
                in_fpt(i) = in_fpt(i)/2;
                k = k + 1;
            end
            e = 127 + k;
        elseif (in_fpt(i) < 1)
            while (in_fpt(i) < 1)
                in_fpt(i) = in_fpt(i)*2;
                k = k + 1;
            end
            e = 127 - k;
        else
            e = 127;
        end
        f = in_fpt(i)-1;
        f_str = [];
        for j =1:23
            f_t = f - 2^(-j);
            if (f_t >= 0)
                f = f_t;
                f_str = [f_str,'1'];
            else
                f_str = [f_str,'0'];
            end
        end
        e_str = dec2bin(e,8);
        fpt_str = [s,e_str,f_str];
        
    else               %Negative number 
        k = 0;
        s = '1';
        in_fpt(i) = -in_fpt(i);
        if (in_fpt(i) >= 2)
            while (in_fpt(i) >= 2)
                in_fpt(i) = in_fpt(i)/2;
                k = k + 1;
            end
            e = 127 + k;
        elseif (in_fpt(i) < 1)
            while (in_fpt(i) < 1)
                in_fpt(i) = in_fpt(i)*2;
                k = k + 1;
            end
            e = 127 - k;
        else
            e = 127;
        end
        f = in_fpt(i)-1;
        f_str = [];
        for j =1:23
            f_t = f - 2^(-j);
            if (f_t >= 0)
                f = f_t;
                f_str = [f_str,'1'];
            else
                f_str = [f_str,'0'];
            end
        end
        e_str = dec2bin(e,8);
        fpt_str = [s,e_str,f_str];
    end
    bin_kq = fpt_str;
 
bin_temp = [bin_kq(1),bin_kq(2:5),bin_kq(6:9),bin_kq(10:13),bin_kq(14:17),bin_kq(18:21),bin_kq(22:25),bin_kq(26:29),bin_kq(30:32),' '];
bin_result = [bin_result;bin_temp];
end
bin_create = bin_result;

end
%% Display result
function kholi2()
    clc; clear;
    %vectorA()
    %T_hat(1);
    %vectorX();
    Temp()
    A = vectorA();
    B = vectorB();
    A_p = A(1:end,1:end);
    B_p = B(1:end,1);
    X = pinv(A_p)*B_p;
    


end


%% Parameter function
function [Tmax,Tmin,rng,Tavg,N,Data] = par()
    Data =[25.55,24.84,23.75,22.73,22.26,22.59,23.6,25.13,27.09,29.48,32.14,34.66,36.49,37.28,37.05,36.16,34.95,33.53,31.84,29.92,28.07,26.71,26.03,25.82];
    %Temp_Data2 = [25.75, 25.07, 24.03, 23.05, 22.61, 22.91, 23.89, 25.35, 27.22, 29.5,32.04, 34.45,36.2,36.95, 36.74, 35.89, 34.73, 33.37, 31.76, 29.92, 28.15, 26.85, 26.21, 26.01];
    %Temp2 = [24.28,23.57,22.475,21.445,20.98,21.3,22.325,23.86,25.83,28.225,30.895,33.43,35.27,36.06,35.835,34.94,33.725,32.295,30.6,28.665,26.805,25.44,24.765,24.555];
    %Data = [31.21,30.59,29.63,28.73,28.33,28.61,29.5,30.85,32.56,34.66,36.99,39.2,40.81,41.5,41.31,40.53,39.46,38.22,36.73,35.04,33.42,32.23,31.64,31.45];
    Tmax = max(Data);
    Tmin = min(Data);
    Tavg = mean(Data);
    rng = max(Data)-min(Data);
    %rng = 14.12;
    N = 12;
end
%% vector A
function A = vectorA()
    N = 12; % Maximum number of harmonics
    T = 24;
    N0 = 12;
    %[N,N0,M] = var();
    A = zeros(28,28);
    sum_s = 0;
    sum_c = 0;
    sum_sc = 0;
    ao = 1;
    aa = 1;
    Hmax = 1;
    Hmin = 1;
    %rng = 5;
  %row 1-12, col 1-24

      
    for n1=1:N   
        for t1=1:T
            sum_s = 2*sin(n1*pi*t1/12);
            sum_c = 2*cos(n1*pi*t1/12);
            %A(1,1) = A(1,1) +( 2* aa);
            if n1<=1
                A(1,1) = A(1,1) +( 2* aa);
            end
            A(1,n1+1) = A(1,n1+1) + sum_s;
            A(1,13+n1) = A(1,13+n1) + sum_c;
            sum_s = 0;
            sum_c = 0;
        end



    end
%end



%for nd=1:NofD     
    for x = 1:12
        for n=1:N   
            for t=1:T
                sum_s = sin(n*pi*t/12)*(2*sin(x*pi*t/12));
                sum_c = cos(n*pi*t/12)*(2*sin(x*pi*t/12));
                A(x+1,1) = A(x+1,1)+ (ao * (2*sin(x*pi*t/12)));
                A(x+1,n+1) = A(x+1,n+1) + sum_s;
                A(x+1,13+n) = A(x+1,13+n) + sum_c;
                sum_s = 0;
                sum_c = 0;
            end



        end
    end


    for x3 = 1:12
        for n3=1:N   
            for t3=1:T
                sum_s = sin(n3*pi*t3/12)*(2*cos(x3*pi*t3/12));
                sum_c = cos(n3*pi*t3/12)*(2*cos(x3*pi*t3/12));
                A(x3+13,1) = A(x3+13,1)+ (ao * (2*cos(x3*pi*t3/12)));
                A(x3+13,n3+1) = A(x3+13,n3+1) + sum_s;
                A(x3+13,13+n3) = A(x3+13,13+n3) + sum_c;
                sum_s = 0;
                sum_c = 0;
            end



        end
    end


    % lamda1,lamda2,lamda3 (equation 1)
    A(1,N*2+2) = -1;
    A(1,N*2+3) = 1;
    A(1,N*2+4) = 24/23;
    
    for n4=1:N
       % sin (eqn 2-13)
       A(n4+1,N*2+2) = -1*sin(n4*pi*Hmax/12); 
       A(n4+1,N*2+3) =  1*sin(n4*pi*Hmin/12);
       %cos (eqn 14-25)
       A(n4+13,N*2+2) = -1*cos(n4*pi*Hmax/12); 
       A(n4+13,N*2+3) =  1*cos(n4*pi*Hmin/12);
       for t4=1:T
           % sin (eqn 2-13) %lamda 3
           A(n4+1,N*2+4)= A(n4+1,N*2+4)+ (1/23)*sin(n4*pi*t4/12);
           % sin (eqn 2-13)
           A(n4+13,N*2+4)= A(n4+13,N*2+4)+ (1/23)*cos(n4*pi*t4/12);
       end
    end
    
    
    % (equation 26-28, lamdas 1-3)
    A(N*2+2:end,N*2+2) = 0;
    A(N*2+2:end,N*2+3) = 0;
    A(N*2+2:end,N*2+4) = 0;
    
    % (equation 26-28, a0)
    A(N*2+2,1) = 1;
    A(N*2+3,1) = 1;
    A(N*2+4,1) = 24/23;
    
    % aN & bN
    for n5=1:N
        %eqn 26
        A(N*2+2,n5+1) = 1*sin(n5*pi*Hmax/12);
        A(N*2+2,n5+13) = 1*cos(n5*pi*Hmax/12);
        %eqn 27
        A(N*2+3,n5+1) = 1*sin(n5*pi*Hmin/12);
        A(N*2+3,n5+13) = 1*cos(n5*pi*Hmin/12);
        for t5=1:T
           % sin (eqn 27)
           A(N*2+4,n5+1)= A(N*2+4,n5+1) + (1/23)*sin(n5*pi*t5/12);
           % sin (eqn 27)
           A(N*2+4,n5+13)= A(N*2+4,n5+13) + (1/23)*cos(n5*pi*t5/12);
       end
    end
    A_adj = A([1:12,14:end],[1:12,14:end]);
    size(A_adj);
    %A;

end

%% vector B
function B = vectorB()  
    B =zeros(28,1);
    [Tmax,Tmin,rng,Tavg,N,Temp_Data] = par();
    T = 24;
    M =12;%(N before)
    %rng = Tmax - Tmin;
    std_sum = 0;
    sum_s = 0;
    sum_c = 0;
    
    for t=1:T
        %B(1,1) = B(1,1) + 2*Jan(t); 
        B(1,1) = B(1,1) + Temp_Data(t); 
    end
    B(1,1) = 2 * B(1,1);
    
    for n=1:M
        for t=1:T
            sum_s = sum_s + (sin(n*pi*t/12))*Temp_Data(t);
            sum_c = sum_c + (cos(n*pi*t/12))*Temp_Data(t); 
        end
        B(n+1,1) = 2*sum_s;   
        B(n+13,1) = 2* sum_c;
        sum_s = 0;
        sum_c = 0;
        
    end
    
    % eqn 26
    B(M*2+2,1) = rng;
    %eqn 27
    B(M*2+3,1) = 0;
   
    for t2=1:T
        %std_sum = std_sum + (Tmax -Temp_Data(t2));
        std_sum = std_sum + ((Tavg + rng/2) -Temp_Data(t2));
    end
    
    std = std_sum/23;
    
    %eqn 28
    B(M*2+4,1) = (24/23*rng)-std;
    B;
    B_adj = B([1:12,14:end],1);
end

%%
function [ao,an,bn] = vectorX()
    [Tmax,Tmin,rng,Tavg,N,Temp_Data] = par();
    %N = 4;
    A = vectorA();
    B = vectorB();
    X = pinv(A)* B;
    ao = X(1,1);
    %an = X(2:13,1).';
    %bn = X(14:25,1).';
    an = X(2:N+1,1).';
    bn = X(14:13+N,1).';
end

%% Predict

function Temp()
    [ao,an,bn] = vectorX();
    [Tmax,Tmin,rng,Tavg,N,Temp_Data] = par();
    t = 2;
    %N = 12;
    T = 24;
    sum_an = 0;
    sum_bn = 0;
    sum_s = 0;
    sum_c = 0;
    sum_sigma_avg =0;
    total_sigma_avg = 0;
    sum_dt = 0;
    %rng = Tmax - Tmin;
    sum_temp =0;
    dt = zeros(1,24);
    time = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
    predict_temp = zeros(1,24);
    for n=1:N
        sum_an = sum_an + an(n)*sin(n*pi*t/12);
        sum_bn = sum_bn + bn(n)*cos(n*pi*t/12);
    end
    
    Dt = ao + sum_an + sum_bn;
    Tt =  Tmin + Dt;
    
    for t=1:T
        for n=1:N
            sum_s = sum_s + an(n)*sin(n*pi*t/12);
            sum_c = sum_c + bn(n)*cos(n*pi*t/12);
        end
        sum_sigma_avg = sum_sigma_avg + rng -(ao + sum_s + sum_c);
        total_sigma_avg = total_sigma_avg + rng -(ao + sum_s + sum_c);
        sum_dt = sum_dt + ao + sum_s + sum_c;
        dt(1,t) = dt(1,t) + ao + sum_s + sum_c;
        %predict_temp(1,t)  =  predict_temp(1,t) + Tmax - sum_sigma_avg;
        predict_temp(1,t)  =  predict_temp(1,t) + (Tavg + rng/2) - sum_sigma_avg
        sum_s = 0;
        sum_c = 0;
        sum_sigma_avg = 0;
    end
    total_sigma_avg/23;
    predict_temp
    plot(time,dt)
    
%     for t1 = 1:24
%         n = 1;
%         sum_temp = sum_temp + (35.17 - Temp_Data(t1));
%     end
%     sum_temp/23;
    
end
    
    


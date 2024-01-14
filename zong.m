%% Display result
function p_model_step()
    %T_hat(1);
    t = 2;
    range = 15.92;
    To_max = 35.17;
    %Th(t,range,To_max);
    %Y = vectorY()
    vectorX();
    T_hat(2);
    Th();


end

%% variable funtion
function [N,N0,M] = var()
    clc; clear;
    N = 24;
    N0 = 17;
    M = 4;
end

%% Temperature Data
function Temp = tempData()
    %Temp = [26.52, 29.11, 30.55, 29.44, 27.55, 25.97, 25.03, 24.60, 24.93, 25.71, 26.48, 26.09];
    Temp = [22.81,22.07,20.92,19.84,19.35,19.69,20.76,22.37,24.44,26.95,29.75,32.41,34.34,35.17,34.93,33.99,32.72,31.22,29.44,27.41,25.46,24.03,23.32,23.1 ];
    RH = [44.70, 40.10, 48.82, 65.01, 76.53, 81.91, 86.06, 87.40, 85.75, 82.43, 66.71, 52.20];
    Month = [1,2,3,4,5,6,7,8,9,10,11,12];
    %plot(Month,RH);
    %plot(Month,Temp);
end
%% a0_hat function
function a0_hat = A0_p()
    T = tempData();
    [N,N0,M] = var();  %N0 should be less than N
    sum_a0_hat = 0;
    for i=1:N0
        sum_a0_hat = sum_a0_hat + T(i);
    end
    a0_hat = sum_a0_hat/N0;
end




%% vector A
function A = vectorA()
    %M = 12; % Maximum number of harmonics
    %N = 12;
    %N0 = 11;
    [N,N0,M] = var();
    A = zeros(M*2,M*2);
    sum_cc = 0;
    sum_sc = 0;
    sum_cs = 0;
    sum_ss = 0;
    for i=1:M   %row 1-12, col 1-24
        for j=1:M   
            for t=1:N0
                sum_cc =  sum_cc + (cos(j*2*pi*t/N)*cos(i*2*pi*t/N)); %cos cos
                sum_sc =  sum_sc + (sin(j*2*pi*t/N)*cos(i*2*pi*t/N)); %sin cos
                sum_cs =  sum_cs + (cos(j*2*pi*t/N)*sin(i*2*pi*t/N)); %cos sin
                sum_ss =  sum_ss + (sin(j*2*pi*t/N)*sin(i*2*pi*t/N)); %sin sin
            end
            %A(i,j) = A(i,j)+ sum_cc;
            A(i,j) = sum_cc;
            A(i,M+j) = sum_sc;  %A(i,M) = sum_sc;
            A(i+M,j) = sum_cs;
            A(i+M,j+M) = sum_ss;
            sum_cc = 0;
            sum_sc = 0;
            sum_cs = 0;
            sum_ss = 0;
            
            
        end
    end
    
    A;
end

%% vector Y
function Y = vectorY()
    T = tempData();
    [N,N0,M] = var();
    a0_p = A0_p();
    Y =zeros(M,1);
    sum_c = 0;
    sum_s = 0;
    for i=1:M
        for t=1:N0
            sum_c = sum_c + ((T(t)-a0_p)*cos(i*2*pi*t/N));
            sum_s = sum_s + ((T(t)-a0_p)*sin(i*2*pi*t/N));
        end
        Y(i,1) = sum_c;
        Y(i+(M),1) = sum_s;
        sum_c = 0;
        sum_s = 0;
    end
    Y;
end

%% vector X

function X = vectorX()
    A = vectorA();
    Y = vectorY();
    X = pinv(A)*Y
end

%% Dry bulb temperature and Relative Humidity Distribution

function T = T_hat(t)
    %t = 12;  %Month of the year
    [N,N0,M] = var();
    X = vectorX();
    a_hat = (X(1:M,1)).';
    b_hat = (X(M+1:end,1)).';
    a0_hat = A0_p();
    t_sum = 0;
    
    for n=1:M
       %t_sum = t_sum + (a_hat(n)*cos(n*2*pi*t/N) + b_hat(n)*sin(n*2*pi*t/N)); 
       t_sum = t_sum + (a_hat(n)*cos(n*2*pi*t/N) + b_hat(n)*sin(n*2*pi*t/N));
    end
    
    T = a0_hat + t_sum;
    %display(T);
end


%% Prediction Model
function Th()
    [N,N0,M] = var();
    Temp = tempData();
    time = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
    T_predict = zeros(1,24);
    for i = 1:N
        T_predict(1,i) = T_hat(i);
    end
    T_predict

    hold on;
    plot(time,Temp,'bs-','LineWidth',2,'MarkerSize',5);
    plot(time,T_predict,'ro-','LineWidth',2,'MarkerSize',5);
    xlabel('TIME OF THE DAY (HRS)');
    ylabel('MEAN DRY BULB TEMPERATURE (deg CELCIUS)');
    legend('Original_Data','Predicted_Data');

    
end
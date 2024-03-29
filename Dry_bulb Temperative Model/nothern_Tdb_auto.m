%% Display result
function kholi_mod2()
    clc; clear;
    Temp();
    A = vectorA();
    B = vectorB();
    A_p = A(1:end,1:end);
    B_p = B(1:end,1);
    X = pinv(A_p)*B_p;
    MAD();
    [Rng_adj,Tmin] = Estimation_Data();
    


end

%% value of k
function k = kfunc()
    Jan = 1; Feb = 2; Mar = 3; Apr = 4; May = 5; Jun = 6; Jul = 7; Aug = 8; Sep = 9; Oct = 10; Nov = 11; Dec = 12;
    Yola = 1; Bauchi = 2; Maiduguri = 3; Kaduna = 4; Kano = 5; Ilorin = 6; Minna = 7; Jos = 8; Sokoto = 9;
    
    City = Bauchi;  
    Month = Jan; 
    
    k = Month + 12*(City-1);  %This particular program is indepedent of k
    %k = 100;
    
    
end

%% k designation
function [City,Month] = Data_Tag()

    k = kfunc();
    
    if (k==1)
        City = "Yola";
        Month = "January";
    
    elseif (k==2)
        City = "Yola";
        Month = "February";

    elseif (k==3)
    City = "Yola";
    Month = "March";

    elseif (k==4)
    City = "Yola";
    Month = "April";

    elseif (k==5)
    City = "Yola";
    Month = "May";

    elseif (k==6)
    City = "Yola";
    Month = "June";

    elseif (k==7)
    City = "Yola";
    Month = "July";

    elseif (k==8)
    City = "Yola";
    Month = "August";

    elseif (k==9)
    City = "Yola";
    Month = "September";

    elseif (k==10)
    City = "Yola";
    Month = "October";

    elseif (k==11)
    City = "Yola";
    Month = "November";

    elseif (k==12)
    City = "Yola";
    Month = "December";
    
    elseif (k==13)
        City = "Bauchi";
        Month = "January";
    
    elseif (k==14)
        City = "Bauchi";
        Month = "February";

    elseif (k==15)
    City = "Bauchi";
    Month = "March";

    elseif (k==16)
    City = "Bauchi";
    Month = "April";

    elseif (k==17)
    City = "Bauchi";
    Month = "May";

    elseif (k==18)
    City = "Bauchi";
    Month = "June";

    elseif (k==19)
    City = "Bauchi";
    Month = "July";

    elseif (k==20)
    City = "Bauchi";
    Month = "August";

    elseif (k==21)
    City = "Bauchi";
    Month = "September";

    elseif (k==22)
    City = "Bauchi";
    Month = "October";

    elseif (k==23)
    City = "Bauchi";
    Month = "November";

    elseif (k==24)
    City = "Bauchi";
    Month = "December";
    
    elseif (k==25)
    City = "Maiduguri";
    Month = "January";
    
    elseif (k==26)
        City = "Maiduguri";
        Month = "February";

    elseif (k==27)
        City = "Maiduguri";
        Month = "March";

    elseif (k==28)
        City = "Maiduguri";
        Month = "April";

    elseif (k==29)
        City = "Maiduguri";
        Month = "May";

    elseif (k==30)
        City = "Maiduguri";
        Month = "June";

    elseif (k==31)
        City = "Maiduguri";
        Month = "July";

    elseif (k==32)
        City = "Maiduguri";
        Month = "August";

    elseif (k==33)
        City = "Maiduguri";
        Month = "September";

    elseif (k==34)
        City = "Maiduguri";
        Month = "October";

    elseif (k==35)
        City = "Maiduguri";
        Month = "November";

    elseif (k==36)
        City = "Maiduguri";
        Month = "December";
    
    elseif (k==37)
        City = "Kaduna";
        Month = "January";

    elseif (k==38)
        City = "Kaduna";
        Month = "February";

    elseif (k==39)
        City = "Kaduna";
        Month = "March";

    elseif (k==40)
        City = "Kaduna";
        Month = "April";

    elseif (k==41)
        City = "Kaduna";
        Month = "May";

    elseif (k==42)
        City = "Kaduna";
        Month = "June";

    elseif (k==43)
        City = "Kaduna";
        Month = "July";

    elseif (k==44)
        City = "Kaduna";
        Month = "August";

    elseif (k==45)
        City = "Kaduna";
        Month = "September";

    elseif (k==46)
        City = "Kaduna";
        Month = "October";

    elseif (k==47)
        City = "Kaduna";
        Month = "November";

    elseif (k==48)
        City = "Kaduna";
        Month = "December";
    
    elseif (k==49)
        City = "Kano";
        Month = "January";
    
    elseif (k==50)
        City = "Kano";
        Month = "February";

    elseif (k==51)
        City = "Kano";
        Month = "March";

    elseif (k==52)
        City = "Kano";
        Month = "April";

    elseif (k==53)
        City = "Kano";
        Month = "May";

    elseif (k==54)
        City = "Kano";
        Month = "June";

    elseif (k==55)
        City = "Kano";
        Month = "July";

    elseif (k==56)
        City = "Kano";
        Month = "August";

    elseif (k==57)
        City = "Kano";
        Month = "September";

    elseif (k==58)
        City = "Kano";
        Month = "October";

    elseif (k==59)
        City = "Kano";
        Month = "November";

    elseif (k==60)
        City = "Kano";
        Month = "December";
        
    elseif (k==61)
        City = "Ilorin";
        Month = "January";
    
    elseif (k==62)
        City = "Ilorin";
        Month = "February";

    elseif (k==63)
        City = "Ilorin";
        Month = "March";

    elseif (k==64)
        City = "Ilorin";
        Month = "April";

    elseif (k==65)
        City = "Ilorin";
        Month = "May";

    elseif (k==66)
        City = "Ilorin";
        Month = "June";

    elseif (k==67)
        City = "Ilorin";
        Month = "July";

    elseif (k==68)
        City = "Ilorin";
        Month = "August";

    elseif (k==69)
        City = "Ilorin";
        Month = "September";

    elseif (k==70)
        City = "Ilorin";
        Month = "October";

    elseif (k==71)
        City = "Ilorin";
        Month = "November";

    elseif (k==72)
        City = "Ilorin";
        Month = "December";
        
    elseif (k==73)
        City = "Minna";
        Month = "January";
    
    elseif (k==74)
        City = "Minna";
        Month = "February";

    elseif (k==75)
        City = "Minna";
        Month = "March";

    elseif (k==76)
        City = "Minna";
        Month = "April";

    elseif (k==77)
        City = "Minna";
        Month = "May";

    elseif (k==78)
        City = "Minna";
        Month = "June";

    elseif (k==79)
        City = "Minna";
        Month = "July";

    elseif (k==80)
        City = "Minna";
        Month = "August";

    elseif (k==81)
        City = "Minna";
        Month = "September";

    elseif (k==82)
        City = "Minna";
        Month = "October";

    elseif (k==83)
        City = "Minna";
        Month = "November";

    elseif (k==84)
        City = "Minna";
        Month = "December";
        
    elseif (k==85)
        City = "Jos";
        Month = "January";
    
    elseif (k==86)
        City = "Jos";
        Month = "February";

    elseif (k==87)
        City = "Jos";
        Month = "March";

    elseif (k==88)
        City = "Jos";
        Month = "April";

    elseif (k==89)
        City = "Jos";
        Month = "May";

    elseif (k==90)
        City = "Jos";
        Month = "June";

    elseif (k==91)
        City = "Jos";
        Month = "July";

    elseif (k==92)
        City = "Jos";
        Month = "August";

    elseif (k==93)
        City = "Jos";
        Month = "September";

    elseif (k==94)
        City = "Jos";
        Month = "October";

    elseif (k==95)
        City = "Jos";
        Month = "November";

    elseif (k==96)
        City = "Jos";
        Month = "December";
        
    elseif (k==97)
        City = "Sokoto";
        Month = "January";
    
    elseif (k==98)
        City = "Sokoto";
        Month = "February";

    elseif (k==99)
        City = "Sokoto";
        Month = "March";

    elseif (k==100)
        City = "Sokoto";
        Month = "April";

    elseif (k==101)
        City = "Sokoto";
        Month = "May";

    elseif (k==102)
        City = "Sokoto";
        Month = "June";

    elseif (k==103)
        City = "Sokoto";
        Month = "July";

    elseif (k==104)
        City = "Sokoto";
        Month = "August";

    elseif (k==105)
        City = "Sokoto";
        Month = "September";

    elseif (k==106)
        City = "Sokoto";
        Month = "October";

    elseif (k==107)
        City = "Sokoto";
        Month = "November";

    elseif (k==108)
        City = "Sokoto";
        Month = "December";
    end
end

%% Nothern Cities Average Hourly Dry Bulb Temperature Data
function [Temp_Data,Temp_data_raw] = Northern_Data()

    %*****YOLA****
    Jan_Yola = [22.39,21.65,20.51,19.44,18.96,19.3,20.36,21.95,23.99,26.47,29.24,31.87,33.78,34.6,34.36,33.44,32.18,30.7,28.93,26.93,25.0,23.59,22.89,22.67];
    Feb_Yola = [25.55,24.84,23.75,22.73,22.26,22.59,23.6,25.13,27.09,29.48,32.14,34.66,36.49,37.28,37.05,36.16,34.95,33.53,31.84,29.92,28.07,26.71,26.03,25.82];
    Mar_Yola = [29.42,28.77,27.76,26.82,26.39,26.69,27.63,29.04,30.85,33.04,35.5,37.83,39.52,40.24,40.04,39.22,38.1,36.79,35.23,33.45,31.74,30.49,29.87,29.67];
    Apr_Yola = [30.68,30.15,29.32,28.54,28.19,28.44,29.21,30.37,31.85,33.66,35.67,37.59,38.98,39.57,39.4,38.73,37.81,36.73,35.45,33.99,32.59,31.56,31.05,30.89];
    May_Yola = [28.76,28.33,27.66,27.02,26.74,26.94,27.56,28.51,29.71,31.18,32.82,34.37,35.5,35.98,35.84,35.3,34.55,33.68,32.63,31.45,30.31,29.47,29.06,28.93];
    Jun_Yola = [26.5,26.16,25.63,25.14,24.91,25.07,25.56,26.3,27.24,28.4,29.68,30.9,31.79,32.17,32.06,31.63,31.04,30.36,29.54,28.61,27.71,27.06,26.73,26.63];
    Jul_Yola = [25.77,25.47,25.0,24.56,24.36,24.5,24.93,25.59,26.44,27.46,28.61,29.7,30.49,30.83,30.73,30.35,29.83,29.21,28.48,27.65,26.86,26.27,25.98,25.89];
    Aug_Yola = [25.07,24.79,24.37,23.97,23.79,23.92,24.31,24.9,25.66,26.59,27.62,28.6,29.31,29.61,29.52,29.18,28.71,28.16,27.5,26.76,26.04,25.51,25.25,25.17];
    Sep_Yola = [25.13,24.81,24.32,23.85,23.64,23.79,24.25,24.94,25.83,26.9,28.11,29.25,30.08,30.43,30.33,29.93,29.38,28.74,27.97,27.1,26.27,25.65,25.35,25.25];
    Oct_Yola = [25.8,25.39,24.76,24.17,23.91,24.09,24.68,25.56,26.68,28.05,29.58,31.03,32.09,32.54,32.41,31.9,31.2,30.39,29.41,28.31,27.24,26.46,26.07,25.95];
    Nov_Yola = [23.82,23.08,21.95,20.88,20.4,20.74,21.79,23.38,25.42,27.89,30.65,33.27,35.18,35.99,35.76,34.84,33.58,32.1,30.35,28.35,26.43,25.01,24.32,24.1];
    Dec_Yola = [21.81,21.02,19.81,18.67,18.16,18.51,19.65,21.34,23.52,26.16,29.11,31.91,33.95,34.82,34.57,33.58,32.24,30.66,28.78,26.65,24.59,23.09,22.34,22.1];
    
    %*****BAUCHI*****
    Jan_Bauchi = [19.89,19.13,17.97,16.87,16.38,16.72,17.81,19.44,21.53,24.07,26.91,29.6,31.56,32.4,32.16,31.21,29.92,28.4,26.59,24.54,22.57,21.12,20.4,20.17];
    Feb_Bauchi = [23.25,22.53,21.43,20.38,19.91,20.24,21.28,22.83,24.82,27.24,29.95,32.51,34.37,35.17,34.94,34.04,32.81,31.36,29.64,27.69,25.81,24.43,23.75,23.53];
    Mar_Bauchi = [26.49,25.81,24.76,23.78,23.33,23.64,24.62,26.09,27.96,30.25,32.8,35.22,36.98,37.73,37.51,36.66,35.5,34.14,32.51,30.67,28.9,27.59,26.95,26.74];
    Apr_Bauchi = [28.47,27.87,26.95,26.08,25.68,25.96,26.82,28.12,29.78,31.8,34.06,36.2,37.76,38.42,38.23,37.48,36.45,35.25,33.81,32.18,30.61,29.45,28.88,28.7];
    May_Bauchi = [27.68,27.19,26.44,25.73,25.41,25.63,26.34,27.4,28.75,30.4,32.24,33.99,35.26,35.8,35.65,35.03,34.19,33.21,32.04,30.71,29.42,28.48,28.02,27.87];
    Jun_Bauchi = [25.39,24.99,24.36,23.77,23.5,23.69,24.27,25.15,26.28,27.65,29.19,30.64,31.69,32.15,32.02,31.5,30.81,29.99,29.01,27.91,26.84,26.06,25.67,25.55];
    Jul_Bauchi = [24.26,23.94,23.43,22.95,22.74,22.89,23.36,24.07,24.98,26.09,27.33,28.5,29.35,29.72,29.61,29.2,28.64,27.97,27.19,26.29,25.43,24.8,24.49,24.39];
    Aug_Bauchi = [23.44,23.13,22.65,22.21,22.0,22.14,22.59,23.25,24.11,25.15,26.31,27.41,28.2,28.55,28.45,28.06,27.53,26.91,26.18,25.34,24.53,23.94,23.65,23.55];
    Sep_Bauchi = [23.99,23.61,23.02,22.47,22.23,22.4,22.95,23.76,24.81,26.09,27.51,28.86,29.84,30.26,30.14,29.67,29.02,28.26,27.35,26.32,25.33,24.6,24.25,24.13];
    Oct_Bauchi = [24.83,24.31,23.52,22.78,22.44,22.68,23.42,24.53,25.95,27.67,29.6,31.43,32.76,33.33,33.17,32.53,31.65,30.62,29.39,27.99,26.65,25.67,25.18,25.02];
    Nov_Bauchi = [23.0,22.32,21.26,20.27,19.82,20.13,21.12,22.6,24.5,26.8,29.38,31.82,33.6,34.36,34.14,33.28,32.11,30.73,29.09,27.23,25.44,24.12,23.47,23.26];
    Dec_Bauchi = [20.44,19.69,18.52,17.43,16.93,17.28,18.36,20.0,22.09,24.63,27.47,30.17,32.13,32.96,32.72,31.77,30.48,28.96,27.16,25.1,23.13,21.68,20.96,20.73];

    %*****MAIDUGURI****
    Jan_Maiduguri = [17.99,17.13,15.81,14.58,14.02,14.41,15.64,17.48,19.84,22.72,25.93,28.97,31.19,32.13,31.86,30.79,29.33,27.61,25.57,23.25,21.02,19.38,18.57,18.31];
    Feb_Maiduguri = [21.05,20.21,18.92,17.71,17.16,17.54,18.75,20.55,22.87,25.69,28.84,31.82,33.99,34.92,34.65,33.6,32.17,30.49,28.49,26.21,24.02,22.41,21.62,21.37];
    Mar_Maiduguri = [25.29,24.45,23.17,21.96,21.42,21.8,23.0,24.8,27.1,29.91,33.04,36.01,38.17,39.09,38.83,37.78,36.36,34.68,32.69,30.43,28.25,26.65,25.86,25.6];
    Apr_Maiduguri = [29.46,28.75,27.65,26.62,26.15,26.48,27.5,29.04,31.0,33.4,36.07,38.6,40.44,41.23,41.01,40.11,38.9,37.47,35.77,33.84,31.98,30.62,29.94,29.73];
    May_Maiduguri = [29.95,29.36,28.44,27.58,27.19,27.46,28.32,29.6,31.24,33.24,35.47,37.59,39.12,39.78,39.59,38.85,37.83,36.64,35.22,33.61,32.06,30.92,30.36,30.18];
    Jun_Maiduguri = [28.04,27.58,26.87,26.2,25.9,26.11,26.78,27.77,29.05,30.6,32.34,33.99,35.18,35.69,35.55,34.97,34.18,33.25,32.15,30.89,29.68,28.8,28.36,28.22];
    Jul_Maiduguri = [26.21,25.84,25.27,24.73,24.49,24.65,25.19,25.99,27.01,28.26,29.66,30.98,31.94,32.35,32.23,31.77,31.13,30.39,29.5,28.49,27.52,26.81,26.46,26.35];
    Aug_Maiduguri = [24.92,24.61,24.14,23.69,23.48,23.63,24.07,24.74,25.6,26.64,27.8,28.9,29.7,30.05,29.95,29.56,29.03,28.41,27.67,26.83,26.02,25.43,25.13,25.04];
    Sep_Maiduguri = [25.54,25.14,24.51,23.92,23.66,23.84,24.43,25.3,26.43,27.79,29.32,30.76,31.81,32.26,32.13,31.63,30.93,30.12,29.15,28.04,26.98,26.21,25.82,25.7];
    Oct_Maiduguri = [25.3,24.67,23.69,22.78,22.36,22.65,23.56,24.93,26.68,28.8,31.18,33.43,35.07,35.77,35.56,34.77,33.69,32.42,30.91,29.19,27.54,26.33,25.73,25.54];
    Nov_Maiduguri = [21.66,20.81,19.51,18.28,17.73,18.11,19.33,21.16,23.5,26.35,29.54,32.56,34.75,35.69,35.42,34.36,32.91,31.21,29.18,26.88,24.67,23.04,22.24,21.98];
    Dec_Maiduguri = [18.49,17.6,16.23,14.94,14.36,14.77,16.04,17.96,20.42,23.41,26.75,29.91,32.21,33.2,32.91,31.8,30.28,28.5,26.37,23.96,21.64,19.93,19.09,18.82];
    
    %****KADUNA****
    Jan_Kaduna = [19.61,18.91,17.82,16.8,16.34,16.66,17.67,19.2,21.15,23.52,26.18,28.69,30.52,31.3,31.08,30.19,28.98,27.57,25.88,23.96,22.12,20.76,20.1,19.88];
    Feb_Kaduna = [22.2,21.5,20.42,19.41,18.95,19.27,20.27,21.78,23.72,26.07,28.7,31.19,33.0,33.77,33.55,32.67,31.48,30.08,28.4,26.51,24.68,23.34,22.67,22.46];
    Mar_Kaduna = [25.49,24.83,23.8,22.84,22.41,22.71,23.67,25.1,26.93,29.16,31.65,34.01,35.73,36.46,36.25,35.42,34.29,32.96,31.37,29.57,27.84,26.57,25.94,25.74];
    Apr_Kaduna = [26.36,25.82,24.99,24.2,23.85,24.1,24.87,26.04,27.53,29.35,31.38,33.3,34.7,35.3,35.13,34.45,33.53,32.44,31.15,29.69,28.27,27.24,26.73,26.56];
    May_Kaduna = [24.91,24.46,23.76,23.1,22.81,23.01,23.67,24.64,25.9,27.42,29.12,30.73,31.91,32.41,32.26,31.7,30.92,30.01,28.93,27.7,26.52,25.65,25.22,25.08];
    Jun_Kaduna = [23.33,22.95,22.36,21.8,21.55,21.73,22.28,23.1,24.16,25.45,26.89,28.25,29.24,29.66,29.54,29.06,28.41,27.64,26.73,25.69,24.69,23.95,23.59,23.48];
    Jul_Kaduna = [22.73,22.4,21.9,21.43,21.22,21.37,21.84,22.54,23.43,24.53,25.74,26.9,27.74,28.1,28.0,27.59,27.04,26.38,25.61,24.73,23.88,23.26,22.95,22.85];
    Aug_Kaduna = [22.3,22.01,21.55,21.13,20.93,21.07,21.49,22.13,22.95,23.94,25.05,26.1,26.86,27.19,27.1,26.73,26.22,25.63,24.92,24.12,23.35,22.78,22.51,22.42];
    Sep_Kaduna = [22.45,22.06,21.47,20.91,20.66,20.84,21.39,22.22,23.29,24.58,26.03,27.4,28.39,28.82,28.7,28.21,27.56,26.78,25.86,24.82,23.81,23.08,22.71,22.59];
    Oct_Kaduna = [22.84,22.37,21.65,20.96,20.65,20.87,21.55,22.57,23.87,25.46,27.23,28.92,30.14,30.66,30.51,29.92,29.11,28.16,27.04,25.75,24.52,23.61,23.17,23.02];
    Nov_Kaduna = [21.15,20.44,19.33,18.29,17.82,18.15,19.18,20.73,22.72,25.13,27.83,30.39,32.24,33.04,32.81,31.91,30.68,29.24,27.53,25.58,23.7,22.32,21.64,21.43];
    Dec_Kaduna = [19.62,18.89,17.75,16.69,16.21,16.54,17.6,19.19,21.22,23.69,26.45,29.07,30.97,31.78,31.55,30.63,29.37,27.9,26.14,24.15,22.23,20.82,20.12,19.9];
    
    %****KANO****
    Jan_Kano = [18.1,17.37,16.25,15.19,14.71,15.05,16.1,17.67,19.68,22.14,24.88,27.47,29.36,30.17,29.93,29.02,27.78,26.31,24.57,22.59,20.69,19.29,18.6,18.38];
    Feb_Kano = [20.91,20.18,19.04,17.98,17.5,17.83,18.89,20.48,22.51,24.98,27.74,30.36,32.26,33.07,32.84,31.92,30.66,29.19,27.43,25.44,23.52,22.11,21.41,21.19];
    Mar_Kano = [25.09,24.37,23.25,22.2,21.73,22.06,23.1,24.67,26.67,29.1,31.82,34.4,36.28,37.08,36.85,35.94,34.71,33.25,31.52,29.55,27.66,26.27,25.59,25.37];
    Apr_Kano = [28.64,27.99,27.0,26.06,25.64,25.93,26.86,28.26,30.05,32.22,34.65,36.95,38.63,39.34,39.14,38.33,37.22,35.92,34.38,32.62,30.93,29.69,29.08,28.88];
    May_Kano = [28.54,27.98,27.13,26.32,25.96,26.21,27.01,28.21,29.75,31.62,33.71,35.69,37.13,37.74,37.57,36.87,35.92,34.8,33.47,31.96,30.51,29.44,28.92,28.75];
    Jun_Kano = [26.43,25.97,25.26,24.6,24.29,24.5,25.17,26.16,27.43,28.98,30.71,32.35,33.55,34.06,33.91,33.33,32.55,31.62,30.52,29.27,28.07,27.18,26.75,26.61];
    Jul_Kano = [24.72,24.35,23.77,23.23,22.99,23.16,23.7,24.5,25.53,26.78,28.18,29.51,30.47,30.89,30.77,30.3,29.66,28.92,28.03,27.01,26.04,25.33,24.97,24.86];
    Aug_Kano = [23.96,23.63,23.11,22.63,22.41,22.56,23.04,23.77,24.69,25.82,27.08,28.28,29.15,29.52,29.41,28.99,28.42,27.74,26.94,26.03,25.15,24.51,24.19,24.09];
    Sep_Kano = [24.68,24.27,23.64,23.05,22.78,22.96,23.55,24.44,25.57,26.95,28.5,29.96,31.02,31.47,31.34,30.83,30.13,29.3,28.32,27.21,26.14,25.35,24.96,24.84];
    Oct_Kano = [24.91,24.34,23.47,22.64,22.27,22.53,23.35,24.57,26.14,28.05,30.19,32.21,33.68,34.31,34.13,33.42,32.45,31.31,29.95,28.41,26.92,25.83,25.3,25.12];
    Nov_Kano = [21.69,20.95,19.82,18.76,18.28,18.61,19.67,21.25,23.28,25.75,28.51,31.12,33.02,33.83,33.6,32.68,31.43,29.95,28.2,26.21,24.29,22.88,22.19,21.96];
    Dec_Kano = [18.56,17.8,16.64,15.55,15.06,15.4,16.49,18.11,20.2,22.73,25.56,28.24,30.2,31.03,30.79,29.85,28.56,27.05,25.25,23.2,21.23,19.79,19.07,18.84];
    
    %****ILORIN****
    Jan_Ilorin = [23.72,23.09,22.12,21.2,20.79,21.08,21.98,23.35,25.09,27.21,29.58,31.83,33.46,34.16,33.96,33.17,32.09,30.82,29.32,27.6,25.96,24.74,24.15,23.96];
    Feb_Ilorin = [25.7,25.09,24.15,23.26,22.86,23.14,24.02,25.34,27.03,29.08,31.37,33.55,35.13,35.81,35.61,34.85,33.8,32.58,31.12,29.46,27.86,26.69,26.11,25.93];
    Mar_Ilorin = [26.73,26.18,25.34,24.55,24.19,24.44,25.23,26.41,27.92,29.75,31.8,33.75,35.16,35.77,35.59,34.91,33.98,32.88,31.57,30.09,28.67,27.62,27.1,26.94];
    Apr_Ilorin = [26.19,25.75,25.08,24.44,24.16,24.36,24.99,25.93,27.15,28.62,30.27,31.83,32.97,33.45,33.31,32.76,32.01,31.13,30.08,28.89,27.75,26.91,26.49,26.36];
    May_Ilorin = [25.25,24.88,24.31,23.77,23.53,23.7,24.23,25.03,26.06,27.3,28.69,30.01,30.97,31.38,31.26,30.8,30.17,29.42,28.54,27.53,26.56,25.85,25.5,25.39];
    Jun_Ilorin = [24.16,23.83,23.33,22.86,22.65,22.8,23.26,23.96,24.86,25.96,27.18,28.33,29.18,29.54,29.43,29.02,28.47,27.82,27.04,26.16,25.31,24.69,24.38,24.28];
    Jul_Ilorin = [23.58,23.31,22.89,22.49,22.31,22.43,22.83,23.42,24.18,25.1,26.12,27.1,27.81,28.11,28.02,27.68,27.21,26.66,26.01,25.27,24.55,24.03,23.77,23.68];
    Aug_Ilorin = [23.23,22.97,22.57,22.19,22.02,22.14,22.51,23.08,23.8,24.68,25.66,26.59,27.27,27.56,27.47,27.15,26.7,26.18,25.55,24.84,24.16,23.66,23.41,23.33];
    Sep_Ilorin = [23.46,23.15,22.67,22.22,22.02,22.16,22.61,23.28,24.14,25.19,26.36,27.47,28.27,28.62,28.52,28.13,27.6,26.97,26.23,25.38,24.57,23.97,23.68,23.58];
    Oct_Ilorin = [23.99,23.62,23.04,22.5,22.26,22.43,22.97,23.77,24.81,26.06,27.47,28.8,29.77,30.18,30.06,29.59,28.96,28.2,27.31,26.3,25.32,24.6,24.25,24.14];
    Nov_Ilorin = [24.41,23.88,23.06,22.29,21.94,22.18,22.95,24.1,25.57,27.35,29.35,31.24,32.62,33.21,33.04,32.37,31.47,30.4,29.13,27.68,26.3,25.28,24.77,24.61];
    Dec_Ilorin = [23.61,22.98,22.01,21.1,20.68,20.97,21.88,23.24,24.98,27.1,29.47,31.71,33.35,34.04,33.84,33.05,31.98,30.71,29.2,27.49,25.85,24.64,24.04,23.85];
    
    
    %****MINNA****
    Jan_Minna = [24.55,23.95,23.02,22.14,21.75,22.03,22.89,24.19,25.86,27.88,30.15,32.29,33.85,34.52,34.33,33.57,32.54,31.34,29.9,28.26,26.69,25.53,24.96,24.78];
    Feb_Minna = [27.08,26.49,25.56,24.69,24.3,24.58,25.44,26.73,28.39,30.4,32.65,34.79,36.34,37.0,36.81,36.06,35.04,33.83,32.4,30.78,29.21,28.06,27.49,27.31];
    Mar_Minna = [29.28,28.74,27.91,27.13,26.77,27.02,27.8,28.97,30.46,32.28,34.31,36.24,37.64,38.23,38.06,37.38,36.46,35.38,34.08,32.62,31.2,30.17,29.65,29.49];
    Apr_Minna = [28.65,28.18,27.46,26.77,26.46,26.68,27.36,28.38,29.68,31.27,33.04,34.73,35.95,36.47,36.32,35.73,34.92,33.97,32.85,31.56,30.33,29.42,28.98,28.83];
    May_Minna = [26.67,26.29,25.71,25.15,24.9,25.08,25.63,26.45,27.51,28.79,30.22,31.58,32.57,32.99,32.87,32.39,31.74,30.97,30.06,29.03,28.03,27.3,26.94,26.82];
    Jun_Minna = [24.97,24.64,24.12,23.63,23.41,23.56,24.05,24.77,25.7,26.83,28.09,29.29,30.16,30.53,30.43,30.0,29.43,28.76,27.95,27.04,26.16,25.52,25.2,25.1];
    Jul_Minna = [24.2,23.91,23.48,23.07,22.88,23.01,23.42,24.03,24.81,25.76,26.82,27.83,28.56,28.88,28.79,28.43,27.95,27.38,26.71,25.94,25.2,24.66,24.39,24.3];
    Aug_Minna = [23.89,23.64,23.26,22.9,22.74,22.85,23.21,23.74,24.42,25.25,26.18,27.05,27.69,27.96,27.89,27.58,27.16,26.66,26.07,25.4,24.76,24.29,24.05,23.98];
    Sep_Minna = [24.03,23.72,23.25,22.81,22.61,22.75,23.19,23.85,24.69,25.71,26.86,27.95,28.73,29.07,28.98,28.59,28.07,27.46,26.73,25.9,25.11,24.52,24.23,24.14];
    Oct_Minna = [24.54,24.14,23.52,22.94,22.68,22.86,23.44,24.31,25.42,26.77,28.28,29.71,30.75,31.19,31.07,30.56,29.88,29.07,28.11,27.02,25.97,25.2,24.82,24.7];
    Nov_Minna = [24.62,24.0,23.05,22.15,21.75,22.03,22.92,24.25,25.96,28.04,30.37,32.57,34.17,34.85,34.66,33.88,32.83,31.58,30.11,28.43,26.81,25.63,25.04,24.85];
    Dec_Minna = [24.1,23.43,22.4,21.43,20.99,21.3,22.26,23.7,25.55,27.8,30.32,32.7,34.43,35.17,34.96,34.12,32.98,31.63,30.04,28.22,26.47,25.19,24.55,24.35];

    %****JOS****
    Jan_Jos = [15.64,14.92,13.8,12.75,12.27,12.6,13.65,15.22,17.22,19.67,22.39,24.98,26.86,27.66,27.43,26.52,25.28,23.82,22.09,20.12,18.22,16.83,16.14,15.92];
    Feb_Jos = [18.1,17.38,16.28,15.24,14.77,15.1,16.13,17.68,19.66,22.07,24.76,27.31,29.17,29.96,29.73,28.84,27.61,26.17,24.46,22.52,20.64,19.27,18.59,18.37];
    Mar_Jos = [20.85,20.21,19.23,18.3,17.88,18.17,19.09,20.47,22.24,24.4,26.8,29.08,30.74,31.44,31.24,30.44,29.35,28.06,26.53,24.79,23.12,21.89,21.29,21.09];
    Apr_Jos = [21.38,20.89,20.14,19.43,19.11,19.33,20.04,21.1,22.45,24.1,25.94,27.69,28.96,29.5,29.35,28.73,27.89,26.91,25.74,24.41,23.12,22.18,21.72,21.57];
    May_Jos = [20.5,20.13,19.55,19.01,18.77,18.94,19.47,20.28,21.31,22.56,23.96,25.29,26.26,26.67,26.55,26.08,25.45,24.7,23.81,22.8,21.82,21.11,20.75,20.64];
    Jun_Jos = [19.59,19.27,18.78,18.32,18.11,18.25,18.71,19.4,20.28,21.35,22.55,23.68,24.51,24.86,24.76,24.36,23.82,23.18,22.42,21.55,20.72,20.11,19.8,19.71];
    Jul_Jos = [18.57,18.32,17.95,17.6,17.44,17.55,17.9,18.42,19.1,19.92,20.84,21.7,22.34,22.61,22.53,22.22,21.81,21.32,20.73,20.07,19.43,18.97,18.73,18.66];
    Aug_Jos = [18.51,18.27,17.89,17.53,17.37,17.48,17.84,18.37,19.04,19.87,20.79,21.67,22.3,22.57,22.5,22.19,21.77,21.28,20.69,20.02,19.38,18.91,18.68,18.6];
    Sep_Jos = [18.88,18.54,18.02,17.53,17.31,17.46,17.95,18.68,19.62,20.75,22.02,23.23,24.11,24.48,24.37,23.95,23.37,22.69,21.88,20.96,20.08,19.43,19.11,19.01];
    Oct_Jos = [19.03,18.57,17.86,17.19,16.88,17.09,17.76,18.76,20.03,21.58,23.32,24.96,26.16,26.67,26.52,25.94,25.15,24.23,23.12,21.87,20.67,19.78,19.34,19.2];
    Nov_Jos = [17.17,16.52,15.53,14.6,14.18,14.47,15.4,16.79,18.57,20.74,23.16,25.46,27.13,27.84,27.63,26.83,25.73,24.43,22.89,21.14,19.46,18.22,17.61,17.41];
    Dec_Jos = [16.03,15.31,14.21,13.17,12.7,13.03,14.06,15.61,17.59,20.0,22.69,25.25,27.1,27.9,27.67,26.77,25.55,24.11,22.39,20.45,18.57,17.2,16.52,16.3];
    
    %****SOKOTO****
    Jan_Sokoto = [21.04,20.33,19.23,18.19,17.73,18.05,19.08,20.62,22.6,25.0,27.68,30.22,32.07,32.86,32.64,31.74,30.52,29.09,27.38,25.44,23.58,22.21,21.53,21.31];
    Feb_Sokoto = [23.56,22.83,21.71,20.66,20.18,20.51,21.56,23.13,25.15,27.6,30.33,32.93,34.81,35.62,35.39,34.48,33.23,31.77,30.03,28.05,26.15,24.75,24.06,23.84];
    Mar_Sokoto = [27.73,27.03,25.94,24.92,24.46,24.78,25.79,27.32,29.27,31.64,34.3,36.81,38.64,39.42,39.2,38.31,37.1,35.69,34.0,32.08,30.24,28.88,28.22,28.0];
    Apr_Sokoto = [31.21,30.59,29.63,28.73,28.33,28.61,29.5,30.85,32.56,34.66,36.99,39.2,40.81,41.5,41.31,40.53,39.46,38.22,36.73,35.04,33.42,32.23,31.64,31.45];
    May_Sokoto = [30.62,30.1,29.3,28.54,28.2,28.44,29.19,30.32,31.76,33.51,35.47,37.32,38.67,39.25,39.08,38.43,37.54,36.49,35.25,33.83,32.47,31.47,30.98,30.82];
    Jun_Sokoto = [28.65,28.19,27.49,26.83,26.53,26.74,27.4,28.38,29.64,31.17,32.88,34.5,35.68,36.19,36.04,35.47,34.69,33.78,32.69,31.45,30.26,29.39,28.96,28.82];
    Jul_Sokoto = [26.43,26.06,25.48,24.94,24.7,24.87,25.41,26.21,27.24,28.49,29.89,31.22,32.18,32.6,32.48,32.01,31.37,30.63,29.74,28.72,27.75,27.04,26.68,26.57];
    Aug_Sokoto = [25.19,24.86,24.36,23.88,23.67,23.82,24.29,24.99,25.9,27.0,28.22,29.39,30.24,30.6,30.49,30.08,29.53,28.87,28.09,27.2,26.35,25.72,25.41,25.31];
    Sep_Sokoto = [25.85,25.45,24.82,24.24,23.97,24.16,24.74,25.61,26.73,28.09,29.61,31.05,32.1,32.55,32.42,31.91,31.22,30.41,29.44,28.34,27.29,26.51,26.13,26.0];
    Oct_Sokoto = [27.17,26.57,25.65,24.78,24.39,24.66,25.52,26.81,28.47,30.48,32.73,34.86,36.41,37.08,36.89,36.14,35.11,33.91,32.48,30.85,29.29,28.14,27.58,27.39];
    Nov_Sokoto = [25.42,24.67,23.52,22.44,21.95,22.29,23.37,24.98,27.04,29.55,32.35,35.0,36.94,37.76,37.52,36.59,35.32,33.82,32.04,30.01,28.06,26.63,25.92,25.7];
    Dec_Sokoto = [22.5,21.75,20.61,19.53,19.05,19.39,20.45,22.06,24.11,26.62,29.41,32.06,33.98,34.8,34.57,33.63,32.37,30.87,29.1,27.08,25.14,23.71,23.01,22.78];
        
    
    
    
    Temp_data_raw = vertcat(Jan_Yola,Feb_Yola,Mar_Yola,Apr_Yola,May_Yola,Jun_Yola,Jul_Yola,Aug_Yola,Sep_Yola,Oct_Yola,Nov_Yola,Dec_Yola,Jan_Bauchi,Feb_Bauchi,Mar_Bauchi,Apr_Bauchi,May_Bauchi,Jun_Bauchi,Jul_Bauchi,Aug_Bauchi,Sep_Bauchi,Oct_Bauchi,Nov_Bauchi,Dec_Bauchi,Jan_Maiduguri,Feb_Maiduguri,Mar_Maiduguri,Apr_Maiduguri,May_Maiduguri,Jun_Maiduguri,Jul_Maiduguri,Aug_Maiduguri,Sep_Maiduguri,Oct_Maiduguri,Nov_Maiduguri,Dec_Maiduguri,Jan_Kaduna,Feb_Kaduna,Mar_Kaduna,Apr_Kaduna,May_Kaduna,Jun_Kaduna,Jul_Kaduna,Aug_Kaduna,Sep_Kaduna,Oct_Kaduna,Nov_Kaduna,Dec_Kaduna,Jan_Kano,Feb_Kano,Mar_Kano,Apr_Kano,May_Kano,Jun_Kano,Jul_Kano,Aug_Kano,Sep_Kano,Oct_Kano,Nov_Kano,Dec_Kano,Jan_Ilorin,Feb_Ilorin,Mar_Ilorin,Apr_Ilorin,May_Ilorin,Jun_Ilorin,Jul_Ilorin,Aug_Ilorin,Sep_Ilorin,Oct_Ilorin,Nov_Ilorin,Dec_Ilorin,Jan_Minna,Feb_Minna,Mar_Minna,Apr_Minna,May_Minna,Jun_Minna,Jul_Minna,Aug_Minna,Sep_Minna,Oct_Minna,Nov_Minna,Dec_Minna,Jan_Jos,Feb_Jos,Mar_Jos,Apr_Jos,May_Jos,Jun_Jos,Jul_Jos,Aug_Jos,Sep_Jos,Oct_Jos,Nov_Jos,Dec_Jos,Jan_Sokoto,Feb_Sokoto,Mar_Sokoto,Apr_Sokoto,May_Sokoto,Jun_Sokoto,Jul_Sokoto,Aug_Sokoto,Sep_Sokoto,Oct_Sokoto,Nov_Sokoto,Dec_Sokoto);
    
    %%Standardized Temperature
    row = size(Temp_data_raw,1);
    col = size(Temp_data_raw,2);
    Temp_Data = zeros(size(Temp_data_raw));
    for r=1:row
        for c=1:col
        Temp_Data(r,c) = (Temp_data_raw(r,c)-mean(Temp_data_raw(r,:)))/(max(Temp_data_raw(r,:))-min(Temp_data_raw(r,:)));
        end
    end
    
    Temp_Data;
    %min(Feb_Sokoto(:))

end

%% Raw Temperature Data
function [Rng_adj,Tmin] = Estimation_Data()
    %%*****YOLA****
    Rng = [15.39,14.78,13.64,11.20,9.10,7.14,6.37,5.73,6.68,8.50,15.35,16.40;
           15.77,15.02,14.17,12.54,10.23,8.51,6.87,6.44,7.91,10.72,14.31,15.78;
           17.83,17.48,17.40,14.84,12.39,9.64,7.74,6.46,8.47,13.19,17.68,18.54;
           14.73,14.59,13.83,11.27,9.45,7.98,6.77,6.16,8.03,9.85,14.98,15.33;
           15.21,15.33,15.11,13.49,11.60,9.61,7.77,7.00,8.56,11.85,15.31,15.72;
           13.16,12.74,11.39,9.15,7.73,6.78,5.71,5.45,6.50,7.80,11.09,13.15;
           12.57,12.50,11.28,9.85,7.96,7.01,5.90,5.14,6.36,8.38,12.90,13.96;
           15.15,14.95,13.35,10.23,7.78,6.65,5.09,5.12,7.06,9.63,13.45,14.96;
           14.90,15.20,14.73,12.97,10.87,9.50,7.77,6.82,8.44,12.49,15.56,15.51;
           ]';
      
    Rng_adj = Rng(:);
       
    Avg = [26.05,29.07,32.67,33.35,30.93,28.20,27.29,26.43,26.72,27.82,27.47,25.71;
           23.64,26.83,29.86,31.46,30.12,27.42,25.90,24.97,25.87,27.38,26.41,24.20;
           22.23,25.21,29.43,32.99,32.90,30.34,28.05,26.46,27.56,28.44,25.87,22.90;
           23.12,25.67,28.78,29.04,27.16,25.23,24.34,23.77,24.36,25.19,24.72,23.27;
           21.72,24.56,28.69,31.85,31.30,28.72,26.57,25.63,26.72,27.73,25.33,22.30;
           26.85,28.73,29.44,28.37,27.09,25.77,24.94,24.53,25.01,25.85,27.05,26.74;
           27.54,30.06,31.97,31.00,28.57,26.64,25.60,25.11,25.54,26.54,27.69,27.42;
           19.25,21.66,24.03,23.82,22.35,21.17,19.78,19.73,20.56,21.32,20.37,19.59;
           24.59,27.18,31.24,34.30,33.21,30.91,28.28,26.81,27.86,30.14,29.12,26.19;
           ]';
    
    Avg_adj = Avg(:);
    Tmin = (Avg_adj - Rng_adj./2);
    
    size(Rng_adj);

end

%% Parameter function
function [rng,Tavg,N,Temp_Data,Temp_data_raw] = par()
    [Temp_Data,Temp_data_raw] = Northern_Data();
    k = kfunc();
    N = 12;
    Tavg = mean(Temp_data_raw(k,:));
    rng = max(Temp_data_raw(k,:))-min(Temp_data_raw(k,:));
end

%% G value
function G = NoM()
    [~,~,~,Temp_Data,~] = par();
    G = size(Temp_Data,1);
end
%% vector A
function A = vectorA()
    G = NoM();
    N = 12; % Maximum number of harmonics
    T = 24;
    A = zeros(26,26);
    ao = 1;
    aa = 1;
  %row 1-12, col 1-24

 for m=1:G     
    for n1=1:N   
        for t1=1:T
            sum_s = 2*sin(n1*pi*t1/12);
            sum_c = 2*cos(n1*pi*t1/12);
            if n1<=1
                A(1,1) = A(1,1) +( 2* aa);
            end
            A(1,n1+1) = A(1,n1+1) + sum_s;
            A(1,13+n1) = A(1,13+n1) + sum_c;
        end



    end
end



%for nd=1:NofD 
for m2=1:G
    for x = 1:12
        for n=1:N   
            for t=1:T
                sum_s = sin(n*pi*t/12)*(2*sin(x*pi*t/12));
                sum_c = cos(n*pi*t/12)*(2*sin(x*pi*t/12));
                A(x+1,1) = A(x+1,1)+ (ao * (2*sin(x*pi*t/12)));
                A(x+1,n+1) = A(x+1,n+1) + sum_s;
                A(x+1,13+n) = A(x+1,13+n) + sum_c;
            end



        end
    end
end


for m3=1:G
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
end


    % lamda3 (equation 1)
    A(1,N*2+2) = 24/23;
    
    for n4=1:N
       for t4=1:T
           % sin (eqn 2-13) lamda3
           A(n4+1,N*2+2)= A(n4+1,N*2+2)+ (1/23)*sin(n4*pi*t4/12);
           % sin (eqn 2-13)
           A(n4+13,N*2+2)= A(n4+13,N*2+2)+ (1/23)*cos(n4*pi*t4/12);
       end
    end
    
    
    % (equation 26, lamdas 3)
    A(N*2+2:end,N*2+2) = 0;
    
    % (equation 26, a0)
    A(N*2+2,1) = 24/23;  %changed
    
    % aN & bN
    for n5=1:N
        for t5=1:T
           % sin (eqn 27)
           A(N*2+2,n5+1)= A(N*2+2,n5+1) + (1/23)*sin(n5*pi*t5/12);
           % sin (eqn 27)
           A(N*2+2,n5+13)= A(N*2+2,n5+13) + (1/23)*cos(n5*pi*t5/12);
       end
    end
    A_adj = A([1:12,14:end],[1:12,14:end]);
    size(A_adj);
    A;

end

%% vector B
function B = vectorB()  
    G = NoM();
    B =zeros(26,1);
    [~,~,~,Temp_Data,~] = par();

    T = 24;
    M =12;%(N before)
    std_sum = 0;
    sum_s = 0;
    sum_c = 0;
    std_summa = 0;
    range_sum = 0;
    
    for m=1:G
        for t=1:T 
            B(1,1) = B(1,1) + Temp_Data(m,t); 
        end
    end
    B(1,1) = 2 * B(1,1);
    
    for n=1:M
        for m1=1:G
            for t=1:T
                sum_s = sum_s + (sin(n*pi*t/12))*Temp_Data(m1,t);
                sum_c = sum_c + (cos(n*pi*t/12))*Temp_Data(m1,t); 
            end
        end
        B(n+1,1) = 2*sum_s;   
        B(n+13,1) = 2* sum_c;
        sum_s = 0;
        sum_c = 0;
        
    end
    
    % eqn 28
    for m2=1:G
        for t2=1:T
            std_sum = std_sum + (max(Temp_Data(m2,:)) -Temp_Data(m2,t2));       
        end
        range_sum = range_sum + (max(Temp_Data(m2,:)) - min(Temp_Data(m2,:)));
        std = std_sum/(23);
        std_summa = std_summa + std;
        std_sum = 0;
    end
    std_avg = std_summa / G;
    range_avg = range_sum/G;
    %range_avg = max(Temp_Data(:)) - min(Temp_Data(:))
      %change
    
    %eqn 26
    B(M*2+2,1) =((24/23*range_avg)-std_avg);
    B;
end

%%
function [ao,an,bn] = vectorX()
    [~,~,N,~,~] = par();
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

function [actual_data,predicted_temp_adj] = Temp()
    [ao,an,bn] = vectorX();
    %[rng,Tavg,N,Temp_Data,Temp_data_raw] = par();
    [~,~,N,Temp_Data,Temp_data_raw] = par();
    %[rnge,Tmine] = Estimation_Data();
    %rng= 8.51;
    %Tmin = 23.165;
    
    G = NoM();
    MAPE = zeros(G,1);
    MFE = zeros(G,1);
    MAD = zeros(G,1);
    h = zeros(G,1);
    T_min = zeros(G,1);
    T_rng = zeros(G,1);
    predict_t = zeros(G,24);
    predict_t_daily = zeros(G,1);
    predict_rng_daily = zeros(G,1);
    incr = 1;
    k = kfunc();
    
    for g=1:G
        %k = kfunc();
        t = 2;
        %N = 12;
        T = 24;
        Tmin = min(Temp_data_raw(g,:));
        rng =  max(Temp_data_raw(g,:)) -  min(Temp_data_raw(g,:));
        T_min(g,1) = Tmin;
        T_rng (g,1) = max(Temp_data_raw(g,:)) -  min(Temp_data_raw(g,:));
        %rng = rnge(g);
        %Tmin = Tmine(g);
        
        sum_an = 0;
        sum_bn = 0;
        sum_s = 0;
        sum_c = 0;
        dt = zeros(1,24);
        time = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
        predicted_temp = zeros(1,24);
        for n=1:N
            sum_an = sum_an + an(n)*sin(n*pi*t/12);
            sum_bn = sum_bn + bn(n)*cos(n*pi*t/12);
        end


        for t=1:T
            for n=1:N
                sum_s = sum_s + an(n)*sin(n*pi*t/12);
                sum_c = sum_c + bn(n)*cos(n*pi*t/12);
            end
            dt(1,t) = ao + sum_s + sum_c;
            predicted_temp(1,t) = dt(1,t);
            sum_s = 0;
            sum_c = 0;
        end
        [City,Month] = Data_Tag();

        predicted_temp;
        predicted_temp_adj = (predicted_temp * rng)+Tmin;
        %plot(time,predicted_temp_adj,time,Temp_data_raw(g,:));
        %title({City,Month}, 'FontSize',12);
        %legend('predicted','actual data');

        actual_data = Temp_data_raw(g,:);
        rng;
        %Tavg;
        %actual_data = Temp_data_raw(1,:)
        
        %To evaluate error
        n = 24;
        abs_dev = abs(actual_data - predicted_temp_adj);
        dev = actual_data - predicted_temp_adj;
        percent_error= abs(actual_data - predicted_temp_adj)./abs(actual_data);
        mean_abs_dev = sum(abs_dev)/n;
        mean_dev = sum(dev)/n;
        mean_percent_error = (sum(percent_error)/n)*100;
        
        MAPE(g,1) = mean_percent_error; %MAPE has 108 data for 12 mon * 9 loc
        MFE(g,1) = mean_dev;
        MAD(g,1) = mean_abs_dev;
        predict_t(g,:) = predicted_temp_adj;
        predict_t_daily(g,1) = mean(predicted_temp_adj);
        predict_rng_daily(g,1) = max(predicted_temp_adj)-min(predicted_temp_adj);
        
        if(mean_dev < 0)
            h(incr,1) = g;
            incr = incr +1;
        end
    end
    
    l=1;  %l helps to extract 12 month data for each location
    %w = T_min(12*(l-1)+1:12*(l-1)+12);  %To generate Tmin for our model
    w = T_rng(12*(l-1)+1:12*(l-1)+12);  %To generate Trng for our model
    %w_adj = round(w',2)
    %fprintf('%.2f\t',w)
    x = MAPE(12*(l-1)+1:12*(l-1)+12);
    y = MFE(12*(l-1)+1:12*(l-1)+12);
    z = MAD(12*(l-1)+1:12*(l-1)+12);
    %fprintf('%.4f\t',x);
    
    
    size(x);
    size(predict_rng_daily);
    size(predict_t_daily);
    predict_t;
    [max_value_MAPE,max_index_MAPE] = max(MAPE);
    [min_value_MAPE,min_index_MAPE] = min(MAPE);
    [max_value_MFE,max_index_MFE] = max(MFE);
    [min_value_MFE,min_index_MFE] = min(MFE);
    [max_value_MAD,max_index_MAD] = max(MAD);
    [min_value_MAD,min_index_MAD] = min(MAD);
    [max_value_T_avg,max_index_T_avg] = max(predict_t_daily)
    [min_value_T_rng,min_index_T_rng] = min(predict_rng_daily);
    
    %[sorted_Tavg, indices] = sort(predict_t_daily,'descend');
    %highestFour = sorted_Tavg(1:4)
    %highestFourInd = indices(1:4)
    
    h;
    count = sum(h(:)>0);
    
    
end

%% Mean Absolute Deviation

function mean_abs_dev = MAD()
    [actual_data,predicted_temp_adj] = Temp();
    n = 24;
    sum_dev = 0;
    abs_dev = abs(actual_data - predicted_temp_adj);
    dev = actual_data - predicted_temp_adj;
    percent_error= abs(actual_data - predicted_temp_adj)./abs(actual_data);
    mean_abs_dev = sum(abs_dev)/n;
    mean_dev = sum(dev)/n;
    mean_percent_error = (sum(percent_error)/n)*100;
    
    
    
    
end
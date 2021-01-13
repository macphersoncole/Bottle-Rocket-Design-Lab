%No variables in 301-4

%housekeeping
clear;
clc;
close all;

name = [];

data_303 = [];
for k = 1:8
    filename = sprintf('LA_Demo_303_%d',k);
    str = sprintf('303_%d',k);
    [t_thrust, I_sp, max_thrust] = datacleaner(filename);
    
    data_303 = [t_thrust I_sp max_thrust; data_303];
    name = [name; convertCharsToStrings(str)];
end

data_302 = [];
for k = 1:5
    filename = sprintf('LA_Demo_302_%d',k);
    str = sprintf('302_%d',k);
    [t_thrust, I_sp, max_thrust] = datacleaner(filename);
    
    data_302 = [t_thrust I_sp max_thrust; data_302];
    name = [name; convertCharsToStrings(str)];
end

data_301 = [];
for k = 1:4
    if k ~= 2
        filename = sprintf('LA_Demo_301_%d',k);
        str = sprintf('301_%d',k);
        [t_thrust, I_sp, max_thrust] = datacleaner(filename);
    
        data_301 = [t_thrust I_sp max_thrust; data_301];
        name = [name; convertCharsToStrings(str)];
    end
end

data = [data_303; data_302; data_301];
t_thrust = data(:, 1);
ISP = data(:, 2);
max_thrust = data(:, 3);

avg_t = mean(t_thrust);
avg_ISP = mean(ISP);
avg_max = mean(max_thrust);

std_t = std(t_thrust);
std_max = std(max_thrust);
std_isp = std(ISP);

% number of tests of 95% accuracy of 0.1 [s]
z = 1.96;
N_95_01 = ((z*std_isp)/(0.1))^2;

% number of tests of 97.5% accuracy of 0.1 [s]
z = 2.24;
N_975_01 = ((z*std_isp)/(0.1))^2;

% number of tests of 99% accuracy of 0.1 [s]
z = 2.58;
N_99_01 = ((z*std_isp)/(0.1))^2;

% number of tests of 95% accuracy of 0.01 [s]
z = 1.96;
N_95_001 = ((z*std_isp)/(0.01))^2;

% number of tests of 97.5% accuracy of 0.01 [s]
z = 2.24;
N_975_001 = ((z*std_isp)/(0.01))^2;

% number of tests of 99% accuracy of 0.01 [s]
z = 2.58;
N_99_001 = ((z*std_isp)/(0.01))^2;

%SEM STuff
N = [];
std_stuff = [];

for k = 2:numel(data(:,1))
   N = [N; k];
   std_stuff = [std_stuff; std(t_thrust(1:k)) std(ISP(1:k)) std(max_thrust(1:k))]; %t, isp, maxthrusyt
end

SEM_t = [];
SEM_ISP = [];
SEM_max = [];
for k = 1:numel(std_stuff(:,1))
    SEM_t = [SEM_t; std_stuff(k, 1)/sqrt(N(k))];
    SEM_ISP = [SEM_ISP; std_stuff(k, 2)/sqrt(N(k))];
    SEM_max = [SEM_max; std_stuff(k, 3)/sqrt(N(k))];
end

figure
plot(N, SEM_t)
title('SEM for Time of Thrust vs N')
xlabel('N')
ylabel('SEM of Time of Thrust')

figure
plot(N, SEM_ISP)
title('SEM for ISP vs N')
xlabel('N')
ylabel('SEM of ISP')

figure
plot(N, SEM_max)
title('SEM for Peak Thrust vs N')
xlabel('N')
ylabel('SEM of Peak Thrust')

figure
bar(ISP)
set(gca, 'XTick', 1:numel(ISP), 'XTickLabel', name);
title('ISP values from Static Test Stands')

figure
bar(t_thrust)
set(gca, 'XTick', 1:numel(t_thrust), 'XTickLabel', name);
title('Times of Thrust from Static Test Stands')

figure
bar(max_thrust)
set(gca, 'XTick', 1:numel(max_thrust), 'XTickLabel', name);
title('peak Thrusts from Static Test Stands')
function [t_thrust, I_sp, max_thrust] = datacleaner(filename)
    %Does not work for files: 303_9, 301_2
    %Weird diagrams in: 303_8, 301_4

    %setup
    data = load(filename);
    data = data * 4.44822;           %N
    summ = data(:, 3);
    %Sampling rate 1625 Hz
    f = 1625;
    t = linspace(0, numel(summ)/f, numel(summ));

    %index of t = 0.3s and t = 3.5s
    i1 = numel(t) - numel(t(t>0.3));
    i2 = numel(t(t<3.5));

    %Truncate data
    summt = summ(i1:i2);
    t = t(i1:i2);

    %tolerance of 1N
    tol = 25;

    %avg from phase 1 and 2
    avg1 = mean(summ(1:i1));
    avg2 = mean(summ(i2:end));

    tolmet = true;
    index1 = 1;
    while tolmet
        index1 = index1 + 1;
        if abs(summt(index1) - avg1) > tol
            tolmet = false;
        end
    end
    
    tol = 1;
    index2 = index1;
    tolnotmet = true;
    count = 0;
    while tolnotmet
        index2 = index2 + 1;

        %it prefires otherwise
        if abs(summt(index2)-avg2) < tol
            count = count + 1;
        else
            count = 0;
        end

        if abs(summt(index2)-avg2) < tol && count > 20
            tolnotmet = false;
        end
    end

    t_thrust = t(index1:index2);
    sum = summt(index1:index2);

    figure
    hold on
    plot(t, summt)
    yline(0);
    xline(t(index1));
    xline(t(index2));
    xlabel('Time [s]')
    ylabel('Thrust [N]')
    title('Thrust [N] vs Time [s]')

    %mass from the data file
    m = 1;           %kg
    g_0 = 9.81;           %m/s^2
    avg_offset = avg2/2;     %N
    sum = sum - avg_offset;     %adjusted
    I_sp = trapz(t_thrust, sum)/(m*g_0);
    max_thrust = max(sum);
    t_thrust = t_thrust(end) - t_thrust(1);
    %{
        figure
    hold on
    plot(t, summt-avg_offset)
    yline(0);
    xline(t(index1));
    xline(t(index2));
    xlabel('Time [s]')
    ylabel('Thrust [N]')
    title('Thrust [N] vs Time [s]')
    %}
end


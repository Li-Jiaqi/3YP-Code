function [time,segment] = extractSeg(EMG, number, time_stamp, time_window)
start = time_stamp(number);
if number==12
    ending = 215;
else
    ending = time_stamp(number+1);
end
period = find(time_window>start & time_window<ending);
segment = EMG(period);
time = time_window(period);


end
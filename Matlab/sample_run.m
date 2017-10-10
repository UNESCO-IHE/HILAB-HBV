load HourTest

avg_prec = BrueHour.Prec;
temp = BrueHour.Temp;
area = BrueHour.Area;
evap = BrueHour.Evap/24.0;
flow = BrueHour.Flow;

param = [  9.95858456e-01,   1.99174032e+00,   9.97787883e-01,...
         3.93987236e-01,   2.49999155e+02,   8.80942905e-01,...
         4.96144101e-02,   4.52884743e-01,   5.42086567e-03,...
         6.20000000e-04,   4.28494565e-01,   1.14400004e+00,...
         8.48238270e-02,   8.00013673e-01,   7.01014615e-02,...
         4.16941937e-02,   1.26898027e+00,   9.98357889e-01, 0];


% Snow switch is turn off by putting 0 in the last argument of p2, 1 is to turn on
p2 = [1, area]; % No argument implies that snow pack is on
p2 = [1, area, 0]; % 0 implies snow pack is off

ll_temp = mean(temp)*ones(length(avg_prec), 1);
v = [avg_prec(1), temp(1), evap(1), ll_temp(1)];
st = zeros(length(avg_prec), 5);
st(1,:) = [30, 30, 30, 30, 30];
[q, s1, s2, s3, s4, s5] = step_run(param, p2, v, st(1,:));

[qrout, st] = simulate(avg_prec, temp, evap, param, p2);


plot(flow, '.k')
hold on
plot(qrout,'.b');
legend('rec', 'sim')
hold off

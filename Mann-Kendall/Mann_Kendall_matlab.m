tb = csvread('data.csv', 1, 0);
[m, n] = size(tb);
ts1 = [1:m]';
ts2 = tb(:,2);
ts = [ts1,ts2];
[taub tau h sig Z S sigma sen] = ktaub(ts, 0.05)
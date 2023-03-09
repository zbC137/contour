load('data\simData.mat');
genVideo('sim', data);

load('data\lessData.mat');
genVideo('less', data);

load('data\switchData.mat');
genVideo('switch', data);

load('data\faultData.mat');
genVideo('fault', data);

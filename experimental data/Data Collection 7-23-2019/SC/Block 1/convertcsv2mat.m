clc;
clear;
close all;

data_files = dir('*.csv');

N = length(data_files);

for i=1:N
    
    filename = data_files(i).name;
    
    fid = fopen(filename, 'rt');
    for j=1:6
        parameter = textscan(fgets(fid), '%s %f %s','Delimiter', ',');
    end
    duration=parameter{2};
    for j=1:2
        parameter = textscan(fgets(fid), '%s %f %s','Delimiter', ',');
    end
    resangle=parameter{2};
    fclose(fid);
    
    data = dlmread(filename,',',35,0);
    t = data(:,1);
    theta = data(:,2);
    omega = data(:,3);
    alpha = data(:,4);
    pos = data(:,5);
    vel = data(:,6);
    acc = data(:,7);
    Fball = data(:,8);
    Fx = data(:,9);
    Fy = data(:,10);
    Fz = data(:,11);
    newfilename = strrep(data_files(i).name,'csv','mat');
    save(newfilename,'t','theta','omega','alpha','pos','vel','acc','Fball','Fx','Fy','Fz','duration','resangle');
    
end
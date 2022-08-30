load('peak times/EFb3.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 5
    peakArray{1,1}(:,3:4) = [];
end

if length(peakArray{8,1}(1,:)) == 4
    peakArray{8,1}(:,4) = [];
end
    
if length(peakArray{15,1}(1,:)) == 5
    peakArray{15,1}(:,4:5) = [];
end

if length(peakArray{16,1}(1,:)) == 5
    peakArray{16,1}(:,4:5) = [];
end

if length(peakArray{18,1}(1,:)) == 5
    peakArray{18,1}(:,4:5) = [];
end

if length(peakArray{23,1}(1,:)) == 6
    peakArray{23,1}(:,4:6) = [];
end

if length(peakArray{26,1}(1,:)) == 6
    peakArray{26,1}(:,4:6) = [];
end

if length(peakArray{27,1}(1,:)) == 5
    peakArray{27,1}(:,4:5) = [];
end

if length(peakArray{47,1}(1,:)) == 4
    peakArray{47,1}(:,4) = [];
end

if length(peakArray{50,1}(1,:)) == 4
    peakArray{50,1}(:,4) = [];
end

save('peak times/EFb3.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/EFb4.mat', 'peakArray')

if length(peakArray{9,1}(1,:)) == 4
    peakArray{9,1}(:,4) = [];
end

if length(peakArray{11,1}(1,:)) == 4
    peakArray{11,1}(:,4) = [];
end

if length(peakArray{19,1}(1,:)) == 8
    peakArray{19,1}(:,4:8) = [];
end

if length(peakArray{25,1}(1,:)) == 4
    peakArray{25,1}(:,4) = [];
end

if length(peakArray{28,1}(1,:)) == 4
    peakArray{28,1}(:,4) = [];
end

if length(peakArray{33,1}(1,:)) == 7
    peakArray{33,1}(:,4:7) = [];
end

if length(peakArray{36,1}(1,:)) == 5
    peakArray{36,1}(:,4:5) = [];
end

if length(peakArray{37,1}(1,:)) == 4
    peakArray{37,1}(:,4) = [];
end

if length(peakArray{38,1}(1,:)) == 4
    peakArray{38,1}(:,4) = [];
end

if length(peakArray{40,1}(1,:)) == 5
    peakArray{40,1}(:,4:5) = [];
end

if length(peakArray{49,1}(1,:)) == 5
    peakArray{49,1}(:,1:2) = [];
end

save('peak times/EFb4.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/NOb3.mat', 'peakArray')

if length(peakArray{2,1}(1,:)) == 4
    peakArray{2,1}(:,4) = [];
end

% Need to add trial 9 here

if length(peakArray{13,1}(1,:)) == 4
    peakArray{13,1}(:,4) = [];
end

if length(peakArray{19,1}(1,:)) == 4
    peakArray{19,1}(:,4) = [];
end

if length(peakArray{26,1}(1,:)) == 5
    peakArray{26,1}(:,4:5) = [];
end

if length(peakArray{32,1}(1,:)) == 4
    peakArray{32,1}(:,4) = [];
end

save('peak times/NOb3.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/NOb4.mat', 'peakArray')

if length(peakArray{10,1}(1,:)) == 5
    peakArray{10,1}(:,4:5) = [];
end

if length(peakArray{13,1}(1,:)) == 4
    peakArray{13,1}(:,4) = [];
end

if length(peakArray{21,1}(1,:)) == 4
    peakArray{21,1}(:,4) = [];
end

if length(peakArray{24,1}(1,:)) == 4
    peakArray{24,1}(:,4) = [];
end

if length(peakArray{25,1}(1,:)) == 4
    peakArray{25,1}(:,4) = [];
end

if length(peakArray{35,1}(1,:)) == 4
    peakArray{35,1}(:,4) = [];
end

if length(peakArray{36,1}(1,:)) == 4
    peakArray{36,1}(:,4) = [];
end

if length(peakArray{40,1}(1,:)) == 4
    peakArray{40,1}(:,4) = [];
end

if length(peakArray{41,1}(1,:)) == 5
    peakArray{41,1}(:,4:5) = [];
end

if length(peakArray{43,1}(1,:)) == 4
    peakArray{43,1}(:,4) = [];
end

if length(peakArray{50,1}(1,:)) == 4
    peakArray{50,1}(:,4) = [];
end

save('peak times/NOb4.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/MHHb3.mat', 'peakArray')

if length(peakArray{45,1}(1,:)) == 5
    peakArray{45,1}(:,4:5) = [];
end

save('peak times/MHHb3.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/MHHb4.mat', 'peakArray')

if length(peakArray{46,1}(1,:)) == 5
    peakArray{46,1}(:,4:5) = [];
end

save('peak times/MHHb4.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/CMb3.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 8
    peakArray{1,1}(:,4:8) = [];
end

if length(peakArray{2,1}(1,:)) == 4
    peakArray{2,1}(:,4) = [];
end

if length(peakArray{14,1}(1,:)) == 4
    peakArray{14,1}(:,4) = [];
end

if length(peakArray{24,1}(1,:)) == 4
    peakArray{24,1}(:,4) = [];
end

if length(peakArray{33,1}(1,:)) == 8
    peakArray{33,1}(:,2:3) = [];
    peakArray{33,1}(:,3:4) = [];
    peakArray{33,1}(:,4) = [];
end

if length(peakArray{43,1}(1,:)) == 5
    peakArray{43,1}(:,4:5) = [];
end

if length(peakArray{50,1}(1,:)) == 4
    peakArray{50,1}(:,4) = [];
end

save('peak times/CMb3.mat', 'peakArray')

clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/CMb4.mat', 'peakArray')

if length(peakArray{42,1}(1,:)) == 5
    peakArray{42,1}(:,4:5) = [];
end

save('peak times/CMb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/APb3.mat', 'peakArray')

if length(peakArray{16,1}(1,:)) == 5
    peakArray{16,1}(:,4:5) = [];
end

if length(peakArray{18,1}(1,:)) == 4
    peakArray{18,1}(:,4) = [];
end

if length(peakArray{41,1}(1,:)) == 4
    peakArray{41,1}(:,4) = [];
end

save('peak times/APb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/APb4.mat', 'peakArray')

if length(peakArray{28,1}(1,:)) == 4
    peakArray{28,1}(:,4) = [];
end

if length(peakArray{38,1}(1,:)) == 4
    peakArray{38,1}(:,4) = [];
end

if length(peakArray{40,1}(1,:)) == 4
    peakArray{40,1}(:,4) = [];
end

if length(peakArray{41,1}(1,:)) == 4
    peakArray{41,1}(:,4) = [];
end

save('peak times/APb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/DDb3.mat', 'peakArray')

if length(peakArray{3,1}(1,:)) == 5
    peakArray{3,1}(:,4:5) = [];
end

if length(peakArray{5,1}(1,:)) == 4
    peakArray{5,1}(:,4) = [];
end

if length(peakArray{6,1}(1,:)) == 5
    peakArray{6,1}(:,4:5) = [];
end

if length(peakArray{7,1}(1,:)) == 5
    peakArray{7,1}(:,4:5) = [];
end

if length(peakArray{8,1}(1,:)) == 5
    peakArray{8,1}(:,4:5) = [];
end

if length(peakArray{9,1}(1,:)) == 5
    peakArray{9,1}(:,4:5) = [];
end

if length(peakArray{11,1}(1,:)) == 4
    peakArray{11,1}(:,4) = [];
end

if length(peakArray{12,1}(1,:)) == 5
    peakArray{12,1}(:,4:5) = [];
end

if length(peakArray{18,1}(1,:)) == 5
    peakArray{18,1}(:,4:5) = [];
end

if length(peakArray{23,1}(1,:)) == 4
    peakArray{23,1}(:,4) = [];
end

if length(peakArray{25,1}(1,:)) == 5
    peakArray{25,1}(:,4:5) = [];
end

if length(peakArray{27,1}(1,:)) == 4
    peakArray{27,1}(:,4) = [];
end

if length(peakArray{29,1}(1,:)) == 4
    peakArray{29,1}(:,4) = [];
end

if length(peakArray{34,1}(1,:)) == 4
    peakArray{34,1}(:,4) = [];
end

if length(peakArray{42,1}(1,:)) == 4
    peakArray{42,1}(:,4) = [];
end

save('peak times/DDb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/DDb4.mat', 'peakArray')

if length(peakArray{6,1}(1,:)) == 4
    peakArray{6,1}(:,4) = [];
end

if length(peakArray{37,1}(1,:)) == 4
    peakArray{37,1}(:,4) = [];
end

if length(peakArray{39,1}(1,:)) == 4
    peakArray{39,1}(:,4) = [];
end

save('peak times/DDb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/SCb3.mat', 'peakArray')

if length(peakArray{2,1}(1,:)) == 5
    peakArray{2,1}(:,4:5) = [];
end

if length(peakArray{9,1}(1,:)) == 4
    peakArray{9,1}(:,4) = [];
end

if length(peakArray{14,1}(1,:)) == 5
    peakArray{14,1}(:,4:5) = [];
end

if length(peakArray{26,1}(1,:)) == 4
    peakArray{26,1}(:,4) = [];
end

if length(peakArray{28,1}(1,:)) == 4
    peakArray{28,1}(:,4) = [];
end

if length(peakArray{32,1}(1,:)) == 5
    peakArray{32,1}(:,4:5) = [];
end

if length(peakArray{39,1}(1,:)) == 4
    peakArray{39,1}(:,4) = [];
end

if length(peakArray{44,1}(1,:)) == 5
    peakArray{44,1}(:,4:5) = [];
end

if length(peakArray{48,1}(1,:)) == 4
    peakArray{48,1}(:,4) = [];
end

save('peak times/SCb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/SCb4.mat', 'peakArray')

if length(peakArray{23,1}(1,:)) == 5
    peakArray{23,1}(:,4:5) = [];
end

if length(peakArray{31,1}(1,:)) == 4
    peakArray{31,1}(:,4) = [];
end

if length(peakArray{47,1}(1,:)) == 4
    peakArray{47,1}(:,4) = [];
end

if length(peakArray{48,1}(1,:)) == 4
    peakArray{48,1}(:,4) = [];
end

save('peak times/SCb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/AKb3.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 5
    peakArray{1,1}(:,4:5) = [];
end

if length(peakArray{4,1}(1,:)) == 5
    peakArray{4,1}(:,4:5) = [];
end

if length(peakArray{16,1}(1,:)) == 7
    peakArray{16,1}(:,3:4) = [];
    peakArray{16,1}(:,4:5) = [];
end

if length(peakArray{19,1}(1,:)) == 5
    peakArray{19,1}(:,4:5) = [];
end

if length(peakArray{22,1}(1,:)) == 4
    peakArray{22,1}(:,4) = [];
end

if length(peakArray{29,1}(1,:)) == 4
    peakArray{29,1}(:,4) = [];
end

if length(peakArray{44,1}(1,:)) == 4
    peakArray{44,1}(:,4) = [];
end

if length(peakArray{50,1}(1,:)) == 5
    peakArray{50,1}(:,4:5) = [];
end

save('peak times/AKb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/AKb4.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 5
    peakArray{1,1}(:,4:5) = [];
end

if length(peakArray{19,1}(1,:)) == 4
    peakArray{19,1}(:,4) = [];
end

if length(peakArray{22,1}(1,:)) == 5
    peakArray{22,1}(:,4:5) = [];
end

if length(peakArray{26,1}(1,:)) == 4
    peakArray{26,1}(:,4) = [];
end

if length(peakArray{45,1}(1,:)) == 5
    peakArray{45,1}(:,4:5) = [];
end

save('peak times/AKb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/GBb3.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 5
    peakArray{1,1}(:,4:5) = [];
end

if length(peakArray{27,1}(1,:)) == 4
    peakArray{27,1}(:,4) = [];
end

if length(peakArray{43,1}(1,:)) == 5
    peakArray{43,1}(:,4:5) = [];
end

if length(peakArray{45,1}(1,:)) == 5
    peakArray{45,1}(:,4:5) = [];
end

if length(peakArray{48,1}(1,:)) == 7
    peakArray{48,1}(:,4:7) = [];
end

save('peak times/GBb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/GBb4.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 4
    peakArray{1,1}(:,4) = [];
end

if length(peakArray{4,1}(1,:)) == 4
    peakArray{4,1}(:,4) = [];
end

if length(peakArray{8,1}(1,:)) == 5
    peakArray{8,1}(:,4:5) = [];
end

if length(peakArray{14,1}(1,:)) == 4
    peakArray{14,1}(:,4) = [];
end

if length(peakArray{39,1}(1,:)) == 5
    peakArray{39,1}(:,4:5) = [];
end

if length(peakArray{40,1}(1,:)) == 4
    peakArray{40,1}(:,4) = [];
end

if length(peakArray{50,1}(1,:)) == 7
    peakArray{50,1}(:,2:3) = [];
    peakArray{50,1}(:,4:5) = [];
end

save('peak times/GBb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/HCb3.mat', 'peakArray')

if length(peakArray{16,1}(1,:)) == 5
    peakArray{16,1}(:,4:5) = [];
end

if length(peakArray{18,1}(1,:)) == 4
    peakArray{18,1}(:,4) = [];
end

if length(peakArray{24,1}(1,:)) == 5
    peakArray{24,1}(:,4:5) = [];
end

if length(peakArray{48,1}(1,:)) == 4
    peakArray{48,1}(:,4) = [];
end

save('peak times/HCb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/HCb4.mat', 'peakArray')

if length(peakArray{2,1}(1,:)) > 3
    peakArray{2,1}(:,3:end) = [];
end

if length(peakArray{8,1}(1,:)) == 4
    peakArray{8,1}(:,4) = [];
end

if length(peakArray{14,1}(1,:)) == 4
    peakArray{14,1}(:,4) = [];
end

if length(peakArray{39,1}(1,:)) == 4
    peakArray{39,1}(:,4) = [];
end

if length(peakArray{43,1}(1,:)) == 4
    peakArray{43,1}(:,4) = [];
end

if length(peakArray{50,1}(1,:)) == 5
    peakArray{50,1}(:,4:5) = [];
end

save('peak times/HCb4.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/JMb3.mat', 'peakArray')

if length(peakArray{1,1}(1,:)) == 4
    peakArray{1,1}(:,4) = [];
end

if length(peakArray{5,1}(1,:)) == 7
    peakArray{5,1}(:,2:3) = [];
    peakArray{5,1}(:,4:5) = [];
end

if length(peakArray{14,1}(1,:)) == 4
    peakArray{14,1}(:,4) = [];
end

if length(peakArray{20,1}(1,:)) == 4
    peakArray{20,1}(:,4) = [];
end

if length(peakArray{31,1}(1,:)) == 4
    peakArray{31,1}(:,4) = [];
end

if length(peakArray{43,1}(1,:)) == 4
    peakArray{43,1}(:,4) = [];
end

if length(peakArray{48,1}(1,:)) == 4
    peakArray{48,1}(:,4) = [];
end
if length(peakArray{49,1}(1,:)) == 4
    peakArray{49,1}(:,4) = [];
end

save('peak times/JMb3.mat', 'peakArray')
clear peakArray

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('peak times/JMb4.mat', 'peakArray')

if length(peakArray{11,1}(1,:)) == 5
    peakArray{11,1}(:,4:5) = [];
end

if length(peakArray{45,1}(1,:)) == 5
    peakArray{45,1}(:,4:5) = [];
end

if length(peakArray{46,1}(1,:)) == 5
    peakArray{46,1}(:,4:5) = [];
end

if length(peakArray{50,1}(1,:)) == 4
    peakArray{50,1}(:,4) = [];
end

save('peak times/JMb4.mat', 'peakArray')
clear peakArray

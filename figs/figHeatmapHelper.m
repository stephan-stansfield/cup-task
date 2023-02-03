
meanK = readmatrix('best fit simulation/best fit impedance.xlsx','Sheet','meanK','Range','B2:L15');
stdK = readmatrix('best fit simulation/best fit impedance.xlsx','Sheet','stdK','Range','B2:L15');
meanB = readmatrix('best fit simulation/best fit impedance.xlsx','Sheet','meanB','Range','B2:L15');
stdB = readmatrix('best fit simulation/best fit impedance.xlsx','Sheet','stdB','Range','B2:L15');

meanK = [nan(1,11); meanK];
stdK = [nan(1,11); stdK];
meanB = [nan(1,11); meanB];
stdB = [nan(1,11); stdB];

save('best fit simulation/best fit impedance.mat','meanK','stdK','meanB','stdB')

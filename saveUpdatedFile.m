% Save updated workspace file. This was specifically made while manually
% fixing misidentified velocity peaks and valleys in the experimental data;
% user beware!

saveFolder = strcat("best fit simulation/Experimental Trials/Peak and Duration Data/");
saveFile = strcat(subjNum, "_B", block3or4, ".mat");
fileName = strcat(saveFolder,saveFile);
save(fileName, 'velAmpArray', 'velTimeArray', 'durationArray');
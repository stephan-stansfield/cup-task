function block3or4 = blockNumFunc(blockNum)
% BLOCKNUMFUNC
% Returns a string indicating whether current block is 3 or 4

if mod(blockNum,2)
   block3or4 = '3';
else
   block3or4 = '4';
end
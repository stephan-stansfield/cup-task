function y = mydigits(z)
% MYDIGITS
% Returns order of magnitude of decimal number less than one
    z = abs(z);
    y = 0;
    while (floor(z)~=z)
        y = y+1;
        z = z*10;
    end
end
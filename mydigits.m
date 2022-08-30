% Returns order of magnitude of decimal number less than one
function y = mydigits(z)
    z = abs(z);
    y = 0;
    while (floor(z)~=z)
        y = y+1;
        z = z*10;
    end
end
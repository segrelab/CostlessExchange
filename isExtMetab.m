% Detect and split metabolite from compartment
%
% Alan R. Pacheco 10/24/16

function answer = isExtMetab(met,delimiter)

answer = 0;
spl = strsplit(met,delimiter);
c = spl{2};
compartment = lower(c(1));
if strcmp(compartment,'e')
    answer = 1;
end
function plotLinModelFitsWrap(timebase, avgTS, stdTS, MSE, modelTS, idCell, startTimesSorted_A, stimValuesSorted_A, startTimesSorted_B, stimValuesSorted_B)

for i = 1:size(avgTS,1)
   if strfind(char(idCell(i)),' A')
        display('A');
   elseif strfind(char(idCell(i)),' B')
        display('B');
   else
      error('plotLinModelFitsWrap: unknown stimulus type');
   end
end

end
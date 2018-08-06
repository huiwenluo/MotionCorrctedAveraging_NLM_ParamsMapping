function ff = computeFF( outParams )

        
curf = outParams.species(2).amps;
curw = outParams.species(1).amps;


ff1 = 100*abs(curf)./(abs(curf) + abs(curw));

fatregions = ff1>50;
watregions = ff1<=50;
all_ff = 100*abs(curf)./(abs(curf) + abs(curw));
all_ff(watregions) = 100 - 100*abs(curw(watregions))./abs(curf(watregions) + curw(watregions));
all_ff(fatregions) = 100*abs(curf(fatregions))./abs(curf(fatregions) + curw(fatregions));


ff = all_ff;
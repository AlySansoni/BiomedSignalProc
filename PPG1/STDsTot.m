%% Function that prints the plots in the right format
% plot_lb is the boolean parameters to the optional print of some plots

function []= STDsTot(STDhr,STDamp,STDsd,STDwav,label,plot_lb)
   % STD multiplied for 7s (a segment)
   STDhrTot = zeros(1,210);
   STDampTot = zeros(1,210);
   STDsdTot = zeros(1,210);
   STDwavTot = zeros(1,210);
   LabelTot = zeros(1,210);
   
    n=0;
   %Assigning the computed value to the whole segment duration
   for j = 1:29
    STDhrTot(1,(1+n*7):(7*(n+1))) = STDhr(j);
    STDampTot(1,(1+n*7):(7*(n+1))) = STDamp(j);
    STDsdTot(1,(1+n*7):(7*(n+1))) = STDsd(j);
    STDwavTot(1,(1+n*7):(7*(n+1))) = STDwav(j);
    LabelTot(1,(1+n*7):(7*(n+1))) = label(j);
    n=n+1;
   end
   
  if plot_lb 
   figure();
   subplot(5,1,1);
   plot(STDhrTot);
   title('STDhr');
   xlabel('Seconds');
   ylabel('STDhr [bpm]');

   subplot(5,1,2);
   plot(STDampTot);
   title('STDamp');
   xlabel('Seconds');
   ylabel('STDamp [au]');
   subplot(5,1,3);
   plot(STDsdTot);
   title('STDsd');
   xlabel('Seconds');
   ylabel('STDsd [s/s]');

   subplot(5,1,4);
   plot(STDwavTot);
   title('STDwav');
   xlabel('Seconds');
   ylabel('STDwav [au]');
   
   subplot(5,1,5);
   plot(LabelTot);
   title('Label');
   xlabel('Seconds');
   ylabel('Decision');
  end
end
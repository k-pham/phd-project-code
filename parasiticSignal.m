for fignum = 1:30
figure(fignum)
hold on
title(['slice_y = ',num2str(fignum)])
filename = ['D:\PROJECT\figures\coatingCharac2D\180424 parasitic signal BA59\parasiticSignal_Yslice',num2str(fignum),'.png'];
saveas(gcf,filename)
end
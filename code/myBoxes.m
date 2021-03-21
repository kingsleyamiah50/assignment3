function [ ] = myBoxes(BoxBnd1x, BoxBnd1y, BoxBnd2x, BoxBnd2y, BoxBnd3x, BoxBnd3y, BoxBnd4x, BoxBnd4y)

plot([BoxBnd1x BoxBnd2x], [BoxBnd1y BoxBnd2y],'k');
hold on
plot([BoxBnd2x BoxBnd3x], [BoxBnd2y BoxBnd3y],'k');
hold on
plot([BoxBnd3x BoxBnd4x], [BoxBnd3y BoxBnd4y],'k');
hold on
plot([BoxBnd4x BoxBnd1x], [BoxBnd4y BoxBnd1y],'k');
hold on


end




global NB Radius wall_L wall_R wall_B wall_T
figure;
hold on

cmArray = zeros(2*NB, 1);
for i = 1:NB
    cmArray(1+2*(i-1)) = RenderInfo(1, 1+2*(i-1));
    cmArray(2+2*(i-1)) = RenderInfo(1, 2+2*(i-1));
end

radArray = zeros(NB,1);
for i = 1:NB
    radArray(i) = Radius(i);
end

th = 0:pi/50:2*pi;
drawSphere = zeros(length(th), 2*NB);
for i = 1:NB
    drawSphere(:,1+2*(i-1)) = radArray(i)*cos(th) + cmArray(1+2*(i-1));
    drawSphere(:,2+2*(i-1)) = radArray(i)*sin(th) + cmArray(2+2*(i-1));
end


handleArray = cell(NB,1);
for i = 1:NB
    handleArray{i} = plot(drawSphere(:, 1+2*(i-1)), drawSphere(:, 2+2*(i-1)));
end

axis([wall_L wall_R wall_B wall_T]);
axis square

for k = 2:length(RenderInfo);
    
    for i = 1:NB
        handleArray{i}.XData = handleArray{i}.XData + RenderInfo(k,2*i-1) - RenderInfo(k-1,2*i-1);
        handleArray{i}.YData = handleArray{i}.YData + RenderInfo(k,2*i) - RenderInfo(k-1,2*i);
    end
    
    pause(0.01)
    
end
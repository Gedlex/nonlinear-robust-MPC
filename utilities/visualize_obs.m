function visualize_obs(obs)
    r = obs(:,3);
    pos = [obs(:,[1,2])-r,2*r,2*r];
    for i = 1:size(obs,1)
        rectangle('Position',pos(i,:),'Curvature',[1 1],'FaceColor','k')
    end
    axis equal
end
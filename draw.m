function draw(X)
    h = 25;
    w = 18;
    r = 1;

    for i = 1:size(X,2)
        x = X(1,i);
        y = X(2,i);
        z = X(3,i);
        
        column_tip1 = [0,-w/2,h];
        column_tip2 = [0,w/2,h];

        plot([column_tip1(2),y - r],[column_tip1(3),z]); hold on
        plot([y + r,column_tip2(2)],[z,column_tip2(3)]); hold off
        axis([-10 10 -5 40]);
        daspect([1 1 1]);
        drawnow
        %pause(0.1);
    end
end
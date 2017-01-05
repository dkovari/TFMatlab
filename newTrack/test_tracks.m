
clear P
[xx,yy] = meshgrid(1:30,1:30);

P{1} = [reshape(xx,[],1), reshape(yy,[],1)];
for f=2:100
    P{f} = P{f-1}+ 0.4*(rand(size(P{1}))-.5);
end

%tracks = find_tracks(P,'MinLength',98,'MaxDisp',1,'ExpandingSearch',false);
tracks = find_displacement(P,1,'MaxDisp',1);

figure(1);
clf;
gca;
hold on;

% for f=1:numel(P)
%     plot(P{f}(:,1),P{f}(:,2),'.','markersize',20);
% end

%%
hold on
for t=1:numel(tracks)
    plot(tracks{t}(:,1),tracks{t}(:,2),'-');
end
plot(xx(:),yy(:),'+k');

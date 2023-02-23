function show_match_sift(im1,im2,matches,f1,f2,inlierIdx,tform)
% matches: [index_in_f1; index_in_f2]
% inlierIdx are the indices in matches that are good matches
Nm = size(matches,2);
figure
subplot(1,3,1), imagesc(im1), axis image, colormap('gray'), hold on
h11 = vl_plotframe(f1(:,matches(1,:)));
h21 = vl_plotframe(f1(:,matches(1,:)));
set(h11,'color','k','linewidth',3) ;
set(h21,'color','y','linewidth',2) ;
h12 = vl_plotframe(f2(:,matches(2,:)));
h22 = vl_plotframe(f2(:,matches(2,:)));
set(h12,'color','m','linewidth',3) ;
set(h22,'color','r','linewidth',2) ;
X = f1(1,matches(1,:)); Y = f1(2,matches(1,:));
X2 = f2(1,matches(2,:)); Y2 = f2(2,matches(2,:));
U = X2 - X; V = Y2 - Y;
bad_matches = ones(Nm,1,'logical');
bad_matches(inlierIdx) = false;
quiver(X(inlierIdx),Y(inlierIdx),U(inlierIdx),V(inlierIdx),0,'b')
quiver(X(bad_matches),Y(bad_matches),U(bad_matches),V(bad_matches),0,'r')
[expected_u, expected_v] = transformPointsInverse(tform,f1(1,matches(1,:)),f1(2,matches(1,:)));
plot(expected_u, expected_v,'go')

subplot(1,3,2), imagesc(im1), axis image, colormap('gray'), hold on
h11 = vl_plotframe(f1(:,matches(1,:)));
h21 = vl_plotframe(f1(:,matches(1,:)));
set(h11,'color','k','linewidth',3) ;
set(h21,'color','y','linewidth',2) ;
subplot(1,3,3), imagesc(im2), axis image, colormap('gray'), hold on
h12 = vl_plotframe(f2(:,matches(2,:)));
h22 = vl_plotframe(f2(:,matches(2,:)));
set(h12,'color','k','linewidth',3) ;
set(h22,'color','y','linewidth',2) ;
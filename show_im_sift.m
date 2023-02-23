function show_im_sift(im, f, pts, d)
Nf = size(f,2);
if ~exist('pts', 'var') || isempty(pts)
   pts = Nf;
end
figure
imagesc(im), axis image, colormap('gray'), hold on
perm = randperm(Nf) ;
sel = perm(1:pts) ;
h1 = vl_plotframe(f(:,sel)) ;
h2 = vl_plotframe(f(:,sel)) ;
set(h1,'color','k','linewidth',3) ;
set(h2,'color','y','linewidth',2) ;
if exist('d','var') || ~isempty(pts)
   h3 = vl_plotsiftdescriptor(d(:,sel),f(:,sel)) ;
   set(h3,'color','g') ;
end
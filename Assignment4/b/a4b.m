% 2D/2D image registration

global I1 I2 original_I2 trans2 rot2 I1_title I2_title

I1_filename = 'ct_1.png';
I2_filename = 'mr_1.png';

I1_title = 'CT';
I2_title = 'MR'; 

interactive = true;

% Get the files

I1 = rgb2gray( imread( I1_filename ) );
I2 = rgb2gray( imread( I2_filename ) );

original_I2 = I2; % keep a copy, as I2 will be modified


% ---------------- INTERACTIVE VERSION ----------------

if interactive

  % Interactive transformations to apply to I2

  trans2 = [0 0];
  rot2 = 0;

  % transform I2

  I2 = imtranslate( imrotate( original_I2, rot2, 'bilinear', 'crop' ), trans2 );

  % Compute similarity measures

  [MI, joint_histo] = compute_MI( I1, I2 );

  RMS = compute_RMS( I1, I2 );

  NCC = compute_NCC( I1, I2 ); 

  disp( sprintf( 'trans %3d,%3d, rot %3d: MI %.4f, RMS %6.0f, NCC %6.0f', trans2(1), trans2(2), rot2, MI, RMS, NCC ) );

  % Draw everything and wait for key presses

  draw_all_MI( I1_title, I2_title, joint_histo, MI );
  return;

end

% ---------------- GRAPHING VERSION ----------------


% Build 2D plots of similarity measure with varying offsets of image I2 from (0,0).
%
% These cases are done:
%
%   rotation =  0, CT image vs MR image, three similarity measures
%   rotation = 45, CT image vs MR image, three similarity measures
%   rotation = 90, CT image vs MR image, three similarity measures
%
%   rotation =  0, CT image vs CT image, three similarity measures
%   rotation = 45, CT image vs CT image, three similarity measures
%   rotation = 90, CT image vs CT image, three similarity measures
%
% All results are displayed in a two image windows.
%
% In your initial testing, set 'step' to 40.  For your final images, set
% it to 4 and save the images at the biggest size you can.

step = 40;

% CT/MR measures

figure(1);
tiledlayout(3,3);
axis equal;
for rotation = [0 45 90]
  [MI, RMS, NCC] = collect_measures( I1, I2, rotation, step );  % note the 'I1, I2' arguments
  draw_measure( MI,  'MI',  rotation );
  draw_measure( RMS, 'RMS', rotation );
  draw_measure( NCC, 'NCC', rotation );
end
sgtitle( sprintf( '%s / %s measures', I1_title, I2_title ) );


% CT/T1 measures

figure(2);
tiledlayout(3,3);
axis equal;
for rotation = [0 45 90]
  [MI, RMS, NCC] = collect_measures( I1, I1, rotation, step );  % note the 'I1, I1' arguments
  draw_measure( MI,  'MI',  rotation );
  draw_measure( RMS, 'RMS', rotation );
  draw_measure( NCC, 'NCC', rotation );
end
sgtitle( sprintf( '%s / %s measures', I1_title, I1_title ) );

return;

% ---------------- END OF MAIN SCRIPT ----------------


% For two images, I1 and I2, collect the three similarity measures of
% Mutual Information, RMS Error, and Normalized Cross Correlation with
% various offsets of I2 with respect to I1.
%
% If the 'step' parameter is 10, for example, the (row,column) offsets (i,j) would be
%
%   for i = -nrows/4 ... -40 -30 -20 -10 0 10 20 30 40 ... +3/4*nrows
%       for j = -ncols/4 ... -40 -30 -20 -10 0 10 20 30 40 ... +3/4*ncols
%
% For each (i,j) pair, first rotate I2 FROM ITS ORIGINAL POSITION by
% the 'rotation' parameter, then offset I2 by (i,j) and compute values
% for the three similarity measures (MI, RMS, NCC).  Store those
% values in matrices.  Note that compute_MI, compute_RMS, and compute_NCC
% are already provided.
%
% For the image transformations, use 'imrotate' (with 'bilinear' and
% 'crop') and use 'imtranslate'.
%
% I1 and I2 must have the same dimensions.


function [MI, RMS, NCC] = collect_measures( I1, I2, rotation, step )

  nrows = size(I1,1);
  ncols = size(I1,2);

  rowRange = round(1.5*nrows/step);  % number of row offsets to iterate over (some of these are positive offsets and some are negative offsets).
  colRange = round(1.5*ncols/step);  % number of column offsets ...

  MI  = zeros( rowRange, colRange );  % store similarity measures here
  RMS = zeros( rowRange, colRange );
  NCC = zeros( rowRange, colRange );

  % YOUR CODE HERE

end



% Draw three surfaces in the next three tiles

function draw_measure( graph, name, rotation )

  nexttile;
    
  if size(graph,1) < 100
    surf( graph, 'edgecolor', [0.5 0.5 0.5] )
  else
    surf( graph, 'edgecolor','none' );
  end

  title( sprintf( '%s - %d degrees', name, rotation ) );
end


% Compute the mutual information of two images

function [MI, joint_histo] = compute_MI( I1, I2 )
 
  % Compute joint histogram

  x_edges = [-0.5:255.5];
  y_edges = [-0.5:255.5];

  joint_histo = histcounts2( I1, I2, x_edges, y_edges );

  % Compute the joint probability, P(x,y)

  N = sum( joint_histo, 'all' );

  Pxy = joint_histo/N;

  % Compute the marginal probabilities, P(x) and P(y)

  Py = sum( Pxy, 1 ); % column sums = marginal probilities for one image
  Px = sum( Pxy, 2 ); % row sums = marginal probilities for other image

  % Compute MI

  MI = 0;
  for i = 1:length(Px)
      for j = 1:length(Py)
          denom = Px(i) * Py(j);
          if denom ~= 0 && Pxy(i,j) ~= 0
              MI = MI + Pxy(i,j) * log( Pxy(i,j) / denom );
          end
      end
  end
end


% Compute the RMS of two images

function RMS = compute_RMS( I1, I2 )
    
  RMS = - sqrt( sum( (double(I1) - double(I2)) .^ 2, 'all' ) ); 
end


% Compute the Normalize Cross Correlation of two images

function NCC = compute_NCC( I1, I2 )

  I1d = double(I1);
  I2d = double(I2);

  mu1 = mean( I1d, 'all' );
  mu2 = mean( I2d, 'all' );

  stdev1 = std( I1d, 0, 'all' );
  stdev2 = std( I2d, 0, 'all' );

  NCC = sum( (I1d - mu1) .* (I2d - mu2), 'all' ) / (stdev1 * stdev2);
end


% For interactive use: Draw two images and compute similarity measures

function fig = draw_all_MI( title1, title2, joint_histo, MI )

  global I1 I2;

  % Draw two images

  fig = tiledlayout(2,3).Parent;
  set( fig, 'keypressfcn', @keyPressHandler );

  nexttile(1);
  imshow(I1);
  title( title1, 'FontSize', 11, 'FontWeight', 'normal' );

  nexttile(4);
  imshow(I2);
  title( title2, 'FontSize', 11, 'FontWeight', 'normal' );

  % ---- Draw joint histogram ----

  % Zero the dark areas, which have high counts and correspond to background pixels
  
  joint_histo(1:10,1:10) = 0;

  % log-scale the histogram to show more detail

  log_joint_histo = log( joint_histo + 1 );

  % Draw joint histogram

  x_edges = [-0.5:255.5];
  y_edges = [-0.5:255.5];

  nexttile(2,[2 2]);
  histogram2( 'xbinedges', x_edges, 'ybinedges', y_edges, 'bincounts', log_joint_histo, 'facecolor', 'flat', 'displaystyle', 'tile' );
  
  title( sprintf( 'MI = %.4f', MI ), 'FontSize', 16, 'FontWeight', 'normal' );
end



function fig = draw_all_RMS_NCC( title1, title2, joint_histo, RMS, NCC )

  global I1 I2;

  % Draw two images

  t = tiledlayout(1,2);
  fig = t.Parent;
  txt = title( t, sprintf( 'RMS = %.0f, NCC = %.0f', RMS, NCC ), 'Color', 'w' );

  set( fig, 'keypressfcn', @keyPressHandler );
  set(gcf,'color','black');

  nexttile(1);
  imshow(I1);
  title( title1, 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'w' );

  nexttile(2);
  imshow(I2);
  title( title2, 'FontSize', 11, 'FontWeight', 'normal', 'Color', 'w' );
end



% For interactive use: Key press handler

function keyPressHandler( src, event )

  global I1 I2 original_I2 trans2 rot2 I1_title I2_title;

  transInc = 2; % in pixels
  rotInc = 1;   % in degrees

  % Translate or rotate with arrows and < >

  if strcmp( event.Key, 'rightarrow' )
    trans2 = trans2 + [transInc 0];
  elseif strcmp( event.Key, 'leftarrow' )
    trans2 = trans2 + [-transInc 0];
  elseif strcmp( event.Key, 'uparrow' )
    trans2 = trans2 - [0 transInc];       % use '-' because image origin is TOP-left corner
  elseif strcmp( event.Key, 'downarrow' )
    trans2 = trans2 - [0 -transInc];
  elseif strcmp( event.Key, 'comma' )
    rot2 = rot2 + rotInc;
  elseif strcmp( event.Key, 'period' )
    rot2 = rot2 - rotInc;
  end

  % Second image title includes description of translation and rotation

  if trans2 == [0 0] & rot2 == 0
    title2 = sprintf( '%s\n', I2_title );
  else
    title2 = sprintf( '%s\ntrans %d,%d pixels, rot %d degrees', I2_title, trans2(1), -trans2(2), rot2 );
  end
  
  % Apply transformation to I2

  I2 = imtranslate( imrotate( original_I2, rot2, 'bilinear', 'crop' ), trans2 );

  % Recompute joint histogram and similarity measures

  [MI, joint_histo] = compute_MI( I1, I2 );
  RMS = compute_RMS( I1, I2 );
  NCC = compute_NCC( I1, I2 ); 

  disp( sprintf( 'trans %3d,%3d, rot %3d: MI %.4f, RMS %6.0f, NCC %6.0f', trans2(1), trans2(2), rot2, MI, RMS, NCC ) );

  % Draw everything

  fig = draw_all_MI( I1_title, title2, joint_histo, MI );
  %fig = draw_all_RMS_NCC( I1_title, title2, joint_histo, RMS, NCC );

  % Save the figure if 's' is pressed

  if strcmp( event.Key, 's' )
    disp( 'Figure saved in fig.png' );
    saveas( fig, 'fig.png' );
  end
end

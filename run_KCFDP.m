%   This is the visual object tracker KCFDP presented in:
%   "Enable Scale and Aspect Ratio Adaptability in Visual Tracking with Detection Proposals" BMVC, 2015
%   
%   Dafei Huang, Lei Luo, Mei Wen, Zhaoyun Chen and Chunyuan Zhang
%
%   Utilization of EdgeBoxes to enable scale and aspect ratio adaptibility

%  "Exploiting the Circulant Structure of Tracking-by-detection with Kernels" ECCV, 2012
%
%   J. F. Henriques, R. Caseiro, P. Martins and J. Batista
%
%   Original CSK tracker implementation

%  "High-Speed Tracking with Kernelized Correlation Filters"  TPAMI, 2014
%   http://www.isr.uc.pt/~henriques/circulant/
%
%   J. F. Henriques, R. Caseiro, P. Martins, J. Batista
%
%   Original KCF and DCF tracker implementation

%  "Online Object Tracking: A Benchmark", CVPR, 2013.
%   http://visual-tracking.net/
% 
%   Y. Wu, J. Lim, M.-H. Yang
%
%   Benchmark sequences and toolkit

%  "Piotr's Image and Video Matlab Toolbox (PMT)"
%   http://vision.ucsd.edu/~pdollar/toolbox/doc/index.html
% 
%   P. Dollar
%
%   Tools untilized

%  "Adaptive Color Attributes for Real-Time Visual Tracking" CVPR, 2014.
%
%   Martin Danelljan, Fahad Shahbaz Khan, Michael Felsberg and Joost van de Weijer
%
%   Extended CSK tracker implementation with color attributes, new updating scheme and adaptive dimensionality reduction 

%  "Structured Forests for Fast Edge Detection", ICCV 2013
% 
%   P. Dollar and C. Zitnick
%
%   Very fast edge detector (up to 60 fps depending on parameter settings) that achieves excellent accuracy

%  "Edge Boxes: Locating Object Proposals from Edges", ECCV 2014.
%
%   C. Zitnick and P. Dollar
%
%   Edge Boxes object proposal generation

%   Codes above are integrated and modified by Dafei Huang

function results=run_KCFDP(seq, res_path, bSaveImage)
% 'res=run_' t.name '(subS, rp, bSaveImage);'  sequence info, result path: ./tmp/evalType/sequence name_SCKCF_1.mat, whether to save result image

addpath(genpath('/Users/dafei/Projects/toolbox-master/'));
addpath(('../../rstEval'));

close all

s_frames = seq.s_frames;

% parameters according to the paper
padding = 1.5;  %extra area surrounding the target
lambda = 1e-4;  %regularization
% output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
output_sigma_factor = 0.06;  %spatial bandwidth (proportional to target)
% interp_factor = 0.02;
interp_factor = 0.01;
sigma = 0.5;
hog_orientations = 9;
cell_size = 4;


target_sz = [seq.init_rect(1,4), seq.init_rect(1,3)];
pos = floor([seq.init_rect(1,2), seq.init_rect(1,1)]) + floor(target_sz/2);

% general parameters
params.visualization = 1;
params.init_pos = pos;
params.wsize = floor(target_sz);
params.video_path = seq.path;
params.s_frames = s_frames;
params.bSaveImage = bSaveImage;
params.res_path = res_path;


% load pre-trained edge detection model and set opts
model = load('./models/forest/modelBsds'); 
model = model.model;
model.opts.multiscale = 0;
model.opts.sharpen = 0;
model.opts.nThreads = 4;

% set up parameters for edgeBoxes (see edgeBoxes.m)
opts = edgeBoxes;
opts.alpha = .65;      % step size of sliding window search
opts.beta = .75;       % nms threshold for object proposals
%opts.maxBoxes = 1e4;  % max number of boxes to detect
opts.maxBoxes = 1e3;   % don't need that many proposals (default is 1e4)
% opts.minBoxArea = 200;
% opts.maxAspectRatio = 3.5;
% opts.minScore = .01;  % min score of boxes to detect
opts.minScore = 0.0005;
opts.kappa = 1.4;     % 1.5 as default, can be changed for larger overlapping
opts.minBoxArea = 200;
opts.edgeMinMag = 0.1;


% set up parameters for using edgeBoxes in scale detection
scale_params.scale_detect_window_factor = 1.4;
scale_params.proposal_num_limit = 200;
scale_params.pos_shift_damping = 0.7;
scale_params.rescale_damping = 0.7;


[rect_position, fps] = tracker(params, ...
			padding, sigma, lambda, output_sigma_factor, interp_factor, ...
			cell_size, hog_orientations, ...
            model, opts, scale_params);


disp(['fps: ' num2str(fps)])

results.type = 'rect';
results.res = rect_position;%each row is a rectangle
results.fps = fps;

end

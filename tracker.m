%   Original code is from Kernelized/Dual Correlation Filter (KCF/DCF)
%   by Joao F. Henriques, 2014

%   Integrated and modified by Dafei Huang

function [rect_position, fps] = tracker(params, ...
    padding, sigma, lambda, output_sigma_factor, interp_factor, ...
    cell_size, hog_orientations, ...
    model, opts, scale_params)

% general parameters
bSaveImage = params.bSaveImage;
res_path = params.res_path;
video_path = params.video_path;
s_frames = params.s_frames;
pos = floor(params.init_pos);
target_sz = floor(params.wsize);
visualization = params.visualization;

% scale detection parameters
scale_detect_window_factor = scale_params.scale_detect_window_factor;
proposal_num_limit = scale_params.proposal_num_limit;
pos_shift_damping = scale_params.pos_shift_damping;
rescale_damping = scale_params.rescale_damping;


num_frames = numel(s_frames);

%if the target is large, lower the resolution, we don't need that much
%detail
resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
if resize_image,
    pos = floor(pos / 2);
    target_sz = floor(target_sz / 2);
end
org_target_sz = target_sz;


%window size, taking padding into account
window_sz = floor(target_sz * (1 + padding));
org_window_sz = window_sz;

% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);


%create regression labels, gaussian shaped, with a bandwidth
%proportional to target size
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));

%store pre-computed cosine window
cos_window = hann(size(yf,1)) * hann(size(yf,2))';

rect_position = zeros(numel(s_frames), 4);

temp = load('w2crs');
w2c = temp.w2crs;

%note: variables ending with 'f' are in the Fourier domain.

time = 0;  %to calculate FPS

total_proposal_num = 0;
total_flitered_propsal_num = 0;

for frame = 1:numel(s_frames),
    %load image
    im = imread(s_frames{frame});
    org_im = imread(s_frames{frame});
    if size(im,3) > 1,
        im = rgb2gray(im);
    end
    if resize_image,
        im = imresize(im, 0.5);
        org_im = imresize(org_im, 0.5, 'bilinear');
    end
    
    tic()
    
    if frame > 1,
        %obtain a subwindow for detection at the position from last
        %frame, and convert to Fourier domain (its size is unchanged)
        patch = get_subwindow(im, pos, window_sz);
        patch = imresize(patch, org_window_sz, 'bilinear');  % scale its size to original size
        org_patch = get_subwindow(org_im, pos, window_sz);
        org_patch = imresize(org_patch, org_window_sz, 'bilinear');
        zf = fft2(get_features(patch, org_patch, hog_orientations, cell_size, cos_window, w2c));
        
        %calculate response of the classifier at all shifts
        kzf = gaussian_correlation(zf, model_xf, sigma);
        response = real(ifft2(model_alphaf .* kzf));  %equation for fast detection
        
        %target location is at the maximum response. we must take into
        %account the fact that, if the target doesn't move, the peak
        %will appear at the top-left corner, not at the center (this is
        %discussed in the paper). the responses wrap around cyclically.
        max_response = max(response(:));
        [vert_delta, horiz_delta] = find(response == max_response, 1);
        if vert_delta > size(zf,1) / 2,  %wrap around to negative half-space of vertical axis
            vert_delta = vert_delta - size(zf,1);
        end
        if horiz_delta > size(zf,2) / 2,  %same for horizontal axis
            horiz_delta = horiz_delta - size(zf,2);
        end
        
        scale = window_sz ./ org_window_sz;
        pos = pos + ( (cell_size * [vert_delta - 1, horiz_delta - 1]) .* scale );
        
        pre_target_rect = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        
        % begin scale detection
        detect_sz = floor(target_sz * scale_detect_window_factor); % window size for scale detection
        mid_pt = detect_sz * 0.5;   % center position in the window for scale detection
        
        % get the window for scale detection
        edgeBoxes_window = get_subwindow(org_im, pos, detect_sz);
        if size(org_im,3) == 1 % for gray sequences
            edgeBoxes_window = single(edgeBoxes_window / 255);
            edgeBoxes_window = cat(3, edgeBoxes_window, edgeBoxes_window, edgeBoxes_window);
        end
        
        
        % dynamically adjust edgeBoxes parameters
        opts.maxAspectRatio = max([target_sz(1)/target_sz(2), target_sz(2)/target_sz(1)]) * 1.5;
        opts.minBoxArea = floor( prod(target_sz) * 0.3);
        
            
        % edgeBoxes proposals
        edgeBoxes_window_proposals = edgeBoxes(edgeBoxes_window, model, opts);
        cur_total_prop_num = size(edgeBoxes_window_proposals,1);
        
        % choose candidate proposals
        num_of_proposals = 0;
		proposals = zeros(proposal_num_limit,4); % center_y, center_x, rows, cols
        proposals_xywh = zeros(proposal_num_limit,4);
        
        target_in_edgeBoxes_window = [mid_pt([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        
        % find candidate proposals among the top proposal_num_limit proposals detected by edgeBoxes
        for i = 1 : min([size(edgeBoxes_window_proposals,1) proposal_num_limit])
            if calcRectInt(edgeBoxes_window_proposals(i,[1:4]), target_in_edgeBoxes_window) > 0.6 && ...
                    calcRectInt(edgeBoxes_window_proposals(i,[1:4]), target_in_edgeBoxes_window) < 0.9
                proposal_sz = [edgeBoxes_window_proposals(i,4) edgeBoxes_window_proposals(i,3)];
                proposal_pos = [edgeBoxes_window_proposals(i,2) edgeBoxes_window_proposals(i,1)] + floor(proposal_sz/2);
                num_of_proposals = num_of_proposals + 1;
                proposals(num_of_proposals,:) = [pos+proposal_pos-mid_pt, proposal_sz];
                proposals_xywh(num_of_proposals,:) = [proposals(num_of_proposals,[2,1]) - proposal_sz([2,1])/2, proposal_sz([2,1])];
            end
        end
        
%         disp(['proposal num before rejection: ' num2str(cur_total_prop_num) ', after rejection: ' num2str(num_of_proposals)]);
        total_proposal_num = total_proposal_num + cur_total_prop_num;
        total_flitered_propsal_num = total_flitered_propsal_num + num_of_proposals;
        % evaluate all the candidate proposals using kernel correlation
        model_alpha = ifft2(model_alphaf);
        max_proposal_response = max_response;

        new_pos = pos;
        new_target_sz = target_sz;
        for j = 1 : num_of_proposals
            proposal_patch = get_subwindow( im, proposals(j,1:2), proposals(j,3:4)*(1 + padding) );
            proposal_patch = imresize(proposal_patch, org_window_sz, 'bilinear');
            proposal_org_patch = get_subwindow( org_im, proposals(j,1:2), proposals(j,3:4)*(1 + padding) );
            proposal_org_patch = imresize(proposal_org_patch, org_window_sz, 'bilinear');
           
            proposal_zf = fft2(get_features(proposal_patch, proposal_org_patch, hog_orientations, cell_size, cos_window, w2c));
            proposal_kz = gaussian_correlation_nofft(proposal_zf, model_xf, sigma); % no fft needed here
            % calculate the response of the classifier without considering the cyclic shifts
            proposal_response = model_alpha(:)' * proposal_kz(:);
            if proposal_response > max_proposal_response
                max_proposal_response = proposal_response;
                new_pos = proposals(j,1:2);
                new_target_sz = proposals(j,3:4);
            end
        end
        
        chosen_proposal = [new_pos([2,1]) - new_target_sz([2,1])/2, new_target_sz([2,1])];
        
        pos = floor( ( 1 - pos_shift_damping ) * pos + pos_shift_damping * new_pos );
        target_sz = floor( ( 1 - rescale_damping ) * target_sz + rescale_damping * new_target_sz );
        
        window_sz = floor( target_sz * (1 + padding) );
    end
    
    %obtain a subwindow for training at newly estimated target position
    patch = get_subwindow(im, pos, window_sz);
    org_patch = get_subwindow(org_im, pos, window_sz);
    if frame > 1
        patch = imresize(patch, org_window_sz, 'bilinear');  %resize to original window size then train the new model
        org_patch = imresize(org_patch, org_window_sz, 'bilinear');
    end
    xf = fft2(get_features(patch, org_patch, hog_orientations, cell_size, cos_window, w2c));
    
    %Kernel Ridge Regression, calculate alphas (in Fourier domain)
    kf = gaussian_correlation(xf, xf, sigma);
%     alphaf = yf ./ (kf + lambda);   %equation for fast training
    %utilize the updating scheme proposed in ACT
    new_alphaf_num = yf .* kf;
    new_alphaf_den = kf .* (kf + lambda);
    
    
    if frame == 1,  %first frame, train with a single image
%         model_alphaf = alphaf;
        %utilize the updating scheme proposed in ACT
        alphaf_num = new_alphaf_num;
        alphaf_den = new_alphaf_den;
        model_xf = xf;
    else
        %subsequent frames, interpolate model
%         model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
        %utilize the updating scheme proposed in ACT
        alphaf_num = (1 - interp_factor) * alphaf_num + interp_factor * new_alphaf_num;
        alphaf_den = (1 - interp_factor) * alphaf_den + interp_factor * new_alphaf_den;
        model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
    end
    model_alphaf = alphaf_num ./ alphaf_den;
    
    %save position and timing
    rect_position(frame,:) = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    time = time + toc();
    
    %visualization (uncomment the commented code below for detailed visual effects)
    if visualization == 1      
        if frame == 1,  %first frame, create GUI
            figure('Number','off', 'Name',['Tracker - ' video_path]);
            im_handle = imshow(uint8(org_im), 'Border','tight', 'InitialMag', 100 + 100 * (length(im) < 500));
            rect_handle = rectangle('Position',rect_position(frame,:), 'EdgeColor','g', 'LineWidth', 5);
            % top scored proposals
%             for j = 1 : 4
%                 proposal_handle(j) = rectangle('Position',[0,0,1,1], 'EdgeColor', 'r', 'LineWidth', 3);
%             end
            % preniminary location and previous size
%             pre_target_handle = rectangle('Position', [0,0,1,1], 'EdgeColor', [0.3,0.3,1], 'LineWidth', 4, 'lineStyle', '--');
            % the most promising proposal
%             chosen_proposal_handle = rectangle('Position', [0,0,1,1], 'EdgeColor', [1,1,0], 'LineWidth', 3);
            % frame index
%             text_handle = text(30, 30, ['# ' int2str(frame)]);
%             set(text_handle, 'color', [0 1 1], 'FontSize', 45, 'FontWeight', 'Bold');
        else
            try  %subsequent frames, update GUI
                set(im_handle, 'CData', org_im)
                set(rect_handle, 'Position', rect_position(frame,:));
%                 for j = 1 : min([num_of_proposals 4])
%                     set(proposal_handle(j), 'Position', proposals_xywh(j,:));
%                 end
%                 for j = (min([num_of_proposals 4])+1) : 4
%                     set(proposal_handle(j), 'Position', [0,0,1,1]);
%                 end
%                 set(pre_target_handle, 'Position', pre_target_rect);
%                 set(chosen_proposal_handle, 'Position', chosen_proposal);
%                 set(text_handle, 'string', ['# ' int2str(frame)]);
            catch
                return
            end
        end     
        drawnow
%         pause
    end
    
    if bSaveImage == 1
        imwrite(frame2im(getframe(gcf)),['../.' res_path num2str(frame) '.jpg']);
    end

    
end

if resize_image,
    rect_position = rect_position * 2;
end

fps = num_frames/time;

% Statistics for analyzing the effect of proposal rejection stage:

% avg_proposal_num = total_proposal_num / num_frames;
% avg_proposal_num_after_rejetion = total_flitered_propsal_num / num_frames;
% 
% disp([num2str(avg_proposal_num) ' vs ' num2str(avg_proposal_num_after_rejetion)]);
% 
% dlmwrite(['avg_proposal_num_per_seq.txt'], avg_proposal_num, '-append');
% dlmwrite(['avg_proposal_num_after_rejetion_per_seq.txt'], avg_proposal_num_after_rejetion, '-append');
% 

end


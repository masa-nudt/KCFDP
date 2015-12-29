function x = get_features(im, org_im, hog_orientations, cell_size, cos_window, w2c)
%GET_FEATURES
%   Extracts dense features from image.
%
%   Extracts features specified in struct FEATURES, from image IM. The
%   features should be densely sampled, in cells or intervals of CELL_SIZE.
%   The output has size [height in cells, width in cells, features].
%
%   To specify HOG features, set field 'hog' to true, and
%   'hog_orientations' to the number of bins.
%
%   To experiment with other features simply add them to this function
%   and include any needed parameters in the FEATURES struct. To allow
%   combinations of features, stack them with x = cat(3, x, new_feat).
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/


%HOG features, from Piotr's Toolbox
hog = double(fhog(single(im) / 255, cell_size, hog_orientations));
hog(:,:,end) = [];  %remove all-zeros channel ("truncation feature")

resized_im = imresize(org_im, floor(size(im)/cell_size), 'bilinear');

[gray, cn] = get_patch_feature(resized_im, 'gray', 'cn', w2c);
	
x = cat(3, hog, gray, cn);

%process with cosine window if needed
if ~isempty(cos_window),
    x = bsxfun(@times, x, cos_window);
end

end

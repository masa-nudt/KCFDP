function [out_npca, out_pca] = get_patch_feature(im_patch, non_pca_features, pca_features, w2c)

    size_3d = size(im_patch);
    sz = size_3d(1:2);
    
    % compute non-pca feature map
    if ~isempty(non_pca_features)
        out_npca = get_feature_map(im_patch, non_pca_features, w2c);
    else
        out_npca = [];
    end

    % compute pca feature map
    if ~isempty(pca_features)
%         temp_pca = get_feature_map(im_patch, pca_features, w2c);
          out_pca = get_feature_map(im_patch, pca_features, w2c); % no need to reshape
%         out_pca = reshape(temp_pca, [prod(sz), size(temp_pca, 3)]);
    else
        out_pca = [];
    end

end


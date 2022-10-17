function snr = calc_snr(imageRecon,imageRef)
% imageRecon - reconstructed image
% imageRef - reference image
% 
% Returns the average SNR over the image (or set of images, if multiple
% partitions and/or time frames are passed in).
% ---------------------------------------------------
% MIMOSA code repository
% Author: Jesse Hamilton
% Last updated: Dec 2013
% ---------------------------------------------------

imageRecon = mat2gray(abs(imageRecon));
imageRef = mat2gray(abs(imageRef));

switch ndims(imageRecon)
    case 2
        snr = snr_2d(imageRecon,imageRef);
        
    case 3
        nMeas = size(imageRecon,3);
        snr = zeros(1,nMeas);
        for i = 1:nMeas
            snr(i) = snr_2d(imageRecon(:,:,i),imageRef(:,:,i));
        end
        snr = mean(snr);
        
    case 4
        [~,~,nPartition,nMeas] = size(imageRecon);
        snr = zeros(nPartition,nMeas);
        for i = 1:nPartition
            for j = 1:nMeas
                snr(i,j) = snr_2d(imageRecon(:,:,i,j),imageRef(:,:,i,j));
            end
        end
        snr = mean(mean( snr ));
        
    otherwise
        error('Images have wrong dimensions');
end

end
% END function calc_snr()

function snr = snr_2d(imageRecon,imageRef)
varI = var(imageRef(:),1);
mse = mean( (imageRecon(:)-imageRef(:)).^2 );
snr = 10*log10(varI/mse);
end
% END function snr_2d()
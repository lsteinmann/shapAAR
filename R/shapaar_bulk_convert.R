#' Extracting the conture from a bulk of scanned images in one folder
#'
#' @param p the path to the directory with the images (results will be placed alongside originals)
#' @param ext the file extention to be used; at the moment only PNG is supported)
#' @param thresh threshold value for the comparison of left and right side
#' @param v the velocity of the active contour (default = 1)
#' @param dt the time interval that is used for the development of the contour (default = 1)
#' @param n_iter maximum iteration length of the algorithm (default = 1000)
#' @param buffer buffer of the contour (default = 2)
#' @param vectorlength desired length of the resulting contour vector (default = 100)
#' @param side_of_profile side of the drawing to be used as contour ('both' (mean of l+r), 'left', 'right' defaulting to 'left')
#'
#' @export


shapaar_bulk_convert <- function(p,
                        ext = 'png',
                        thresh = 0.02,
                        v = 1,
                        dt = 1,
                        n_iter = 1000,
                        buffer = 2,
                        vectorlength = 100,
                        side_of_profile = c('both', 'left', 'right')) {

  # as doParallel and parallel are only used with function I put them in Suggests and checking here:
  if (!requireNamespace('doParallel', quietly = TRUE)) {
    stop('doParallel needed for this function to work. Please install it.',
         call. = FALSE)
  }
  if (!requireNamespace('parallel', quietly = TRUE)) {
    stop('parallel needed for this function to work. Please install it.',
         call. = FALSE)
  }

  cl <- parallel::makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  f <- list.files(path = p, pattern = paste('\\.', ext, '$', sep = ''))

  if (length(side_of_profile) != 1) {
    side_of_profile <- 'left'
  }

  # loop over all images
  a <- foreach(i = 1:length(f),
    .packages = c('shapAAR',
    'EBImage')) %dopar% {

    # see https://github.com/ISAAKiel/shapAAR/tree/master/vignettes
    # load image and prepare for further processing
    img <- EBImage::readImage(paste(p, f[i], sep = ''))
    img <- EBImage::normalize(img)^5
    # resize to vectorlength * 4 to avoid large images while keeping them
    # large enough for the resulting vector, multiplier should fit most
    # cases, but consider cropping large whitespace beforehand
    img <- EBImage::resize(img, h = (vectorlength * 5), antialias = T)

    # posterization
    img <- shapAAR::add_canvas(img, 10, 10, center = T)
    img <- EBImage::channel(img, 'grey')
    img <- img > EBImage::otsu(img)
    img <- EBImage::gblur(img, 1)
    img <- round(1 + img * 255)
    img <- img - mean(img)

    # contours
    g <- shapAAR::stopping_fun(img)
    phi <- shapAAR::default_phi(g)
    phi_out <- shapAAR::active_contour(phi,
                                       g,
                                       n_iter = n_iter,
                                       v = v,
                                       dt = dt,
                                       show = FALSE,
                                       buffer = buffer)

    # two-color image
    inner <- EBImage::Image(phi_out <= 0)
    inner <- EBImage::fillHull(inner)
    for (j in 1:buffer) {
      inner <- EBImage::erode(inner)
    }

    # extraction of largest object
    labelled_img <- EBImage::bwlabel(inner)
    features <- EBImage::computeFeatures.shape(labelled_img)
    all_objects <- 1:max(labelled_img)
    biggest_object <- which.max(features[, 1])
    img_biggest_only <- EBImage::rmObjects(labelled_img,
                                           all_objects[all_objects != biggest_object])
    img_biggest_only <- shapAAR::img_crop_background(img_biggest_only)
    fg_points <- which(img_biggest_only != 0, arr.ind = T)
    minbbox <- shapAAR::getMinBBox(fg_points)

    # rotate image to get horizontal orifice line
    moment <- 90 - minbbox$angle
    if (abs(moment) > 45) {moment = moment - sign(moment) * 90}
    img_rect <- EBImage::rotate(img_biggest_only, moment, bg.col='black')
    img_rect <- img_rect > EBImage::otsu(img_rect)

    # crop the background and top and bottom artefacts
    img_rect <- shapAAR::img_crop_background(img_rect)
    img_crop <- shapAAR::crop_bottom(shapAAR::crop_top(img_rect))

    # devide images
    dim_x <- nrow(img_crop)
    splitaxis <- dim_x %/% 2
    forward_ind <- 1:splitaxis
    backward_ind <- sort((dim_x + 1) - forward_ind)
    left_img <- img_crop[forward_ind,]
    right_img <- img_crop[backward_ind,]

    # scale to desired height
    right_img <- EBImage::resize(right_img,
                                 h = vectorlength,
                                 antialias = T)
    # scale to desired height
    left_img <- EBImage::flop(left_img)
    left_img <- EBImage::resize(left_img,
                                h = vectorlength,
                                antialias = T)

    # produce contour vectors
    c_left <- shapAAR::contour_vector_bg(left_img)
    c_right <- shapAAR::contour_vector_bg(right_img)

    # Thumbnail for later use
    img_crop <- EBImage::resize(img_crop,
                                h = vectorlength,
                                antialias = T)
    # invert so thumb is black
    img_crop <- max(img_crop) - img_crop

    # save thumb
    EBImage::writeImage(img_crop,
                        paste(p,
                              f[i],
                              '_thumb.png',
                              sep = ''))

    # saving the contour vector according to specifications
    deviance <- sd(c_left - c_right) / length(c_left)

    if (side_of_profile == 'both' || side_of_profile == 'mean') {
      if (deviance < thresh) {
        c_mean <- round(apply(cbind(c_left, c_right), 1, mean))
        write.table(c_mean, paste(p, f[i], '.txt', sep = ''), sep=';')
        write(paste(f[i], ' - contour produced from mean. Left/Right Deviance of ', deviance, '.'),
              paste(p, 'shapAAR.log', sep = ''), append = TRUE)
        return()
      } else {
        write.table(c_left, paste(p, f[i], '.txt', sep = ''), sep=';')
        write(paste(f[i], 'ERROR - deviance > threshhold of ', thresh,
                    ', used left profile instead of mean. Left/Right Deviance of ', deviance, '.'),
              paste(p, 'shapAAR.log', sep = ''), append = TRUE)
        return()
      }
    } else {
      if (side_of_profile == 'left') {
        write.table(c_left, paste(p, f[i], '.txt', sep = ''), sep=';')
        write(paste(f[i], ' - contour produced from left profile. Left/Right Deviance of ', deviance, '.'),
              paste(p, 'shapAAR.log', sep = ''), append = TRUE)
        return()
      }
      if (side_of_profile == 'right') {
        write.table(c_right, paste(p, f[i], '.txt', sep = ''), sep=';')
        write(paste(f[i], ' - contour produced from right profile. Left/Right Deviance of ', deviance, '.'),
              paste(p, 'shapAAR.log', sep = ''), append = TRUE)
        return()
      }

  }
}
}

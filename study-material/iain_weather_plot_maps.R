library(tidyverse)
library(ggmap)
#library(raster)
library(fda)
library(fda.usc)


# ********************** Formatting Weather Data ***************************** #

index_wrap <- function(index, acceptable_values) {
  if(index == acceptable_values[1] - 1) 
    return(acceptable_values[length(acceptable_values)])
  if(index == acceptable_values[length(acceptable_values)] + 1)
    return(acceptable_values[1])
  return(index)
} # End of index_wrap

fix_na <- function(wdf, cc, max_tol = 300) {
  # for each row
  rows_fixed <- rep(T, nrow(wdf))
  col_nums <- grep(cc, names(wdf))
  for(i in 1:nrow(wdf)) {
    # get only those that are na's
    curr_col_na <- col_nums[is.na(wdf[i,col_nums])]
    if(any(is.na(wdf[i,col_nums])) && length(curr_col_na) <= max_tol) {
      for(j in curr_col_na) {
        ls <- index_wrap(j-1, col_nums)
        rs <- index_wrap(j+1, col_nums)
        if(!is.na(wdf[i, ls]) && !is.na(wdf[i, rs])) {
          wdf[i, j] <- (wdf[i, ls] + wdf[i, rs]) / 2
        }
      }
    } else {
      if(any(is.na(wdf[i,col_nums]))) rows_fixed[i] <- F
    }
    
  }
  
  return(wdf[rows_fixed,])
} # End of fix_na

format_weather_fd <- function(wdf, cc, nsplines, valid_years = 1880:2021, 
                              areas = unique(wdf$area), plain = F) {
  wfds <- wdf %>% dplyr::select(contains(cc, ignore.case = T, vars = NULL))
  k_rows <- c()
  wfds <- fix_na(wfds, cc)
  for(i in 1:nrow(wfds)) if (!any(is.na(wfds[i,]))) k_rows <- append(k_rows, i)
  m <- data.matrix(t(wfds[k_rows,]))
  if(plain) {
    sbd = fda.usc::fdata(t(m))
  } else {
    basis <- create.fourier.basis(c(0,365), nbasis = nsplines)
    sbd <- smooth.basis(argvals = seq(0,365, length.out = nrow(m)), y = m, 
                        fdParobj = basis)$fd
  }
  return(list(fd = sbd, rows = k_rows, labels = wdf$area[k_rows], 
              years = wdf$Year[k_rows], c_labels = wdf$city[k_rows], 
              p_label = wdf$prov[k_rows], k_rows = k_rows))
} # End of format_weather_fd


# plot just all of the data
# plot_weather_station_map_simple()

# plot just a certain area and restrict it to the data shown
# plot_weather_station_map_simple(ar = c("PAC"), use_bb=T)

# *********************** Plotting Weather Data ****************************** #
plot_weather_curves <- function(by_area=F, show_outliers=T, by_group = F,
                                weather_file="weather_2021.csv") {
  
  group.color <- c(PAC = "cyan", SUB = "magenta", 
                   ATL = "orange", CON = "green", 
                   OUT = "red")
  
  w_data <- read_csv(weather_file)
  wfd <- format_weather_fd(w_data, "tPrecip", 65, plain=F)
  
  if(!by_group) {
    if(show_outliers) w_data$area[w_data$outlier == 1] <- "OUT"
    if(!show_outliers) plot(wfd$fd, col=group.color[w_data$area], ylim = c(-10,70),
                            main="Weather Groups Plot", ylab="Total Precipitation",
                            xlab="Day")
    else {
      plot(wfd$fd[w_data$outlier != 1], 
           col=group.color[w_data$area[w_data$outlier != 1]],ylim = c(-10,70),
           main="Weather Groups Plot", ylab="Total Precipitation",
           xlab="Day")
      lines(wfd$fd[w_data$outlier == 1], col=group.color["OUT"])
    }
    legend(x="topleft", legend=unique(w_data$area), lty = 1,
           col=group.color[unique(w_data$area)])
  } else {
    par(mfrow=c(2,2))
    plot(wfd$fd[w_data$area=="CON" & w_data$outlier != 1], 
         col = group.color["CON"], main = "Continental Weather", ylab = "Total Precipitation", 
         xlab = "Day", ylim = c(-10,70))
    lines(wfd$fd[w_data$area=="CON" & w_data$outlier == 1],
         col = group.color["OUT"])
    
    plot(wfd$fd[w_data$area=="PAC" & w_data$outlier != 1], 
         col = group.color["PAC"], main = "Pacific Weather", ylab = "Total Precipitation", 
         xlab = "Day", ylim = c(-10,70))
    lines(wfd$fd[w_data$area=="PAC" & w_data$outlier == 1],
          col = group.color["OUT"])
    
    plot(wfd$fd[w_data$area=="ATL" & w_data$outlier != 1], 
         col = group.color["ATL"], main = "Atlantic Weather", ylab = "Total Precipitation", 
         xlab = "Day", ylim = c(-10,70))
    lines(wfd$fd[w_data$area=="ATL" & w_data$outlier == 1],
          col = group.color["OUT"])

    plot(wfd$fd[w_data$area=="SUB" & w_data$outlier != 1],
         col = group.color["SUB"], main = "Subarctic Weather", ylab = "Total Precipitation",
         xlab = "Day", ylim = c(-10,70))
    lines(wfd$fd[w_data$area=="SUB" & w_data$outlier == 1],
          col = group.color["OUT"])
    
    
  }
}


plot_weather_station_map_simple <- function(ar = NULL, use_bb=F, 
                                            weather_file="weather_2021.csv") {
  # use constant colors
  # include in the final weather data the outliers that are present in a diff color (black?)
  group.color <- c(PAC = "cyan", SUB = "magenta", 
                   ATL = "orange", CON = "green", 
                   OUT = "red")
  
  
  if(is.null(ar)) {  # get all thee data and area 
    title_text <- "Weather Station Locations"
    w_data <- read_csv(weather_file)
    w_data$area[w_data$outlier==1] <- "OUT"
    ca.provinces <- raster::getData("GADM", country="CAN", level = 1)
  } else {  # get only the represented provinces
    title_text <- "Weather Station Locations"
    for(a in ar) {   # add the area names
      title_text <- paste0(title_text," ",a)
    }
    w_data <- read_csv(weather_file) %>%
      filter(area %in% ar)
    # maybe edit this
    w_data$area[w_data$outlier==1] <- "OUT"
    # translate the abbreviations to full name
    p_trans <- c(AB="Alberta", BC="British Columbia", MB="Manitoba",
                 NB="New Brunswick", NL="Newfoundland and Labrador",
                 NT="Northwest Territories", NS="Nova Scotia",
                 NU="Nunavut", ON="Ontario", PE="Prince Edward Island",
                 QC="Quebec", SK="Saskatchewan", YT="Yukon Territoty")
    ca.provinces <- raster::getData("GADM", country="CAN", level = 1)
    ca.provinces <- ca.provinces[ca.provinces$NAME_1 %in% 
                                   p_trans[unique(w_data$prov)],]
  }
  
  # plot the provinces and data
  plt <- ggplot(ca.provinces, aes(x=long, y=lat,group=group)) + 
    geom_path() + 
    coord_map() + 
    geom_point(data=w_data, aes(x=longitude, y=latitude, group=1, fill=area),
               colour = "black", pch=21, size=2) +
    scale_fill_manual(values = group.color) +
    ggtitle(title_text) +
    xlab("Longitude") +
    ylab("Latitude")
  
  if(use_bb) {  # plot using a bounding box
    w_bbox <- make_bbox(lat = latitude, lon = longitude, data = w_data)
    plt <- plt + xlim(w_bbox["left"], w_bbox["right"]) + 
      ylim(w_bbox["bottom"], w_bbox["top"])
  }
  
  # plot the data and area
  plt
  
} # End of plot_weather_station_map_simple

# res3 <- funHDDC::funHDDC(wfd$fd, K=4, model = models, threshold = 0.1, itermax = 200, nb.rep=50, init="kmeans")


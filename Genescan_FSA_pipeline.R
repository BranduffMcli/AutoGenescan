#####Load packages#####
library(tidyverse)
library(data.table)
library(Fragman)

get.scores2<-function (my.scores, mark = "mark") 
{
  y <- lapply(my.scores, function(x) {
    x$hei
  })
  n <- lapply(my.scores, function(x) {
    length(x$hei)
  })
  nn <- max(unlist(n))
  da <- data.frame(matrix(NA, byrow = T, ncol = nn, nrow = length(my.scores)))
  rownames(da) <- names(my.scores)
  for (i in 1:dim(da)[1]) {
    da[i, 1:(length(y[[i]]))] <- y[[i]]
  }
  col.na <- vector(mode = "character")
  for (k in 1:nn) {
    col.na[k] <- paste(mark, "A.", k, sep = "")
  }
  names(da) <- col.na
  return(da)
}

score.markers2<-function (my.inds, channel = 1, n.inds = NULL, panel = NULL, 
                          shift = 0.8, ladder, channel.ladder = NULL, ploidy = 2, left.cond = c(0.6, 
                                                                                                3), right.cond = 0.35, warn = FALSE, window = 0.5, init.thresh = 200, 
                          ladd.init.thresh = 200, method = "iter2", env = parent.frame(), 
                          my.palette = NULL, plotting = TRUE, electro = FALSE, pref = 3) 
{
  oldw <- getOption("warn")
  options(warn = -1)
  ci.upp = 1.96
  ci.low = 1.96
  dev = 50
  thresh = NULL
  if (length(n.inds) > length(my.inds)) {
    print(paste("Hey! you are trying to examine more individuals than the ones you actually read? You selected in 'n.inds' argument", 
                length(n.inds), "individuals but you only provided", 
                length(my.inds), " individuals. Please select a number of individuals smaller or same size than the ones contained in 'my.inds' argument"))
    stop
  }
  else {
    cat(paste("\n1) You have used a shift of", shift, 
              "base pairs. All peaks at that distance from the tallest peak will be ignored and be considered noise. \n2) In addition the window used is", 
              window, ". Which means that all peaks closer by that distance to panel peaks will be accounted as peaks. \n3) Remember using the get.scores() function to extract the results from this output as a dataframe. \n\n"))
  }
  if (method == "ci") {
    print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
  }
  if (is.null(channel.ladder)) {
    channel.ladder <- dim(my.inds[[1]])[2]
  }
  else {
    channel.ladder <- channel.ladder
  }
  if (dim(my.inds[[1]])[2] < channel.ladder) {
    print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
    stop
  }
  if (is.null(n.inds)) {
    n.inds <- c(1:length(my.inds))
  }
  else {
    n.inds <- n.inds
  }
  if (is.null(thresh)) {
    thresh <- rep(list(c(1, 1, 1, 1, 1)), length(my.inds))
  }
  else {
    thresh <- thresh
  }
  count <- 0
  tot <- length(n.inds)
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  my.inds2 <- list(NA)
  thresh2 <- list(NA)
  for (i in 1:length(n.inds)) {
    count <- count + 1
    v1 <- n.inds[i]
    my.inds2[[i]] <- my.inds[[v1]]
    names(my.inds2)[i] <- names(my.inds)[i]
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  ncfp <- c("COLOR 1", "COLOR 2", "COLOR 3", 
            "COLOR 4", "COLOR 5", "COLOR 6")
  if (!is.null(my.palette)) {
    cfp <- rep(my.palette, 100)
  }
  else {
    cfp <- c("cornflowerblue", "chartreuse4", 
             "gold2", "red", "orange", "purple")
  }
  col.list <- list(NA)
  att1 <- numeric()
  list.data <- list(NA)
  if (exists("list.data.covarrubias")) {
    list.data <- env$list.data.covarrubias
  }
  else {
    list.ladders <- lapply(my.inds2, function(x) {
      y <- x[, channel.ladder]
      return(y)
    })
    list.data <- lapply(list.ladders, find.ladder, ladder = ladder, 
                        ci.upp = ci.upp, ci.low = ci.low, draw = F, dev = dev, 
                        warn = warn, method = method, init.thresh = ladd.init.thresh)
  }
  list.models <- lapply(list.data, function(da) {
    y <- da[[3]]
    x <- da[[1]]
    mod <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5), 
              data = da)
    return(mod)
  })
  list.models.inv <- lapply(list.data, function(da) {
    x <- da[[3]]
    y <- da[[1]]
    mod <- lm(y ~ x, data = da)
    return(mod)
  })
  xx <- lapply(my.inds2, function(x, channel) {
    1:length(x[, channel])
  }, channel = channel)
  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)
  for (h in 1:length(xx)) {
    h1 <- n.inds[h]
    count <- count + 1
    newxx <- as.vector(predict(list.models[[h1]], newdata = data.frame(x = xx[[h]])))
    newyy <- my.inds2[[h]][, channel]
    new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  top <- max(unlist(lapply(new.whole.data, function(x) {
    max(x$yy)
  })))
  bott <- min(unlist(lapply(new.whole.data, function(x) {
    min(x$yy)
  })))
  list.weis <- list(NA)
  lower.bounds <- numeric()
  for (k in 1:length(my.inds2)) {
    newtt <- init.thresh
    lower.bounds[k] <- newtt
    plant <- big.peaks.col(new.whole.data[[k]]$yy, newtt)
    plant$wei <- new.whole.data[[k]]$xx[plant$pos]
    plant <- separate(plant, shift, type = "bp")
    list.weis[[k]] <- plant
    if (plotting == TRUE) {
      count <- count + 1
    }
    else {
      count <- count + 2
    }
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  list.weis <- lapply(list.weis, function(x) {
    x$wei <- round(x$wei, digits = 4)
    return(x)
  })
  names(list.weis) <- names(my.inds2)
  if (length(panel) > 0) {
    list.weis <- lapply(list.weis, reals, panel = panel, 
                        shi = shift, ploidy = ploidy, left.cond = left.cond, 
                        right.cond = right.cond, window = window)
    list.weis2 <- lapply(list.weis, FUN = homo.panel, panel = panel, 
                         window = window)
  }
  else {
    list.weis2 <- list.weis
  }
  if (plotting == TRUE) {
    layout(matrix(1:pref, pref, 1))
    if (length(panel) > 0) {
      xm <- round(min(panel, na.rm = TRUE) - 10, digits = 0)
      xl <- round(max(panel, na.rm = TRUE) + 10, digits = 0)
    }
    else {
      xm <- 0
      xl <- max(ladder)
    }
    for (g in 1:length(n.inds)) {
      hh4 <- n.inds[g]
      if (length(which(new.whole.data[[g]]$xx > xm & new.whole.data[[g]]$xx < 
                       xl)) > 0) {
        mylim <- max(new.whole.data[[g]]$yy[which(new.whole.data[[g]]$xx > 
                                                    xm & new.whole.data[[g]]$xx < xl)], na.rm = TRUE) + 
          100
      }
      else {
        mylim = 1000
      }
      if (is.infinite(mylim)) {
        mylim = 1000
      }
      plot(new.whole.data[[g]]$xx, new.whole.data[[g]]$yy, 
           type = "l", col = cfp[channel], xaxt = "n", 
           xlim = c(xm, xl), ylim = c(-200, mylim), ylab = "Intensity", 
           main = paste(ncfp[channel], "electrogram", hh4), 
           xlab = names(list.models)[hh4], lwd = 2, las = 2)
      axis(1, at = c(xm:xl), labels = xm:xl, cex.axis = 0.8)
      abline(v = list.weis[[g]]$wei, lty = 3, col = "red", 
             cex = 0.5)
      abline(h = lower.bounds[g], lty = 2, col = "chocolate", 
             cex = 0.5)
      legend("topright", legend = c("Peak found","Minimum Detected"), 
             col = c("red", "chocolate"), bty = "n", 
             lty = c(3, 3), lwd = c(1, 1), cex = 0.75)
      count <- count + 1
      setTxtProgressBar(pb, (count/tot) * 0.25)
    }
    if (electro == TRUE) {
      if (length(n.inds) > 1) {
        layout(matrix(1, 1, 1))
        forjet <- lapply(new.whole.data, function(x) {
          x$yy
        })
        forjet2 <- matrix(unlist(forjet), ncol = length(new.whole.data), 
                          byrow = F)
        forjet2[which(forjet2 < 0)] <- 0
        jet.colors <- colorRampPalette(c("#00007F", 
                                         "blue", "#007FFF", "cyan", 
                                         "#7FFF7F", "yellow", "#FF7F00", 
                                         "red"))
        palette <- jet.colors(25)
        image(forjet2, col = palette, xaxt = "n", 
              yaxt = "n", main = paste("Electrogram for plants scored in color", 
                                       channel))
        labb <- paste(rep("Plant", (dim(forjet2)[2])), 
                      1:(dim(forjet2)[2]))
        axis(side = 2, at = seq(1/(dim(forjet2)[2]), 
                                1, by = 1/(dim(forjet2)[2])), labels = labb, 
             las = 2, cex.axis = 0.5)
      }
    }
  }
  close(pb)
  options(warn = oldw)
  return(list.weis2)
}

score.markers3<-function (my.inds, channel = 1, n.inds = NULL, panel = NULL, 
                          shift = 0.8, ladder, channel.ladder = NULL, ploidy = 2, left.cond = c(0.6, 
                                                                                                3), right.cond = 0.35, warn = FALSE, window = 0.5, init.thresh = 200, 
                          ladd.init.thresh = 200, method = "iter2", env = parent.frame(), 
                          my.palette = NULL, plotting = TRUE, electro = FALSE, pref = 3) 
{
  oldw <- getOption("warn")
  options(warn = -1)
  ci.upp = 1.96
  ci.low = 1.96
  dev = 50
  thresh = NULL
  if (length(n.inds) > length(my.inds)) {
    print(paste("Hey! you are trying to examine more individuals than the ones you actually read? You selected in 'n.inds' argument", 
                length(n.inds), "individuals but you only provided", 
                length(my.inds), " individuals. Please select a number of individuals smaller or same size than the ones contained in 'my.inds' argument"))
    stop
  }
  else {
    cat(paste("\n1) You have used a shift of", shift, 
              "base pairs. All peaks at that distance from the tallest peak will be ignored and be considered noise. \n2) In addition the window used is", 
              window, ". Which means that all peaks closer by that distance to panel peaks will be accounted as peaks. \n3) Remember using the get.scores() function to extract the results from this output as a dataframe. \n\n"))
  }
  if (method == "ci") {
    print(paste("Please make sure you have used the same 'dev' value you found convenient for your ladder detection or probably your call will not match"))
  }
  if (is.null(channel.ladder)) {
    channel.ladder <- dim(my.inds[[1]])[2]
  }
  else {
    channel.ladder <- channel.ladder
  }
  if (dim(my.inds[[1]])[2] < channel.ladder) {
    print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
    stop
  }
  if (is.null(n.inds)) {
    n.inds <- c(1:length(my.inds))
  }
  else {
    n.inds <- n.inds
  }
  if (is.null(thresh)) {
    thresh <- rep(list(c(1, 1, 1, 1, 1)), length(my.inds))
  }
  else {
    thresh <- thresh
  }
  count <- 0
  tot <- length(n.inds)
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  my.inds2 <- list(NA)
  thresh2 <- list(NA)
  for (i in 1:length(n.inds)) {
    count <- count + 1
    v1 <- n.inds[i]
    my.inds2[[i]] <- my.inds[[v1]]
    names(my.inds2)[i] <- names(my.inds)[i]
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  ncfp <- c("COLOR 1", "COLOR 2", "COLOR 3", 
            "COLOR 4", "COLOR 5", "COLOR 6")
  if (!is.null(my.palette)) {
    cfp <- rep(my.palette, 100)
  }
  else {
    cfp <- c("cornflowerblue", "chartreuse4", 
             "gold2", "red", "orange", "purple")
  }
  col.list <- list(NA)
  att1 <- numeric()
  list.data <- list(NA)
  if (exists("list.data.covarrubias")) {
    list.data <- env$list.data.covarrubias
  }
  else {
    list.ladders <- lapply(my.inds2, function(x) {
      y <- x[, channel.ladder]
      return(y)
    })
    list.data <- lapply(list.ladders, find.ladder, ladder = ladder, 
                        ci.upp = ci.upp, ci.low = ci.low, draw = F, dev = dev, 
                        warn = warn, method = method, init.thresh = ladd.init.thresh)
  }
  list.models <- lapply(list.data, function(da) {
    y <- da[[3]]
    x <- da[[1]]
    mod <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5), 
              data = da)
    return(mod)
  })
  list.models.inv <- lapply(list.data, function(da) {
    x <- da[[3]]
    y <- da[[1]]
    mod <- lm(y ~ x, data = da)
    return(mod)
  })
  xx <- lapply(my.inds2, function(x, channel) {
    1:length(x[, channel])
  }, channel = channel)
  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)
  for (h in 1:length(xx)) {
    h1 <- n.inds[h]
    count <- count + 1
    newxx <- as.vector(predict(list.models[[h1]], newdata = data.frame(x = xx[[h]])))
    newyy <- my.inds2[[h]][, channel]
    new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  top <- max(unlist(lapply(new.whole.data, function(x) {
    max(x$yy)
  })))
  bott <- min(unlist(lapply(new.whole.data, function(x) {
    min(x$yy)
  })))
  list.weis <- list(NA)
  lower.bounds <- numeric()
  for (k in 1:length(my.inds2)) {
    newtt <- init.thresh
    lower.bounds[k] <- newtt
    plant <- big.peaks.col(new.whole.data[[k]]$yy, newtt)
    plant$wei <- new.whole.data[[k]]$xx[plant$pos]
    plant <- separate(plant, shift, type = "bp")
    list.weis[[k]] <- plant
    if (plotting == TRUE) {
      count <- count + 1
    }
    else {
      count <- count + 2
    }
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }
  list.weis <- lapply(list.weis, function(x) {
    x$wei <- round(x$wei, digits = 4)
    return(x)
  })
  names(list.weis) <- names(my.inds2)
  if (length(panel) > 0) {
    list.weis <- lapply(list.weis, reals, panel = panel, 
                        shi = shift, ploidy = ploidy, left.cond = left.cond, 
                        right.cond = right.cond, window = window)
    list.weis2 <- lapply(list.weis, FUN = homo.panel, panel = panel, 
                         window = window)
  }
  else {
    list.weis2 <- list.weis
  }
  if (plotting == TRUE) {
    layout(matrix(1:pref, pref, 1))
    if (length(panel) > 0) {
      xm <- round(min(panel, na.rm = TRUE) - 10, digits = 0)
      xl <- round(max(panel, na.rm = TRUE) + 10, digits = 0)
    }
    else {
      xm <- 0
      xl <- max(ladder)
    }
    for (g in 1:length(n.inds)) {
      hh4 <- n.inds[g]
      if (length(which(new.whole.data[[g]]$xx > xm & new.whole.data[[g]]$xx < 
                       xl)) > 0) {
        mylim <- max(new.whole.data[[g]]$yy[which(new.whole.data[[g]]$xx > 
                                                    xm & new.whole.data[[g]]$xx < xl)], na.rm = TRUE) + 
          100
      }
      else {
        mylim = 1000
      }
      if (is.infinite(mylim)) {
        mylim = 1000
      }
      plot(new.whole.data[[g]]$xx, new.whole.data[[g]]$yy, 
           type = "l", col = cfp[channel], xaxt = "n", 
           xlim = c(xm, xl), ylim = c(-200, mylim), ylab = "Intensity", 
           main = paste(ncfp[channel], "electrogram", hh4), 
           xlab = names(list.models)[hh4], lwd = 2, las = 2)
      axis(1, at = c(xm:xl), labels = xm:xl, cex.axis = 0.8)
      count <- count + 1
      setTxtProgressBar(pb, (count/tot) * 0.25)
    }
    if (electro == TRUE) {
      if (length(n.inds) > 1) {
        layout(matrix(1, 1, 1))
        forjet <- lapply(new.whole.data, function(x) {
          x$yy
        })
        forjet2 <- matrix(unlist(forjet), ncol = length(new.whole.data), 
                          byrow = F)
        forjet2[which(forjet2 < 0)] <- 0
        jet.colors <- colorRampPalette(c("#00007F", 
                                         "blue", "#007FFF", "cyan", 
                                         "#7FFF7F", "yellow", "#FF7F00", 
                                         "red"))
        palette <- jet.colors(25)
        image(forjet2, col = palette, xaxt = "n", 
              yaxt = "n", main = paste("Electrogram for plants scored in color", 
                                       channel))
        labb <- paste(rep("Plant", (dim(forjet2)[2])), 
                      1:(dim(forjet2)[2]))
        axis(side = 2, at = seq(1/(dim(forjet2)[2]), 
                                1, by = 1/(dim(forjet2)[2])), labels = labb, 
             las = 2, cex.axis = 0.5)
      }
    }
  }
  close(pb)
  options(warn = oldw)
  return(list.weis2)
}

#####PARAMETERS TO CHANGE#####
#Folder
folder<-"PUT PATHWAY TO FOLDER WITH FSA SCRIPTS HERE"

#Position to look for peaks (adjust depending on CAG length)
BPselection<-c(400:600)

#Utlised Ladder
ladder_Liz<-c(20,40,60,80,100,114,120,140,160,180,200,214,220,240,250,260,280,300,314,
                 320,340,360,380,400,414,420,440,460,480,500,514,520,540,560,580,600)

#ladder_Liz<-c(20,30,40,60,80,100,114,120,140,160,180,200,214,220,240,250,260,280,300,
              #314,320,340,360,380,400,414,420,440,460,480,500,514,520,540,560,580,600,
              #614,620,640,660,680,700,714,720,740,760,780,800,820,840,850,860,880,900,
              #920,940,960,980,1000,1020,1040,1080,1100,1120,1160,1200)

#Remove peaks manually; Example: removemepls<-c(1,2,3) would remove the first 3 peaks
removemepls<-c()

#####Pipeline part 1: Obtain peaks#####
FSA_list<-storing.inds(folder,channels=5,rawPlot = F,fourier=T)

pdf(paste(folder,'Ladder.pdf',sep=''))
ladder.info.attach(stored=FSA_list,ladder=ladder_Liz,ladd.init.thresh=400,method='iter2',attempt=10)
dev.off()

#ladder.corrector(FSA_list,to.correct = '5F__D_WK3__1_F03_LJ_LOQUS_POLQ_PLATE1A_011.fsa',ladder=ladder_Liz)

pdf(paste(folder,'Mainplots.pdf',sep=''))
res<-score.markers2(my.inds = FSA_list, channel = 1, panel=BPselection,ladder=ladder_Liz,init.thresh = 0,
                    ploidy = 60,left.cond = c(0.1,1),right.cond=c(0.1),window = 1,shift=1.5)
dev.off()

pdf(paste(folder,'UnmarkedMainplots.pdf',sep=''))
res<-score.markers3(my.inds = FSA_list, channel = 1, panel=BPselection,ladder=ladder_Liz,init.thresh = 0,
                    ploidy = 60,left.cond = c(0.1,1),right.cond=c(0.1),window = 1,shift=1.5)
dev.off()

#####Pipeline part 2: Extract data#####
Oregon_final<-left_join(gather(setDT(get.scores(res), keep.rownames = TRUE)[],
                               condition,measurement,2:ncol(setDT(get.scores(res), keep.rownames = TRUE)[])),
                        gather(setDT(get.scores2(res), keep.rownames = TRUE)[],
                               condition,measurement,2:ncol(setDT(get.scores2(res), keep.rownames = TRUE)[])),
                        by=c('rn','condition'))
ifelse(is.null(removemepls), print("No peaks skipped"),
       Oregon_final<-Oregon_final[-c(removemepls), ])
Oregon_final<-Oregon_final %>% tidyr::separate(condition,c('Dye','Dye_no')) 
Oregon_final$Dye_no<-as.numeric(Oregon_final$Dye_no)
Oregon_final<- na.omit(Oregon_final)

#####Pipeline part 3: Determine modal, II, EI & CI#####
##Find modal
Modal_C<-Oregon_final[order(Oregon_final$Dye_no),]
Modal_C<-Modal_C[order(Modal_C$measurement.y,decreasing=T),]
Modal_C<-Modal_C[!duplicated(Modal_C[,c('rn')]),]
Modal_C<-Modal_C %>% dplyr::select(rn,Dye_no)

Oregon_final<-left_join(Oregon_final,Modal_C,by='rn')
Oregon_final$Rel_PSN<-Oregon_final$Dye_no.x - Oregon_final$Dye_no.y

#####Pipeline part 4: Export Modal output#####
Summarytable_c<-Oregon_final
Summarytable_c<-subset(Summarytable_c,Rel_PSN==0)
Summarytable_c$CAG_modal<-(Summarytable_c$measurement.x - 80) / 3
Summarytable_c<-Summarytable_c %>%
  select(rn,CAG_modal)
write.csv(Summarytable_c,paste(folder,'Summary_modalonly.csv',sep=''))

Summarytable_b<-Oregon_final %>%
  select(rn,Dye_no.x,Rel_PSN,measurement.x,measurement.y) %>%
  rename(BPSize = measurement.x,Intensity = measurement.y,Peak_no = Dye_no.x)
write.csv(Summarytable_b,paste(folder,'All_peak_report.csv',sep=''))

files <- list.files(path=folder,pattern = "\\.fsa$")
files<-data.frame(files)
files$startcag=NA
files<-left_join(files,Summarytable_c,by=c('files'='rn')) %>%
  rename(foundcag = CAG_modal)
write.csv(files,paste(folder,'filelist.csv',sep=''))

#####STOP!!! MAKE SURE THE MANIFEST FILE IS UPDATED#####
#####Pipeline part 5: Calculate II/EI/CI/SI#####
manifest<-read.csv(paste(folder,'manifest.csv',sep=''))

##Get max value (modal CAG)
get_max<-Oregon_final %>%
  group_by(rn) %>%
  summarise(max = max(measurement.y))
Oregon_final<-left_join(Oregon_final,get_max,by='rn')

Oregon_final<-left_join(Oregon_final,manifest,by = c('rn'='files'))
Oregon_final$dasCAG<-(Oregon_final$measurement.x - 80) / 3
Oregon_final$CAGchange<- Oregon_final$dasCAG - Oregon_final$startcag

##Normalised Instability/Exp/Contraction indices
BlackRangeMtn<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  summarise(allsum = sum(measurement.y))

Oregon_final<-left_join(Oregon_final,BlackRangeMtn,by='rn')
Oregon_final$NPH<-Oregon_final$measurement.y / Oregon_final$allsum
Oregon_final$N_II<-Oregon_final$NPH * Oregon_final$CAGchange
Oregon_final$N_II2<- Oregon_final$NPH * Oregon_final$Rel_PSN

BlackRangeMtn2<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  summarise(InstabIndex = sum(N_II))

BlackRangeMtn3<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  filter(Rel_PSN > 0) %>%
  summarise(ExpIndex = sum(N_II))

BlackRangeMtn4<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  filter(Rel_PSN < 0) %>%
  summarise(ConIndex = sum(N_II))

BlackRangeMtn5<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  filter(Rel_PSN > 0) %>%
  summarise(NeueExpIndex = sum(N_II2))

##Somatic mosacism
BlueRangeMtn<-Oregon_final %>%
  group_by(rn) %>%
  filter(measurement.y > (0.1*max(measurement.y))) %>%
  filter(Rel_PSN > 0) %>%
  summarise(EIndex_sum = sum(measurement.y))

BlueRangeMtn<-left_join(get_max,BlueRangeMtn,by='rn')
BlueRangeMtn$SomMos<-(BlueRangeMtn$EIndex_sum - BlueRangeMtn$max) / BlueRangeMtn$max

Summarytable_a<-left_join(BlueRangeMtn,Oregon_final,by='rn')
Summarytable_a<-left_join(Summarytable_a,BlackRangeMtn2,by='rn')
Summarytable_a<-left_join(Summarytable_a,BlackRangeMtn3,by='rn')
Summarytable_a<-left_join(Summarytable_a,BlackRangeMtn4,by='rn')
Summarytable_a<-left_join(Summarytable_a,BlackRangeMtn5,by='rn')
Summarytable_a<-subset(Summarytable_a,Rel_PSN==0)
Summarytable_a$CAG_modal<-(Summarytable_a$measurement.x - 80) / 3
Summarytable_a<-Summarytable_a %>%
  select(rn,CAG_modal,CAGchange,InstabIndex,ExpIndex,ConIndex,NeueExpIndex,SomMos)
write.csv(Summarytable_a,paste(folder,'Summary.csv',sep=''))

readdepfiles.adj.func <- function(file,acc, adjfactor) {
  tmpdf <- read.table(file=file, sep="\t", header = F)
  tmpdf$sample <- acc
  colnames(tmpdf) <- c("chr", "loc", "depth", "sample")
  subdf <- tmpdf
  ttnt <- sum(tmpdf$depth)
  subdf$adjdepth <- subdf$depth / ttnt * adjfactor
  
  # names(depl) <- names
  return(subdf)
}

extract.norm.depth <- function(dep, chr, start,end) {
  # print(chr)
  depf <- dep[dep$chr==chr,]
  # print(unique(depf$chr))
  depf <- subset(depf, loc >= start & loc <= end )
  return(depf)
}

set.curvilear <- function(start,end,strand,len,arrow,r1=1.5,r2=2,ws=0.5) {
  if(strand == "+") {
    startloc <- start
    endloc <- end
    arrowl=round(arrow*len/360)
    
  }else{
    startloc <- end
    endloc <- start
    ws <- -ws
    arrowl <- -round(arrow*len/360)
    
  }
  
  if( (end-round(arrow*len/360)) > start) {
    spf <- seq(startloc,endloc-arrowl,ws)
    spb <- seq(endloc-arrowl,startloc,-ws)
    h=(r2-r1)/4
    df <- data.frame(r=c(rep(r2-h,length(spf)),r2,(r1+r2)/2,r1,rep(r1+h,length(spb)),r2-h),
                     theta=c(spf,endloc-arrowl,endloc,endloc-arrowl,spb,startloc))
  }else{
    print("less")
    h=(r2-r1)/4
    df <- data.frame(r=c(r2,(r1+r2)/2,r1,r2),
                     theta=c(startloc,endloc,startloc,startloc))
  }
  return(df)
}

set.polygon <- function(start,end,strand,len,arrow=5,r1=190,r2=210) {
  if(strand == "+") {
    startloc <- start
    endloc <- end
    arrowl=round(arrow*len/360)
    
  }else{
    startloc <- end
    endloc <- start
    arrowl <- -round(arrow*len/360)
    
  }
  if( (end-round(arrow*len/360)) > start) {
    h=(r2-r1)/4
    df <- data.frame(y=c(r2-h,r2-h,r2,(r1+r2)/2,r1,r1+h,r1+h,r2-h),
                     x=c(startloc,endloc-arrowl,endloc-arrowl,endloc,endloc-arrowl,endloc-arrowl,startloc,startloc))
  }else{
    h=(r2-r1)/4
    df <- data.frame(y=c(r2,(r1+r2)/2,r1,r2),
                     x=c(startloc,endloc,startloc,startloc))
  }
  return(df)
}

plot.circular <- function(features,depthfile,chr, select=NULL, source="cir_plot") {
  viztab <- extract.norm.depth(depthfile,chr,1,features$end[features$type=="source"])
  
  maxv <-  max(viztab$adjdepth)
  relative <- maxv/16
  fig <- plot_ly(
    source = source,
    type = 'scatterpolar',
    mode = 'lines'
  ) %>% add_trace(r=maxv+2*relative,theta=seq(1,features$end[features$type=="source"],1),
                  hovertemplate = paste('<i>Locus</i>: %{theta}'), line=list(color="grey"),name="Plasmid") %>%
    add_trace(r=maxv+2*relative,theta=features$end[features$type=="source"],mode = 'markers', marker = list(
      color = 'white',
      size = 15,
      line = list(
        color = 'grey',
        width = 2
      )
    ),name="Source")
  
  newviz <- data.frame(theta=viztab$loc,r=viztab$adjdepth)
  nummax <- features$end[features$type=="source"]
  temp <- data.frame(theta=(1:nummax)[!1:nummax %in% viztab$loc],r=0)
  newviz <- rbind(newviz,temp)
  newviz <- newviz[order(newviz$theta),]
  
  fig <-  fig %>%
    add_trace(
      r = c(newviz$r,newviz$r[1]),
      theta =c(newviz$theta,newviz$theta[1]),
      name="Coverage Depth",
      fill = 'toself'
    )
  
  if(length(features$type[features$type %in% c("rep_origin", "CDS")]) > 0) {
    
    
    cdstab <- features[features$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.curvilear(start=cdstab$start[i],end=cdstab$end[i],
                              strand =cdstab$strand[i] ,
                              len=features$end[features$type=="source"],arrow=5,ws=1,r1=maxv+relative,r2=maxv+3*relative)
      
      fig <-  fig %>%
        add_trace(
          r = testdf$r,
          theta =testdf$theta,
          name=cdstab$label[i],
          fill = 'toself')
      
    }
  }
  
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.curvilear(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],len=features$end[features$type=="source"],
                              ws=1,arrow=5,r1=maxv+4*relative,r2=maxv+6*relative)
      
      fig <-  fig %>%
        add_trace(
          r = testdf$r,
          theta =testdf$theta,
          name=selected$label[i],
          fill = 'toself')
    }
  }
  
  
  fig <- fig %>% layout(
    polar=list(
      angularaxis = list(
        direction="clockwise",
        type="category",
        tickfont=list(size=10),
        nticks=12,
        #  period=1
        #  tick0=0,
        #  dtick=4033
        color="grey"
        # categoryarray = c("a", "b", "c", "other","CDS")
      ),
      radialaxis=list(
        color="blue",
        # type="category",
        visible=F
        #    categoryarray = c("a", "b", "c","d","e","f", "other","CDS")
      )
    )
  )
  
  return(fig) 
}

plot.circular.in <- function(features,depthfile,chr, select=NULL, source="cir_plot") {
  viztab <- extract.norm.depth(depthfile,chr,1,features$end[features$type=="source"])
  
  maxv <-  max(viztab$adjdepth) 
  relative <- maxv/16
  addv <- maxv+ 4*relative
  fig <- plot_ly(
    source = source,
    type = 'scatterpolar',
    mode = 'lines'
  ) %>% add_trace(r=addv-3*relative,theta=seq(1,features$end[features$type=="source"],1),
                  hovertemplate = paste('<i>Locus</i>: %{theta}'), line=list(color="grey"),name="Plasmid") %>%
    add_trace(r=addv-3*relative,theta=features$end[features$type=="source"],mode = 'markers', marker = list(
      color = 'white',
      size = 15,
      line = list(
        color = 'grey',
        width = 2
      )
    ),name="Source")
  
  newviz <- data.frame(theta=viztab$loc,r=viztab$adjdepth)
  nummax <- features$end[features$type=="source"]
  temp <- data.frame(theta=(1:nummax)[!1:nummax %in% viztab$loc],r=0)
  newviz <- rbind(newviz,temp)
  newviz <- newviz[order(newviz$theta),]
  
  fig <-  fig %>%
    add_trace(
      r = c(newviz$r + addv,rep(addv,length(newviz$r))),
      theta =c(newviz$theta,rev(newviz$theta)),
      name="Coverage Depth",
      fill = 'toself'
    ) 
  
  if(length(features$type[features$type %in% c("rep_origin", "CDS")]) > 0) {
    
    
    cdstab <- features[features$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.curvilear(start=cdstab$start[i],end=cdstab$end[i],
                              strand =cdstab$strand[i] ,
                              len=features$end[features$type=="source"],arrow=5,
                              ws=1,r1=addv-4*relative,r2=addv-2*relative)
      
      fig <-  fig %>%
        add_trace(
          r = testdf$r,
          theta =testdf$theta,
          name=cdstab$label[i],
          fill = 'toself')
      
    }
  }
  
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.curvilear(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],
                              len=features$end[features$type=="source"],
                              ws=1,arrow=5,r1=addv-7*relative,r2=addv-5*relative)
      
      fig <-  fig %>%
        add_trace(
          r = testdf$r,
          theta =testdf$theta,
          name=selected$label[i],
          fill = 'toself')
    }
  }
  
  
  fig <- fig %>% layout(
    polar=list(
      angularaxis = list(
        direction="clockwise",
        type="category",
        tickfont=list(size=10),
        nticks=12,
        #  period=1
        #  tick0=0,
        #  dtick=4033
        color="grey"
        # categoryarray = c("a", "b", "c", "other","CDS")
      ),
      radialaxis=list(
        color="blue",
        # type="category",
        visible=F
        #    categoryarray = c("a", "b", "c","d","e","f", "other","CDS")
      )
    )
  )
  
  return(fig) 
}

plot.linear <- function(features,depthfile,chr, select=NULL, source="lin_plot") {
  viztab <- extract.norm.depth(depthfile,chr,1,features$end[features$type=="source"])
  
  maxv <-  max(viztab$adjdepth)
  
  relative <- maxv/16
  
  fig <- plot_ly(
    source = source,
    mode = 'lines'
  ) %>% add_trace(x=1:features$end[features$type=="source"],
                  y=-2*relative, 
                  hovertemplate = paste('<i>Locus</i>: %{y}'),
                  line=list(color="grey"),name="Plasmid") %>%
    add_trace(x=0,
              y=-2*relative,
              mode = 'markers', marker = list(
                color = 'white',
                size = 15,
                line = list(
                  color = 'grey',
                  width = 2
                )
              ),name="Source")
  
  fig <-  fig %>%
    add_trace(
      y = viztab$adjdepth,
      x =viztab$loc,
      hovertemplate = paste('Locus: %{x}',
                            '<br>Depth: %{y}<br>'),
      name="Coverage Depth",
      fill = 'tozeroy'
    )
  
  if(length(features$type[features$type %in% c("rep_origin", "CDS")]) > 0) {
    
    cdstab <- features[features$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.polygon(start=cdstab$start[i],
                            end=cdstab$end[i],
                            strand =cdstab$strand[i] ,
                            len=features$end[features$type=="source"],arrow=10,r1=-relative,r2=-3*relative)
      
      
      fig <-  fig %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=cdstab$label[i],
          fill = 'toself')
      
    }
  }
  
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.polygon(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],len=features$end[features$type=="source"],arrow=10,r1=-4*relative,r2=-6*relative)
      
      fig <-  fig %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=selected$label[i],
          fill = 'toself')
    }
  }
  
  return(fig)
} 

plot.boxes <- function(features,depthfile,chr, select=NULL, source="box_plot") {
  
  
  fig <- plot_ly(
    source = source) 
  
  if(length(features$type[features$type=="CDS"]) > 0) {
    cdstab <- features[features$type=="CDS",]
    
    for(i in 1:length(cdstab$type)) {
      
      values <- extract.norm.depth(depthfile,chr,cdstab$start[i],cdstab$end[i])$adjdepth
      
      fig <-  fig %>%
        add_boxplot(
          x = cdstab$label[i],
          y =values,
          name=cdstab$label[i])
    }
  }
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      sel <- extract.norm.depth(depthfile,chr,selected$start[i],selected$end[i])$adjdepth
      
      fig <-  fig %>%
        add_boxplot(
          x = selected$label[i],
          y =sel,
          name=selected$label[i])
    }
    
  }
  
  
  
  return(fig)
} 

plot.linear.compare <- function(features,features2,depthfile,depthfile2,chr,chr2, select=NULL, select2=NULL,source="lin_plot_comp") {
  viztab <- extract.norm.depth(depthfile,chr,1,features$end[features$type=="source"])
  viztab2 <- extract.norm.depth(depthfile2,chr2,1,features2$end[features2$type=="source"])
  
  maxv <-  max(viztab$adjdepth,viztab2$adjdepth)
  relative <- maxv/16
  fig <- plot_ly(
    source = source,
    mode = 'lines'
  ) 
  
  fig1 <-  fig %>%
    add_trace(
      y = viztab$adjdepth,
      x =viztab$loc,
      hovertemplate = paste(
        '<br>Locus: %{x}<br>',
        '<br>Depth: %{y}<br>'),
      name=paste("Coverage Depth:",unique(depthfile$sample)),
      fill = 'tozeroy'
    ) %>%
    add_trace(
      y = -viztab2$adjdepth,
      x =viztab2$loc,
      customdata = viztab2$adjdepth,
      hovertemplate = paste(
        '<br>Locus: %{x}<br>',
        '<br>Depth: %{customdata}<br>'),
      name=paste("Coverage Depth:",unique(depthfile2$sample)),
      fill = 'tozeroy'
    ) %>% layout(yaxis=list(showticklabels=F))
  
  fig2 <- plot_ly(mode = 'lines') %>% 
    add_lines(x=1:features$end[features$type=="source"],
              y=0, 
              hovertemplate = paste('<i>Locus</i>: %{y}'),
              line=list(color="grey"),name="Plasmid") %>%
    add_trace(x=0,
              y=0,
              mode = 'markers', marker = list(
                color = 'white',
                size = 15,
                line = list(
                  color = 'grey',
                  width = 2
                )
              ),name="Source") 
  
  
  
  if(length(features$type[features$type %in% c("rep_origin", "CDS")]) > 0) {
    
    cdstab <- features[features$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.polygon(start=cdstab$start[i],
                            end=cdstab$end[i],
                            strand =cdstab$strand[i] ,
                            len=features$end[features$type=="source"],arrow=10,r1=-1*relative,r2=1*relative)
      
      
      fig2 <-  fig2 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=cdstab$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
      
    }
  }
  
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.polygon(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],
                            len=features$end[features$type=="source"],
                            arrow=10,r1=-2*relative,r2=-4*relative)
      
      fig2 <-  fig2 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=selected$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
    }
  }
  
  fig3 <- plot_ly(mode = 'lines') %>% 
    add_lines(x=1:features$end[features$type=="source"],
              y=0, 
              hovertemplate = paste('<i>Locus</i>: %{y}'),
              line=list(color="grey"),name="Plasmid") %>%
    add_trace(x=0,
              y=0,
              mode = 'markers', marker = list(
                color = 'white',
                size = 15,
                line = list(
                  color = 'grey',
                  width = 2
                )
              ),name="Source")
  
  if(length(features2$type[features2$type %in% c("rep_origin", "CDS")]) > 0) {
    
    cdstab <- features2[features2$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.polygon(start=cdstab$start[i],
                            end=cdstab$end[i],
                            strand =cdstab$strand[i] ,
                            len=features2$end[features2$type=="source"],arrow=10,r1=-1*relative,r2=1*relative)
      
      
      fig3 <-  fig3 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=cdstab$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
      
    }
  }
  
  if(length(select2) >0) {
    selected <- features2[features2$label %in% select2, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.polygon(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],len=features$end[features$type=="source"],
                            arrow=10,r1=-2*relative,r2=-4*relative)
      
      fig3 <-  fig3 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=selected$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
    }
  }
  
  return(subplot(fig2, fig1,fig3,nrows = 3,shareX = T,heights=c(0.15,0.7,0.15)))
}

plot.lin.comp.multi <- function(features,depthfiles,chrs, select=NULL, source="lin_align_plot") {
  viztab <- lapply(depthfiles, function(depthfile) {
    extract.norm.depth(depthfile,chrs,1,features$end[features$type=="source"])
  })
  names(viztab) <- names(depthfiles)
  
  fig2 <- plot_ly(mode = 'lines') %>% 
    add_lines(x=1:features$end[features$type=="source"],
              y=0, 
              hovertemplate = paste('<i>Locus</i>: %{y}'),
              line=list(color="grey"),name="Plasmid") %>%
    add_trace(x=0,
              y=0,
              mode = 'markers', marker = list(
                color = 'white',
                size = 15,
                line = list(
                  color = 'grey',
                  width = 2
                )
              ),name="Source") 
  
  
  if(length(features$type[features$type %in% c("rep_origin", "CDS")]) > 0) {
    
    cdstab <- features[features$type %in% c("rep_origin", "CDS"),]
    for(i in 1:length(cdstab$type)) {
      
      
      testdf <- set.polygon(start=cdstab$start[i],
                            end=cdstab$end[i],
                            strand =cdstab$strand[i] ,
                            len=features$end[features$type=="source"],arrow=10,r1=-10,r2=10)
      
      
      fig2 <-  fig2 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=cdstab$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
      
    }
  }
  
  if(length(select) >0) {
    selected <- features[features$label %in% select, ]
    
    for(i in 1:length(selected$start)) {
      testdf <- set.polygon(start=selected$start[i],end=selected$end[i],strand=selected$strand[i],len=features$end[features$type=="source"],arrow=10,r1=-20,r2=-40)
      
      fig2 <-  fig2 %>%
        add_trace(
          x = testdf$x,
          y =testdf$y,
          name=selected$label[i],
          fill = 'toself') %>% layout(yaxis=list(showticklabels=F))
    }
  }
  
  figs <- list(fig2)
  
  for(i in 1:length(names(viztab))) {
    fig <- plot_ly(
      source = source,
      mode = 'lines'
    ) %>%
      add_trace(
        y = viztab[[i]]$adjdepth,
        x =viztab[[i]]$loc,
        hovertemplate = paste('Locus: %{x}',
                              '<br>Depth: %{y}<br>'),
        name=paste("Coverage Depth:", names(viztab)[i]),
        fill = 'tozeroy'
      )
    figs <- c(figs, list(fig))
    
  }
  
  numpl <- length(names(viztab))+1
  figs <- c(figs,list(nrows=numpl,shareX=T,heights=c(1/numpl/2, rep(1/numpl, length(names(viztab))))))
  
  
  return(do.call("subplot",figs))
}
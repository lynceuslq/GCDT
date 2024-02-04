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

gcl.sankey <- function(sanktab,linkvals=1, source="sankey_plot") {
  
  if(length(sanktab$GC_Units$items) >0) {
    nodetab <- foreach::foreach(a=1:length(sanktab$GC_Units$items),.combine = "rbind") %do% data.frame(items=sanktab$GC_Units$items[[a]],unit=paste0(sanktab$GC_Units$type[[a]],":",sanktab$GC_Units$label[[a]]),color=sanktab$GC_Units$color[[a]])
  }
  
  if(length(sanktab$Link_Info$donor) >0) {
    linklist = list(
      source = match(sanktab$Link_Info$donor,nodetab$items)-1,
      target = match(sanktab$Link_Info$recept, nodetab$items)-1,
      value =  rep(1,length(sanktab$Link_Info$donor)),
      label =  as.character(sanktab$Link_Info$label)
    )
  }
  
  if(length(sanktab$GC_Units$items) >0 & length(sanktab$Link_Info$donor) >0) {
    fig <- plot_ly(
      source=source,
      type = "sankey",
      orientation = "h",
      
      
      node = list(
        label = nodetab$items,
        color = nodetab$color,
        customdata=nodetab$unit,
        hovertemplate=paste("<b>%{customdata}</b><br><br>"),
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      
      
      link = linklist
    )
    fig <- fig %>% layout(
      title = "Flowchart for Genetic Circuit Logic",
      font = list(
        size = 10
      )
    )
    
    return(fig)
  }
  
  
}

# Genetic elements

semicir_x <- function(center=c(0,0), diameter=1, npoints=10){
  tt <- seq(0, pi, length.out=npoints)
  return(center[1] + diameter / 2 * cos(tt))
  
}

semicir_y <- function(center=c(0,0), diameter=1, npoints=10){
  tt <- seq(0, pi, length.out=npoints)
  return(center[2] + diameter / 2 * sin(tt))
  
}

set.point <- function(post,l,h=1) {
  return(
    data.frame(
      x=c(post-l/2,post,post+l/2,post-l/2),
      y=c(-3^(1/2)/2*l*h,0,-3^(1/2)/2*l*h,-3^(1/2)/2*l*h)
    )
  )
}

cir_x <- function(center=c(0,0), diameter=1, npoints=20){
  tt <- c(seq(-pi/2, 3/2*pi, length.out=npoints))
  
  
  return(c(center[1] + diameter / 2 * cos(tt)))
  
}

cir_y <- function(center=c(0,0), diameter=1, npoints=20){
  tt <- c(seq(-pi/2, 3/2*pi, length.out=npoints))
  
  return(c(center[2] + diameter / 2 * sin(tt)))
  
}

ele.func.suite <- list(
  promoter=function(post,len,h=1,arrow=0.3) {
    return(
      data.frame(
        x=c(post,post,post+len-arrow*len,post+len-arrow*len,post+len,post+len-arrow*len,post+len-arrow*len,post,post),
        y=c(0,len*h,len*h,len*h+arrow*len*h*0.4,len*h,len*h-arrow*len*h*0.4,len*h,len*h,0)
      )
    )
  },
  promoter.rc=function(post,len,h=1,arrow=0.3) {
    return(
      data.frame(
        x=c(post,post,post-len+arrow*len,post-len+arrow*len,post-len,post-len+arrow*len,post-len+arrow*len,post,post),
        y=c(0,-len*h,-len*h,-len*h-arrow*len*h*0.4,-len*h,-len*h+arrow*h*len*0.4,-len*h,-len*h,0)
      )
    )
  },
  RBS=function(post,d) {
    return(
      data.frame(
        x=c(semicir_x(center=c(post,0),diameter=d,npoints=20),post+d/2),
        y=c(semicir_y(center=c(post,0),diameter=d,npoints=20),0)
      )
    )
  },
  RBS.rc = function(post,d) {
    return(
      data.frame(
        x=c(semicir_x(center=c(post,0),diameter=d,npoints=20),post+d/2),
        y=-c(semicir_y(center=c(post,0),diameter=d,npoints=20),0)
      )
    )
  },
  terminator=function(post,len,h=1) {
    return(
      data.frame(
        x=c(post-len/40,post-len/40,post-len/40-len/2,post-len/40-len/2,post+len/40+len/2,post+len/40+len/2,post+len/40,post+len/40,post-len/40),
        y=c(0,h*len-h*len/20,h*len-h*len/20,h*len,h*len,h*len-h*len/20,h*len-h*len/20,0,0)
      )
    )
  },
  terminator.rc=function(post,len,h=1) {
    return(
      data.frame(
        x=c(post-len/40,post-len/40,post-len/40-len/2,post-len/40-len/2,post+len/40+len/2,post+len/40+len/2,post+len/40,post+len/40,post-len/40),
        y=-c(0,h*len-h*len/20,h*len-h*len/20,h*len,h*len,h*len-h*len/20,h*len-h*len/20,0,0)
      )
    )
  },
  insulator=function(post,h,d=0.5) {
    return(
      data.frame(
        x=c(post,cir_x(center=c(post,h),diameter=d*h),post),
        y=c(0,cir_y(center=c(post,h),diameter=d*h),0)
      ))
  },
  insulator.rc=function(post,h,d=0.5) {
    return(
      data.frame(
        x=c(post,cir_x(center=c(post,h),diameter=d*h),post),
        y=-c(0,cir_y(center=c(post,h),diameter=d*h),0)
      ))
  },
  CDS=function(start,end,arrow=0.2,h=1) {
    len=end-start
    return(data.frame(
      x=c(start,start+(1-arrow)*len,start+(1-arrow)*len,end,start+(1-arrow)*len,start+(1-arrow)*len,start,start),
      y=c(h*len/10/2,h*len/10/2,h*len/10,0,-h*len/10,-h*len/10/2,-h*len/10/2,h*len/10/2)
    ))
  },
  CDS.rc=function(start,end,arrow=0.2,h=1) {
    len=end-start
    return(data.frame(
      x=c(end,end-(1-arrow)*len,end-(1-arrow)*len,start,end-(1-arrow)*len,end-(1-arrow)*len,end,end),
      y=c(h*len/10/2,h*len/10/2,h*len/10,0,-h*len/10,-h*len/10/2,-h*len/10/2,h*len/10/2)
    ))
  }
)


plot.with.ele.height <- function(profile,chr,featuretab) {
  leng <- featuretab$end-featuretab$start 
  if(leng <= 80) {
    start <- round(featuretab$start - 30 )
    end <- round(featuretab$end + 30)
  }else{
    start <- round(featuretab$start - leng*0.4 )
    end <- round(featuretab$end + leng*0.4)
  }
  
  depdata <- extract.norm.depth(profile,chr,start,end)
  height <- max(depdata$adjdepth)
  return(height)
}

plot.with.ele <- function(profile,chr,featuretab,pointstart, pointend) {
  name <- featuretab$label
  leng <- featuretab$end-featuretab$start 
  if(leng <= 80) {
    start <- round(featuretab$start - 30 )
    end <- round(featuretab$end + 30)
  }else{
    start <- round(featuretab$start - leng*0.4 )
    end <- round(featuretab$end + leng*0.4)
  }
  
  type <- featuretab$type
  strand <- featuretab$strand
  
  depdata <- extract.norm.depth(profile,chr,start,end)
  height <- max(depdata$adjdepth)
  ratio <- height/leng
  
  pointdata <- set.point(pointstart,ratio*leng/12,h=1)
  pointdata2 <- set.point(pointend,ratio*leng/12,h=1)
  
  fig <- plot_ly() %>%
    add_trace(x=depdata$loc,y=depdata$adjdepth,mode="lines",fill="tozeroy",
              name="RNAP Flux",
              hovertemplate = paste(
                '<br>Locus: %{x}<br>',
                '<br>Depth: %{y}<br>')) %>%
    add_trace(x=pointdata$x,y=pointdata$y+depdata$adjdepth[depdata$loc== pointstart],mode="lines",line = list(
      width = 2,color="black" ),name="start", fill="toself",fillcolor="white") %>%
    add_trace(x=pointdata2$x,y=pointdata2$y+depdata$adjdepth[depdata$loc== pointend],mode="lines",line = list(
      width = 2,color="black"  ),name="end", fill="toself",fillcolor="white")
  
  if(strand=="-") {
    eletype <- paste0(type,".rc")
  }else{
    eletype <- type
  }
  
  if(eletype %in% names(ele.func.suite)) {
    myfunc <-  ele.func.suite[[eletype]]
    
    fig <- fig  %>%
      add_trace(x=depdata$loc,y=rep(height*1.1,length(depdata$loc)),name="Plasmid",mode="lines",line = list(
        width = 3,
        color="grey"
      ))
    
    color="grey"
    fillcolor ="white"
    
    
    if(eletype %in% c("promoter","promoter.rc")){
      eletab <- myfunc(pointstart,leng/6,h=2*ratio,arrow=0.3)
      
      fig <- fig %>%
        add_trace(x=eletab$x,y=height*1.1+eletab$y,mode="lines",line = list(
          color=color,
          width = 3
        ),names=name, fill="toself",fillcolor=fillcolor)
    }
    
    if(eletype %in% c("RBS","RBS.rc")){
      eletab <- myfunc(pointstart,height/8)
      
      fig <- fig %>%
        add_trace(x=eletab$x,y=height*1.1+eletab$y,mode="lines",line = list(
          color=color,
          width = 3
        ),names=name, fill="toself",fillcolor=fillcolor)
    }
    
    if(eletype %in% c("terminator","terminator.rc")){
      eletab <- myfunc(pointstart,leng/6,h=2*ratio)
      
      fig <- fig %>%
        add_trace(x=eletab$x,y=height*1.1+eletab$y,mode="lines",line = list(
          color=color,
          width = 3
        ),names=name, fill="toself",fillcolor=fillcolor)
    }
    
    if(eletype %in% c("insulator","insulator.rc")){
      eletab <- myfunc(pointstart,leng/4)
      
      fig <- fig %>%
        add_trace(x=eletab$x,y=height*1.1+eletab$y,mode="lines",line = list(
          color=color,
          dash="dot",
          width = 3
        ),names=name, fill="toself",fillcolor=fillcolor)
    }
    
    if(eletype %in% c("CDS","CDS.rc")){
      eletab <- myfunc(pointstart,pointend,h=0.4*ratio)
      
      fig <- fig %>%
        add_trace(x=eletab$x,y=height*1.1+eletab$y,mode="lines",line = list(
          color=color,
          width = 3
        ),names=name, fill="toself",fillcolor=fillcolor)
    }
  }
  
  
  
  return(fig)
}

plot.pl.reg <- function(chr,features,logicfile,adj=100, source="regulation_plot") {
  fig <- plot_ly(source=source)
  
  if(length(logicfile$GC_Units$items) >0)
  {
    
    elelist <- foreach::foreach(a=logicfile$GC_Units$items,.combine = "c") %do% a
    
    newf <- features[na.omit(match(elelist,features$label)),]
    
    lenpl <- round(max(newf$end)/adj - min(newf$start)/adj)
    startpl <- round(min(newf$start)/adj - lenpl/10)
    endpl <- round(max(newf$end)/adj + lenpl/10)
    
    fig <- fig %>%
      add_trace(x=startpl:endpl,
                y=rep(0,length(startpl:endpl)),mode="lines",line = list(
                  width = 3,
                  color="grey"
                ),name=paste("Plasmid",chr))
    
    
    if(length(elelist) >0) {
      for( x in elelist )  {
        
        if(x %in% features$label ) {
          print(x)
          featuretab <- na.omit(features[features$label == x,])
          for(i in 1:length(featuretab$label)) {
            name <- featuretab$label[i]
            type <- featuretab$type[i]
            strand <- featuretab$strand[i]
            
            if(strand=="-") {
              eletype <- paste0(type,".rc")}else{
                eletype <- type
              }
            
            if(eletype %in% names(ele.func.suite)) {
              myfunc <-  ele.func.suite[[eletype]]
              
              len=1
              hei <- 2/len
              
              if(eletype %in% c("promoter","promoter.rc")){
                eletab <- myfunc(featuretab$start[i]/adj,len,h=hei,arrow=0.3)
                fig <- fig   %>%
                  add_trace(x=eletab$x,y=eletab$y,mode="lines",
                            name=name,line = list(
                              width = 2
                            ), fill="toself")
              }
              
              if(eletype %in% c("RBS","RBS.rc")){
                eletab <- myfunc(featuretab$start[i]/adj,d=len/2)
                fig <- fig   %>%
                  add_trace(x=eletab$x,y=eletab$y,mode="lines",
                            name=name,line = list(
                              width = 2
                            ), fill="toself")
              }
              
              if(eletype %in% c("terminator","terminator.rc")){
                eletab <- myfunc(featuretab$start[i]/adj,len,h=hei)
                fig <- fig   %>%
                  add_trace(x=eletab$x,y=eletab$y,mode="lines",
                            name=name,line = list(
                              width = 2
                            ), fill="toself")
              }
              
              if(eletype %in% c("insulator","insulator.rc")){
                eletab <- myfunc(featuretab$start[i]/adj,h=1,d=len/4)
                fig <- fig   %>%
                  add_trace(x=eletab$x,y=eletab$y,mode="lines",
                            name=name,line = list(
                              width = 2,
                              dash="dot"
                            ), fill="toself")
              }
              
              if(eletype %in% c("CDS","CDS.rc")){
                eletab <- myfunc(featuretab$start[i]/adj,featuretab$end[i]/adj,h=10/(featuretab$end[i]/adj-featuretab$start[i]/adj))
                fig <- fig   %>%
                  add_trace(x=eletab$x,y=eletab$y,mode="lines",
                            name=name,line = list(
                              width = 2
                            ), fill="toself")
              }
              
              
              
            }
          }
          
        }
        
      }
    }
    
  }
  fig <- fig   %>%
    layout(xaxis=list(showticklabels=F),yaxis=list(showticklabels=F,range=c(-5,5)))
  
  return(fig)
}

# RNAP flux calculations
rnap.func <- function(x, pars) {
  #str(pars)
  Jx <- as.numeric(pars["cpratio"] * pars["ratio1"] * pars["convert1"] * 10^(pars["convert2"] * log10(pars["gamma"] * x)))
  return(Jx)
}

jskew.func <- function(start, end, profile, pars) {
  jxlist <- c()
  acc=start:(end+1)
  for (i in 1:(length(acc)-1)) {
    jxlist[i] <- rnap.func(profile$adjdepth[profile$loc==acc[i+1]], pars) / rnap.func(profile$adjdepth[profile$loc==acc[i]], pars)
  }
  return(data.frame(loc=start:end,value=jxlist))
}

sumjx.func <- function(loc, window, profile, pars) {
  jxlist <- c()
  acc <- loc:(loc+window)
  for(i in 1:length(acc)) {
    jxlist[i] <- rnap.func(profile$adjdepth[acc[i]], pars)
  }
  return(sum(jxlist))
}

aws.func <- function(start, end, window=5, profile, pars) {
  awslist <- c()
  if(start < end ) {
    acc <- start:(end+1)
    for(i in 1:(length(acc)-1)){
      if(acc[i]-window <= 0 | acc[i] + window >= end) {
        awslist[i] <- 0
      }else {
        awslist[i] <- sumjx.func(acc[i]-window, window, profile=profile, pars=pars) / sumjx.func(acc[i], window, profile=profile, pars=pars)
      }
    }
    return(data.frame(loc=start:end,value=awslist))
  }else{
    print("invalid input! please check your start and end loci")
  }
}

findpeak.func <- function(vec, threshold) {
  post <- (1:length(vec))[vec >= threshold]
  #print(paste("site at positions beyond threshold", post, sep=" "))
  peak <- c()
  post <- na.omit(post)
  for(i in post) {
    if(length(vec[(i-1):(i+1)][is.na(vec[(i-1):(i+1)])]) == 0) {
      if(vec[i] > vec[i+1] & vec[i] > vec[i-1]){
        peak <- c(peak, i)
      }
    }else{print(i)}
  }
  return(peak)
} 

termi.str.func <- function(locus, gap, window, profile, pars) {
  if(locus - gap - window > 0 ) {
    ter <- sumjx.func(locus - gap - window,  window, profile=profile, pars=pars) /  sumjx.func(locus + gap, window, profile=profile, pars=pars)
    return(ter)
  }else{
    print("invalid input! please check your locus")
  }
}

prom.str.func <- function(locus, gap, window, profile, pars) {
  if(locus - gap - window > 0 ) {
    prom <-  sumjx.func(locus + gap, window, profile=profile, pars=pars) - sumjx.func(locus - gap - window,  window, profile=profile, pars=pars)
    gapsum <-  sumjx.func(locus - gap , 2*gap, profile=profile, pars=pars)
    prom <- prom /gapsum
    return(prom)
  }else{
    print("invalid input! please check your locus")
  }
}

ce.func <- function(jx1, jx2) {
  ce <- (jx2 - jx1) / jx2
  return(ce)
}

diffj.func <- function(locus, gap, window, profile, pars) {
  if(locus -gap - window > 0 ) {
    dj <- sumjx.func(locus + gap, window, profile=profile, pars=pars) / window -  sumjx.func(locus - gap -window, window, profile=profile, pars=pars) / window
    return(dj)
  }else{
    print("invalid input! please check your locus")
  }
}

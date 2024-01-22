#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(bslib)
library(rlang)
library(plotly)
library(plyr)
library(foreach)
library(dplyr)
library(readxl)
library(xlsx)
library(BiocManager)
#options(repos = BiocManager::repositories())
library(genbankr)
library(htmltools)
library(DT)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


source("global.R")

# Define server logic required to draw a histogram
function(input, output, session) {
  
  rv <- reactiveValues(depinaccs=list(),pfinaccs=list())
  
  rv$depth <- reactive({
    req(rv$project)
    file <- input$uploaddep
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    shiny::validate(need(ext %in% c("depth"), "Please upload a file in the right format"))
    
    readdepfiles.adj.func(file=file$datapath,acc=rv$project,adjfactor=input$adjfactor)
    
  })
  
  rv$logic <- reactive({
    file <- input$uploadlogic
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    shiny::validate(need(ext %in% c("json"), "Please upload a file in the right format"))
    
    jsonlite::fromJSON(jsonlite::read_json(file$datapath)[[1]])
    
  })
  
  
  output$cirfiles <- renderUI({
    req(rv$project)
    
    if(is.null(rv[[rv$project]]$depfile)) {
      mydiv <- div()
    }else{
      depf <- rv[[rv$project]]$depfile
      rv[[rv$project]]$chrs <- unique(depf$chr)
       print(rv[[rv$project]]$chrs)
      mydiv <- div(foreach::foreach(a=rv[[rv$project]]$chrs) %do% do.call("fileInput", list(inputId=paste0("file_",a), 
                                                                              label=paste0("Plasmid file for ", a,": "), 
                                                                              multiple = FALSE,
                                                                              accept = c(".gb",".gbk"))))
    }
    mydiv
  })
  
  observeEvent(input$acc,{
    rv$project <-  paste(rv$analysis,input$acc,sep="_")
    print(rv$project)
  })
  
  observeEvent(input$uploaddep,{
    req(rv$project)
    rv[[rv$project]]$depfile <-  rv$depth()
    
    print(head(rv$depth()))
  })
  
  observeEvent(input$uploadlogic,{
    req(rv$project)
    rv[[rv$project]]$logicfile <-  rv$logic()
  })
  
  observeEvent(input$godepth,{
    req(rv$project)
    rv$depinaccs <-  c(rv$depinaccs,list(rv$project))
    print(rv$depinaccs)
  })
  
  
  observeEvent(list(input$godepth,input$gopf), {
    
    req(rv$project)
    req(rv[[rv$project]]$chrs)
    cir = list()
    
    for(nr in 1:length(rv[[rv$project]]$chrs)){
      file <- input[[paste0("file_",rv[[rv$project]]$chrs[nr])]]
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      shiny::validate(need(ext %in% c("gb","gbk"), "Please upload a file in genbank format"))
      
      if(ext %in% c("gb","gbk")) {
        gb <- genbankr::readGenBank(file$datapath)
        
        cir[[nr]] <- bind_rows(as.data.frame(sources(gb)),as.data.frame(cds(gb)),as.data.frame(otherFeatures(gb)))[,c("start","end","strand","type","label")]
        
        #  bind_rows(as.data.frame(cds(gb)),as.data.frame(otherFeatures(gb)))
        
      }
    }
    
    names(cir) <- rv[[rv$project]]$chrs
    
    print(names(cir))
    rv[[rv$project]]$circuit <- cir
    
  })
  
  
  observeEvent(input$rnadep,{
    
    rv$analysis <- "CovDep"
    
    removeModal()
  })
  
  observeEvent(input$rnapf,{
    
    rv$analysis <- "RnapFlux"
    
    removeModal()
  })
  
  observeEvent(input$rnade,{
    rv$analysis <- "Deseq"
    removeModal()
  })
  
  observeEvent(input$gopf,{
    req(rv[[rv$project]]$depfile)
    rv[[rv$project]]$params <-  c(cpratio = input$cpratio, #adjust plasmid copy number of different backbones
                                               ratio1 = input$ratio1, #1RPU =0.019 RNAP/s, conversion factor estimated by single-molecule measurement
                                               convert1 = input$intercept, #intercept of RPU conversion 
                                               convert2 = input$coefficient, #coefficient of RPU conversion
                                               gamma = input$gamma #degradation rate
    )
    
    rv[[rv$project]]$rnapflux <- foreach::foreach(a=rv[[rv$project]]$chr,.combine = "rbind") %do% {
      viztab <- rv[[rv$project]]$depfile
      viztab <- viztab[viztab$chr == a, ]
      
      print(paste(rv$project, a, length(viztab$loc)))
      newviz <- data.frame(loc=viztab$loc,adjdepth=viztab$adjdepth,chr=a)
      nummax <- max(newviz$loc)
      temp <- data.frame(loc=(1:nummax)[!1:nummax %in% viztab$loc],adjdepth=0,chr=a)
      newviz <- rbind(newviz,temp)
      newviz <- newviz[order(newviz$loc),]
      print(paste(rv$project, a, length(newviz$loc)))
      
      rnapf <- data.frame(chr=a,loc=newviz$loc,adjdepth=rnap.func(newviz$adjdepth,rv[[rv$project]]$params))
      return(rnapf)
    }
    
    print(head(rv[[rv$project]]$rnapf))
    print(tail(rv[[rv$project]]$rnapf))
    
  })
   
  output$sampselect <- renderUI({
    if(length(rv$depinaccs) > 0) {
     
      selectInput("selectsamp","Please select samples for comparison from:",
                  choices = unlist(rv$depinaccs),multiple = T)
     
    }else{
      div()
    }
  })
    
  output$sampselectchrs <- renderUI({
    req(input$selectsamp)
    if(input$compall) {
      
      choosefrom <- foreach::foreach(a=input$selectsamp,.combine = "c") %do% rv[[a]]$chrs
      div(selectInput("selectchrin", "Please select plasmids to compare:",choices = unique(choosefrom)
                      #,multiple = T
                      ))
    }else{
    if(length(input$selectsamp) ==2) {
      div(foreach::foreach(a=input$selectsamp) %do% do.call("selectInput",
                                                            list(inputId=paste("selectsampchr",a,sep = "_"),
                                                                label=paste0("Please select plasmid to compare for",a,": "),
                                                                choices=rv[[a]]$chrs)))
    
      
    }else{
      div()
    }}
  })
  
  output$selectcompele <- renderUI({
    req(input$selectsamp)
    if(input$compall) {
      
     samplewithchr <- input$selectsamp[foreach::foreach(a=input$selectsamp,.combine = "c") %do% input$selectchrin %in% rv[[a]]$chrs]
     print(samplewithchr)
     i <- match(input$selectchrin, rv[[samplewithchr[1]]]$chrs)
     features=rv[[samplewithchr[1]]]$circuit[[i]]
      
     div(selectInput("selectelein","Please select elements to compare:",
                     choices = na.omit(features$label[!features$type %in% c("rep_origin","CDS")]),
                     multiple = T)
         )
    }else{
      
      if(length(input$selectsamp) ==2) {
        
        
        i <- match(input[[paste("selectsampchr",input$selectsamp[1],sep = "_")]], rv[[input$selectsamp[1]]]$chrs)
        k <- match(input[[paste("selectsampchr",input$selectsamp[2],sep = "_")]], rv[[input$selectsamp[2]]]$chrs)
        
        features=list(rv[[input$selectsamp[1]]]$circuit[[i]],rv[[input$selectsamp[2]]]$circuit[[k]])
        
        div(foreach::foreach(a=1:2) %do% do.call("selectInput",
                                                 list(inputId=paste("selectsampele",input$selectsamp[a],sep = "_"),
                                                      label=paste0("Please select elements to compare for",input$selectsamp[a],": "),
                                                      multiple = T,
                                                      choices=na.omit(features[[a]]$label[!features[[a]]$type %in% c("rep_origin","CDS")]))
                                                 )
            )
        
      }else{
        div()
      }
    }

  })
  
  output$lincompare <- renderPlotly({
    req(input$selectsamp)
    if(input$compall) {
      samplewithchr <- input$selectsamp[foreach::foreach(a=input$selectsamp,.combine = "c") %do% input$selectchrin %in% rv[[a]]$chrs]
    #  print(samplewithchr)
      i <- match(input$selectchrin, rv[[samplewithchr[1]]]$chrs)
      features=rv[[samplewithchr[1]]]$circuit[[i]]
      
      plotlist <- foreach::foreach(a=samplewithchr) %do% rv[[a]]$depfile
      names(plotlist) <- samplewithchr
     plot.lin.comp.multi(features,plotlist,chrs = input$selectchrin, select=input$selectelein)
      
      
    }else{
      if(length(input$selectsamp) ==2) {
        i <- match(input[[paste("selectsampchr",input$selectsamp[1],sep = "_")]], rv[[input$selectsamp[1]]]$chrs)
        k <- match(input[[paste("selectsampchr",input$selectsamp[2],sep = "_")]], rv[[input$selectsamp[2]]]$chrs)
        print(rv[[input$selectsamp[1]]]$chrs[i])
        print(rv[[input$selectsamp[2]]]$chrs[k])
      plot.linear.compare(features=rv[[input$selectsamp[1]]]$circuit[[i]],
                            features2=rv[[input$selectsamp[2]]]$circuit[[k]],
                            depthfile=rv[[input$selectsamp[1]]]$depfile,
                            depthfile2=rv[[input$selectsamp[2]]]$depfile,
                            chr=rv[[input$selectsamp[1]]]$chrs[i],
                            chr2=rv[[input$selectsamp[2]]]$chrs[k],
                              select=input[[paste("selectsampele",input$selectsamp[1],sep = "_")]],
                               select2=input[[paste("selectsampele",input$selectsamp[2],sep = "_")]]
        )
      }
    }

  })
  
  
  observeEvent(input$godepth, {
    id <- paste0("dep",rv$project, input$prepend, "p")
    print(id)
    insertTab(inputId = "main",
              # target = "Wellcome",
              navbarMenu(paste0("Coverage Depth Analysis for ", rv$project),
                         
                         tabPanel("Circular Coverage Distribution", 
                                  card(card_header("Coverage Plot"),
                                       layout_sidebar(
                                         fillable = TRUE,
                                         sidebar = sidebar(uiOutput(paste("cirsidebar",rv$project,sep="_")),
                                                           uiOutput(paste("selectele",rv$project,sep="_")),
                                                           checkboxInput(paste("outvals",rv$project,sep="_"),"Check to change style")
                                         ),
                                         
                                         plotlyOutput(paste("cirdepplot",rv$project,sep="_"))
                                       ),
                                       full_screen = T),
                                  layout_column_wrap(
                                    width = 1/2,
                                    card(card_header("Element Distribution"),
                                         plotlyOutput(paste("cirdist",rv$project,sep="_"))
                                    ),
                                    card(card_header("Total Element Distribution"),
                                         plotlyOutput(paste("boxes",rv$project,sep="_")),
                                         full_screen = T))
                                  
                         ),
                         tabPanel("Linear Coverage Map", 
                                  card(card_header("Linear Coverage Depth"),
                                       layout_sidebar(
                                         fillable = TRUE,
                                         sidebar = sidebar(uiOutput(paste("linsidebar",rv$project,sep="_")),
                                                           uiOutput(paste("linselectele",rv$project,sep="_"))),
                                         plotlyOutput(paste("lincov",rv$project,sep="_"))
                                       ),full_screen = T)
                         ),
                         "------",
                         "Sample Info",
                         paste0("Depth File: ", input$uploaddep$name),
                         paste0("Circuit Files: \n"),
                         do.call("paste0",foreach::foreach(a=rv[[rv$project]]$chrs) %do% paste0(input[[paste0("file_",a)]]$name, "\n"))
              )
    )
    removeModal()
  })
  
  
  observe({
  
  lapply(rv$depinaccs, function(x) {
    output[[paste("cirsidebar",x,sep="_")]]  <- renderUI({
      req(rv[[x]]$circuit)
      
      mydiv <- div(
        selectInput(paste("sidecir",x,sep="_"),"Select plasmid:", choices = rv[[x]]$chrs,selected =rv[[x]]$chrs[1])
      )
    })
    
    output[[paste("selectele",x,sep="_")]] <- renderUI({
      req(rv[[x]]$circuit)
      req(input[[paste("sidecir",x,sep="_")]])
      
      i <- match(input[[paste("sidecir",x,sep="_")]], rv[[x]]$chrs)
      features <- rv[[x]]$circuit[[i]]
      choice <- na.omit(features$label[!features$type %in% c("rep_origin","CDS")])
      
      mydiv <- div(
        selectInput(paste("selectele",x,sep="_"),"Select additional elements:", choices = choice,
                    #  selected =choice[1],
                    multiple = T)
      )
    })
    
    output[[paste("linsidebar",x,sep="_")]]  <- renderUI({
      req(rv[[x]]$circuit)
      
      mydiv <- div(
        selectInput(paste("linsidecir",x,sep="_"),"Select plasmid:", choices = rv[[x]]$chrs,selected =rv[[x]]$chrs[1])
      )
    })
    
    output[[paste("linselectele",x,sep="_")]] <- renderUI({
      req(rv[[x]]$circuit)
      req(input[[paste("linsidecir",x,sep="_")]])
      
      i <- match(input[[paste("linsidecir",x,sep="_")]], rv[[x]]$chrs)
      features <- rv[[x]]$circuit[[i]]
      choice <- na.omit(features$label[!features$type %in% c("rep_origin","CDS")])
      
      mydiv <- div(
        selectInput(paste("linselectele",x,sep="_"),"Select additional elements:", choices = choice,
                    #  selected =choice[1],
                    multiple = T)
      )
    })
    
    output[[paste("cirdepplot",x,sep="_")]] <- renderPlotly({
      #req(input$godepth)
      req(rv[[x]]$depfile)
      req(rv[[x]]$circuit)
      req(input[[paste("sidecir",x,sep="_")]])
      
      depf <- rv[[x]]$depfile
      
      i <- match(input[[paste("sidecir",x,sep="_")]], rv[[x]]$chrs)
      print(paste("seleted",rv[[x]]$chrs[i]))
      
      if(input[[paste("outvals",x,sep="_")]]) {
        fig <- plot.circular(rv[[x]]$circuit[[i]],depf,
                                chr=rv[[x]]$chrs[i],
                                select =input[[paste("selectele",x,sep="_")]], 
                                source=paste("cir_plot",x,sep = "_"))
      }else{
        fig <- plot.circular.in(rv[[x]]$circuit[[i]],depf,
                             chr=rv[[x]]$chrs[i],
                             select =input[[paste("selectele",x,sep="_")]], 
                             source=paste("cir_plot",x,sep = "_"))
      }
      
      
      fig
      
    })
    
    output[[paste("cirdist",x,sep="_")]] <- renderPlotly({
      req(rv[[x]]$depfile)
      req(rv[[x]]$circuit)
      req(input[[paste("sidecir",x,sep="_")]])
      clickData <- event_data("plotly_click", source = paste("cir_plot",x,sep = "_"))
      
      if (is.null(clickData)) return(NULL)
      
      print(clickData)
      
      depf <- rv[[x]]$depfile
      
      i <- match(input[[paste("sidecir",x,sep="_")]], rv[[x]]$chrs)
      print(paste("seleted",rv[[x]]$chrs[i]))
      
      features <- rv[[x]]$circuit[[i]]
      if(length(features$label[features$type %in% c("rep_origin","CDS")]) >0) {
        mylabs <-  features$label[features$type %in% c("rep_origin","CDS")]
      }
      
      if(length(input[[paste("selectele",x,sep="_")]]) >0) {
        selected <- features$label[features$label %in% input[[paste("selectele",x,sep="_")]]]
        mylabs <- c(mylabs,selected)
      }
      
      
      if(clickData[["curveNumber"]] > 3) {
        k <- clickData[["curveNumber"]]
        print(paste("clicked",mylabs[k-3]))
        print(features[features$label  %in% mylabs[k-3], ])
        newfend <- features$end[features$label  %in% mylabs[k-3]]
        newfstart <- features$start[features$label  %in% mylabs[k-3]]
        vari <- extract.norm.depth(depf,rv[[x]]$chrs[i],newfstart,newfend)$adjdepth
        print(vari)
        fig1 <- plot_ly() %>% add_histogram(x=vari,name=paste0("selected=",mylabs[k-3]))
        fig2 <- plot_ly() %>% add_boxplot(x=vari,name=paste0("selected=",mylabs[k-3]))
        
        fig <- subplot(fig1,fig2,nrows = 2,shareX = T)
        
      }else{
        fig <- plot_ly() 
      }
      
      fig
    })
    
    output[[paste("boxes",x,sep="_")]] <- renderPlotly({
      req(rv[[x]]$depfile)
      req(rv[[x]]$circuit)
      req(input[[paste("sidecir",x,sep="_")]])
      
      depf <- rv[[x]]$depfile
      
      i <- match(input[[paste("sidecir",x,sep="_")]], rv[[x]]$chrs)
      print(paste("seleted",rv[[x]]$chrs[i]))
      
      fig <- plot.boxes(rv[[x]]$circuit[[i]],depf,chr=rv[[x]]$chrs[i],select =input[[paste("selectele",x,sep="_")]], source=paste("box_plot",x,sep = "_"))
      
      fig
    })
    
    output[[paste("lincov",x,sep="_")]] <- renderPlotly({
      req(rv[[x]]$depfile)
      req(rv[[x]]$circuit)
      req(input[[paste("linsidecir",x,sep="_")]])
      
      depf <- rv[[x]]$depfile
      
      i <- match(input[[paste("linsidecir",x,sep="_")]], rv[[x]]$chrs)
      print(paste("seleted",rv[[x]]$chrs[i]))
      
      fig <- plot.linear(rv[[x]]$circuit[[i]],depf,chr=rv[[x]]$chrs[i],select =input[[paste("linselectele",x,sep="_")]], source=paste("lin_plot",x,sep = "_"))
      
      fig
    })
  }) 

  })
  
  observeEvent(input$gopf,{
    req(rv$project)
    rv$pfinaccs <-  c(rv$pfinaccs,list(rv$project))
    print(rv$pfinaccs)
  })
  observeEvent(input$gopf, {
    id <- paste0("pf",rv$project, input$prepend, "p")
    print(id)
    insertTab(inputId = "main",
              # target = "Wellcome",
              navbarMenu(paste0("RNAP Flux Analysis for ", rv$project),
                         
                         tabPanel("Logic Flowchart",
                                  
                                  layout_column_wrap(
                                    width = 1/2,
                                    card(card_header("Logic Chart by GraphViz"),
                                         grVizOutput(paste("gvgraph",rv$project,sep="_")),
                                         full_screen = T),
                                    card(card_header("Sankey Plot of Genetic Circuit Logic"),
                                         plotlyOutput(paste("logicflow",rv$project,sep="_")),
                                        # 
                                         full_screen = T)
                                        
                                    ),
                                  layout_column_wrap(
                                    width = 1/2,
                                    card(card_header("Logic Elements on Plasmids"),
                                           plotlyOutput(paste("rflogic",rv$project,sep="_")),
                                         full_screen = T),
                                    card(card_header("RNAP Flux of selected elements"),
                                         layout_sidebar(
                                           sidebar = sidebar(
                                             checkboxInput(paste("cal_flux",rv$project,sep="_"),
                                                           label = "Check to display RNAP Flux",
                                                           value = F)
                                           ),
                                           fillable = TRUE,
                                           plotlyOutput(paste("rfele",rv$project,sep="_"))
                                         ),
                                         full_screen = T)
                                  )
                                  
                                  
                         ),
                         "------",
                         "Sample Info",
                         paste0("Depth File: ", input$uploaddep$name),
                         paste0("Circuit Files: \n"),
                         do.call("paste0",foreach::foreach(a=rv[[rv$project]]$chrs) %do% paste0(input[[paste0("file_",a)]]$name, "\n"))
              )
    )
    removeModal()
  })
  observe({
    
    lapply(rv$pfinaccs, function(x) {
      output[[paste("logicflow",x,sep="_")]] <- renderPlotly({
        print(rv[[x]]$logicfile)
       fig <- gcl.sankey(rv[[x]]$logicfile)
       
       fig
      })
      
      output[[paste("gvgraph",x,sep="_")]] <- renderGrViz({
        req(rv[[x]]$logicfile)
        
        d <- DiagrammeR::grViz(rv[[x]]$logicfile$Dot_Lang)
        d
      })
      
      output[[paste("rflogic",x,sep="_")]] <- renderPlotly({
        req(rv[[x]]$logicfile)
        req(rv[[x]]$circuit)
        
        pllist <- foreach::foreach(a=rv[[x]]$chr) %do% {
          return(plot.pl.reg(a,rv[[x]]$circuit[[a]],rv[[x]]$logicfile, source=paste("reg_plot",x,a,sep = "_")))
        }
        
       plot <-  subplot(pllist,nrows = length(pllist))
        
       plot
      })
      
      
      output[[paste("rfele",x,sep="_")]] <- renderPlotly({
        req(rv[[x]]$logicfile)
        req(rv[[x]]$depfile)
        req(rv[[x]]$circuit)
        req(rv[[x]]$chr)
        
        clickData <- event_data("plotly_click", source=paste("reg_plot",x,rv[[x]]$chr[[length(rv[[x]]$chr)]],sep = "_"))
        
        
        if (is.null(clickData)) return(NULL)
        
        print(clickData)
        
        num <-  clickData[["curveNumber"]]
        print(num)
        
        eledf <- foreach::foreach(a=rv[[x]]$chr,.combine = "rbind") %do% {
          if(length(rv[[x]]$logicfile$GC_Units$items) >0){
            elelist <- foreach::foreach(b=rv[[x]]$logicfile$GC_Units$items,.combine = "c") %do% b
            cir <- rv[[x]]$circuit[[a]]
            newf <- cir[na.omit(match(elelist,cir$label)),]
            if(length(newf$label)>0) {
              newf <- newf[newf$type %in% names(ele.func.suite),]
              newf <- rbind(newf,data.frame(start=0,end=0,strand="line",type="line",label="line"))
              newf$chr <- a
              return(newf)
            }
            
          }
        }
        
        print(eledf[num,])
        plot <- plot_ly()
        
        if(clickData[["curveNumber"]] >0) {
          
          if(input[[paste("cal_flux",x,sep="_")]] == F) {
            plot <- plot.with.ele(rv[[x]]$depfile,eledf[num,]$chr,eledf[num,])
            print("Cov!")
          }else{
            plot <- plot.with.ele(rv[[x]]$rnapflux,eledf[num,]$chr,eledf[num,])
            print("RNAP!")
          }
          
        }
        
        plot
      })
      
    })
  })
  
  dataModal <- function(failed = FALSE) {
    modalDialog(size="l",title="Genetic Circuit Analysis Options",
                helpText("Please choose one of the analyses below:"
                ),
                layout_column_wrap(
                  width = 1/3,
                  fill = FALSE,
                  value_box(
                    "RNA-seq Data Coverage Depth Analysis",
                    #  uiOutput("total_tippers", container = h2),
                    #  textInput("depname","Plsease enter the name of analysis"),
                    actionBttn("rnadep","GO",color = "success"),
                    showcase = bsicons::bs_icon("bar-chart-steps")
                  ),
                  value_box(
                    "RNA Polymerase Flux Analysis",
                    theme_color = "warning",
                    actionBttn("rnapf","GO",color = "success"),
                    showcase = bsicons::bs_icon("chevron-double-right")
                  ),
                  value_box(
                    "RNA-seq Data Differentiation Expression Analysis",
                    theme_color = "info",
                    #  uiOutput("total_tippers", container = h2),
                    #     textInput("dename","Plsease enter the name of analysis"),
                    actionBttn("rnade","GO",color = "success"),
                    showcase = bsicons::bs_icon("activity")
                  )
                ),
                footer = tagList(
                  modalButton("Cancel")
                  #   actionButton("gomcmc", "Confirm")
                )
    )
  }
  
  depthmodal <- function(failed = FALSE) {
    modalDialog(size="m",title="Genetic Circuit Coverage Depth Analysis",
                helpText("Please upload files as required below:"
                ),
                
                textInput("acc", "Sample Accession:"),
                #  fileInput("uploadc", NULL, buttonLabel = "Upload plasmids...", multiple = TRUE,accept = c(".gb")),
                fileInput("uploaddep", NULL, buttonLabel = "Upload one coverage depth file...", multiple = FALSE,accept = c(".depth")),
                numericInputIcon("adjfactor","Adjust Factor Value:",value = 1e6, min=1, max = 1e9),
                uiOutput("cirfiles"),
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("godepth", "Confirm")
                )
    )
  }
  
  pfmodal <- function(failed = FALSE) {
    modalDialog(size="m",title="Genetic Circuit RNA Polymerase Flux Analysis",
                helpText("Please upload files as required below:"
                ),
                
                textInput("acc", "Sample Accession:"),
                #  fileInput("uploadc", NULL, buttonLabel = "Upload plasmids...", multiple = TRUE,accept = c(".gb")),
                fileInput("uploaddep", NULL, buttonLabel = "Upload one coverage depth file...", multiple = FALSE,accept = c(".depth")),
                fileInput("uploadlogic", NULL, buttonLabel = "Upload circuit logic file...", multiple = FALSE,accept = c(".json")),
                uiOutput("cirfiles"),
                numericInputIcon("adjfactor","Adjust Factor Value:",value = 1e6, min=1, max = 1e9),
                dropdownButton(
                  inputId = "mydropdown",
                  label = "Set RNAP Flux Paramters",
                  icon = icon("sliders"),
                  circle = T,
                  status = "primary",
                
                numericInputIcon("cpratio","Adjust Plasmid Copy Number:",value=1, min=1, max = 1e5),
                numericInputIcon("ratio1","RPU-RNAP Conversion Factor:",value=0.019, min=0, max = 10),
                numericInputIcon("intercept","RPU-RNAP Conversion Intercept:",value=5.05e-5, min=0, max = 100),
                numericInputIcon("coefficient","RPU-RNAP Conversion Coefficient:",value=1.64, min=0, max = 100),
                numericInputIcon("gamma","RNA Degradation Rate:",value=0.0067, min=0, max = 10)
                ),
               
                footer = tagList(
                  modalButton("Cancel"),
                  actionButton("gopf", "Confirm")
                )
    )
  }
  

  
  observeEvent(input$rnadep, {
    showModal(depthmodal())
  })
  
  observeEvent(input$rnapf, {
    showModal(pfmodal())
  })
  
  # Show modal when button is clicked.
  observeEvent(input$start, {
    showModal(dataModal())
  })
  
  
}

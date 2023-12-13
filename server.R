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
library(MASS)
library(shinybusy)
library(genbankr)
library(htmltools)
library(DT)


source("global.R")

# Define server logic required to draw a histogram
function(input, output, session) {
  
  rv <- reactiveValues(inaccs=list())
  
  rv$depth <- reactive({
    file <- input$uploaddep
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext %in% c("depth"), "Please upload a file in the right format"))
    
    readdepfiles.adj.func(file=file$datapath,acc=input$acc,adjfactor=input$adjfactor)
    
  })
  
  
  output$cirfiles <- renderUI({
    req(input$acc)
    if(is.null(rv$depth())) {
      mydiv <- div()
    }else{
      rv[[input$acc]]$chrs <- unique(rv$depth()$chr)
      # print(rv[[input$acc]]$chrs)
      mydiv <- div(foreach::foreach(a=rv[[input$acc]]$chrs) %do% do.call("fileInput", list(inputId=paste0("file_",a), 
                                                                              label=paste0("Plasmid file for ", a,": "), 
                                                                              multiple = FALSE,
                                                                              accept = c(".gb",".gbk"))))
    }
    mydiv
  })
  
  observeEvent(input$uploaddep,{
    req(input$acc)
    rv[[input$acc]]$depfile <-  rv$depth()
  })
  
  observeEvent(input$uploaddep,{
    req(input$acc)
    rv$inaccs <-  c(rv$inaccs,list(input$acc))
    print(rv$inaccs)
  })
  
  observeEvent(input$godepth, {
    
    cir = list()
    
    for(nr in 1:length(rv[[input$acc]]$chrs)){
      file <- input[[paste0("file_",rv[[input$acc]]$chrs[nr])]]
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext %in% c("gb","gbk"), "Please upload a file in genbank format"))
      
      if(ext %in% c("gb","gbk")) {
        gb <- genbankr::readGenBank(file$datapath)
        
        cir[[nr]] <- bind_rows(as.data.frame(sources(gb)),as.data.frame(cds(gb)),as.data.frame(otherFeatures(gb)))[,c("start","end","strand","type","label")]
        
        #  bind_rows(as.data.frame(cds(gb)),as.data.frame(otherFeatures(gb)))
        
      }
    }
    
    names(cir) <- rv[[input$acc]]$chrs
    
    print(names(cir))
    rv[[input$acc]]$circuit <- cir
    
  })
  
  observeEvent(input$rnadep,{
    
    rv$analysis <- "covdep"
    
    removeModal()
  })
  
  observeEvent(input$rnade,{
    rv$analysis <- "deseq"
    removeModal()
  })
  
  output$sampselect <- renderUI({
    if(length(rv$inaccs) > 0) {
     
      selectInput("selectsamp","Please select samples for comparison from:",
                  choices = unlist(rv$inaccs),multiple = T)
     
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
    id <- paste0("dep",input$acc, input$prepend, "p")
    print(id)
    insertTab(inputId = "main",
              # target = "Wellcome",
              navbarMenu(paste0("Coverage Depth Analysis for ", input$acc),
                         
                         tabPanel("Circular Coverage Distribution", 
                                  card(card_header("Coverage Plot"),
                                       layout_sidebar(
                                         fillable = TRUE,
                                         sidebar = sidebar(uiOutput(paste("cirsidebar",input$acc,sep="_")),
                                                           uiOutput(paste("selectele",input$acc,sep="_")),
                                                           checkboxInput(paste("outvals",input$acc,sep="_"),"Check to change style")
                                         ),
                                         
                                         plotlyOutput(paste("cirdepplot",input$acc,sep="_"))
                                       ),
                                       full_screen = T),
                                  layout_column_wrap(
                                    width = 1/2,
                                    card(card_header("Element Distribution"),
                                         plotlyOutput(paste("cirdist",input$acc,sep="_"))
                                    ),
                                    card(card_header("Total Element Distribution"),
                                         plotlyOutput(paste("boxes",input$acc,sep="_")),
                                         full_screen = T))
                                  
                         ),
                         tabPanel("Linear Coverage Map", 
                                  card(card_header("Linear Coverage Depth"),
                                       layout_sidebar(
                                         fillable = TRUE,
                                         sidebar = sidebar(uiOutput(paste("linsidebar",input$acc,sep="_")),
                                                           uiOutput(paste("linselectele",input$acc,sep="_"))),
                                         plotlyOutput(paste("lincov",input$acc,sep="_"))
                                       ),full_screen = T)
                         ),
                         "------",
                         "Sample Info",
                         paste0("Depth File: ", input$uploaddep$name),
                         paste0("Circuit Files: \n"),
                         do.call("paste0",foreach::foreach(a=rv[[input$acc]]$chrs) %do% paste0(input[[paste0("file_",a)]]$name, "\n"))
              )
    )
    removeModal()
  })
  
  observe({
  
  lapply(rv$inaccs, function(x) {
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
                
                textInput("acc", "Sample Accesstion:"),
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
  

  
  observeEvent(input$rnadep, {
    showModal(depthmodal())
  })
  
  # Show modal when button is clicked.
  observeEvent(input$start, {
    showModal(dataModal())
  })
  
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white',
         xlab = 'Waiting time to next eruption (in mins)',
         main = 'Histogram of waiting times')
    
  })
  
}

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(bslib)
library(shinyBS)
library(shinycssloaders)
library(shinyjs)
library(shinyWidgets)
library(ggpubr)
library(survminer)
library(tidyverse)
library(rlang)
library(plyr)
library(tippy)
library(plotly)
library(kableExtra)
library(reactable)

# Define UI for application that draws a histogram
page_navbar(
  theme = bs_theme(version = 5, bootswatch = "lumen"),
  id="main",
  
  # Application title
  title="Genetic Circuit Debugging Toolkits",
  fillable="Dashboard",
  fillable_mobile=T,
  collapsible = T,
  nav_panel("Home",
            layout_column_wrap(
              width = 1,
              fill = FALSE,
              h1(),
              h1("Welcome to Genetic Circuit Debugging Toolkits!"),
              h4("Please click on the button and choose you analysis options to start your analysis."),
              layout_column_wrap(width=1/3,actionBttn("start", "New Project"))
            )
  ),
  nav_spacer(),
  
  nav_panel("Compare Samples",
            
            card(card_header("Compare uploaded samples"),
                 layout_sidebar(
                   fillable = TRUE,
                   sidebar = sidebar(uiOutput("selectcompele")),
                   plotlyOutput('lincompare')),
                 full_screen = T)
            
  ),
  sidebar = sidebar(position = "right",
                    title = "Samples",
                    open = "open",
                    uiOutput("sampselect"),
                    checkboxInput("compall","Compare all selected samples"),
                    uiOutput("sampselectchrs")
                    )
                    
  
)

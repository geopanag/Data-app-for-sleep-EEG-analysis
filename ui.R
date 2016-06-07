options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)

shinyUI(fluidPage(
    titlePanel("EEG Sleep Study Spectral Analysis"),
    
    sidebarLayout(
        sidebarPanel(
            width = 3,
            h4("Upload Subject's Recording"),
            fileInput(
                'file', '',
                accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
            ),
            
            
            h4("Define Analysis Parameters"),
            selectInput(
                "freq", "Choose frequency bands",
                choices = c(
                    "Delta (0.5 - 3.5 Hz)",
                    "Theta (4 - 8 Hz)",
                    "Alpha (8.5 - 12 Hz)",
                    "Sigma (12.5 - 16 Hz)",
                    "Beta (16.5 - 30 Hz)",
                    "Gamma (30.5 - 60 Hz)"
                )
            ),
            checkboxInput("ica", "Run Independent Component Analysis", value = FALSE),
            selectInput("channel", "Select channel to plot",
                        choices = c(""))
            
        ),
        
        
        mainPanel(
            column(
                6,
                h3("Mean Spectral Density Before Sleep"),
                webGLOutput("ScalpPlotBefore")
            ),
            column(
                6,
                h3("Mean Spectral Density After Sleep"),
                webGLOutput("ScalpPlotAfter")
            ),
            fluidRow(
                height = 3,
                h3("Filtered Signal by Epoch"),
                tabsetPanel(
                    tabPanel("E1",plotOutput("ChannelE1")),
                    tabPanel("E2",plotOutput("ChannelE2")),
                    tabPanel("E3",plotOutput("ChannelE3")),
                    tabPanel("E4",plotOutput("ChannelE4")),
                    tabPanel("E5",plotOutput("ChannelE5")),
                    tabPanel("E6",plotOutput("ChannelE6")),
                    tabPanel("E7",plotOutput("ChannelE7")),
                    tabPanel("E8",plotOutput("ChannelE8")),
                    tabPanel("E9",plotOutput("ChannelE9")),
                    tabPanel("E10",plotOutput("ChannelE10")),
                    tabPanel("E11",plotOutput("ChannelE11")),
                    tabPanel("E12",plotOutput("ChannelE12"))
                )
            ),
            fluidRow(
                h3("Mean Spectral Density by Epoch-Channel"),
                tableOutput('fullSpec')
            )
        )
        
    )
))

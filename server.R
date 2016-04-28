require(rCharts)

shinyServer(function(input, output, session){
    observe({
        if (input$start) {
            startstop <<- TRUE
            go        <<- TRUE
        }
    })
  
    observe({
        if (input$stp) {
            startstop <<- FALSE
            go        <<- FALSE
        }
    })
  
    observe({
        if (input$add > 0) {
            add.points <<- isolate(input$num)
            go         <<- TRUE
            startstop  <<- TRUE
        }
    })
  
    observe({
        if (input$go > 0) {
            newSession <<- TRUE
        }
    })
    
    observe({
        input$go
        if (newSession) {           
            family.real     <<- isolate({input$input_type})
            family.pa       <<- isolate({input$dynamic})
            family.boot     <<- isolate({input$dynamic_boot}) 
            pattern.control <<- isolate({input$pattern_control})
            InitNewSession(family.real, family.pa, family.boot, pattern.control)
        }
    })
  
    output$Splash = renderUI({ 
        input$go
        if (input$go < 1) {   
            HTML("
              <br>
              <br>
              <br>
              <p align=\"right\"><font size=5><i>All models are wrong, but some are useful.</i></font> </p>      
              <br><p align=\"right\"><i>George E.P.Box</i>
              <br>
              <br>
              <br>
              <br>
              <br>
              <br>
              <br><p align=\"center\">Please select parameters and press &quot;Start demo&quot;.")  
        }
    })
  
    observe({
        input$realTime
        RealTime <<- TRUE
    })
  
    observe({
        if (input$timeAlign > 0) {
            RealTime <<- FALSE
        }
    })
  
    observe({ 
        input$timeAlign
        input$realTime
        input$go
        input$start
        input$add
        input$resetSession
        if (go == TRUE) {
            if (!is.na(add.points)) {
              if (add.points <= 1) {
                add.points <<- NA
                startstop  <<- FALSE
              } else {
                add.points <<- add.points - 1
              }
            }
            if (startstop) {  
                invalidateLater(as.integer(700 / input$slider1 + 300), session)
            }
            i     <<- i + 1
            i.col <<- findInterval(i, as.integer(cp.info['start', ]))
            observed.data[1 : (nrow(observed.data) - 1), ] <<- observed.data[2 : nrow(observed.data), ]
            observed.data$value[nrow(observed.data)]       <<- GenerateDataOnline(family.real, 
                                                                                  as.numeric(cp.info['MeanState', i.col]), 
                                                                                  sd = as.numeric(cp.info['VarState', i.col])) 
            observed.data$MeanValue[nrow(observed.data)]   <<- as.numeric(cp.info['MeanState', i.col])
            observed.data$color[nrow(observed.data)]       <<- cp.info['colour', i.col]
            #change point flag
            if (sum(which(i == as.integer(cp.info['start', ]))) > 0) {
                observed.data$cpFlag[nrow(observed.data)] <<- 1
            } else {
                observed.data$cpFlag[nrow(observed.data)] <<- 0
            }
            LRTestStat[i, 1:3]   <<- sqrt( 2 * abs(lrTestStat(H, observed.data$value[(nrow(observed.data) - 2 * max(H) + 1) : nrow(observed.data)], 
                                                           current.point = i,family.pa)))
            LRTestStat$cpFlag[i] <<- observed.data$cpFlag[nrow(observed.data)]
            if (pattern.control == 'Least squares') {
                LRTestStatTriang$cpFlag[i] <<- observed.data$cpFlag[nrow(observed.data)]
                for (l in seq_along(H)) {
                    dat      = LRTestStat[(i - H[l] + 1): i, l]
                    triangle = lsfit((1:H[l]), dat)
                    if (triangle$coef[2] < 0) {
                        triangle = rep(triangle$coef[1], H[l]) 
                    } else {
                        triangle = triangle$coef[2]*(1:H[l]) + triangle$coef[1]  
                    }     
                    LRTestStatTriang[i, l] <<- sum(dat * triangle / sum(triangle))
                }
            } else if (pattern.control == 'Triangle') {
                LRTestStatTriang$cpFlag[i] <<- observed.data$cpFlag[nrow(observed.data)]
                for (l in seq_along(H)) {
                    dat                    = LRTestStat[(i - H[l] + 1): i, l]
                    triangle               = 2 * (1 : H[l]) / (H[l] * (H[l] + 1))
                    LRTestStatTriang[i, l] <<- sum(dat * triangle)
                }
            }
            if (i >= N + 2) {
                x    <<- (i - theta + 1) : i
                xMAX <<- i
                cps  <<- sum(observed.data[(nrow(observed.data) - length(x) + 1) : nrow(observed.data),]$cpFlag) #number of change points inside the window
                
                output$realTimeData = renderChart({ 
                    p1 = PlotData(input.data = observed.data[(nrow(observed.data) - length(x) + 1) : nrow(observed.data),], 
                                  x = x, cps = cps)
                    p1$plotOptions(scatter = list(marker = list(symbol = 'circle')), 
                                  line = list(marker = list(enabled = FALSE), lineWidth = 1))
                    p1$set(dom = 'realTimeData')
                    return(p1)
                })
                
                output$realTimeStat1 = renderChart({
                    input$realTime
                    input$timeAlign
                    if (RealTime == TRUE) {
                        shift = FALSE
                    }
                    if (RealTime == FALSE) {
                        shift = TRUE
                    }
                    if (pattern.control == 'Off') {
                        z = LRTestStat[x ,c(1,4)]  
                    } else {
                        z = LRTestStatTriang[x ,c(1,4)]
                    }
                    p1 = PlotStatistics(input.data = z, x = x, h = H[1], boot.value = BootStrapLine[1], cps = cps,
                                         pattern = pattern.control, shift = shift)
                    p1$plotOptions(scatter = list(marker = list(symbol = 'circle')), 
                                   line = list(marker = list(enabled = FALSE), lineWidth = 2))
                    p1$set(dom = 'realTimeStat1')
                    return(p1)
                }) 
                
                output$realTimeStat2 <- renderChart({
                    input$realTime
                    input$timeAlign
                    if (RealTime == TRUE) {
                        shift = FALSE
                    }
                    if (RealTime == FALSE) {
                        shift = TRUE
                    }
                    if (pattern.control == 'Off') {
                        z = LRTestStat[x, c(2,4)]  
                    } else {
                        z = LRTestStatTriang[x, c(2,4)]
                    }
                    p1 = PlotStatistics(input.data = z, x = x, h = H[2], boot.value = BootStrapLine[2], 
                                         cps = cps, pattern = pattern.control, shift = shift)
                    p1$plotOptions(scatter = list(marker = list(symbol = 'circle')), 
                                   line = list(marker = list(enabled = FALSE), lineWidth = 2))
                    p1$set(dom = 'realTimeStat2')
                    return(p1)
                }) 
                
                output$realTimeStat3 <- renderChart({
                    input$realTime
                    input$timeAlign
                    if (RealTime == TRUE) {
                      shift = FALSE
                    }
                    if (RealTime == FALSE) {
                      shift = TRUE
                    }
                    if (pattern.control == 'Off') {
                      z = LRTestStat[x, c(3,4)]  
                    } else {
                      z = LRTestStatTriang[x, c(3,4)]
                    }
                    p1 = PlotStatistics(input.data = z, x = x, h = H[3], boot.value = BootStrapLine[3], 
                                         cps = cps, pattern = pattern.control, shift = shift)
                    p1$plotOptions(scatter = list(marker = list(symbol = 'circle')), 
                                   line = list(marker = list(enabled = FALSE), lineWidth = 2))
                    p1$set(dom = 'realTimeStat3')
                    return(p1)
                })
            }   
        }  
    })
  
    output$ui = renderUI({
        if (is.null(input$input_type))
            return()
        switch(input$input_type,
            "Gaussian" = selectInput("dynamic", h5("Parametric model:"),
                                    choices = c("Gaussian" = "Gaussian"),
                                    selected = "Gaussian"
            ),
            "Poisson" = selectInput("dynamic", h5("Parametric model:"),
                                   choices = c("Gaussian" = "Gaussian",
                                               "Poisson" = "Poisson"),
                                   selected = "Poisson"
            ),
            "Bernoulli" = selectInput("dynamic", h5("Parametric model:"),
                                     choices = c("Bernoulli" = "Bernoulli"
                                     ),
                                     selected = "Bernoulli"
            )
        )
    })
  
    output$ui_boot = renderUI({
        if (is.null(input$dynamic))
        return()
        switch(input$dynamic,
            "Gaussian" = selectInput("dynamic_boot", h5("Bootstrap:"),
                                    choices = c("Gaussian" = "Gaussian",
                                                "Poisson" = "Poisson"),
                                    selected = "Gaussian"
            ),
            "Poisson" = selectInput("dynamic_boot", h5("Bootstrap:"),
                                   choices = c("Poisson" = "Poisson"),
                                   selected = "Poisson"
            ),
            "Bernoulli" = selectInput("dynamic_boot", h5("Bootstrap:"),
                                     choices = c("Gaussian" = "Gaussian"),
                                     selected = "Gaussian"
            )
        )
    })
})

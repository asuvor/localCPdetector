require(rCharts)
require(Matrix)
require(markovchain)

plotitnow       = FALSE
c.shifted       = 2
newSession      <<- FALSE
list.of.colours = list('darkcyan'      = '#008B8B',
                       'bluebell'      = '#A2A2D0',
                       'dartmouthgrey' = '#00703C',
                       'chamoisee'     = '#A0785A', 
                       'fuchsiarose'   = '#C74375')
H                = c(15, 30, 40) #set of scales
N.cps            = 100 #the total number of change points
theta            = 120 #the width of main windows
M                = 60
N                = 400 #size of the bootstrap sample
go               = FALSE
startstop        = FALSE
add.points       = NA
dummy            = 10000
i                = 1

InitializeDataArray = function (n, family, info) {
    out           = as.data.frame(matrix(NA, n, 6))
    colnames(out) = c('value', 'color', 'cpFlag', 'MeanValue', 'typeScatter', 'typeLine')
    for (l in 1 : n) {
        clmn                  = findInterval(l, as.integer(info['start', ]))
        out[l, 'value']       = GenerateDataOnline(family, as.numeric(info['MeanState', clmn]), 
                                            sd = as.integer(info['VarState', clmn]))
        out[l, 'typeScatter'] = 'scatter'
        out[l, 'typeLine']    = 'line'
        out[l, 'MeanValue']   = as.integer(info['MeanState', clmn])
        out[l, 'color']       = info['colour', clmn]
        if (sum(which(l == as.numeric(info['start', ]))) > 0) {
            out[l, 'cpFlag'] = 1
        } else {
            out[l, 'cpFlag'] = 0
        }
    }
    return(out)
}

InitNewSession = function(family.real, family.pa, family.boot, pattern.control) {
    cp.info           <<- MCDataGeneration(family.real, palette = list.of.colours, N.cps) #Generate a markov chain
    observed.data     <<- InitializeDataArray(N, family.real, cp.info) #Observed data
    LRTestStat        <<- as.data.frame(matrix(NA, dummy, length(H) + 1)) #Generalized LRT
    names(LRTestStat) <<- c('h1', 'h2', 'h3', 'cpFlag')
    BootStrapLine     <<- rep(0, length(H)) #Bootstrap treshold
    if (pattern.control != 'Off') {
        LRTestStatTriang        <<- as.data.frame(matrix(NA, dummy, length(H) + 1)) #Convolution with pattern
        names(LRTestStatTriang) <<- c('h1', 'h2', 'h3', 'cpFlag')    
    }
    GenerateBootstrapLine(start = 1, n = N, family.real, family.pa, family.boot, 
                          data = observed.data, scales = H, pattern = pattern.control, M, session)
    i          <<- N #The current step of the procedure
    startstop  <<- TRUE
    go         <<- TRUE
    newSession <<- FALSE
    return(NULL)
}

hplot1 = function(..., radius = 3, title = NULL, subtitle = NULL, group.na = NULL) {
    rChart = Highcharts$new()
    d      = getLayer(...)
    data   = data.frame(x = d$data[[d$x]], y = d$data[[d$y]])
    if (!is.null(d$group)) {
        data$group = as.character(d$data[[d$group]])
        if (!is.null(group.na)) {
            data$group[is.na(data$group)] = group.na
        }
    }
    if (!is.null(d$size)) 
        data$size = d$data[[d$size]]
    nrows = nrow(data)
    data = na.omit(data)
    if (nrows != nrow(data)) 
      #warning("Observations with NA has been removed")
        data = data[order(data$x, data$y), ]
    if ("bubble" %in% d$type && is.null(data$size)) 
        stop("'size' is missing")
    if (!is.null(d$group)) {
        groups = sort(unique(data$group))
        types  = rep(d$type, length(groups))
        colrs = rep(d$colrs, length(groups))
        plyr::ddply(data, .(group), function(x) {
            g = unique(x$group)
            i = which(groups == g)
            x$group = NULL
            rChart$series(data = toJSONArray2(x, json = F, names = F), color = colrs[[i]], 
                        name = g, type = types[[i]], marker = list(radius = radius),
                        animation = d$animation, enableMouseTracking = FALSE, 
                        line = list(marker = FALSE))
            return(NULL)
      })
    } else {
        rChart$series(data = toJSONArray2(data, json = F, names = F), 
                      type = d$type[[1]], marker = list(radius = radius))
    }
    if (is.categorical(data$x)) {
        rChart$xAxis(title = list(text = ''), categories = unique(as.character(data$x)), 
                   replace = T)
    } else {
        rChart$xAxis(title = list(text = ''), replace = T)
    }
    if (is.categorical(data$y)) {
        rChart$yAxis(title = list(text = ''), categories = unique(as.character(data$y)), 
                   replace = T)
    } else {
        rChart$yAxis(title = list(text = ''), replace = T)
    }
    rChart$title(text = title, replace = T)
    rChart$subtitle(text = subtitle, replace = T)
    rChart$legend(enabled = FALSE)
    return(rChart$copy())
}

environment(hplot1) = environment(hPlot)

PlotData = function(input.data, x, cps){
    rownames(input.data) = NULL
    MAX                  = max(input.data$value)
    MIN                  = min(input.data$value)
    if (family.real == 'Bernoulli') {
        delta = 0.5
    } else {
        delta = 15
    }
    if (cps > 0) {
        x.cps = x[which(input.data$cpFlag > 0)] #find cp coordiantes
        if (x.cps[1] == x[1]) {
            number.of.segments = length(x.cps)
        } else {
            number.of.segments = length(x.cps) + 1
        }
        group = findInterval(x, x.cps) + 1
    } else {
        group              = rep(1, length(x))
        number.of.segments = 1
    }
    y         = c(input.data$value, input.data$MeanValue)
    groups    = c(group, letters[group]) 
    data.plot = data.frame(x = c(x, x), 
                           y = y,
                           group = groups)
    ind       = match(unique(group), group)
    clrs      = c(input.data$color[ind], input.data$color[ind])
    types     = c(rep('scatter', number.of.segments), rep('line', number.of.segments))
    p1        = hplot1(x = 'x', y = 'y', colrs = clrs, data = data.plot, 
                      type = types, group = "group" , size = 2, animation = FALSE)
    p1$xAxis(min = x[1], max = x[1] + theta)
    p1$yAxis(max = MAX + delta, min = MIN - delta, startOnTick = FALSE, endOnTick = FALSE,
             title = list(text = family.real))
    p1$chart(height = 140, width = 650)
    return(p1)
}

PlotStatistics = function (input.data, x, h, boot.value, cps, pattern, shift = FALSE) {
    rownames(input.data) = NULL
    yMAX                 = max(input.data[, 1])
    yMIN                 = min(input.data[, 1])
    boot.line.x          = c(x[1], x[length(x)] + 2*h)
    boot.line.y          = c(boot.value, boot.value)
    boot.line.group      = c('boot.line', 'boot.line')
    data.group           = rep('data', length(x))
    if (cps == 0) {
        group = rep(1, length(x))
    } else {
      x.cps    = x[which(input.data$cpFlag > 0)]
      cp.x     = rep(x.cps, each = 2)
      cp.y     = rep(c(yMIN - 15, yMAX + 15), length(x.cps))
      cp.group = rep((1 : cps), each = 2)
    } 
    if (shift == TRUE){
        if (pattern == 'Triangle'){
          delta = ceiling(3 * h / 2)
        } else {
          delta = h
        }
        input.data = data.frame(y = c(input.data[delta : length(x), 1], rep(NA, (delta - 1))),
                               cpFlag = c(input.data$cpFlag[delta : length(x)], rep(0, (delta - 1))))
    }
    if (cps == 0) {
        data.plot = data.frame(x = c(x, boot.line.x),
                               y = c(input.data[, 1], boot.line.y),
                               group = c(data.group, boot.line.group))
    } else {
        data.plot = data.frame(x = c(x, boot.line.x, cp.x),
                               y = c(input.data[, 1], boot.line.y, cp.y),
                               group = c(data.group, boot.line.group, cp.group)) 
    }
    clrs  = c(rep('black', cps), '#76FF7A', '#126180')
    types = c(rep('line', cps), 'line', 'scatter')
    p1    = hplot1(x = 'x', y = 'y', colrs = clrs, data = data.plot, 
                 type = types, group = 'group', size = 1.5, animation = FALSE)
    p1$yAxis(max = yMAX + 15, min = yMIN - 15, startOnTick = FALSE, endOnTick = FALSE,
             title = list(text = paste('T(', 2*h, ', t)', sep = '')))
    p1$xAxis(min = x[1], max = x[1] + theta)
    p1$chart(height = 140, width = 650)
    return(p1) 
}

#Generate random data
GenerateDataOnline = function (family = c('Poisson', 'Gaussian', 'Bernoulli'), m, sd) {
    if (family == 'Poisson'){
      out = rpois(1, m)
    }else if (family == 'Gaussian'){
      out = rnorm(1, m, sd)
    }else if (family == 'Bernoulli'){
      out = rbinom(1, 1, m)  
    }  
    return(out)
}

#generate MC
MCDataGeneration = function (true.family, palette, N.cps) {
    list.of.params = list('Poisson'   = list('a' = 10, 'b' = 12, 'c' = 15, 'd' = 17, 'e' =  21),
                           'Gaussian'  = list('a' = c(10, 8), 'b' = c(12, 8.5), 'c' = c(18, 9),
                                              'd' = c(23, 8), 'e' =  c(27, 10)),
                           'Bernoulli' = list('a' = 0.1, 'b' = 0.3, 'c' = 0.5, 
                                              'd' = 0.75, 'e' =  0.8),                         
                           'colours'   = list('a' = palette[['darkcyan']], 
                                              'b' = palette[['bluebell']], 
                                              'c' = palette[['dartmouthgrey']],
                                              'd' = palette[['chamoisee']],
                                              'e' = palette[['fuchsiarose']]))
    
    #generate seq of states
    statesNames = c ('a','b', 'c', 'd', 'e')
    mcA         = new("markovchain", 
                       transitionMatrix = matrix(c(0   , 0.25, 0.25, 0.25, 0.25,
                                                   0.25, 0   , 0.25, 0.25, 0.25,
                                                   0.25, 0.25, 0   , 0.25, 0.25,
                                                   0.25, 0.25, 0.25, 0   , 0.25,
                                                   0.25, 0.25, 0.25, 0.25, 0), 
                                                 byrow = TRUE, nrow = 5, dimnames = list(statesNames, statesNames)))
    list.of.states        = rmarkovchain(n = N.cps, mcA)
    list.of.lengths       = rpois(N.cps, 80)
    #one change point inside the training set
    list.of.lengths[1:2]  = c(N / 2, N / 2) 
    #respectively small one
    list.of.states[1 : 2] = c('a', 'b')
    
    data.colours          = rmarkovchain(n = N.cps, mcA)
    out                   = matrix(NA, 6, N.cps)
    rownames(out)         = c('MeanState','VarState','length', 'start', 'end', 'colour')
    for (i in (1 : N.cps)) {
        out['MeanState', i]  = list.of.params[[true.family]][[list.of.states[i]]][1] 
        out['length', i]     = list.of.lengths[[i]]  
        out['colour', i]     = list.of.params[['colours' ]][[data.colours[i]]]
        if (true.family == 'Gaussian'){
            out['VarState', i] = list.of.params[[true.family]][[list.of.states[i]]][2] 
        }
        if (i == 1) {
            out['start', i] = 1
            out['end', i]   = as.integer(list.of.lengths[[i]])
        } else {
            out['start', i] = as.integer(out['end', (i-1)]) + 1
            out['end', i]   = as.integer(out['start', i]) + as.integer(out['length', i]) - 1  
        }
    }
    return(out)
}

#likelihood
Likelihood = function (data, family = c('Poisson', 'Gaussian', 'Bernoulli')) {
    theta = mean(data)
    if (family == 'Poisson'){
        if (theta == 0) {
            theta = 0.001
        }
        out = log(theta) * sum(data) - length(data) * theta
    }else if (family == 'Gaussian') {
        out = -.5 * sum((data - theta) ^ 2)
    }else if (family == 'Bernoulli') {
        if (theta == 0) {
            theta = 0.001
        }else if(theta == 1){
            theta = 0.999
        }
        out = sum(data) * (log(theta) - log(1 - theta)) + length(data) * log(1 - theta)
    }
    return(out)
}

#Likelihood-ratio test
LRTest = function (data.left, data.right, family = c('Poisson', 'Gaussian', 'Bernoulli')) {
    LR.left  = Likelihood(data.left, family)
    LR.right = Likelihood(data.right, family)
    LR       = Likelihood(c(data.left, data.right), family)
    LRT      = LR.left + LR.right - LR
    return(LRT)
}

#LRTStatistics 
lrTestStat = function(H, data, current.point, family.pa){
    lrtest.stat = rep(0, length(H))
    for (h in seq_along(H)){
        if (current.point >= 2 * H[h]){
            lrtest.stat[h] = LRTest(data[(length(data) - 2 * H[h] + 1) : (length(data) - H[h])], 
                                   data[(length(data) - H[h] + 1) : length(data)], family.pa)
        }
    }
    return(lrtest.stat)
}

GenerateBootstrapLine = function(start, n, family.real, family.pa, family.boot, data, scales, pattern, M, session){
  #i - current iteration
  #N num of observations for bootstrap
    if (start == 1) {
        for (j in (max(scales) * 2 + 1 ):n) {
            LRTestStat[j, 1:3] <<- sqrt(2 * abs(lrTestStat(scales, data$value[1:j], current.point = j, family.pa)  ))
            LRTestStat$cpFlag[j] <<- data$cpFlag[j]
            if (pattern == 'Least squares') {
                if (j >= 3 * max(scales) + 1) {
                    for (l in seq_along(scales)) {
                        dat      = LRTestStat[(j - scales[l] + 1): j, l]
                        triangle = lsfit((1:scales[l]), dat)
                        if (triangle$coef[2] < 0) {
                            triangle = rep(triangle$coef[1], scales[l]) 
                        } else {
                            triangle = triangle$coef[2]*(1:scales[l]) + triangle$coef[1]  
                        }
                        LRTestStatTriang[j, l] <<- sum(dat * triangle / sum(triangle))
                    }
                    LRTestStatTriang$cpFlag[j] <<- data$cpFlag[j]
                }          
            } else if (pattern == 'Triangle') {
                if (j >= 3 * max(scales) + 1) {
                    for (l in seq_along(scales)) {
                        dat                    = LRTestStat[(j - scales[l] + 1): j, l]
                        triangle               = 2 * (1 : H[l]) / (H[l] * (H[l] + 1))
                        LRTestStatTriang[j, l] <<- sum(dat * triangle)
                    }
                    LRTestStatTriang$cpFlag[j] <<- data$cpFlag[j]
                }          
            }
        }
    }
    BootStrapLine <<- BootstrapValues(data$value, family.pa, family.boot, scales, pattern = pattern, M, 
                                    alpha = 0.05, session)
}

BootLikelihood = function(data, coeff,family = c('Poisson', 'Gaussian', 'Bernoulli')){
    rm.flag = 0
    if (sum(coeff) == 0){
        rm.flag = 1
    } else {
      theta = sum(data * coeff) / sum(coeff)
      if (family == 'Poisson') {
          if (theta == 0) {
              rm.flag = 1
          } else {
            data = data * log(theta) - theta  
          }
      }
      if (family == 'Gaussian'){
          data = -.5 * (data - theta) * (data - theta)
      }
      if (family == 'Bernoulli') {
          if ((theta <= 0) | (theta >= 1)) {
              rm.flag = 1
          } else {
            data = data * (log(theta) - log(1 - theta)) + log(1 - theta) 
          }
      }
    }
    if (rm.flag == 1){
        out = NA
    } else {
      out = sum(data * coeff)   
    }
    return(out)
}

BootLikelihoodSm <- function(data, coeff, family = c('Poisson', 'Gaussian', 'Bernoulli')){
    rm.flag = 0
    r  = length(data) / 2
    sm = mean(data[1:r]) - mean(data[(r + 1): (2 * r)])  #bias correction
    if (sum(coeff) == 0){
        rm.flag = 1
    } else {
        theta = sum(data * coeff)  / sum(coeff)
        if (family == 'Gaussian') {
            data1 = -.5 * (data[1:r] - theta - sm) * (data[1:r] - theta - sm)
            data2 = -.5 * (data[(r + 1): (2 * r)] - theta) * (data[(r + 1): (2 * r)] - theta)
        }
        if (family == 'Poisson') {
            if ((theta + sm <= 0) | (theta <= 0)) {
                rm.flag = 1
            } else {
                data1 = log(theta + sm) * data[1:r] - (theta + sm)
                data2 = log(theta) * data[(r + 1):(2 * r)] - theta
            }
        }
        if (family == 'Bernoulli') {
            if ((theta >= 1) | (theta <= 0) | (theta + sm <= 0) | (theta + sm >= 1)) {
                rm.flag = 1
            } else {
                data1 = data[1:r] * (log(theta + sm) - log(1 - (theta + sm))) + log(1 - (theta + sm))
                data2 = data[(r + 1):(2 * r)] * (log(theta) - log(1 - theta)) + log(1 - theta)              
            }
        }
    }
    if (rm.flag == 1){
      out = NA
    } else {
      out = sum(c(data1, data2) * coeff)  
    }
    return(out)
}

BootLRT = function(data.left, data.right, coeff.left, coeff.right, family){
    LL.left  = BootLikelihood(data.left, coeff.left, family)
    LL.right = BootLikelihood(data.right, coeff.right, family)
    LL.lr    = BootLikelihoodSm(c(data.left, data.right), c(coeff.left, coeff.right), family)
    return(LL.left + LL.right - LL.lr)
}

BootstrapValues = function(data, family.pa, family.boot, H, pattern, M, alpha, session){
    stat.triang = matrix(0, 2 * length(H), M)
    stat        = matrix(NA, 2 * length(H), M * length(data))
    out         = rep(0, length(H)) 
    iter        = rep(0, length(H))
    withProgress(message = 'Applying bootstrap procedure', value = 0.01, {
        for (m in 1:M) {
            stat.max = matrix(NA, length(H), length(data))
            #generate bootstrap coefficients
            if (family.boot == 'Poisson') {
                coeff = rpois(length(data), 1)
            }
            if (family.boot == 'Gaussian') {
                coeff = rnorm(length(data), 1, 1)
            }
            for (h in seq_along(H)) {
                for (i in (2 * H[h] + 2) : length(data) ) {
                    stat.max[h, i] = sqrt(2* abs(BootLRT(data.left = data[(i - 2 * H[h] + 1):(i - H[h])], 
                                                    data.right = data[(i - H[h] + 1 ):i],  
                                                    coeff.left = coeff[(i - 2 * H[h] + 1):(i - H[h])], 
                                                    coeff.right = coeff[(i - H[h] + 1 ):i],
                                                    family = family.pa)))
                    if (is.na(stat.max[h,i])){
                        stat.max[h, i] = (stat.max[h, i - 2] + stat.max[h, i - 1]) / 2
                    }

                    if (pattern == 'Least squares') {
                        if (i >= 3 * H[h] + 1) {
                            iter[h]  = iter[h] + 1
                            triangle = lsfit((1:H[h]), stat.max[h,(i- H[h] + 1):i])
                        if (triangle$coef[2] < 0) {
                            triangle = rep(triangle$coef[1], H[h]) 
                        } else {
                            triangle = triangle$coef[2]*(1:H[h]) + triangle$coef[1]  
                        }
                        stat.triang[h, m] = max(stat.triang[h, m], sum(stat.max[h,(i- length(triangle) + 1):i] * triangle / sum(triangle)) )
                        stat[h, iter[h]]  = sum(stat.max[h,(i- length(triangle) + 1):i] * triangle / sum(triangle))
                      }            
                    }else if (pattern == 'Off') {
                        iter[h] <- iter[h] + 1
                        stat.triang[h, m] = max(stat.triang[h, m], stat.max[h,i] )
                        stat[h, iter[h]]  = stat.max[h,i]
                    } else {
                        if (i >= 3 * H[h] + 1) {
                            iter[h]           = iter[h] + 1
                            triangle          = 2 * (1 : H[h]) / (H[h] * (H[h] + 1))
                            stat.triang[h, m] = max(stat.triang[h, m], sum(stat.max[h,(i- length(triangle) + 1):i] * triangle) )
                            stat[h, iter[h]]  = sum(stat.max[h,(i- length(triangle) + 1):i] * triangle)
                        }
                    }
                }
            }
            incProgress(1/ M)
            Sys.sleep(0.1)
        }
    })
    ololo <<- stat.max
    alphas = seq(alpha,0 ,-0.0001)
    thlds  = matrix(0,length(H), length(alphas))
    for (h in seq_along(H)) {
        thlds[h, ] = quantile(stat[h,!is.na(stat[h, ]) ], probs = 1 - alphas)
    }
    withProgress(message = 'Computing critical values', value = 0.01 ,{
        for (i in seq_along(alphas)) {
            stat.normalized = stat.triang - thlds[, i]
            stat.normalized = apply(stat.normalized, 2, function(x) max(x))
            Fn = ecdf(stat.normalized)
            if ((1 - Fn(0)) <= alpha) {
                out = thlds[, i]
                break
            }
            incProgress(1/ length(alphas))
        }
    })
  return(out)
} 

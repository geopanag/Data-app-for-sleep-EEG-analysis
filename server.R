options(shiny.maxRequestSize=30*1024^2) 
options(rgl.useNULL=TRUE)
#options(warn=-1)

pkgTest <- function(x)
{
    if (!require(x,character.only = TRUE))
    {
        install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
}

libs<-c("shiny","ggplot2","data.table","eegkit","knitr","caret","signal","shinyRGL","fastICA","dplyr")
#suppressPackageStartupMessages(lapply(libs, require, character.only = TRUE))
sapply(libs,pkgTest)

aggregateEpoch=function(epoch,band){
    epoch = epoch %>% dplyr::filter(freq>band[1],freq<band[2])
    epoch = epoch[-ncol(epoch)]
    epoch=as.numeric(round(sapply(epoch,mean),digits=2))
    c(epoch,mean(epoch))
}

deriveSpec=function(data){
    s = spectrum(data,plot=FALSE)
    EEGspectrum= data.frame(s$spec)
    names(EEGspectrum) = names(data)
    EEGspectrum$freq = s$freq
    EEGspectrum
}

filterChannel =function(data,band,channel){
    filt = butter(3,c(band[1],band[2]),"pass")
    signal::filtfilt(filt, data[,channel])
}


shinyServer(function(input, output, session) {
    
    EEG=reactive({
        file = input$file
        validate(
            need(!is.null(file),"Blank File")
        )
        
        EEG=read.csv(file$datapath)
        
        if(input$ica){
            withProgress(message = 'Deriving ICA', value = 0, {
                preprocess=preProcess(EEG,method="ica",n.comp=ncol(EEG))
                EEGICA=predict(preprocess, EEG)
                names(EEGICA)=names(EEG)
                EEG = EEGICA
            })
        }
        
        EEG
    })
    
    
    observe({
        updateSelectInput(session, "channel", choices = c(names(EEG())))
    })
    
    
    band = reactive({
        switch(input$freq,
               "Delta (0.5 - 3.5 Hz)"=c(0.5,3.5)/100,
               "Theta (4 - 8 Hz)"=c(4,8)/100,
               "Alpha (8.5 - 12 Hz)"=c(8.5,12)/100,
               "Sigma (12.5 - 16 Hz)"=c(12.5,16)/100, 
               "Beta (16.5 - 30 Hz)"=c(16.5,30)/100,
               "Gamma (30.5 - 60 Hz)"=c(30.5,60)/100)##100= nyquist = sampling rate/2
    })
    
    
    E1spec=reactive({
        deriveSpec(EEG()[1:6000,])
    })
    
    E2spec=reactive({
        deriveSpec(EEG()[6001:12000,])
    })
    
    E3spec=reactive({
        deriveSpec(EEG()[12001:18000,])
    })
    
    E4spec=reactive({
        deriveSpec(EEG()[18001:24000,])
    })
    
    E5spec=reactive({
        deriveSpec(EEG()[24001:30000,])
    })
    
    E6spec=reactive({
        deriveSpec(EEG()[30001:36000,])
    })
    
    E7spec=reactive({
        deriveSpec(EEG()[36001:42000,])
    })
    
    E8spec=reactive({
        deriveSpec(EEG()[42001:48000,])
    })
    
    E9spec=reactive({
        deriveSpec(EEG()[48001:54000,])
    })
    
    E10spec=reactive({
        deriveSpec(EEG()[54001:60000,])
    })
    
    
    E11spec=reactive({
        deriveSpec(EEG()[60001:66000,])
    })
    
    E12spec=reactive({
        deriveSpec(EEG()[66001:72000,])
    })
    
    SpecBefore = reactive({
        band = band()
        sp = (E1spec()+E2spec())/2
        names(sp) = names(E1spec())
        sp = sp %>% dplyr::filter(freq>band[1],freq<band[2])
        sp = sp[-ncol(sp)]
        round(sapply(sp,mean),digits=2)
    })
    
    SpecAfter = reactive({
        band = band()
        sp = (E3spec()+E4spec()+E5spec()+E6spec()+E7spec()+E8spec()+E9spec()+E10spec()+E11spec()+E12spec())/10
        names(sp) = names(E3spec())
        sp = sp %>% dplyr::filter(freq>band[1],freq<band[2])
        sp = sp[-ncol(sp)]
        round(sapply(sp,mean),digits=2)
    })
    
    
    output$ScalpPlotBefore = renderWebGL({# renderWebGL({
        sp = SpecBefore()
        colfunc = colorRampPalette(c("white","red"))
        sp = sort(sp)
        
        elecs = names(sp)
        eegcap(electrodes=elecs,cex.point=50,col.point=colfunc(length(sp)))
    })
    
    output$ScalpPlotAfter= renderWebGL({
        sp = SpecAfter()
        colfunc = colorRampPalette(c("white","red"))
        sp = sort(sp)
        
        elecs = names(sp)
        #withProgress(message = 'Plotting Spectrum Scalp After Sleep', value = 0, {
        eegcap(electrodes=elecs,cex.point=50,col.point=colfunc(length(sp)))
        #})
    })
    
    output$ChannelE1= renderPlot({
        channel = filterChannel(EEG()[1:6000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE2= renderPlot({
        channel = filterChannel(EEG()[6001:12000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE3= renderPlot({
        channel = filterChannel(EEG()[12001:18000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE4= renderPlot({
        channel = filterChannel(EEG()[18001:24000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE5= renderPlot({
        channel = filterChannel(EEG()[24001:30000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE6= renderPlot({
        channel = filterChannel(EEG()[30001:36000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE7= renderPlot({
        channel = filterChannel(EEG()[36001:42000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE8= renderPlot({
        channel = filterChannel(EEG()[42001:48000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE9= renderPlot({
        channel = filterChannel(EEG()[48001:54000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE10= renderPlot({
        channel = filterChannel(EEG()[54001:60000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE11= renderPlot({
        channel = filterChannel(EEG()[60001:66000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    output$ChannelE12= renderPlot({
        channel = filterChannel(EEG()[66001:72000,],band(),input$channel)
        qplot(y=channel,x=seq_along(channel),geom="line")+labs(y="",x="")
    })
    
    
    output$fullSpec=function(){
        nyq=100
        band = band()
        tmp = E1spec()
        fullSpec = data.frame(matrix(nrow=ncol(tmp),ncol=13))
        names(fullSpec) = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","Mean")
        
        
        tmp = tmp %>% dplyr::filter(freq>band[1]/nyq,freq<band[2]/nyq)
        tmp = tmp[,-ncol(tmp)]
        row.names(fullSpec) = c(names(tmp),"Mean")
        tmp[,1]
        E1=as.numeric(round(sapply(tmp,mean),digits=2))
        
        fullSpec[,1] = c(E1,round(mean(E1),2))
        fullSpec[,2] = aggregateEpoch(E2spec(),band/nyq)
        fullSpec[,3] = aggregateEpoch(E3spec(),band/nyq)
        fullSpec[,4] = aggregateEpoch(E4spec(),band/nyq)
        fullSpec[,5] = aggregateEpoch(E5spec(),band/nyq)
        fullSpec[,6] = aggregateEpoch(E6spec(),band/nyq)
        fullSpec[,7] = aggregateEpoch(E7spec(),band/nyq)
        fullSpec[,8] = aggregateEpoch(E8spec(),band/nyq)
        fullSpec[,9] = aggregateEpoch(E9spec(),band/nyq)
        fullSpec[,10] = aggregateEpoch(E10spec(),band/nyq)
        fullSpec[,11] = aggregateEpoch(E11spec(),band/nyq)
        fullSpec[,12] = aggregateEpoch(E12spec(),band/nyq)
        fullSpec[,13] = round(apply(fullSpec,1,mean,na.rm=T),2)
        
        paste(knitr::kable(
            fullSpec, format = 'html', output = FALSE,
            table.attr="class='data table table-bordered table-condensed'"),
            sep = '\n')
    }
    
    
})



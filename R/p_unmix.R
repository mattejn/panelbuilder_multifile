
menuUnmix <- function(input,output,session,wsp) {
  output$menuUnmix <- shiny::renderUI("Unmixing")
}

serveUnmix <- function(input, output, session, wsp) {
  data <- shiny::reactiveValues(
    inputMtxList = NULL,
    inputNameList = NULL,
    outputNameList = NULL,
    outputMtxList = NULL,
    outputColnames = NULL,
    postCompensationList = NULL,
    postCompHistoryList = list(),
    postCompPtsList = NULL,
    len=NULL,
    currSelection=NULL
  )

  output$uiUnmixLoad <- shiny::renderUI(shiny::tagList(
    shiny::h3("Open a FCS"),
    shiny::fileInput('fileUnmixFCS',
      "FCS for unmixing", accept='.fcs', multiple = T),
  ))

  observeEvent(input$fileUnmixFCS, {
    data$len <- nrow(input$fileUnmixFCS)
    tempList <- vector(mode = "list", length = data$len)
    data$inputMtxList <- tempList
    data$inputNameList <- tempList
    data$outputMtxList <- tempList
    rm(tempList)
    data$outputColnames <- NULL
    shiny::withProgress(tryCatch({
      for (i in 1:data$len){
        m <- flowCore::read.FCS(input$fileUnmixFCS$datapath[i])@exprs
        colnames(m) <- unname(colnames(m))
        data$inputMtxList[[i]] <- m
        data$inputNameList[[i]] <- input$fileUnmixFCS$name[i]
        incProgress(1)
        }
      }, error=function(e) shiny::showNotification(type='error', paste("Loading failed:", e))),
      min=0,max=data$len,value=0,message="Loading FCS file(s)...")
  })

  output$uiUnmixControl <- shiny::renderUI(shiny::tagList(
    shiny::h3("Unmixing control"),
    shiny::uiOutput('uiUnmixLoadedSample'),
    shiny::selectizeInput('unmixMethod', "Unmixing method",
      choices = c(`OLS`='ols',
        `OLS weighted by spectra`='ols-spw',
        `OLS weighted by channel power`='ols-chw',
        `Event-weighted OLS (MAPE-like)`='eols-rw',
        `Positive-weighted gradient descent`='gd-pw'),
      selected="ols", multiple=F),
    shiny::checkboxInput('unmixIncludeFluorochromes', "Include fluorochrome names", value=T),
    shiny::checkboxInput('unmixIncludeOriginals', "Retain original values in the raw channels that were used for unmixing", value=F),
    shiny::checkboxInput('unmixIncludeResiduals', "Include per-channel residuals", value=F),
    shiny::checkboxInput('unmixIncludeRMSE', "Include total unmixing RMSE information", value=T),
    shiny::actionButton('doUnmix', "Run unmixing"),
    shiny::selectInput('PickUnmixedFile', 'Select unmixed sample',data$outputNameList,multiple=F)
  ))
  observeEvent(input$PickUnmixedFile, {
    data$currSelection=which(data$outputNameList==input$PickUnmixedFile)
  })

  output$uiUnmixLoadedSample <- shiny::renderUI(if(length(data$currSelection)==0) "No data loaded." else shiny::tagList(
    shiny::div(shiny::strong("Loaded: "), data$inputNameList[[data$currSelection]], paste0("(", nrow(data$inputMtxList[[data$currSelection]]), " events)")),
    {
      mcs <- matchingChans(data$inputMtxList[[data$currSelection]], getUnmixingInfo(wsp))
      umcs <- nat.sort(colnames(data$inputMtxList[[data$currSelection]])[!colnames(data$inputMtxList[[data$currSelection]]) %in% mcs])
      shiny::tagList(
        do.call(shiny::div, c(
          list(shiny::strong(paste0("Unmixing channels (",length(mcs),"):"))),
          lapply(mcs, function(mc) shiny::span(class="badge", mc)))),
        do.call(shiny::div, c(
          list(shiny::strong(paste0("Other channels (",length(umcs),"):"))),
          lapply(umcs, function(umc) shiny::span(class="badge", umc)))))
    }
  ))

  observeEvent(input$doUnmix,
    if(!is.null(data$inputMtxList))
      shiny::withProgress({
        tryCatch({
          for (i in 1:data$len){
            data$outputMtxList[[i]] <- doUnmix(data$inputMtxList[[i]], getUnmixingInfo(wsp),
            method=input$unmixMethod,
            fcNames=input$unmixIncludeFluorochromes,
            inclOrigs=input$unmixIncludeOriginals,
            inclResiduals=input$unmixIncludeResiduals,
            inclRmse=input$unmixIncludeRMSE)
            incProgress(1)
            data$outputNameList[[i]]=data$inputNameList[[i]]}},
          error=function(e) shiny::showNotification(type='error',
            paste("Unmixing failed:", e)))
          data$outputColnames <- colnames(data$outputMtxList[[1]])
          
      },min=0,max=data$len,value=0, message="Unmixing..."))


  output$uiUnmixPreview <- shiny::renderUI(shiny::tagList(
    shiny::h3("Result preview"),
    shiny::fluidRow(
      shiny::column(4,
        shiny::uiOutput('uiUnmixPlotOpts')),
      shiny::column(4, 
        shiny::h4("Data preview"),
        shiny::uiOutput('uiUnmixPlot'),
        shiny::sliderInput('unmixPlotAlpha', "Point alpha",
          min=0, max=1, step=0.01, value=0.5)),
      shiny::column(4,
        shiny::radioButtons('unmixTool',
          "Tool",
          choices=c(`Leveling tool`='level', `Gating`='gate'),
          selected='level'),
        shiny::uiOutput('uiUnmixTools'),
        shiny::h4("Results"),
        shiny::downloadButton('downloadUnmixFCS', "Download selected unmixed FCS"),
        shiny::downloadButton('downloadAllUnmixFCS', "Download all unmixed FCS(s)")))
  ))

  output$uiUnmixPlotOpts <- shiny::renderUI(if(!is.null(data$outputColnames)) shiny::tagList(
    shiny::selectizeInput('unmixPlotX',
      "Preview column X",
      multiple=F,
      choices=data$outputColnames,
      selected=defaultFSCChannel(data$outputColnames)),
    shiny::checkboxInput('unmixAsinhX',
      "Transform X",
      value=F),
    shiny::sliderInput('unmixCofX',
      "Cofactor X (dB)",,
      min=-10, max=80, step=1, value=30),
    shiny::selectizeInput('unmixPlotY',
      "Preview column Y",
      multiple=F,
      choices=data$outputColnames,
      selected=defaultSSCChannel(data$outputColnames)),
    shiny::checkboxInput('unmixAsinhY',
      "Transform Y",
      value=F),
    shiny::sliderInput('unmixCofY',
      "Cofactor Y (dB)",,
      min=-10, max=80, step=1, value=30),
    shiny::selectizeInput('unmixPlotCol',
      "Preview column color",
      multiple=F,
      choices=c('(Density)', data$outputColnames),
      selected='(Density)'),
    shiny::checkboxInput('unmixAsinhCol',
      "Transform Color",
      value=F),
    shiny::sliderInput('unmixCofCol',
      "Cofactor color (dB)",,
      min=-10, max=80, step=1, value=30)
  ))
  
  output$uiUnmixPlot <- shiny::renderUI(if(length(data$currSelection)>0) shiny::tagList(
    shiny::plotOutput('plotUnmix',
                      width="30em",
                      height="30em",
                      click=if(input$unmixTool=='level') 'clickUnmixPlot' else NULL,
                      brush=if(input$unmixTool=='gate') shiny::brushOpts('brushUnmixPlot')),
    if(!is.null(data$outputMtxList[[data$currSelection]]))
      shiny::div(paste0("(", nrow(data$outputMtxList[[data$currSelection]]), " events)"))
  ))

  getTransFns <- function() list(
    tx = if(input$unmixAsinhX) function(v)asinh(v/db2e(input$unmixCofX)) else identity,
    ty = if(input$unmixAsinhY) function(v)asinh(v/db2e(input$unmixCofY)) else identity,
    tc = if(input$unmixAsinhCol) function(v)asinh(v/db2e(input$unmixCofCol)) else identity,
    itx = if(input$unmixAsinhX) function(v)sinh(v)*db2e(input$unmixCofX) else identity,
    ity = if(input$unmixAsinhY) function(v)sinh(v)*db2e(input$unmixCofY) else identity)

  getCompData <- function(selected=data$currSelection)
    if(is.null(data$postCompensationList[[selected]])){ 
      data$outputMtxList[[selected]]
      }else{
      data$outputMtxList[[selected]] %*% data$postCompensationList[[selected]]
      }

  output$uiUnmixTools <- shiny::renderUI(if(input$unmixTool=='level') shiny::tagList(
    shiny::div("the tool aligns the selected × cross to the level of + cross"),
    ilDiv(
      shiny::actionButton('doUnmixLevelH', "Make horizontal"),
      shiny::actionButton('doUnmixLevelV', "Make vertical"),
      shiny::actionButton('doUnmixLevelUndo', "Undo"),
      shiny::actionButton('doUnmixLevelReset', "Reset")))
    else ilDiv(
      shiny::actionButton('doUnmixGateIn', "Keep only gate"),
      shiny::actionButton('doUnmixGateOut', "Remove gate")))

  observeEvent(input$doUnmixLevelReset, {
    data$postCompensationList[[data$currSelection]] <- diag(1, ncol(data$outputMtxList[[data$currSelection]]))
    colnames(data$postCompensationList[[data$currSelection]]) <- colnames(data$outputMtxList[[data$currSelection]])
    rownames(data$postCompensationList[[data$currSelection]]) <- colnames(data$outputMtxList[[data$currSelection]])
    data$postCompHistoryList[[data$currSelection]] <- list()
    data$postCompPtsList[[data$currSelection]] <- NULL
  })

  observeEvent(input$doUnmixLevelUndo, if(length(data$postCompHistoryList[[1]])>0) {
    data$postCompensation[[1]] <- data$postCompHistory[[data$currSelection]][[1]]
    data$postCompHistory[[data$currSelection]] <- data$postCompHistory[[data$currSelection]][-1]
  })

  observeEvent(data$outputMtxList, {
    if(length(data$currSelection)>0){
    if(is.null(data$outputMtxList[[data$currSelection]]))
      data$postCompensation <- NULL
    else {
      data$postCompensationList[[data$currSelection]] <- diag(1,ncol(data$outputMtxList[[data$currSelection]]))
      colnames(data$postCompensationList[[data$currSelection]]) <- colnames(data$outputMtxList[[data$currSelection]])
      rownames(data$postCompensationList[[data$currSelection]]) <- colnames(data$outputMtxList[[data$currSelection]])
    }
    data$postCompHistory <- list()
    data$postCompPts <- NULL
  }})

  observeEvent(input$unmixPlotXList,
    data$postCompPtsList[[data$currSelection]] <- NULL)

  observeEvent(input$unmixPlotYList,
    data$postCompPtsList[[data$currSelection]] <- NULL)

  observeEvent(input$clickUnmixPlot, {
    ts <- getTransFns()
    data$postCompPtsList[[data$currSelection]] <- rbind(c(ts$itx(input$clickUnmixPlot$x), ts$ity(input$clickUnmixPlot$y)), data$postCompPtsList[[data$currSelection]])
    if(nrow(data$postCompPtsList[[data$currSelection]])>2) data$postCompPtsList[[data$currSelection]] <- data$postCompPtsList[[data$currSelection]][1:2,,drop=F]
  })

  doAlign <- function(tr) {
    data$postCompHistory[[data$currSelection]] <- c(list(data$postCompensationList[[data$currSelection]]), data$postCompHistoryList[[data$currSelection]])
    ds <- c(input$unmixPlotX, input$unmixPlotY)
    data$postCompensationList[[data$currSelection]][ds,ds] <- data$postCompensationList[[data$currSelection]][ds,ds] %*% tr
    data$postCompPtsList[[data$currSelection]] <- NULL
  }

  observeEvent(input$doUnmixLevelH, if(!is.null(data$postCompPtsList[[data$currSelection]]) && nrow(data$postCompPtsList[[data$currSelection]])==2) {
    dstX <- data$postCompPtsList[[data$currSelection]][1,1]
    dstY <- data$postCompPtsList[[data$currSelection]][1,2]
    srcX <- data$postCompPtsList[[data$currSelection]][2,1]
    srcY <- data$postCompPtsList[[data$currSelection]][2,2]
    if(abs(dstX-srcX)<1) shiny::showNotification(type='error', "Source and destination horizontal coordinates too close")
    else doAlign(matrix(c(1,0,(srcY-dstY)/(dstX-srcX),1), 2))
  })

  observeEvent(input$doUnmixLevelV, if(!is.null(data$postCompPtsList[[data$currSelection]]) && nrow(data$postCompPtsList[[data$currSelection]])==2) {
    dstX <- data$postCompPtsList[[data$currSelection]][1,1]
    dstY <- data$postCompPtsList[[data$currSelection]][1,2]
    srcX <- data$postCompPtsList[[data$currSelection]][2,1]
    srcY <- data$postCompPtsList[[data$currSelection]][2,2]
    if(abs(dstY-srcY)<1) shiny::showNotification(type='error', "Source and destination vertical coordinates too close")
    else doAlign(matrix(c(1,(srcX-dstX)/(dstY-srcY),0,1), 2))
  })

  doGate <- function(b, inv) if(!is.null(b)) {
    ts <- getTransFns()
    flt <- xor(inv,
      data$outputMtxList[[data$currSelection]][,input$unmixPlotX] >= ts$itx(b$xmin) &
      data$outputMtxList[[data$currSelection]][,input$unmixPlotX] <= ts$itx(b$xmax) &
      data$outputMtxList[[data$currSelection]][,input$unmixPlotY] >= ts$ity(b$ymin) &
      data$outputMtxList[[data$currSelection]][,input$unmixPlotY] <= ts$ity(b$ymax))
    data$outputMtxList[[data$currSelection]] <- data$outputMtxList[[data$currSelection]][flt,,drop=F]
  }

  observeEvent(input$doUnmixGateIn, doGate(input$brushUnmixPlot, F))
  observeEvent(input$doUnmixGateOut, doGate(input$brushUnmixPlot, T))

 
  
  output$plotUnmix <- shiny::renderPlot({
    ts <- getTransFns()
    d <- getCompData()
    par(mar=c(0,0,0,0))
    if(!is.null(data$outputMtxList[[data$currSelection]]) && input$unmixPlotX!='' && input$unmixPlotY!='') {
      EmbedSOM::PlotEmbed(
        cbind(ts$tx(d[,input$unmixPlotX]), ts$ty(d[,input$unmixPlotY])),
        data=if(input$unmixPlotCol=='(Density)') NULL else cbind(ts$tc(d[,input$unmixPlotCol])),
        val=if(input$unmixPlotCol=='(Density)') 0 else 1,
        alpha=input$unmixPlotAlpha,
        plotf=scattermore::scattermoreplot)
      abline(h=0)
      abline(v=0)
      if(input$unmixTool=='level') points(
          ts$tx(rev(data$postCompPtsList[[data$currSelection]][,1])),
          ts$ty(rev(data$postCompPtsList[[data$currSelection]][,2])),
          cex=4, lwd=4, pch=c(4,3), col='#00cc00')
    }
  })

  output$downloadUnmixFCS <- shiny::downloadHandler(
    filename=function() paste0("pbUnmixed_",data$inputNameList[[data$currSelection]]),
    content=function(conn) flowCore::write.FCS(new('flowFrame', exprs=getCompData()), conn)
  )
  
  output$downloadAllUnmixFCS <- shiny::downloadHandler(
    filename = function(){
      paste0("unmixed_", Sys.Date(), ".zip")
    },
    content = function(file){
      
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      for (i in 1:length(data$outputNameList)){
        file_name <- paste0("pbUnmixed_",data$inputNameList[[i]])
        flowCore::write.FCS(new('flowFrame', exprs=getCompData(selected=i)), file.path(temp_directory, file_name))
      }
      
      
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
      
      
      
    },
    contentType = "application/zip"
  )

  output$uiUnmix <- shiny::renderUI(shiny::tagList(
    shiny::h1("Unmixing"),
    shiny::uiOutput('uiUnmixLoad'),
    shiny::uiOutput('uiUnmixControl'),
    shiny::uiOutput('uiUnmixPreview')))
}

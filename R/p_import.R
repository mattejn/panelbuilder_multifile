
menuImportSpectra <- function(input,output,session,wsp) {
  output$menuImportSpectra <- shiny::renderUI("Import spectra")
}

serveImportSpectra <- function(input, output, session, wsp) {
  
  data <- shiny::reactiveValues(
    fcsMtx = NULL,
    spectrum = NULL
  )

  observeEvent(input$fileImportFCS, {
    data$spectrum <- NULL
    data$fcsMtx <- NULL
    tryCatch({
      m <- flowCore::read.FCS(input$fileImportFCS$datapath)@exprs
      colnames(m) <- unname(colnames(m))
      data$fcsMtx <- m
      doAutofillAG()
      doAutofillFC()
    }, error=function(e) shiny::showNotification(type='error',
      paste("Loading failed:",e)))
  })

  output$uiImportColumnPicker <- shiny::renderUI(
    shiny::selectizeInput('importDataCols',
      "Fluorescent channels",
      multiple=T,
      choices=colnames(data$fcsMtx),
      selected=defaultFluorescenceChannels(nat.sort(colnames(data$fcsMtx))),
      options=list(`actions-box`=TRUE)))

  output$uiImportEliminatePicker <- shiny::renderUI({
    availSpectra <- if(is.null(data$fcsMtx) || length(wsp$spectra)==0) c()
      else {
        flt <- sapply(wsp$spectra,
          function(s) identical(nat.sort(s$spectrum$channels), nat.sort(input$importDataCols))) 
        nms <- sapply(wsp$spectra[flt],
          function(s) paste(
            s$machine,"/",
            s$mconfig,"/",
            s$sample,"/",
            s$antigen,"/",
            s$fluorochrome,"/",
            s$note))
        avs <- seq_len(length(wsp$spectra))[flt]
        names(avs) <- nms
        avs
      }
    shiny::selectizeInput('importEliminateSpectra',
      "Eliminate existing spectrum from signal",
      multiple=T,
      choices=availSpectra,
      selected=NULL,
      options=list(`actions-box`=TRUE))
  })

  output$uiImportGate <- shiny::renderUI(shiny::tagList(
    shiny::selectizeInput('importGateColX',
      "Gate column X",
      multiple=F,
      choices=colnames(data$fcsMtx),
      selected=defaultFSCChannel(colnames(data$fcsMtx))),
    shiny::selectizeInput('importGateColY',
      "Gate column Y",
      multiple=F,
      choices=colnames(data$fcsMtx),
      selected=defaultSSCChannel(colnames(data$fcsMtx)))
  ))

  getGate <- function() {
    if(is.null(data$fcsMtx)) return(NULL)
    if(input$importGateColX=='' || input$importGateColY=='' || is.null(input$importGateBrush))
      return(T)
    b <- input$importGateBrush
    flt <-
      data$fcsMtx[,input$importGateColX] >= b$xmin &
      data$fcsMtx[,input$importGateColX] <= b$xmax &
      data$fcsMtx[,input$importGateColY] >= b$ymin &
      data$fcsMtx[,input$importGateColY] <= b$ymax
    return(flt)
  }

  getGatedData <- function() {
    flt <- getGate()
    if(is.null(flt)) flt <- T
    return(data$fcsMtx[flt,,drop=F])
  }

  output$plotImportGate <- shiny::renderPlot({
    par(mar=c(0,0,0,0))
    if(!is.null(data$fcsMtx) && input$importGateColX!='' && input$importGateColY!='') {
      EmbedSOM::PlotEmbed(
        data$fcsMtx[,c(input$importGateColX, input$importGateColY)],
        fdens=log,
        plotf=scattermore::scattermoreplot)
    }
  })

  getPowerGate <- function() {
    if(is.null(data$fcsMtx)) return(NULL)
    b <- input$importPowerBrush
    if(is.null(b) || is.null(input$importDataCols)) return(T)
    e <- e2db(sqrt(apply(data$fcsMtx[,input$importDataCols,drop=F]^2,1,sum)))
    flt <- e >= b$xmin & e<=b$xmax
    return(flt)
  }

  output$uiImportPower <- shiny::renderUI(shiny::tagList(
    shiny::h4('Fluorescent power distribution'),
    shiny::plotOutput('plotImportPower',
      width="20em",
      height="15em",
      brush=shiny::brushOpts('importPowerBrush', direction='x'))
  ))

  output$plotImportPower <- shiny::renderPlot({
    par(mar=c(2,2,0,0))
    d <- getGatedData()
    if(!is.null(d) && !is.null(input$importDataCols)) {
      d <- e2db(sqrt(apply(d[,input$importDataCols,drop=F]^2,1,sum)))
      xlim <- c(min(d),max(d))
      d <- density(d, from=xlim[1], to=xlim[2], bw=1)
      xs <- c(xlim[1], d$x, xlim[2])
      ys <- c(0, sqrt(d$y), 0)
      if(max(ys)>0) ys <- ys/max(ys)
      plot(type='n', NULL, xlim=xlim, ylim=c(0,1))
      polygon(xs,ys,lwd=2,col='#dddddd',border='#888888')
    }
  })

  observeEvent(input$doImportGetSpectrum, {
    if(is.null(data$fcsMtx)) return()
    if(is.null(input$importDataCols)) return()
    cols <- nat.sort(input$importDataCols)
    data$spectrum <- NULL
    elim <- if(is.null(input$importEliminateSpectra)) NULL
      else matrix(nrow=length(cols),
        sapply(wsp$spectra[as.integer(input$importEliminateSpectra)],
        function(s)s$spectrum$mS[indexin(cols, s$spectrum$channels)]))
    tryCatch(
      data$spectrum <- extractSpectrum(
        data$fcsMtx[,cols,drop=F],
        input$importMethod,
        getGate(), getPowerGate(), elim),
      error=function(e)
        shiny::showNotification(type='error',
          paste("Spectrum import failed:", e)))
  })

  observeEvent(input$doImportSave, {
    if(is.null(data$spectrum)) return()
    n <- c(
      spectrumMetadataFormGather('importForm', input),
      list(spectrum=data$spectrum))
    if(spectrumExistsIn(n, wsp$spectra)) {
      shiny::showNotification(type='error',
        "This spectrum already exists; use a different note to distinguish duplicates.")
      return ()
    }

    wsp$spectra = c(wsp$spectra, list(n))

    shiny::showNotification(type='message',
      "Saved OK.")
  })

  output$uiImportResult <- shiny::renderUI(shiny::tagList(
    shiny::h3("Computed spectrum"),
    if(is.null(data$spectrum)) "No spectrum computed"
    else shiny::tagList(
      paste0("Mean fluorescent intensity: ", dbf(data$spectrum$mI),
             ", signal σ: ", dbf(data$spectrum$sdI,'')),
      shiny::plotOutput('plotImportSpectra', width='100%', height='40ex')
    )
  ))

  output$plotImportSpectra <- shiny::renderPlot(
    if(!is.null(data$spectrum))
      plotSpectrum(data$spectrum$mS, data$spectrum$sdS, data$spectrum$channels)
  )

  output$uiImportSave <- shiny::renderUI({
    wsp$page # reload the information from wsp$spectra after the page changed
    shiny::tagList(
      shiny::h3("Save the spectrum"),
      shiny::checkboxInput('importAutofillAG',
        "Autofill antigen",
        value=T),
      shiny::checkboxInput('importAutofillFC',
        "Autofill fluorochrome",
        value=T),
      spectrumMetadataForm('importForm', isolate(wsp$spectra)),
      shiny::uiOutput('uiImportSaveBtn')
  )})

  output$uiImportSaveBtn <- shiny::renderUI(
    if(is.null(data$spectrum)) "No spectrum to save"
    else shiny::actionButton('doImportSave', "Save")
  )

  doAutofillAG <- function() {
    if(!is.null(input$fileImportFCS) && !is.null(input$importAutofillAG))
      if(input$importAutofillAG) {
      availAntigens <- gatherFormContent('antigen',
        wsp$spectra,
        defaultAntigens)
      updateSelectizeInput(session,
        'importFormAntigen',
        selected=availAntigens[which.min(
          adist(availAntigens, input$fileImportFCS$name,
            ignore.case=T, costs=c(del=3, sub=4)))])
      }
  }
  doAutofillFC <- function() {
    if(!is.null(input$fileImportFCS) && !is.null(input$importAutofillFC))
      if(input$importAutofillFC) {
        availFluorochromes <- gatherFormContent('fluorochrome',
          wsp$spectra,
          defaultFluorochromes)
        updateSelectizeInput(session,
          'importFormFluorochrome',
          selected=availFluorochromes[which.min(
            adist(availFluorochromes, input$fileImportFCS$name,
              ignore.case=T, costs=c(del=3, sub=4)))])
      }
  }

  shiny::observeEvent('importAutofillAG', {
    doAutofillAG()
  })
  shiny::observeEvent('importAutofillFC', {
    doAutofillFC()
  })

  output$uiImportSpectra <- shiny::renderUI(shiny::tagList(
    shiny::h1("Import spectral profiles"),
    shiny::fluidRow(
      shiny::column(4,style='min-width:20em',
        shiny::fileInput('fileImportFCS', "Upload an FCS file", accept='.fcs'),
        shiny::uiOutput('uiImportColumnPicker'),
        shiny::uiOutput('uiImportEliminatePicker'),
        shiny::selectizeInput('importMethod', "Spectrum computation method",
          choices=c(`OLS`='ols', `Approximated Theil-Sen estimator`='theil-sen'),
          selected='ols',
          multiple=F),
        shiny::actionButton('doImportGetSpectrum', "Compute spectrum")),
      shiny::column(4,style='min-width:20em',
        shiny::uiOutput('uiImportGate'),
        shiny::plotOutput('plotImportGate',
          width="20em",
          height="20em",
          brush=shiny::brushOpts('importGateBrush')),
        shiny::uiOutput('uiImportPower')),
      shiny::column(4,style='min-width:10em',
        shiny::uiOutput('uiImportSave'))),
    shiny::uiOutput('uiImportResult')
  ))
}

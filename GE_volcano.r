#run in python:
#import subprocess
#subprocess.call (["/usr/bin/Rscript", "--vanilla", "/home/danial/Desktop/makescriptfoit/volcano.r"])


library(ggrepel)
NanoString_NKI <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[,x], na.rm=TRUE),
    max(toptable[,x], na.rm=TRUE)),
  ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 0),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 18,
  title = 'NanoString_NKI',
  caption = paste0('Total = ', nrow(toptable), ' variables'),
  titleLabSize = 18,
  captionLabSize = 14,
  pCutoff = 10e-6,
  pLabellingCutoff = pCutoff,
  FCcutoff = 1.0,
  cutoffLineType = 'longdash',
  cutoffLineCol = 'black',
  cutoffLineWidth = 0.4,
  transcriptPointSize = 0.8,
  transcriptLabSize = 3.0,
  transcriptLabCol = 'black',
  transcriptLabFace = 'plain',
  transcriptLabhjust = 0,
  transcriptLabvjust = 1.5,
  boxedlabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  colCustom = NULL,
  colAlpha = 1/2,
  legend = c("NS","Log2 FC","P","P & Log2 FC"),
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 4.0,
  legendVisible = TRUE,
  shade = NULL,
  shadeLabel = NULL,
  shadeAlpha = 1/2,
  shadeFill = "grey",
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey10',
  hline = NULL,
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = "partial",
  borderWidth = 0.8,
  borderColour = "black")
{
  if(!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call.=FALSE)
  }

  if(!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call.=FALSE)
  }

  if(!is.numeric(toptable[,x])) {
    stop(paste(x, " is not numeric!", sep=""))
  }

  if(!is.numeric(toptable[,y])) {
    stop(paste(y, " is not numeric!", sep=""))
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[,x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[,y] < pCutoff)] <- "P"
  toptable$Sig[(toptable[,y] < pCutoff) &
    (abs(toptable[,x]) > FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig,
    levels=c("NS","FC","P","FC_P"))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  #####
  # New functionality in > v1.2:
  # Now convert to 10^-1 lower than lowest non-zero p-value
  if (min(toptable[,y], na.rm=TRUE) == 0) {
    # <= v1.2
    #warning(paste("One or more P values is 0.",
    #  "Converting to minimum possible value..."),
    #  call. = FALSE)
    #toptable[which(toptable[,y] == 0), y] <- .Machine$double.xmin
    warning(paste("One or more p-values is 0.",
      "Converting to 10^-1 * current",
      "lowest non-zero p-value..."),
      call. = FALSE)
    toptable[which(toptable[,y] == 0), y] <- min(
      toptable[which(toptable[,y] != 0), y],
      na.rm = TRUE) * 10^-1
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[,x]
  toptable$yvals <- toptable[,y]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      # title, subtitle, and caption
      plot.title = element_text(
        angle = 0,
        size = titleLabSize,
        face = 'bold',
        vjust = 1),
      plot.caption = element_text(
        angle = 0,
        size = captionLabSize,
        face = 'plain',
        vjust = 1),

      # axis text
      axis.text.x = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.text.y = element_text(
        angle = 0,
        size = axisLabSize,
        vjust = 1),
      axis.title = element_text(
        size = axisLabSize),

      # legend
      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(
        size = legendLabSize),
      title = element_text(
        size = legendLabSize),
      legend.title = element_blank())

  # Create the plot object differently based on whether colCustom 
  # and shapeCustom are NULL or not. This helps to avoid messing up
  # the legend.
  #
  # 1, both colCustom and shapeCustom are activated
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new shape and colour encodings as aes
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = transcriptPointSize,
        na.rm = TRUE) +

      # specify the colour and shape with the supplied encoding
      scale_color_manual(values = colCustom) +
      scale_shape_manual(values = shapeCustom)

  # 2, only colCustom is activated and 'shape' has just a single value
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate IF shape is also
      # included as aes (it is not, here)
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included, but outside aes
      geom_point(
        aes(
          color = factor(names(colCustom))),
        alpha = colAlpha,
        shape = shape,
        size = transcriptPointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values=colCustom) +

      # 'shape' is not included as aes. Specifying guide = TRUE
      # here will result in legends merging
      scale_shape_manual(guide = TRUE)

  # 3, only colCustom is activated and 'shape' has 4 values
  } else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new colour encodings as aes.
      # 'shape' is included in aes and mapped to 4
      # categories of NS, FC, P, FC_P
      geom_point(
        aes(
          color = factor(names(colCustom)),
          shape = factor(Sig)),
        alpha = colAlpha,
        size = transcriptPointSize,
        na.rm = TRUE) +

      # specify the colour with the supplied encoding
      scale_color_manual(values = colCustom) +

      # as it is included as aes, a separate legend
      # for 'shape' will be drawn. Here, over-ride that
      # legend
      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        labels=c(
          NS=legend[1],
          FC=paste(legend[2], sep=""),
          P=paste(legend[3], sep=""),
          FC_P=paste(legend[4], sep="")),
        guide = TRUE)

  # 4, only shapeCustom is activated
  } else if (is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # guide_legends are separate for colour and shape;
      # so, legends will be drawn separate
      guides(
        colour = guide_legend(
          order = 1,
          override.aes=list(
            size = legendIconSize)),
        shape = guide_legend(
          order = 2,
          override.aes=list(
            size = legendIconSize))) +

      # include new shape encodings as aes.
      # Standard colour for NS, FC, P, FC_P,
      # are added to aes, too.
      geom_point(
        aes(
          color = factor(Sig),
          shape = factor(names(shapeCustom))),
        alpha = colAlpha,
        size = transcriptPointSize,
        na.rm = TRUE) +

      # as it is included as aes, a separate legend
      # for 'colour' will be drawn. Here, over-ride that
      # legend
      scale_color_manual(
        values=c(
          NS=col[1],
          FC=col[2],
          P=col[3],
          FC_P=col[4]),
        labels=c(
          NS=legend[1],
          FC=paste(legend[2], sep=""),
          P=paste(legend[3], sep=""),
          FC_P=paste(legend[4], sep=""))) +

      # specify the shape with the supplied encoding
      scale_shape_manual(values = shapeCustom)

  # 5, both colCustom and shapeCustom are null;
  # only a single shape value specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 1) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # including 'shape' in the colour guide_legend here
      # results in the legends merging
      guides(colour = guide_legend(
        order = 1,
        override.aes=list(
          shape = shape,
          size = legendIconSize))) +

      geom_point(
        aes(color = factor(Sig)),
        alpha = colAlpha,
        shape = shape,
        size = transcriptPointSize,
        na.rm = TRUE,
        show.legend = legendVisible) +

      scale_color_manual(
        values = c(
          NS = col[1],
          FC = col[2],
          P = col[3],
          FC_P = col[4]),
        labels = c(
          NS = legend[1],
          FC = paste(legend[2], sep=""),
          P = paste(legend[3], sep=""),
          FC_P = paste(legend[4], sep="")))

  # 6, both colCustom and shapeCustom are null;
  # four shape values are specified
  } else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 4) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      # over-ride legend icon sizes for colour and shape.
      # including 'shape' in the colour guide_legend here
      # results in the legends merging
      guides(colour = guide_legend(
        order = 1,
        override.aes = list(
          shape = c(
            NS = shape[1],
            FC = shape[2],
            P = shape[3],
            FC_P = shape[4]),
          size = legendIconSize))) +

      geom_point(
        aes(
          color = factor(Sig),
          shape = factor(Sig)),
        alpha = colAlpha,
        size = transcriptPointSize,
        na.rm = TRUE,
        show.legend = legendVisible) +

      scale_color_manual(
        values = c(
          NS = col[1],
          FC = col[2],
          P = col[3],
          FC_P = col[4]),
        labels = c(
          NS = legend[1],
          FC = paste(legend[2], sep=""),
          P = paste(legend[3], sep=""),
          FC_P = paste(legend[4], sep=""))) +

      scale_shape_manual(
        values = c(
          NS = shape[1],
          FC = shape[2],
          P = shape[3],
          FC_P = shape[4]),
        guide = FALSE)
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth) +

    geom_hline(yintercept = -log10(pCutoff),
      linetype = cutoffLineType,
      colour = cutoffLineCol,
      size = cutoffLineWidth)

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline),
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # Border around plot
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(
      colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(
      size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }

  # Gridlines
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # user has specified to draw with geom_text or geom_label?
  if (boxedlabels == FALSE) {
    # For labeling with geom_text/label_repel (connectors) and
    # geom_text/label (.., check_overlap = TRUE), 4 possible
    # scenarios can arise
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_text_repel(
        data=subset(toptable,
          !is.na(toptable[,"lab"])),
        aes(label=subset(toptable,
          !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_text(
        data=subset(toptable,
          !is.na(toptable[,"lab"])),
        aes(
          label=subset(toptable,
            !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        check_overlap = TRUE,
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_text(
        data=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        check_overlap = TRUE,
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    }
  } else {
    # For labeling with geom_text/label_repel (connectors) and
    # geom_text/label (.., check_overlap = TRUE), 4 possible
    # scenarios can arise
    if (drawConnectors == TRUE && is.null(selectLab)) {
      plot <- plot + geom_label_repel(
        data=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[,y]<pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == TRUE && !is.null(selectLab)) {
      plot <- plot + geom_label_repel(
        data=subset(toptable,
          !is.na(toptable[,"lab"])),
        aes(label=subset(toptable,
          !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && !is.null(selectLab)) {
      plot <- plot + geom_label(
        data=subset(toptable,
          !is.na(toptable[,"lab"])),
        aes(
          label=subset(toptable,
            !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        #check_overlap = TRUE,
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    } else if (drawConnectors == FALSE && is.null(selectLab)) {
      plot <- plot + geom_label(
        data=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff),
        aes(label=subset(toptable,
          toptable[,y] < pLabellingCutoff &
            abs(toptable[,x]) > FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        #check_overlap = TRUE,
        hjust = transcriptLabhjust,
        vjust = transcriptLabvjust,
        colour = transcriptLabCol,
        fontface = transcriptLabFace,
        na.rm = TRUE)
    }
  }

  # shading
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = 'polygon',
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE) +

      scale_fill_identity(name = shadeLabel,
        labels = shadeLabel,
        guide = 'legend')
  }

  return(plot)
}



library(edgeR)

target <-read.delim("metadata.txt", header=T)
x <- read.csv('Normalized_data.txt', row.names = 1)
group <- factor(target$group)
design <- model.matrix(~0+group)
comparison <- paste(dimnames(design)[[2]], collapse = "-")
x <- read.csv('Normalized_data.txt', sep = "\t", row.names = 1)
x <- x[, 4:ncol(x)]
dge <- DGEList(counts=x)
v <- voom(dge, design)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(comparison, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, trend=FALSE)
results <- decideTests(fit2, p.value = .1, lfc= log2(2) )
topGenes =topTable(fit2, number = 1e12,sort.by="M" )
write.table(topGenes, "DGE.csv", sep = ",")

pdf("NanoString_Volcano.pdf",paper='special') 
NanoString_NKI(topGenes,
    lab = rownames(topGenes),
    x = 'logFC',
    y = 'P.Value',
    xlim = c(-3, 3),
    pCutoff = 0.05,
    FCcutoff = 0.6,
    cutoffLineType = 'twodash',
    cutoffLineWidth = 0.8,
    transcriptPointSize = 3.0,
    transcriptLabSize = 4.0,
    colAlpha = 1,
    legend=c('NS','Log (base 2) fold-change','P value',
      'P value & Log (base 2) fold-change'),
    legendPosition = 'top',
    legendLabSize = 10,
    legendIconSize = 5.0)
dev.off()


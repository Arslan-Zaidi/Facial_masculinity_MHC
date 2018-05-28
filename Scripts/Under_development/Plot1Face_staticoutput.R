#!/usr/bin/Rscript

#Author: Tomas B. Gonzalez-Zarzar
#Date: 5/8/18
#Edits by JDWhite and AAZaidi - 05/23/2018

args<-commandArgs(TRUE)

#first argument is the individual's index  corresponding to the order in the euro_1233_masc_het_03292018.dat file
col.index<-args[1]

#second argument is the name of output file
out.prefix<-args[2]

require(plotly)
require(data.table)
require(mixOmics)
if (!require("webshot")) install.packages("webshot")

#load template facial mask
RefScan <- read.table("../Dataset/RefScan.obj", sep = "\t", header = F)
RefScan_Vertices <- as.matrix(RefScan[which(RefScan$V1 == "v"), 2:4])
RefScan_Facets <- as.matrix(RefScan[which(RefScan$V1 == "f"), 2:4])

#load FM_ql file
message("Loading Facial masculinity per QL")
qlmasc<-as.matrix(fread("../Dataset/qlmasc_1233_noprop_03292018.txt",header=F))
id1.ql<-qlmasc[col.index,]

#define function to plot
Plot1Face <- function(vertices, facets, colormap=NULL, title=NULL,color.min=NULL,color.max=NULL)
{
  # Plot one face, with the corresponding facets and colormap if given.
  # vertices is a matrix with rows = n landmarks and 3 columns (x, y, z) (v in an obj file)
  # facets is a matrix with rows = n landmarks and 3 columns corresponding to the edge connections between landmarks (f in an obj file)
  # colormap is a numeric vector with entries = n landmarks corresponding to the values to be referenced as colors
  
  ax <- list( #Creating axes
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  
  scene <- list( #Creating scene
    xaxis = ax, 
    yaxis = ax, 
    zaxis = ax, 
    camera = list( #Change the position of camera
      eye = list(
        x = 0, 
        y = -0.1, 
        z = 1.7),
      center = list(
        x = 0,
        y = 0,
        z = 0),
      up = list(
        x = 0,
        y = 1,
        z = 0)
    )
  )
  
  if(is.null(colormap)){
    facecolor1 <- rep("lightgray", nrow(facets))
    myPlot <- plot_ly(
      x = vertices[,1], y = vertices[,2], z = vertices[,3],
      i = facets[, 1]-1, j = facets[, 2]-1, k = facets[, 3]-1,
      facecolor = facecolor1, opacity = 1,
      hoverinfo = "none",
      lighting = list(specular=0.0, ambient=0.3, diffuse=0.8, fresnel=0.1, roughness=0.1),  
      lightposition = list(x = 100, y = 200, z = 0), 
      type = "mesh3d"
    ) %>% 
      layout(scene = scene, title = title)
  } else {
    if(is.null(color.min)){color.min = min(colormap)}else{color.min = color.min} 
    if(is.null(color.max)){color.max = max(colormap)}else{color.max = color.max}
    myPlot <- plot_ly(
      x = vertices[,1], y = vertices[,2], z = vertices[,3],
      i = facets[, 1]-1, j = facets[, 2]-1, k = facets[, 3]-1,
      colors = color.jet(length(colormap), alpha = 1), #Colormap colors
      intensity = colormap, 
      opacity = 1,
      hoverinfo = "none",
      lighting = list(specular=0.0, ambient=0.3, diffuse=0.8, fresnel=0.1, roughness=0.1),  
      lightposition = list(x = 100, y = 200, z = 0),
      cmin=color.min,
      cmax=color.max,
      type = "mesh3d"
    ) %>% 
      layout(scene = scene, title = title)
  }
  return(myPlot)
}


plt<-Plot1Face(vertices = RefScan_Vertices, facets = RefScan_Facets, colormap = id1.ql, title = "FM_ql")
tmpFile <- tempfile(fileext = ".png")
export(plt, tmpFile)
browseURL(tmpFile)

#plotly_IMAGE(plt, format = "png", out_file = paste(out.prefix,".png",sep=""))


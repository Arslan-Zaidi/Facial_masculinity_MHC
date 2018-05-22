#Author: Tomas B. Gonzalez-Zarzar
#Date: 5/8/18

Plot1Face <- function(vertices, facets, colormap=NULL, title=NULL)
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
      lighting = list(specular=0.05, ambient=0.3, diffuse=0.8, fresnel=0.1, roughness=0.1),  
      lightposition = list(x = 100, y = 200, z = 0), 
      type = "mesh3d"
    ) %>% 
      layout(scene = scene, title = title)
  } else {
    
    myPlot <- plot_ly(
      x = vertices[,1], y = vertices[,2], z = vertices[,3],
      i = facets[, 1]-1, j = facets[, 2]-1, k = facets[, 3]-1,
      colors = color.jet(length(colormap), alpha = 1), #Colormap colors
      intensity = colormap, 
      opacity = 1,
      hoverinfo = "none",
      lighting = list(specular=0.05, ambient=0.3, diffuse=0.8, fresnel=0.1, roughness=0.1),  
      lightposition = list(x = 100, y = 200, z = 0), 
      type = "mesh3d"
    ) %>% 
      layout(scene = scene, title = title)
  }
  return(myPlot)
}

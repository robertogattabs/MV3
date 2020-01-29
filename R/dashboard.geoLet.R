#' a dashboard for geoLet
#'
#' @description  a GUI for geoLet
#' @import stringr shiny shinydashboard
dashboard.geoLet <- function() {
  ui <- dashboardPage(
    dashboardHeader(title = "geoLet"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("openDICOMFolder", tabName = "openDICOMFolder", icon = icon("th")),
        menuItem("getImageVoxelCube", tabName = "getImageVoxelCube", icon = icon("th") ),
        menuItem("getROIVoxel", tabName = "getROIVoxel", icon = icon("th") ),
        menuItem("getROIList", tabName = "getROIList", icon = icon("th") )
      )
    ),
    dashboardBody(

           
    )
  )
  
  server <- function(input, output) {
    output$menu <- renderMenu({
      sidebarMenu(
        menuItem("Menu item", icon = icon("calendar"))
      )
    })
  }
  
  shinyApp(ui, server)
}


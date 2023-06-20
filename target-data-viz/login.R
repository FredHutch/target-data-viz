library(shiny)

# dataframe that holds usernames, passwords and other user data
user_base <- tibble::tibble(
  user = "dataviz",
  password = "hutch2025",
  permissions = "standard",
  name = "User 1"
)

ui <- fluidPage(
  # add logout button UI
  div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
  
  # add login panel UI function
  shinyauthr::loginUI(id = "login"),
  
  # setup table output to show user info after login
  tableOutput("user_table")
)

server <- function(input, output, session) {
  
  # call login module supplying data frame, 
  # user and password cols and reactive trigger
  credentials <- shinyauthr::loginServer(
    id = "login",
    data = user_base,
    user_col = user,
    pwd_col = password,
    log_out = reactive(logout_init())
  )
  
  # call the logout module with reactive trigger to hide/show
  logout_init <- shinyauthr::logoutServer(
    id = "logout",
    active = reactive(credentials()$user_auth)
  )
  
  output$user_table <- renderTable({
    # use req to only render results when credentials()$user_auth is TRUE
    req(credentials()$user_auth)
    credentials()$info
  })
}

shinyApp(ui = ui, server = server)
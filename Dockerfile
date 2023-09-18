FROM fredhutch/r-shiny-server-base:4.1.1
RUN apt-get update -y && apt-get install -y pandoc libpq-dev
RUN R -q -e 'install.packages(c("shinyalert", "shinydashboard", "shinyjs", "shinyBS", "shinyWidgets", "shinythemes", "survminer", "gtools", "cmprsk", "ggpubr", "DT", "data.table", "viridis", "viridisLite", "ggplot2", "plotly", "fst", "ComplexHeatmap"))'
RUN R -q -e 'install.packages(c("dplyr"), repos="https://cran.r-project.org")'

RUN rm -rf /srv/shiny-server/
ADD *.R /srv/shiny-server/

# adds each subdirectory to the Docker image
ADD data /srv/shiny-server/data/
ADD www /srv/shiny-server/www/

RUN chown -R shiny:shiny /srv/shiny-server/

EXPOSE 3838

WORKDIR /srv/shiny-server

CMD /usr/bin/shiny-server

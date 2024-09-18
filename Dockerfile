FROM fredhutch/r-shiny-server-base:4.3.0
RUN apt-get update -y && apt-get install -y pandoc libpq-dev python3 python3-pip
RUN R -q -e 'install.packages(c("shinyalert", "shinydashboard", "shinyjs", "shinyBS", "shinyWidgets", "shinythemes", "survminer", "survival", "gtools", "cmprsk", "ggpubr", "DT", "data.table", "viridis", "viridisLite", "ggplot2", "plotly", "fst", "BiocManager", "ggsurvfit", "openxlsx"))'
RUN R -q -e 'install.packages(c("dplyr"), repos="https://cran.r-project.org")'
RUN R -q -e 'BiocManager::install(c("ComplexHeatmap", "biomaRt"), ask=FALSE, update=FALSE)'

RUN pip3 install pybiolib

RUN rm -rf /srv/shiny-server/
ADD *.R /srv/shiny-server/

# adds each subdirectory to the Docker image
ADD data /srv/shiny-server/data/
ADD www /srv/shiny-server/www/
ADD check.R /tmp/

RUN chown -R shiny:shiny /srv/shiny-server/

EXPOSE 3838

WORKDIR /srv/shiny-server

RUN R -f /tmp/check.R --args shinyalert shinydashboard shinyjs shinyBS shinyWidgets shinythemes survminer gtools cmprsk ggpubr DT data.table viridis viridisLite ggplot2 plotly fst BiocManager dplyr ComplexHeatmap survival ggsurvfit openxlsx biomaRt

RUN rm /tmp/check.R

CMD R -q -f start.R

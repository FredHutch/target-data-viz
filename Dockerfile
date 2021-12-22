# build me as dockerimages.fhcrc.org/monitords:latest
FROM fredhutch/r-shiny-server-base:4.1.1
RUN apt-get update -y && apt-get install -y pandoc libpq-dev supervisor nginx
RUN R -q -e 'install.packages(c("shinyalert", "shinydashboard", "shinyjs", "shinyWidgets", "shinythemes", "survminer", "gtools", "cmprsk"))'
RUN R -q -e 'install.packages(c("dplyr"), repos="https://cran.r-project.org")'

RUN rm -rf /srv/shiny-server/
ADD *.R /srv/shiny-server/

ADD data /srv/shiny-server/data/


ADD system/. /home/shiny/system/


EXPOSE 8888

RUN ln -sf /dev/stdout /var/log/nginx/access.log && ln -sf /dev/stderr /var/log/nginx/error.log

WORKDIR /srv/shiny-server

CMD ["/bin/sh", "-c", "/usr/bin/supervisord -c /home/shiny/system/sup.conf"]


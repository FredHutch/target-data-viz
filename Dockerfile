# build me as dockerimages.fhcrc.org/monitords:latest
FROM fredhutch/r-shiny-server-base:latest
RUN apt-get update -y && apt-get install -y pandoc libpq-dev supervisor nginx
RUN R -q -e 'install.packages(c("shinydashboard", "shinyjs", "shinyWidgets", "shinythemes", "survminer", "gtools"))'

RUN rm -rf /srv/shiny-server/
ADD *.R /srv/shiny-server/

ADD data /srv/shiny-server/data/


ADD system/. /home/shiny/system/


EXPOSE 8888

RUN ln -sf /dev/stdout /var/log/nginx/access.log && ln -sf /dev/stderr /var/log/nginx/error.log

CMD ["/bin/sh", "-c", "/usr/bin/supervisord -c /home/shiny/system/sup.conf"]


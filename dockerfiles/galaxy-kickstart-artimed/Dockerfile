FROM artbio/galaxy-kickstart-base:ff5fb361875a1220ae3880412083c5cef4a8f1b1

MAINTAINER Marius van den Beek <m.vandenbeek@gmail.com>

COPY  .  /setup
WORKDIR /setup

RUN ansible-playbook -i extra-files/artimed/artimed_inventory.yml -c local galaxy.yml && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy), 9002 (supvisord web app)
EXPOSE :80
EXPOSE :21
EXPOSE :8800
EXPOSE :9002

CMD ansible-playbook galaxy.yml -c local --tags "persists_galaxy,nginx_config,galaxy_config_files" --skip-tags=skip_supervisor_start_in_docker \
           --extra-vars nginx_galaxy_location=$NGINX_GALAXY_LOCATION \
           --extra-vars galaxy_admin=$GALAXY_CONFIG_ADMIN_USERS -i docker_inventory && \
           /usr/bin/python /usr/bin/supervisord -c /etc/supervisor/supervisord.conf --nodaemon


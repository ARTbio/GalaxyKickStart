FROM artbio/galaxy-kickstart-base

MAINTAINER Marius van den Beek <m.vandenbeek@gmail.com>

COPY  .  /setup
WORKDIR /setup

RUN ansible-playbook -i inventory_files/galaxy-kickstart_workflow1 -c local galaxy.yml && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy), 9002 (supvisord web app)
EXPOSE :80
EXPOSE :21
EXPOSE :8800
EXPOSE :9002

CMD ["/startup.sh", "/setup/inventory_files/galaxy-kickstart_workflow1"]

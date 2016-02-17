FROM ubuntu:14.04

MAINTAINER Marius van den Beek <m.vandenbeek@gmail.com>

RUN DEBIAN_FRONTEND=noninteractive  apt-get update  && \
    \
    \
    echo "===> Allow start of services"  && \
    echo "exit 0" > /usr/sbin/policy-rc.d  && \
    \
    apt-get install -y --no-install-recommends \
    apt-transport-https software-properties-common python-dev python-pip build-essential

ONBUILD  RUN  DEBIAN_FRONTEND=noninteractive  apt-get update   && \
              echo "===> Updating TLS certificates..."         && \
              apt-get install -y openssl ca-certificates

RUN pip install --upgrade pip && pip install ansible

COPY  .  /setup
WORKDIR /setup

RUN ansible-playbook -i docker_inventory -c local galaxy.yml && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# work around for AUFS bug
# as per https://github.com/docker/docker/issues/783#issuecomment-56013588
RUN mkdir /etc/ssl/private-copy; mv /etc/ssl/private/* /etc/ssl/private-copy/; rm -r /etc/ssl/private; mv /etc/ssl/private-copy /etc/ssl/private; chmod -R 0700 /etc/ssl/private; chown -R postgres /etc/ssl/private


ENV NGINX_GALAXY_LOCATION=""

ONBUILD  WORKDIR  /setup
ONBUILD  COPY  .  /setup
ONBUILD  RUN  \
              echo "===> Diagnosis: host information..."  && \
              ansible -c local -m setup all

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy), 9002 (supvisord web app)
EXPOSE :80
EXPOSE :21
EXPOSE :8800
EXPOSE :9002

CMD ansible-playbook galaxy.yml -c local --tags "persists_galaxy,nginx_config,galaxy_config_files" --skip-tags=skip_supervisor_start_in_docker \
           --extra-vars nginx_galaxy_location=$NGINX_GALAXY_LOCATION -i docker_inventory && \
           /usr/bin/python /usr/bin/supervisord -c /etc/supervisor/supervisord.conf --nodaemon

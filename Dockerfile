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

COPY  .  /tmp
WORKDIR /tmp

ONBUILD  WORKDIR  /tmp
ONBUILD  COPY  .  /tmp
ONBUILD  RUN  \
              echo "===> Diagnosis: host information..."  && \
              ansible -c local -m setup all

RUN ansible-playbook -c local -i docker_inventory -l "travis_bioblend" --skip-tags=persists_galaxy galaxy.yml

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy), 9002 (supvisord web app)
EXPOSE :80
EXPOSE :21
EXPOSE :8800
EXPOSE :9002

# Mark folders as imported from the host.
VOLUME ["/export", "/var/lib/docker"]

CMD ansible-playbook galaxy.yml -c local --tags persists_galaxy --skip-tags=skip_supervisor_start_in_docker -i docker_inventory && \
           /usr/bin/python /usr/bin/supervisord -c /etc/supervisor/supervisord.conf --nodaemon

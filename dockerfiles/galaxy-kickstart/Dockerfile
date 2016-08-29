FROM artbio/ansible-galaxy-os

MAINTAINER Marius van den Beek <m.vandenbeek@gmail.com>

RUN DEBIAN_FRONTEND=noninteractive  apt-get update  && \
    \
    \
    echo "===> Allow start of services"  && \
    echo "exit 0" > /usr/sbin/policy-rc.d  && \
    \
    apt-get install -y --no-install-recommends \
    apt-transport-https software-properties-common

ONBUILD  RUN  DEBIAN_FRONTEND=noninteractive  apt-get update   && \
              echo "===> Updating TLS certificates..."         && \
              apt-get install -y openssl ca-certificates

# TODO: install as a system package after v2.1.2.0 is released
RUN git clone --recursive -b stable-2.1 https://github.com/ansible/ansible && \
    pip install ansible/

COPY  .  /setup
WORKDIR /setup

ENV LC_ALL=en_US.UTF-8 \
LANG=en_US.UTF-8

RUN locale-gen en_US.UTF-8 && dpkg-reconfigure locales
# work around for AUFS bug
# as per https://github.com/docker/docker/issues/783#issuecomment-56013588
RUN mkdir /var/lib/postgresql-copy; mv /var/lib/postgresql/* /var/lib/postgresql-copy; rm -R /var/lib/postgresql; mv /var/lib/postgresql-copy /var/lib/postgresql; \
    chmod 0700 /var/lib/postgresql; chown -R 1550:1550 /var/lib/postgresql
RUN ansible-playbook -i inventory_files/docker -c local galaxy.yml && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# work around for AUFS bug
RUN mkdir /etc/ssl/private-copy; mv /etc/ssl/private/* /etc/ssl/private-copy/; rm -r /etc/ssl/private; mv /etc/ssl/private-copy /etc/ssl/private; chmod -R 0700 /etc/ssl/private; chown -R postgres /etc/ssl/private


ENV NGINX_GALAXY_LOCATION="" \
GALAXY_CONFIG_ADMIN_USERS="admin@galaxy.org" \
NAT_MASQUERADE=false

ONBUILD  WORKDIR  /setup
ONBUILD  COPY  .  /setup

ADD startup.sh /startup.sh
RUN chmod +x /startup.sh

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy), 9002 (supvisord web app)
EXPOSE :80
EXPOSE :21
EXPOSE :8800
EXPOSE :9002

CMD ["/startup.sh", "/setup/inventory_files/docker"]

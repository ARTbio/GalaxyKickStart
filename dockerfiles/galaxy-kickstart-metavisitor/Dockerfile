FROM artbio/galaxykickstart

LABEL maintainer="Christophe Antoniewski <drosofff@gmail.com>"

COPY  .  /setup
WORKDIR /setup

RUN ansible-playbook -i inventory_files/metavisitor -c local galaxy.yml && \
    apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

CMD ["/startup.sh", "/setup/inventory_files/metavisitor"]

#!/usr/bin/env bash
set -e
sudo /etc/init.d/postgresql stop
sudo apt-get -y --purge remove postgresql libpq-dev libpq5 postgresql-client-common postgresql-common
sudo rm -rf /var/lib/postgresql

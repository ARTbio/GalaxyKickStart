#!/usr/bin/env bash
set -e
apt update -y
apt install -y ansible
ansible --version
git clone https://github.com/ARTbio/galaxykickstart -b Analyse_genomes
cd galaxykickstart/
ansible-galaxy install -r requirements_roles.yml -p roles/ -f
ansible-playbook -i inventory_files/Analyse_genomes galaxy.yml
sleep 15
nginx -s reload
sleep 15
ansible-playbook -i inventory_files/Analyse_genomes galaxy_tool_install.yml
echo "end of deployment\n"

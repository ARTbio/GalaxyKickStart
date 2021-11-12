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
cd /home/galaxy/tool_dependencies
rm -rf _conda
#wget https://mydeepseqbucket.s3.amazonaws.com/_conda_tar.gz
#apt install -y pigz
#pigz -dc _conda_tar.gz > _conda
#chown -R galaxy:galaxy _conda
#rm -rf _conda_tar.gz
cd ~/galaxykickstart
ansible-playbook -i inventory_files/Analyse_genomes galaxy_tool_install.yml
echo "end of deployment\n"

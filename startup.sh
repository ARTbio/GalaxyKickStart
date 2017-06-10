#!/bin/bash
# write out environmental variables to /etc/default/supervisor,
# then source this file (to get both user-specified and default environemntal variables)
# and run ansible to adjust runtime settings. The argument to startup.sh is the location of
# the inventory file to be used
ansible-playbook galaxy.yml -c local \
    --tags "persists_galaxy,nginx_config,galaxy_config_files,galaxy_extras_job_conf,env_vars" --skip-tags=skip_supervisor_start_in_docker \
    --extra-vars nginx_galaxy_location=$NGINX_GALAXY_LOCATION \
    --extra-vars galaxy_admin=$GALAXY_CONFIG_ADMIN_USERS \
    --extra-vars ftp_upload_site=$IP_ADDRESS \
    --extra-vars nat_masquerade=$NAT_MAQUERADE \
    -i "$1" && \
    source /etc/default/supervisor && \
    exec /usr/bin/supervisord -c /etc/supervisor/supervisord.conf --nodaemon

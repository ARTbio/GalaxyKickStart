#!/bin/bash

export DOCKER_INSTANCE=`docker run -d -p 80:80 -p 21:21 -p 8800:8800 \
          --privileged=true \
          -e GALAXY_CONFIG_ALLOW_USER_DATASET_PURGE=True \
          -e GALAXY_CONFIG_ALLOW_LIBRARY_PATH_PASTE=True \
          -e GALAXY_CONFIG_ENABLE_USER_DELETION=True \
          -e GALAXY_CONFIG_ENABLE_BETA_WORKFLOW_MODULES=True \
          -v /mnt/galaxy_tmp:/tmp \
          -v /mnt/galaxy_export:/export \
          artbio/analyse_genomes:2021`

echo "Warning: export of data folder may take more than 15 min"
echo "Ctl-C to interrupt logging of the docker instance"
docker logs -f $DOCKER_INSTANCE

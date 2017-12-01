#!/usr/bin/env bash
set -e
sleep 60s
docker logs $STANDARD
docker ps
curl --fail $BIOBLEND_GALAXY_URL/api/version
docker exec -it $STANDARD supervisorctl status
echo "Starting Bioblend tests\n"
sudo -E su $GALAXY_TRAVIS_USER -c "export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
  cd $GALAXY_HOME && \
  bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy*.py"
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://localhost:8021 --user $GALAXY_USER:$GALAXY_USER_PASSWD
docker stop $STANDARD && docker rm $STANDARD
cd $TRAVIS_BUILD_DIR

#!/usr/bin/env bash
set -e
sleep 60s
docker logs $STANDARD
curl --fail $BIOBLEND_GALAXY_URL/api/version
#curl --fail http://localhost:8181/subdir/api/version| grep version_major
docker exec -it $STANDARD supervisorctl status
#docker exec -it $CUSTOM supervisorctl status #| grep proftpd | grep RUNNING
printenv
#sudo -E su $GALAXY_TRAVIS_USER -c "export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
#  cd $GALAXY_HOME &&
#  bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy*.py"
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://localhost:21 --user $GALAXY_USER:$GALAXY_USER_PASSWD
docker stop $STANDARD && docker rm $STANDARD
#CID3=`docker run -d --privileged=true -p 8181:80 -e NAT_MASQUERADE=true -v /export2:/export galaxy_kickstart` && sleep 60s
#docker logs $CID3
#curl http://localhost:8181/api/version| grep version_major
cd $TRAVIS_BUILD_DIR

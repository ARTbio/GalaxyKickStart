#!/usr/bin/env bash
set -e
echo -e "sleeping 120s, zzzzzz"
sleep 120s
# docker logs $CID2
#export $BIOBLEND_GALAXY_URL="http://localhost:8080/subdir"
curl http://localhost:8181/subdir/api/version| grep version_major
sudo -E su $GALAXY_TRAVIS_USER -c "export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
  export $BIOBLEND_GALAXY_URL=http://localhost:8080/subdir &&
  cd $GALAXY_HOME &&
  bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy*.py"
curl --fail $BIOBLEND_GALAXY_URL/api/version
date > $HOME/date.txt && curl --fail -T $HOME/date.txt ftp://localhost:8021 --user $GALAXY_USER:$GALAXY_USER_PASSWD
curl --fail ftp://localhost:8021 --user $GALAXY_USER:$GALAXY_USER_PASSWD
docker exec -it $CID1 supervisorctl status | grep proftpd | grep RUNNING
docker stop $CID1 $CID2 && docker rm $CID1 $CID2
CID3=`docker run -d --privileged=true -p 8181:80 -e NAT_MASQUERADE=true -v /export2:/export galaxy_kickstart` && sleep 120s
docker logs $CID3
curl http://localhost:8181/api/version| grep version_major
cd $TRAVIS_BUILD_DIR

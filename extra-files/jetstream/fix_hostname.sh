#!/bin/bash
service hostname restart
HN=$(cat /etc/hostname)

cat<<EOF > /etc/hosts
# This file was automatically created by fix_hostname.sh script
127.0.0.1 localhost $HN

# The following lines are desirable for IPv6 capable hosts
::1 ip6-localhost ip6-loopback
fe00::0 ip6-localnet
ff00::0 ip6-mcastprefix
ff02::1 ip6-allnodes
ff02::2 ip6-allrouters
ff02::3 ip6-allhosts
EOF

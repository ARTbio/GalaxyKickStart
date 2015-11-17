#!/bin/bash
echo "To be sure that the current user can ssh to the localhost without password" 
if cat /dev/zero | ssh-keygen -q -N ""
then
	echo "new key created"
	cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
else
	echo "Old key exists"
        if [ -r ~/.ssh/authorized_keys ]
        then
           echo "Authorized key exists"
	   k=`grep -c -f ~/.ssh/id_rsa.pub ~/.ssh/authorized_keys`
        else
           k=0
        fi
	if [ "$k" == "0" ]
	then
		echo "Old key added"
		cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
	fi
fi

if ssh localhost "echo ''"; 
then 
	echo "This user has permission to ssh directly to the localhost without any user input."; 
else 
	echo "This user does not have permission to ssh directly to the localhost without any user input, so galaxy install will not work.";
	exit 1
fi

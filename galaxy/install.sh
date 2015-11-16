#!/bin/bash
#Parsing parameters 
if [ "$1" == "-h" ]
then
	echo "Just 4 optional parameters: the git repository, galaxy user, galaxy port and ftp port."
	exit 0
fi

ftp_port=$4
galaxy_port=$3
galaxy_user=$2
artimed_git_repo=$1
if [ "$artimed_git_repo" == "" ]
then
	artimed_git_repo="--branch supervisor https://github.com/fabiorjvieira/ansible-artimed.git"
fi

OS=`head -n1 /etc/issue | cut -d " " -f 1`
if [ "$OS" == "Debian" ]
then
	IM="apt"
elif [ "$OS" == "Ubuntu" ]
then
	IM="apt"
else
	echo "Not a Debian like system. Cannot determine the Operation System, so aborting."
fi

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
	echo "This user does not have permission to ssh directly to the localhost without any user input, so this script cannot continue.";
	exit 1
fi

#install Dependencies
echo "Minimum requirements: OpenSSH client, Ansible >=1.8 and git (this script will try to install these)."
if [ "$IM" == "apt" ]
then
	sudo apt-get install software-properties-common -y
	sudo apt-add-repository ppa:ansible/ansible -y
	sudo apt-get update -y
	sudo apt-get install git openssh-client openssh-server ansible -y
	sudo service ssh start
fi

#Installation
if git clone --recursive $artimed_git_repo
then
	cd ansible-artimed/galaxy/
	hostIP=`hostname -I | cut -d " " -f 1`
	FTP_PORT=$ftp_port GALAXY_USER=$galaxy_user GALAXY_PORT=$galaxy_user INSTALL_HOSTNAME=$hostIP ansible-playbook -i "localhost," galaxy.yml -vvvv
	echo "Wait until galaxy provide the web service on http://localhost"
	echo -n "Press control+c to stop here or enter key to install Galaxy tools..."
	read
	GALAXY_USER=$galaxy_user GALAXY_PORT=$galaxy_port ansible-playbook -i "localhost," tools.yml -vvvv
fi

When you are done with your own Metavisitor Galaxy instance installation either using [ansible](metavisitor_ansible.md) or [docker](metavisitor_docker), there are a few basic things to know for web access and basic server admin operations

## 1. Connect to your Metavisitor Galaxy instance with your web browser


## 2. Log to the Galaxy server using the admin credentials:

login: admin@galaxy.org
password: admin

Note: you may (should) change the password at any time as soon as you are logged as admin@galaxy.org


But before doing this, we have to perform a low level Galaxy admin task: restart the Galaxy server instance so that the dbkey and the fasta genomes that we've just created for the server are registered and seen by the tools.

Depending on your skill level, you have the simple, dirty way:

- reboot the machine where the galaxy server is running !

or the clean, freaking (for non-Geek normal biologists) way:

- Connect to the server where the Galaxy instance has been installed either through an ssh connection or using your local terminal. and type

`sudo supervisorctl restart galaxy:`

If everything went fine you should see in your terminal
```
# supervisorctl restart galaxy:
galaxy_web: stopped
handler0: stopped
handler1: stopped
handler0: started
handler1: started
galaxy_web: started
```
 freaking insn't it ?


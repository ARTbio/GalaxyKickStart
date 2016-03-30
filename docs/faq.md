What is the username and password of the galaxy admin account ?
----

Username and password of the galaxy account are controlled by the variables `galaxy_admin` and `galaxy_admin_pw` and
default to `admin@galaxy.org` and `admin` (Defaults are defined in group_vars/all). This should be changed in the group or host variables for the host you are working on.
If you have a host in the `mygroup` group, you can edit group_vars/my_group and set
```
galaxy_admin: new_admin@email.com
galaxy_admin_pw: new_password
```

As with each change, run the playbook again.
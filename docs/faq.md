Why does the playbook fail?
----

Make sure that you are on ansible version >=2.1.
You can check your ansible version by typing:

```
ansible --version
```

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


How can I set up GalaxyKickStart behind a proxy?
----

Many commandline utilities can be configured to use a proxy by setting the
`http_proxy` and `https_proxy` environment variables. Tasks launched by ansible
will only see these environment variables if ansible sets these variables for
the task. We have included a global `proxy_env` variable in the galaxy.yml playbook.
You can set the content of this variable in your inventory or group variables 
(See [Customizing the playbook](customizations.md) for details on how to define variable).
To use the proxy at http://proxy.bos.example.com:8080 define the variable `proxy_env` like so:

```
proxy_env:
  http_proxy: http://proxy.bos.example.com:8080
  https_proxy: http://proxy.bos.example.com:8080
  no_proxy: localhost,127.0.0.0,127.0.1.1,127.0.1.1,local.home
```

Adresses that should not be contacted through a proxy should be listed in the `no_proxy` variable.
An example can be found in group_vars/proxy.

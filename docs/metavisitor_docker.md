# Installing Metavisitor with Docker
----

We distribute a docker image of Metavisitor, which can thus be used to run a Metavisitor docker container. For a quick start, go directly to the last section "Persisting to disk".

## Requirements
You need to have docker installed and configured for your user.

## Running images from the dockerhub
You can search for pre-built docker images from the dockerhub by typing in the terminal of the machine where you want to run the docker container:

```
docker search metavisitor
```

Then, to get the docker image, type:

```
docker pull artbio/metavisitor-1.2
```
In this documentation, we recommend to use the `artbio/metavisitor-1.2` which better corresponds to the environment described in our [Metavisitor preprint](http://dx.doi.org/10.1101/048983)

When this pull is done (may take a few minutes depending on your connection speed to the dockerhub), you can start the container by typing:

```
docker run -d -p 80:80 artbio/metavisitor-1.2
```

This command starts a container in daemon mode (`-d`) from the image and serve it on port 80 of the local machine in the standard docker way.

`-p 80:80` forwards requests to nginx inside the container running on port 80. If you want to access the machine hosting the running container through another port (for instance 8080), just change `-p 80:80` to `-p 8080:80`

## Runtime changes to pre-built docker images

If you wish to reach the container on a subdirectory, add `-e NGINX_GALAXY_LOCATION="/my-subdirectory"` to the docker call.

For instance,
```
docker run -d -e NGINX_GALAXY_LOCATION="/my-subdirectory" -p 80:80 artbio/metavisitor-1.2
```

will get the metavisitor docker container serving at `http://127.0.0.1:80/my-subdirectory`.

We recommend also changing the default admin user as well, so the command becomes:
```
docker run -d -e NGINX_GALAXY_LOCATION="/my-subdirectory" -e GALAXY_CONFIG_ADMIN_USERS=admin@artbio.fr -p 80:80 artbio/galaxy-kickstart-base
```
Note that is you do not make this latest change, the admin login for the metavisitor container is by default `admin@galaxy.org` and the password is `admin`.

## Persisting to disk

All changes made to a docker container are by default ephemeral; if you remove the container, the changes are gone.
To persist data (this includes the postgresql database, galaxy's config files and your user data), mount a Volume into
the containers /export folder.
Due to the persistance mechanism (we use bind-mounts inside the container), you need to privilege the container.
Thus, assuming you would like to mount your local `/my/data` folder and persist you Galaxy data in this folder, run
```
docker run -d --privileged -v /my/data:/export -p 80:80 artbio/metavisitor-1.2
```
This will run through the persistence tags of the galaxy.yml and export the required files to /export (now on your machine's /my/data).
From the new location the files are being bind-mounted back into their original location.

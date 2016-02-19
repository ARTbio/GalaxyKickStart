
#Building and deploying galaxy-kickstart in docker
----

## Requirements
You need to have docker installed and configured for your user.

The repository comes with various Dockerfiles that can be used to configure a deployment using Docker,
or you can start with a pre-built docker image.

## Running images from the dockerhub
You can obtain a pre-built docker image from the dockerhub:
```
docker pull artbio/galaxy-kickstart-base
```

Start the image and serve it on port 8080 of your local machine in the standard docker way:
```
CID=`docker run -d -p 8080:80 artbio/galaxy-kickstart-base`
```
`-p 8080:80` will forward requests to nginx inside the container running on port 80.

## Runtime changes to pre-built docker images

If you wish to reach the container on a subdirectory, add `-e NGINX_GALAXY_LOCATION="/my-subdirectory"` to the docker call 
and galaxy will be served at `http://127.0.0.1:8080/my-subdirectory`.

We recommend changing the default admin user as well, so the command becomes:
```
CID=`docker run -d -e NGINX_GALAXY_LOCATION="/my-subdirectory" -e GALAXY_CONFIG_ADMIN_USERS=admin@artbio.fr -p 8080:80 artbio/galaxy-kickstart-base`
```

## Commit changed containers to new images

As with standard docker containers, you can change, tag and commit running containers when you have configured them to your liking:
Commit the changes to my-new-image
```
docker commit $CID my-new-image

```
Stop and remove the original container:
```
docker stop $CID && docker rm $CID
```
Start the new container:
```
CID=`docker run -d -e NGINX_GALAXY_LOCATION="/my-subdirectory" -e GALAXY_CONFIG_ADMIN_USERS=admin@artbio.fr -p 8080:80 my-new-image`
```

## Persisting to disk

All changes made to the container are by default ephemeral; if you remove the container, the changes are gone.
To persist data (this includes the postgresql database, galaxy's config files and your user data), mount a Volume into
the containers /export folder.
Due to the persistance mechanism (we use bind-mounts inside the container), you need to privilege the container.
Assuming you would like to mount your local `/data` folder, run
```
CID=`docker run -d --privileged -v /data:/export -p 8080:80 my-new-image`
```
This will run through the persistence tags of the galaxy.yml and export the required files to /export (now on your machine's /data).
From the new location the files are being bind-mounted back into their original location.

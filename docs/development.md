Fixing submodule issues
----

This project is using roles from different places, and these are included as git submodules. 
In a development cycle it may occur that the url of a repo changes, or that we have removed certain commits (We aim to avoid this). 
If at any point you find yourself in the situation where a git checkout of submodules does not work, you can follow this procedure: 
(galaxyprojectdotorg.galaxy-extras as an example, replace this with the problematic role)

Skip the git clone if you already cloned the repo.

```
git clone --recursive https://github.com/ARTbio/GalaxyKickStart.git
cd GalaxyKickStart/
```

Enter the directory of the problematic submodule

```
cd roles/galaxyprojectdotorg.galaxy-extras/
```

Verify that you have the correct remote.
You can find the correct remote in GalaxyKickStart/.gitmodules
```
git remote -v
```

If this is not the correct remote, change back to the GalaxyKickStart directory and do

```
git submodule sync
```

Get the latest changesets (change back to the submodule directory if you had to sync the submodules)
```
git fetch
```

```
Check which branch you are on. You want to checkout the master branch.

```
git checkout master
git add galaxyprojectdotorg.galaxy-extras
git commit -m "Fix galaxy-extras submodule revision"
```

If you push this work back to github and you compare your repo with GalaxyKickStart/master, you should not see any changes
that look like `Subproject commit 8265bdceaf0eca5ea4daeda06eb0d28583754ef6` (unless you have intentionally updated a submodule).

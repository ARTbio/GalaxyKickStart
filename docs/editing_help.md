To view the help, type
```bash
mkdocs serve
```
This will serve the documentation on http://127.0.0.1:8000.

To add new files, create a markdown document in docs/, and reference it in mkdocs.yml.  
If you want to add a document called `editing_help.md` add the following line.
```
 - Editing the readme: editing_help.md
```

so that [mkdocs.yml](https://github.com/ARTbio/ansible-artimed/blob/master/mkdocs.yml) looks like this.
```
site_name: GalaxyKickStart
pages:
 - Home: index.md
 - What is GalaxyKickStart: about.md
 - Getting started: getting_started.md
 - Customizations: customizations.md
 - Examples: examples.md
 - Available roles: available_roles.md
 - Available variables: available_variables.md
 - Editing the readme: editing_help.md
theme: readthedocs
```

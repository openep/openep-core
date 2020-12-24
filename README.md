# Carto Data Reader

CartoDataReader supports the import of Carto3 studies into Matlab

## Instructions

### Importing a Carto dataset

First, export the data from Carto using File>Export !!! fix this

Then, unzip the resulting file on your local machine, and at the Matlab command prompt, type:

```
userdata = importcarto_mem
```
This will open a dialog asking you to select the .xml file in the unzipped directory. The result will be given in the variable userdata and the option will be given to save the data to file. This is always recommended.

### Adding Visitag data

To import just the Visitag data from a Carto export, type:

```
visitag = importvisitag()
```
This will open a dialog asking you to select the ablationsites.txt file that is located within the Wisetag subdirectory.

It is also possible ammend the Visitag data to the userdata structure as follows:
```
userdata = importvisitag('openfile', userdata)
```
This will return a new data structure with the Visitag data correctly projected onto the shell.

### Plotting activation maps
```
drawMap(userdata, 'type', 'act')
```
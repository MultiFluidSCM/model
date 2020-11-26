# Multi-Fluid Single Column Model
Two-Fluid Single Column Model created by John Thuburn and based on the Met Office's EndGAME dynamical core.

## Installation
- Ensure git is installed on your machine: [Git installation guide.](https://github.com/git-guides/install-git)
- Create a MultiFluidSCM folder
- In the MultiFluidSCM folder, get this and the other required repositories using:
```
git clone "https://github.com/MultiFluidSCM/model/"
git clone "https://github.com/MultiFluidSCM/plots/"
git clone "https://github.com/MultiFluidSCM/test_cases/"
```

## Usage of the model
- Run the [model](https://github.com/MultiFluidSCM/model/) using ```convection.m```
- Rename and copy the SCM data ```SCM_results.mat``` into your personal [test_cases](https://github.com/MultiFluidSCM/test_cases/) folder under ```scm_data```
- In ```get_settings.m```, edit ```settings.scm_data``` to be the filename of your SCM data set
- Adjust your plotting preferences in ```get_settings.m``` (more options to be added in the future)
- Run ```plot_profiles.m```
- If the relevant options are selected, ```.fig``` files will be saved to the ```figures/``` folder in your test case directory and rendered images in the ```images/``` folder.

## Version control using git
NOTE: Version control should currently only be used for the [plots](https://github.com/MultiFluidSCM/test_cases/) and [test_cases](https://github.com/MultiFluidSCM/test_cases/) repositories.

### Download the latest version of the code
Get the latest versions of a repository by navigating to your local repository (e.g. /path/to/MultiFluidSCM/plots/) and using
```
git pull
```

### Upload your code to merge with online version
If you have made changes to the scripts which you would like to merge with the master code, do the following:
- Add all modified files in the current directory using ```git add .```
- Commit the changes with a description of the changes using ```git commit -m "Description of changes to code"```
- Upload and merge the changes using ```git push -u origin main```

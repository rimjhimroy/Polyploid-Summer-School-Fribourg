# Polyploid Summer School 2023
10-20 July
University of Fribourg

This is the Renku project for the [Polyploid Summer School 2023](https://events.unifr.ch/summerschool_polyploidy/en/). Renku project is basically a git repository with some bells and whistles. You'll find we have already created some useful things like `data` and `notebooks` directories and a `Dockerfile`.


## Quick Renku instructions:

1. Go to [https://renkulab.io/](https://renkulab.io/). Click on the 'Login or Sign Up' button and make yourself an account if you do not have one, either via 'Register', a 'SWITCH edu-ID' or via a GitHub account.
2. Once logged in, go to the Polyploid Summer School "project":
[https://renkulab.io/projects/rimjhim.choudhury/polyploid-summer-school-2023](https://renkulab.io/projects/rimjhim.choudhury/polyploid-summer-school-2023)
3. Near the top right, click on the button to fork the project; this should bring up a box where you can modify the name (you don't need to) and you can click on the 'Fork' button, putting you on the landing page of your (new) project.
4. Click on the environments tab, then 'New' to open a new environment. Hopefully, it shows "Docker image available"; if not, perhaps it says "Docker Image building", which is also ok (it might take a couple of minutes to build and then turn to "Docker image available"). For Python, command line work, and for people who are fan of VSCode, select /lab; for R work in Rstudio, select, /rstudio (you can always switch back and forth). Leave the number of CPUs at 0.25 and Memory at 1G. Click 'Start environment'. This may take a couple minutes to boot up. 

Note: you can always go back and forth between the `/rstudio`, `/lab`, and `/vscode` environments by modifying the end of web link, e.g., if you started in /rstudio, but want to move to /lab, you can remove "rstudio", type "lab" in its place and then press enter.

Note: you can upload files from the local computer into your environment via the left tab of the /lab/ environment.


## Working with Renku

Please see [the documentation](https://renku.readthedocs.io/en/latest/) for 
more details about Renku.


## Working with the project

The simplest way to start your project is right from the Renku
platform - just click on the `Sessions` tab and start a new session.
This will start an interactive environment right in your browser.

To work with the project anywhere outside the Renku platform,
click the `Settings` tab where you will find the
git repo URLs - use `git` to clone the project on whichever machine you want.

### Changing interactive session dependencies

Initially we install a very minimal set of packages to keep the images small.
However, you can add python and conda packages in `requirements.txt` and
`environment.yml`, and R packages to `install.R` (listed as, for example,
`install.packages("ggplot2")`), to your heart's content. If you need more fine-grained
control over your environment, please see [the documentation](https://renku.readthedocs.io/en/stable/topic-guides/customizing-sessions.html).

## Project configuration

Project options can be found in `.renku/renku.ini`. In this
project there is currently only one option, which specifies
the default type of environment to open, in this case `/rstudio`.

## Moving forward

Once you feel at home with your project, we recommend that you replace
this README file with your own project documentation! Happy data wrangling!

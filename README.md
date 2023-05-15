# BOO 2023

In this short guide, we will introduce the Noteable environment and show how to install the [QSPRpred](https://github.com/CDDLeiden/QSPRPred) package.

## Noteable

You can access the Noteable environment from the [Brightspace module](https://brightspace.universiteitleiden.nl/d2l/le/lessons/190680/topics/2368865) of the course. Open the link and select "Standard (Collaborative session)". This will launch a [JupyterLab](https://jupyter.org/) server instance and redirect you to it. You can look at it as your personal virtual computer running the Linux operating system. JupyterLab is a development environment that enables you to write and execute code on this machine. This is facilitated either with simple scripts or with interactive notebooks (Jupyter Notebooks). The latter are documents that contain both code and text, and allow you to run the code and see the results in the same document. This is what we will use in this course most often.

## Using the Terminal

We will use the shell (terminal) to install the QSPRpred package, and it is also the best way to execute programs and manage our files. You can open a new terminal by clicking on the `+` sign in the left sidebar and selecting "Terminal". You can also open a terminal by clicking on the "File" menu and selecting "New" -> "Terminal". Or you can use the launch menu which you should see when you open the JupyterLab server. Click the black icon and a terminal will open. You can type commands in the terminal. Try these:

```bash
ls # list files in the current directory
ls -la # list all files in the current directory, including hidden files, and show more information
pwd # print the current working directory
```

The output of the `pwd` command should be `/home/jovyan`. This is your home directory. You can create new directories and files here, and you can also delete them. However, you can also do so in the GUI via the file view on the left side of the screen. Try to create a directory to see it appear in the file view.

```bash
mkdir test # create a new directory called "test"
```

You can open and create files in the GUI, but you can also do so in the terminal. Try to create a new file called `test.txt` in the `test` directory.

```bash
cd test # change directory to "test"
touch test.txt # create a new file called "test.txt"
```

If you need to go back to the previous directory, you can use the `cd ..` command. Try to go back to your home directory.

```bash
cd .. # go back to the previous directory
```

You should now be back in your home directory:

```bash
pwd # print the current working directory
```

Your home directory can also be abbreviated as `~`. List the home directory to see the `test` directory you created earlier:

```bash
ls ~ # list the home directory
```

The command line is a very powerful tool, and you can access basically anything on the system from it. For example, try to list hardware information about your computer with the following commands:

```bash
egrep "proc|vendor|model|MHz" /proc/cpuinfo
grep "Mem" /proc/meminfo
```

It is useful to know how to use the command line, see [this cheatsheet for more examples](https://cheatography.com/davechild/cheat-sheets/linux-command-line/), but let's move on to installing software and eventually the QSPRpred package.

## Installing Software with Anaconda

We like to use the [Anaconda](https://www.anaconda.com/) Python distribution to install software. It is a free and open-source distribution of software packages that comes with a package manager called `conda`. It is already available on our machine:

```bash
conda --version # check the version of conda
```

Conda operates with so called 'environments'. In fact, we are already operating within the `base` environment. You can see this in the terminal prompt. You can list currently available environments with the following command:

```bash
conda env list # list available environments
```

It is good practice to create a new environment for each project. This way, you can install different versions of the same software in different environments, and you can also easily share your environment with others. Let's create a new environment called `qsprpred`. Because notable is a bit special and only persists your home directory across virtual machines, we need to specify the path where we want to create the environment:

```bash
conda create python=3.10 --prefix ~/software/envs/qsprpred
```

This will create a new environment called `qsprpred` in the `~/software/envs` directory. You can activate the environment with the following command:

```bash
conda activate ~/software/envs/qsprpred
```

You should see the environment name change in the terminal prompt. We only installed the Python interpreter in this environment, but we can install many new programs this way. For example, try to install the `htop` utility:

```bash
conda install htop
```

You can now run the `htop` program to see a list of running processes and usage of resources on your machine:

```bash
htop
```

Pretty neat, right? You can exit the program with the `q` key. Let's install the QSPRpred package now:

```bash
pip install git+https://github.com/CDDLeiden/QSPRPred.git@v2.0.0.dev1 # installs a version of the code by specifying the tag (v2.0.0.dev1)
```

This will use the pip installer to fetch a version of the QSPRpred package from GitHub and install it in your environment. You can check by importing the package in Python:

```bash
python -c "import qsprpred; print(qsprpred.__version__)" # execute Python code from the command line
```

## Running the QSPRpred Tutorial

You can browse the [source code of QSPRpred on GitHub](https://github.com/CDDLeiden/QSPRPred/tree/v2.0.0.dev1/) and the [tutorial code](https://github.com/CDDLeiden/QSPRPred/tree/v2.0.0.dev1/tutorial) is there as well. You can browse the notebooks in this directory, but you cannot run or edit them on GitHub. That is why we need to create our own copy of the code. We can do this by cloning the repository:

```bash
git clone -b v2.0.0.dev1 https://github.com/CDDLeiden/QSPRPred.git ~/software/QSPRPred
```

However, before we navigate to the tutorial directory, we need to link the `qsprpred` environment to JupyterLab. We can do this by installing the `ipykernel` from this environment to our home directory:

```bash
python -m ipykernel install --prefix=$HOME/.local/ --name 'Python-QSPRpred'
```

**WARNING:** Make sure to execute this with the `python` in the `qsprpred` environment.

This will enable a new kernel in JupyterLab that we can use to run notebooks. You should be able to see the new option after clicking the large `+` sign in the left sidebar. If it does not appear, refresh your browser. You will now have this kernel option every time you open a Noteable session.

Now we can navigate to the tutorial directory with the file browser and run and edit the notebooks with QSPRpred examples.
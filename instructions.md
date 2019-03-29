# Set Up Instructions:

Contact Dr. Gallant (jgallant@msu.edu) with difficulties getting computers set up!

To enable courses to be taught as effectively as possible, you will need to be able to log in to the HPCC resources without difficulty before the class.  If you do not have an HPCC account, please install the appropriate software before the workshop but you will not be able to test that it works before the class. Please try to come 10-15 minutes early to get a temporary account, and test that you can log in.

## Connecting to HPCC
For Mac users:
Please install the latest version of Xquartz (you can follow the directions here: https://wiki.hpcc.msu.edu/display/hpccdocs/Installing+an+X-server++for+Macs)

## To test that it is working properly:
Open a terminal window (if you do not have the terminal program on your taskbar, you can search for terminal), and
Type:

```bash
ssh -XY [classx]@hpcc.msu.edu #replace classx with the ID of your temporary account [classx]
# type your the password provided for the classx account
```
(type your msu password or the password provided for the temporary account)  This should give you a command prompt.  If you type:

```bash
xeyes
```
You should get a separate window that displays a pair of eyes that follow the mouse around.  If this works, that should be all you need to be ready for the course.  If not, please come several minutes early and we will try to help you get online.


## For WINDOWS users:
We recommend installing Moba Xterm "home edition" (http://mobaxterm.mobatek.net/download-home-edition.html).   Please download the "Home Edition, Installer Version unless you are comfortable installing portable applications yourself.  This program provides both an SSH client (command line) and an X-Windows server system for running Graphical User Interface programs remotely (X11).

Once it is installed, you can open it, and try to connect with the hpcc:

To test that it is working properly:
Open a terminal window (if you do not have the mobaxterm shortcut set up, you can search for mobaxterm), and type:

```bash
ssh -XY yournetid@hpcc.msu.edu # replace yournetid with YOUR ACTUAL NETID or the ID of your temporary account [classx]
# type your the password provided for the classx account
```
This should give you a command prompt.  If you type:

```bash
xeyes
```

You should get a separate window that displays a pair of eyes that follow the mouse around.

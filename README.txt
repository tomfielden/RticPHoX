RticPHoX
========

The Arctic Fox is a tool built to investigate the feasibility of combining two or more random variable without sampling.

sudo apt-get install jython
sudo apt-get install libcommons-math-java
sudo apt-get install ant
sudo apt-get install junit
sudo apt-get install latexml
sudo apt-get install libglpk-java
sudo apt-get install gnuplot

## Get the Jave Developer Kit from Oracle ##
tar zxvf jdk-7u17-linux-i586.tar.gz
sudo mv ~/Downloads/jre-6u31-linux-x64.bin /opt/java/64

sudo update-alternatives --install /usr/bin/java java /usr/local/java/jdk1.7.0_17/bin/java 1
sudo update-alternatives --install /usr/bin/javac javac /usr/local/java/jdk1.7.0_17/bin/javac 1
sudo update-alternatives --set java /usr/local/java/jdk1.7.0_17/bin/java
sudo update-alternatives --set javac /usr/local/java/jdk1.7.0_17/bin/javac
java -version
javac -version

## Assuming you're in the working directory,
tar jxvf ./External/JavaPlot-0.4.0.tar.bz2 -C ~/

NB: The build.xml file must have it's own classpath <fileset>'s

## Here's an example .bashrc fragment #########################################################
export JAVA_HOME=/usr/local/java/jdk1.7.0_17
export JYTHON_HOME=/usr/share/jython
export COMMONS_HOME=/usr/share/java
export JAVAPLOT_HOME=/home/tomf/JavaPlot-0.4.0
export PHOX=/home/tomf/work/RticPHoX
export CLASSPATH=.:$JYTHON_HOME/jython.jar:$PHOX/dist/RticPHoX.jar:/usr/share/java/glpk-java.jar:/usr/share/java/junit.jar:$COMMONS_HOME/commons-math-2.2.jar:$JAVAPLOT_HOME/dist/JavaPlot.jar
export PATH=$PATH:$JYTHON_HOME:$PHOX
################################################################################################

## Now, ensure that the shell script "phox" is executable. Should be in your PATH.
chmod +x phox

## Finally invoke the RticPHoX with
phox

## Test out RticPHoX,
N = Normal(2,5)
Plot(N)


#### To build the documentation the following LaTeX is needed. #################################
sudo apt-get install texlive-latex-base
sudo apt-get install texlive-latex-recommended
sudo apt-get install texlive-fonts-recommended


#!/bin/sh

# maketar should substitute these
ALLINEA_SUITE_FULL_NAME="Allinea Forge"
ALLINEA_SUITE_DIRECTORY=forge
ALLINEA_SUITE_TARBALL_NAME=forge

if [ "$(uname)" = "AIX" ] ; then
    TOOLS_RUN_ARCH="$(uname -p)"
else
    TOOLS_RUN_ARCH="$(uname -m)"
fi

INSTALL_USER="$(whoami)"

TOOLS_BUILD_ARCH=x86_64

accept_licence=
if [ "$1" = "--accept-licence" -o "$1" = "--accept-license" ]; then
	accept_licence=1
	shift
fi

if [ $# -gt 1 ]; then
	echo "Usage: $(basename $0) [--accept-licence] [installation-directory]"
	exit 1
fi

DIRECTORY=
if [ -n "$1" ]; then
	DIRECTORY="$1"
fi

if [ "${TOOLS_RUN_ARCH}" != "${TOOLS_BUILD_ARCH}" -a "${TOOLS_BUILD_ARCH}" != k1om ] ; then
    echo
    echo "Your machine uses the $TOOLS_RUN_ARCH architecture but this ${ALLINEA_SUITE_FULL_NAME} installation tarball is for the $TOOLS_BUILD_ARCH architecure."
    echo
    echo "We recommend you download the $TOOLS_RUN_ARCH architecture tarball from the website instead."
    echo
    echo "Do you wish to continue (y/N)?"
    echo
    read answer
    case $answer in
        [Yy]*) ;;
        *) exit ;;
    esac
fi

if(test ${INSTALL_USER} = root)
then
	DEFAULT=/opt/allinea/${ALLINEA_SUITE_DIRECTORY}
else
	DEFAULT=$HOME/allinea/${ALLINEA_SUITE_DIRECTORY}
fi

cd `dirname $0`
INSTALL_FROM_DIR=`pwd`

echo "${ALLINEA_SUITE_FULL_NAME} (c) Allinea Software 2002-2015"
echo ===============================================================
echo
echo Version `cat ./.version`
echo
echo This is the install script for ${ALLINEA_SUITE_FULL_NAME}.
echo
echo By proceeding with this installation you agree to be bound
echo by our licensing conditions which will now be presented.  
echo These may be found in the file README.LICENCE which was
echo extracted from the same archive file as this script. 
echo

if [ -n "$accept_licence" ]; then
	cat ./README.LICENCE
	echo
	echo "Licence accepted..."
else
	echo Press Enter to read the licence.
	read KEY

	more ./README.LICENCE

	while [ "$KEY" != "a" ]; do
	echo "Select (a) to AGREE to terms of this licence"
	echo "Select (q) if you do not agree to the terms of this licence"
	read KEY

	case $KEY in
		'a') ;;
		'q') exit 1
	esac
	done
fi

echo

if [ -z "$DIRECTORY" ]; then
	echo "You now need to select a directory in which to install ${ALLINEA_SUITE_FULL_NAME}."
	echo "This directory must be accessible on all the nodes in your"
	echo "cluster."
	echo

	echo "The default install base of ${ALLINEA_SUITE_FULL_NAME} is $DEFAULT"
	echo "Enter another directory or press Enter to use default: "
 	read DIRECTORY
	if [ "$DIRECTORY" = "" ]; then
	    DIRECTORY=$DEFAULT
	fi
fi

if  (test -d "$DIRECTORY") ;then
    echo 'Directory exists - continuing'
elif (test -e $DIRECTORY) ;  then
    echo $DIRECTORY exists and is not a directory, aborting install
    exit
else mkdir -p "$DIRECTORY" || {
	echo Unable to create directory, aborting install;
	exit
	}
    echo 'Made directory - continuing'
fi

cd "$DIRECTORY"
gunzip -dc "$INSTALL_FROM_DIR/${ALLINEA_SUITE_TARBALL_NAME}.tgz" |
tar --no-same-owner --no-same-permissions -xvf - ||  {
    echo '** An error was encountered, aborting install'
    exit
}

chmod -R a+rx bin
chmod -R a+rx libexec
chmod -R a+r doc
# we do not supply examples with Licence Server
if [ -d examples ]; then
    chmod -R a+r examples
fi
# we do not supply help files with the Xeon Phi build
if [ -d help ]; then
    chmod -R a+r help
fi

if [ -x /bgsys -a -n "$MMCS_SERVER_IP" ]; then
    cp /usr/X11R6/lib/libX11.so.6 lib/
    cp /usr/X11R6/lib/libXext.so.6 lib/
fi

echo 'Finishing installation. This may take a few minutes...'

if [ -f ./map-mime.xml ]; then
    xdg-mime install map-mime.xml 2>/dev/null
    rm ./map-mime.xml
fi

if [ -x ./post-install.sh ]; then
    ./post-install.sh
    rm ./post-install.sh
fi

echo
echo
echo Please read the README file in the
echo $DIRECTORY/doc directory
echo
echo 'IMPORTANT'
echo '========='
echo
echo 'You now need to visit '
echo 'http://www.allinea.com/products/trials'
echo 'to obtain a combined trial licence for Allinea DDT, Allinea MAP and'
echo 'Allinea Performance Reports.'
echo -----------------------------------------------------
echo 'support@allinea.com      +44 (0) 1926 623231'
echo =====================================================

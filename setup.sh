#Make sure 3DSTNeutrons finds plugin libraries at runtime
export LD_LIBRARY_PATH="`pwd`/`dirname $BASH_SOURCE`/../lib:${LD_LIBRARY_PATH}"
#Tell 3DSTNeutrons where default configuration files are installed
export THREEDSTNEUTRONS_CONF_PATH="`pwd`/`dirname $BASH_SOURCE`/../conf/"
#Make sure bash finds NeutronApp by default
export PATH=$PATH:"`pwd`/`dirname $BASH_SOURCE`/../bin"